#!/usr/bin/env python3
"""
Generate and submit HTCondor jobs that produce HLT tau validation DQMIO files
from scratch (GEN,SIM -> DIGI,L1,HLT,VALIDATION), plus a harvesting script
that merges every job's DQMIO output into a single DQM file.

Each Condor job runs step1 (GEN,SIM) followed by step2 (DIGI ... HLT,
VALIDATION:@hltValidation) in the worker's scratch area and copies only the
small step2 DQMIO file back to <outdir>/step2/.  The large GEN-SIM file is
thrown away.  Once all jobs are done, run <outdir>/harvest.sh once: it runs
HARVESTING:@hltValidation over all step2 files together.  Do NOT harvest each
job separately and hadd the results - efficiencies computed at harvesting time
would be summed instead of recomputed.

The cmsDriver configs are generated ONCE on the submit host (with --no_exec),
so any DAS lookup for the pileup input happens there, not on the workers.
Each job gets distinct random seeds (all modules, both steps) and a distinct
luminosity block, derived from its job index, so jobs are reproducible but
statistically independent.

Example:
    cmsenv
    # Copy some MinBias GEN-SIM locally first (see --pileup-input help), then
    submitHltTauValidationCondor.py --outdir /vols/cms/$USER/DYToTauTau_PU200 \\
        --njobs 1000 --events-per-job 10 \\
        --pileup-input filelist:/vols/cms/$USER/minbias/minbias_files.txt

    # Test one job interactively before submitting everything:
    submitHltTauValidationCondor.py --outdir ... --events-per-job 2 --no-submit
    bash <outdir>/run_job.sh 0

    # After all jobs finish:
    bash <outdir>/harvest.sh
"""

import argparse
import os
import shlex
import stat
import subprocess
import sys
import textwrap

DEFAULT_FRAGMENT = "DYToTauTau_M-50_14TeV_pythia8_tauola_cff"
DEFAULT_CONDITIONS = "auto:phase2_realistic_T35"
DEFAULT_GEOMETRY = "ExtendedRun4D128"
DEFAULT_ERA = "Phase2C26I13M9"
DEFAULT_BEAMSPOT = "DBrealisticHLLHC"
DEFAULT_STEP2 = "DIGI:pdigi_valid,L1TrackTrigger,L1,L1P2GT,DIGI2RAW,HLT:@relvalRun4,VALIDATION:@hltValidation"
DEFAULT_STEP2_CUSTOMISE = ",".join([
    "SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000",
    "SimGeneral/TruthGraphAssociatorProducers/customiseTruthGraphAssociators.customiseTruthHltValidation",
])
# The PU200 MinBias sample the release's own Run4D128 relvals use
# (Configuration/PyReleaseValidation/python/relval_steps.py, PUDataSets).
DEFAULT_MINBIAS_DATASET = (
    "/RelValMinBias_14TeV/"
    "CMSSW_20_1_0_pre3-150X_mcRun4_realistic_v1_STD_D128_GSOnly_GSInputs-v1/GEN-SIM"
)

# Python appended to each generated cfg.  JOB_INDEX is exported by run_job.sh.
# A different salt per step keeps step1/step2 seeds independent.
# (Appended directly rather than via --customise_commands, because cmsDriver
# exec()s those line by line and also copies them into the cfg header comment.)
SEED_CUSTOMISATION = """
# ---- Added by submitHltTauValidationCondor.py: per-job random seeds ----
import os as _os, random as _random
from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper as _RNSH
_job = int(_os.environ["JOB_INDEX"])
_helper = _RNSH(process.RandomNumberGeneratorService)
_rng = _random.Random("{salt}-%d-%d" % ({seed_offset}, _job))
_helper.insertSeeds(*[_rng.randint(1, 900000000) for _ in range(_helper.countSeeds())])
print("JOB_INDEX=%d: reseeded %d random-number engines" % (_job, _helper.countSeeds()))
"""

LUMI_CUSTOMISATION = """
process.source.firstLuminosityBlock = cms.untracked.uint32(_job + 1)
"""


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", required=True,
                   help="Output directory (should be on shared storage, e.g. /vols/cms/...). "
                        "Holds cfgs, job scripts, logs, step2 DQMIO files and the harvested output.")
    p.add_argument("--njobs", type=int, default=1000, help="Number of Condor jobs (default: %(default)s)")
    p.add_argument("--events-per-job", type=int, default=10, help="Events per job (default: %(default)s)")
    p.add_argument("--fragment", default=DEFAULT_FRAGMENT, help="Generator fragment (default: %(default)s)")
    p.add_argument("--pileup", default="AVE_200_BX_25ns",
                   help="Mixing scenario, or 'none' for no pileup (default: %(default)s)")
    p.add_argument("--pileup-input", default=None,
                   help="MinBias GEN-SIM input for mixing, in cmsDriver --pileup_input syntax. "
                        "Recommended: 'filelist:/path/to/list.txt' with one 'file:/abs/path.root' per line, "
                        "on storage the workers can read. 'das:<dataset>' also works but needs a valid "
                        "grid proxy when this script runs, and the workers then stream over xrootd. "
                        f"The release's own D128 PU200 sample is das:{DEFAULT_MINBIAS_DATASET}")
    p.add_argument("--nthreads", type=int, default=1,
                   help="cmsRun threads; also sets request_cpus (default: %(default)s)")
    p.add_argument("--memory", default="12000",
                   help="request_memory in MB. PU200 step2 peaks at ~13-16 GB with 4 threads (it scales "
                        "with the number of events in flight, not events per job); jobs exceeding the "
                        "request are held (default: %(default)s)")
    p.add_argument("--max-runtime", type=int, default=3600,
                   help="+MaxRuntime in seconds (default: %(default)s)")
    p.add_argument("--seed-offset", type=int, default=0,
                   help="Change this to get a statistically independent production with the same job "
                        "indices (default: %(default)s)")
    p.add_argument("--conditions", default=DEFAULT_CONDITIONS)
    p.add_argument("--geometry", default=DEFAULT_GEOMETRY)
    p.add_argument("--era", default=DEFAULT_ERA)
    p.add_argument("--beamspot", default=DEFAULT_BEAMSPOT)
    p.add_argument("--step2", default=DEFAULT_STEP2, help="step2 -s argument (default: %(default)s)")
    p.add_argument("--step2-customise", default=DEFAULT_STEP2_CUSTOMISE,
                   help="step2 --customise argument (default: %(default)s)")
    p.add_argument("--no-submit", action="store_true",
                   help="Write everything but don't call condor_submit")
    p.add_argument("--force", action="store_true",
                   help="Allow reusing an outdir that already contains a submit file")
    return p.parse_args()


def run(cmd, cwd):
    print("+", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def common_opts(a):
    return ["--conditions", a.conditions, "--geometry", a.geometry, "--era", a.era]


def make_configs(a, cfgdir):
    seeds1 = SEED_CUSTOMISATION.format(salt="step1", seed_offset=a.seed_offset) + LUMI_CUSTOMISATION
    seeds2 = SEED_CUSTOMISATION.format(salt="step2", seed_offset=a.seed_offset)

    step1 = ["cmsDriver.py", a.fragment, "-s", "GEN,SIM", "-n", str(a.events_per_job),
             *common_opts(a), "--beamspot", a.beamspot,
             "--datatier", "GEN-SIM", "--eventcontent", "FEVTDEBUG",
             "--nThreads", str(a.nthreads),
             "--fileout", "file:step1.root",
             "--python_filename", "step1_cfg.py", "--no_exec"]

    step2 = ["cmsDriver.py", "step2", "-s", a.step2, "-n", "-1",
             *common_opts(a),
             "--datatier", "DQMIO", "--eventcontent", "DQMIO",
             "--nThreads", str(a.nthreads),
             "--filein", "file:step1.root", "--fileout", "file:step2.root",
             "--python_filename", "step2_cfg.py", "--no_exec"]
    if a.step2_customise:
        step2 += ["--customise", a.step2_customise]
    if a.pileup.lower() != "none":
        step2 += ["--pileup", a.pileup, "--pileup_input", a.pileup_input]

    run(step1, cfgdir)
    run(step2, cfgdir)
    for name, snippet in (("step1_cfg.py", seeds1), ("step2_cfg.py", seeds2)):
        with open(os.path.join(cfgdir, name), "a") as f:
            f.write(snippet)


def write_executable(path, content):
    with open(path, "w") as f:
        f.write(content)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_job_script(a, outdir, cfgdir, step2dir, cmssw_base, scram_arch):
    write_executable(os.path.join(outdir, "run_job.sh"), textwrap.dedent(f"""\
        #!/bin/bash
        # Usage: run_job.sh <job index>
        set -euo pipefail
        export JOB_INDEX=$1
        OUT={step2dir}/step2_${{JOB_INDEX}}.root
        if [ -s "$OUT" ]; then
            echo "$OUT already exists, nothing to do"
            exit 0
        fi

        export SCRAM_ARCH={scram_arch}
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        cd {cmssw_base}/src && eval "$(scramv1 runtime -sh)" && cd - > /dev/null

        WORK=$(mktemp -d "${{_CONDOR_SCRATCH_DIR:-${{TMPDIR:-/tmp}}}}/hlttau_${{JOB_INDEX}}_XXXX")
        trap 'rm -rf "$WORK"' EXIT
        cd "$WORK"
        echo "Job $JOB_INDEX on $(hostname) in $WORK, started $(date)"

        cmsRun {cfgdir}/step1_cfg.py
        echo "step1 done $(date)"
        cmsRun {cfgdir}/step2_cfg.py
        echo "step2 done $(date)"

        cp step2.root "$OUT.tmp"
        mv "$OUT.tmp" "$OUT"
        echo "Wrote $OUT"
        """))


def write_submit_file(a, outdir, logdir):
    path = os.path.join(outdir, "submit.sub")
    with open(path, "w") as f:
        f.write(textwrap.dedent(f"""\
            executable     = {outdir}/run_job.sh
            arguments      = $(ProcId)
            output         = {logdir}/job_$(ProcId).out
            error          = {logdir}/job_$(ProcId).err
            log            = {logdir}/condor_$(ClusterId).log

            getenv         = False
            request_cpus   = {a.nthreads}
            request_memory = {a.memory}M
            +MaxRuntime    = {a.max_runtime}

            queue {a.njobs}
            """))
    return path


def write_harvest_script(a, outdir, step2dir, cmssw_base, scram_arch):
    write_executable(os.path.join(outdir, "harvest.sh"), textwrap.dedent(f"""\
        #!/bin/bash
        # Harvest ALL step2 DQMIO files in one go -> single DQM_V0001_*.root in {outdir}/harvest
        set -euo pipefail
        export SCRAM_ARCH={scram_arch}
        source /cvmfs/cms.cern.ch/cmsset_default.sh
        cd {cmssw_base}/src && eval "$(scramv1 runtime -sh)" && cd - > /dev/null

        mkdir -p {outdir}/harvest
        cd {outdir}/harvest
        find {step2dir} -name 'step2_*.root' -size +0 | sort -V | sed 's|^|file:|' > step2_files.txt
        echo "Harvesting $(wc -l < step2_files.txt) / {a.njobs} step2 files"

        cmsDriver.py step3 -s HARVESTING:@hltValidation -n -1 \\
            --conditions {a.conditions} --geometry {a.geometry} --era {a.era} \\
            --mc --scenario pp --filetype DQM \\
            --filein filelist:step2_files.txt --python_filename step3_cfg.py
        ls -lh DQM_V0001_*.root
        """))


def main():
    a = parse_args()

    cmssw_base = os.environ.get("CMSSW_BASE")
    scram_arch = os.environ.get("SCRAM_ARCH")
    if not cmssw_base or not scram_arch:
        sys.exit("ERROR: run cmsenv first.")
    if a.pileup.lower() != "none" and not a.pileup_input:
        sys.exit("ERROR: --pileup-input is required when pileup is enabled (or pass --pileup none).\n"
                 f"  e.g. --pileup-input filelist:/vols/cms/$USER/minbias/minbias_files.txt\n"
                 f"  or   --pileup-input das:{DEFAULT_MINBIAS_DATASET}  (needs a grid proxy)")
    if a.pileup_input and a.pileup_input.startswith("filelist:"):
        flist = a.pileup_input.split(":", 1)[1]
        if not os.path.isabs(flist):
            sys.exit("ERROR: give an absolute path for the pileup filelist (the workers read it via the cfg).")
        if not os.path.isfile(flist):
            sys.exit(f"ERROR: pileup filelist {flist} not found.")

    outdir = os.path.abspath(a.outdir)
    cfgdir = os.path.join(outdir, "cfg")
    logdir = os.path.join(outdir, "logs")
    step2dir = os.path.join(outdir, "step2")
    if os.path.exists(os.path.join(outdir, "submit.sub")) and not a.force:
        sys.exit(f"ERROR: {outdir} already has a submit.sub; use --force to overwrite "
                 "(existing step2 outputs are kept and their jobs skipped).")
    for d in (cfgdir, logdir, step2dir):
        os.makedirs(d, exist_ok=True)

    make_configs(a, cfgdir)
    write_job_script(a, outdir, cfgdir, step2dir, cmssw_base, scram_arch)
    sub = write_submit_file(a, outdir, logdir)
    write_harvest_script(a, outdir, step2dir, cmssw_base, scram_arch)

    print(f"\nWrote cfgs to {cfgdir}, job script, {sub} and {outdir}/harvest.sh")
    print(f"{a.njobs} jobs x {a.events_per_job} events = {a.njobs * a.events_per_job} events")
    if a.no_submit:
        print(f"Not submitting. Test with:  bash {outdir}/run_job.sh 0\n"
              f"Submit with:               condor_submit {sub}")
        return
    run(["condor_submit", sub], outdir)
    print(f"\nWhen all jobs are done:  bash {outdir}/harvest.sh")


if __name__ == "__main__":
    main()
