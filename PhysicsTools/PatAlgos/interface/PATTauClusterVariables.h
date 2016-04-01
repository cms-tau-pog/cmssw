
/*
 * \class PATTauClusterVariables
 * 
 * Cluster variables to be used in MVA tau ID
 * 
 * Adopted from RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h
 * to enable computation of MVA isolation on MiniAOD
 * 
 * \author Alexander Nehrkorn, RWTH Aachen
 */

namespace
{
	reco::CandidatePtrVector getGammas(const pat::Tau& tau, bool signal = true){
		if(signal)
			return tau.signalGammaCands();
		return tau.isolationGammaCands();
	}
	
	bool isInside(float photon_pt, float deta, float dphi){
	    double stripEtaAssociationDistance_0p95_p0 = 0.197077;
	    double stripEtaAssociationDistance_0p95_p1 = 0.658701;

	    double stripPhiAssociationDistance_0p95_p0 = 0.352476;
	    double stripPhiAssociationDistance_0p95_p1 = 0.707716;
	    
	    if(photon_pt == 0)
	    	return false;
	    
	    if(
	       (dphi < std::min(0.3, std::max(0.05, stripPhiAssociationDistance_0p95_p0*std::pow(photon_pt, -(stripPhiAssociationDistance_0p95_p1))))) && \
	       (deta < std::min(0.15, std::max(0.05, stripEtaAssociationDistance_0p95_p0*std::pow(photon_pt, -(stripEtaAssociationDistance_0p95_p1)))))
	       ){
	      return true;
	    }
	    
	    return false;
	}
	
	float tau_Eratio(const pat::Tau& tau){
		float EcalEnInSignalCands = tau.ecalEnergy();
		float HcalEnInSignalCands = tau.hcalEnergy();
		float total = EcalEnInSignalCands + HcalEnInSignalCands;
		if(total == 0) return -1;
		else return EcalEnInSignalCands/total;
	}
	
	float pt_weighted_dx(const pat::Tau& tau, int mode = 0, int var = 0, int decaymode = -1){
	    float sum_pt = 0.;
	    float sum_dx_pt = 0.;
	    float signalrad = std::max(0.05, std::min(0.1, 3./std::max(1., tau.pt())));
	    int is3prong = (decaymode==10);

	    auto cands = getGammas(tau, mode < 2);

	    for (auto& cand : cands) {
	      // only look at electrons/photons with pT > 0.5
	      if (cand->pt() < 0.5)
	        continue;

	      float dr = reco::deltaR(*cand, tau);
	      float deta = std::abs(cand->eta() - tau.eta());
	      float dphi = std::abs(reco::deltaPhi(cand->phi(), tau.phi()));
	      float pt = cand->pt();

	      bool flag = isInside(pt, deta, dphi);

	      if(is3prong==0){
	        if (mode == 2 || (mode == 0 && dr < signalrad) || (mode == 1 && dr > signalrad)) {
	          sum_pt += pt;
	          if (var == 0)
	            sum_dx_pt += pt * dr;
	          else if (var == 1)
	            sum_dx_pt += pt * deta;
	          else if (var == 2)
	            sum_dx_pt += pt * dphi;
	        }
	      }else if(is3prong==1){

	        if( (mode==2 && flag==false) || (mode==1 && flag==true) || mode==0){
	          sum_pt += pt;
	        
	          if (var == 0)
	            sum_dx_pt += pt * dr;
	          else if (var == 1)
	            sum_dx_pt += pt * deta;
	          else if (var == 2)
	            sum_dx_pt += pt * dphi;
	        }
	      }
	    }
	 
	    if (sum_pt > 0.)
	      return sum_dx_pt/sum_pt;  
	    return 0.;
	}
	
	float tau_pt_weighted_dr_signal(const pat::Tau& tau, int dm){
		return pt_weighted_dx(tau, 0, 0, dm);
	}
	
	float tau_pt_weighted_deta_strip(const pat::Tau& tau, int dm){
		if(dm == 10)
			return pt_weighted_dx(tau, 2, 1, dm);
		return pt_weighted_dx(tau, 1, 1, dm);
	}
	
	float tau_pt_weighted_dphi_strip(const pat::Tau& tau, int dm){
		if(dm == 10)
			return pt_weighted_dx(tau, 2, 2, dm);
		return pt_weighted_dx(tau, 1, 2, dm);
	}
	
	float tau_pt_weighted_dr_iso(const pat::Tau& tau, int dm){
		return pt_weighted_dx(tau, 2, 0, dm);
	}
	
	unsigned int tau_n_photons_total(const pat::Tau& tau){
		unsigned int n_photons = 0;
	    for (auto& cand : tau.signalGammaCands()) {
	      if (cand->pt() > 0.5)
	        ++n_photons;
	    }
	    for (auto& cand : tau.isolationGammaCands()) {
	      if (cand->pt() > 0.5)
	        ++n_photons;
	    }
	    return n_photons;
	}
}
