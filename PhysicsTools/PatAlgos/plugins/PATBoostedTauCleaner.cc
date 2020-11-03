//%%%%%%%%%%%%%%%%%%%%%%%%
// A workaround for fixing the cross-cleaning bug of boostedTau reco at miniAOD level
// Abddollah Mohammadi, UW-Madison
//%%%%%%%%%%%%%%%%%%%%%%%
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "FWCore/Utilities/interface/transform.h"
#include "DataFormats/PatCandidates/interface/TauPFSpecific.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <Math/VectorUtil.h>
#include "TLorentzVector.h"

using namespace std;

class PATBoostedTauCleaner : public edm::stream::EDProducer<>
{
public:
    
    explicit PATBoostedTauCleaner(const edm::ParameterSet&);
    ~PATBoostedTauCleaner(){};
    
    void produce(edm::Event&, const edm::EventSetup&);
    
private:
    
    //--- configuration parameters
    edm::EDGetTokenT<pat::TauCollection> src_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pf2pc_;
    edm::EDGetTokenT<reco::VertexCollection> vtxLabel_;
    edm::EDGetTokenT<edm::View<reco::Jet>> jetsCA8Label_;
    bool  removeOverLap_;
    
};

PATBoostedTauCleaner::PATBoostedTauCleaner(const edm::ParameterSet& cfg)
{
    src_ = consumes<pat::TauCollection>(cfg.getParameter<edm::InputTag>("src"));
    pf2pc_ = consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("pfcands"));
    vtxLabel_ = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vtxLabel"));
    jetsCA8Label_= consumes<edm::View<reco::Jet>>(cfg.getParameter<edm::InputTag>("ca8JetSrc"));
    removeOverLap_ = cfg.getParameter<bool>("removeOverLap");
    produces<std::vector<pat::Tau> >();
}


void PATBoostedTauCleaner::produce(edm::Event& evt, const edm::EventSetup& es)
{
    edm::Handle<pat::TauCollection> inputTaus;
    evt.getByToken(src_, inputTaus);
    
    edm::Handle<pat::PackedCandidateCollection> pf2pc;
    evt.getByToken(pf2pc_, pf2pc);
    
    edm::Handle<reco::VertexCollection> vertices;
    evt.getByToken(vtxLabel_, vertices);
    
    edm::Handle<edm::View<reco::Jet> > ca8jetHandle;
    evt.getByToken(jetsCA8Label_, ca8jetHandle);
    
    auto out = std::make_unique<std::vector<pat::Tau>>();
    out->reserve(inputTaus->size());
    
    for (std::vector<pat::Tau>::const_iterator it = inputTaus->begin(), ed = inputTaus->end(); it != ed; ++it) {
        out->push_back(*it);
        pat::Tau &tau = out->back();
        
        float  chargedPtIsoSum = it->tauID("chargedIsoPtSum");
        float  chargedPtIsoSum03 = it->tauID("chargedIsoPtSumdR03");
        float  neutralPtIsoSum = it->tauID("neutralIsoPtSum");
        float  neutralPtIsoSum03  = it->tauID("neutralIsoPtSumdR03");
        
        reco::CandidatePtrVector  isolationChHPtrs, isolationNHPtrs, isolationGammaPtrs;
        
        
        // Tau vertex
        const pat::PackedCandidate* packedLeadChCand = dynamic_cast<const pat::PackedCandidate*>(tau.leadChargedHadrCand().get());
        
        float minDz = 99;
        int tauVertexIdx = 0;
        int idx = 0;
        for (const auto& vertex : *vertices) {//vertices is handle to vertices
            float dz = std::abs(packedLeadChCand->dz(vertex.position()));
            if (dz < minDz) {
                minDz = dz;
                tauVertexIdx = idx;
            }
            idx++;
        }
        
        //############################################################################
        // filling a collection of IsoCandidates which either overlappes with other tau's signal candidates or other subJet contitients
        //############################################################################
        if (removeOverLap_){
            
            reco::CandidatePtrVector OverLappedIsoCand;
            OverLappedIsoCand.clear();
            
            chargedPtIsoSum = 0;
            chargedPtIsoSum03 = 0;
            neutralPtIsoSum = 0;
            neutralPtIsoSum03  = 0;
            
            for (const reco::CandidatePtr &isoCand1 : tau.isolationCands()) {
                
                // Check iso candidate does not overlap with other signal candidates
                auto out2 = std::make_unique<std::vector<pat::Tau>>();
                out2->reserve(inputTaus->size());
                for (std::vector<pat::Tau>::const_iterator it2 = inputTaus->begin(), ed2 = inputTaus->end(); it2 != ed2; ++it2) {
                    
                    if (it2 == it) continue;
                    
                    out2->push_back(*it2);
                    pat::Tau &tau2 = out2->back();
                    
                    if (ROOT::Math::VectorUtil::DeltaR(tau2.p4(), tau.p4()) > 1.0) continue;
                    
                    for (const reco::CandidatePtr &sigCand2 : tau2.signalCands()) {
                        if (ROOT::Math::VectorUtil::DeltaR(isoCand1->p4(), sigCand2->p4()) < 1e-4)
                            OverLappedIsoCand.push_back(isoCand1);
                    }
                }
                
                // Run on CA8 jets
                for (edm::View<reco::Jet>::const_iterator iJet = ca8jetHandle->begin(); iJet != ca8jetHandle->end(); ++iJet) {
                    
                    if (ROOT::Math::VectorUtil::DeltaR(iJet->p4(), tau.p4()) > 1.0) continue;
                    
                    // Find the subjet that seeds taus : closest subjet to tau
                    float dRClosest=1000;
                    float TauSeedSubJetPt=0;
                    for (edm::View<reco::Jet>::const_iterator iJet = ca8jetHandle->begin(); iJet != ca8jetHandle->end(); ++iJet) {
                        if (iJet->pt() < 14 || fabs(iJet->eta()) > 2.4) continue;
                        if (ROOT::Math::VectorUtil::DeltaR(iJet->p4(), tau.p4()) < dRClosest){
                            dRClosest= ROOT::Math::VectorUtil::DeltaR(iJet->p4(), tau.p4());
                            TauSeedSubJetPt = iJet->pt();
                        }
                    }
                    
                    for (edm::View<reco::Jet>::const_iterator jJet = ca8jetHandle->begin(); jJet != ca8jetHandle->end(); ++jJet) {
                        
                        if (ROOT::Math::VectorUtil::DeltaR(jJet->p4(), tau.p4()) > 1.0) continue;
                        if (fabs (TauSeedSubJetPt - jJet->pt()) < 0.001) continue;
                        
                        for (unsigned id = 0; id < jJet->getJetConstituents().size(); id++) {
                            const edm::Ptr<reco::Candidate> daughter = jJet->getJetConstituents().at(id);
                            if (ROOT::Math::VectorUtil::DeltaR(isoCand1->p4(), daughter->p4()) < 1e-4)
                                OverLappedIsoCand.push_back(isoCand1);
                        }
                    }
                }
            }// end of filling the OverLappedIsoCand collection
            //############################################################################
            // Setting IsolationChargedHadrCands and recalculating ChargedIsoPtSums
            //############################################################################
            
            for (const auto &charged : tau.isolationChargedHadrCands()) {
                
                bool hasOverLap=false;
                for (const reco::CandidatePtr &overLapCand : OverLappedIsoCand) {
                    if (ROOT::Math::VectorUtil::DeltaR(charged->p4(), overLapCand->p4()) < 1e-4){
                        hasOverLap=true;
                        break;
                    }
                }
                if (! hasOverLap){
                    isolationChHPtrs.push_back(charged);
                    
                    //q-cuts my selection
                    if (charged->pt() <= 0.5) continue;
                    //                    if (std::abs(tau.dxy((*vertices)[tauVertexIdx].position())) >= 0.03) continue;
                    const reco::Track *track = charged->bestTrack();
                    if (track == nullptr) continue;
                    if (std::abs(track->dxy((*vertices)[tauVertexIdx].position())) >= 0.03) continue;
                    if (track->normalizedChi2() >= 100) continue;
                    if (track->numberOfValidHits() < 3) continue;
                    //                    double dz = std::abs(tau.dz((*vertices)[tauVertexIdx].position()));
                    double dz = std::abs(track->dz((*vertices)[tauVertexIdx].position()));
                    double dR = deltaR(charged->p4(), tau.p4());
                    if (dz < 0.2) {//from tau vertex
                        //iso cone
                        if (dR < 0.5)
                            chargedPtIsoSum += charged->pt();
                        if (dR < 0.3)
                            chargedPtIsoSum03 += charged->pt();
                    }
                    
                }//check has overlap
            } //loop over all isoCandidates
            tau.setIsolationChargedHadrCands(isolationChHPtrs);
            
            //############################################################################
            // Setting IsolationNeutralHadrCands
            //############################################################################
            
            for (const reco::CandidatePtr &neutral : tau.isolationNeutrHadrCands()) {
                bool hasOverLap=false;
                for (const reco::CandidatePtr &overLapCand : OverLappedIsoCand) {
                    if (ROOT::Math::VectorUtil::DeltaR(neutral->p4(), overLapCand->p4()) < 1e-4){
                        hasOverLap=true;
                        break;
                    }
                }
                if (! hasOverLap)
                    isolationNHPtrs.push_back(neutral);
            }
            tau.setIsolationNeutralHadrCands(isolationNHPtrs);
            
            //############################################################################
            // Setting IsolationGammaCands and recalculating neutralIsoPtSum
            //############################################################################
            
            for (const reco::CandidatePtr &gamma : tau.isolationGammaCands()) {
                bool hasOverLap=false;
                for (const reco::CandidatePtr &overLapCand : OverLappedIsoCand) {
                    if (ROOT::Math::VectorUtil::DeltaR(gamma->p4(), overLapCand->p4()) < 1e-4){
                        hasOverLap=true;
                        break;
                    }
                }
                if (! hasOverLap){
                    isolationGammaPtrs.push_back(gamma);
                    
                    // Fill neutralPtIsoSum03 and neutralPtIsoSum
                    //q-cuts
                    if (gamma->pt() <= 1.) continue;
                    //iso cone
                    double dR = deltaR(gamma->p4(), tau.p4());
                    if (dR < 0.5)
                        neutralPtIsoSum += gamma->pt();
                    if (dR < 0.3)
                        neutralPtIsoSum03 += gamma->pt();
                }
            }
            tau.setIsolationGammaCands(isolationGammaPtrs);
            
            //############################################################################
        }// check if overLap removal is needed
        //############################################################################
        
        size_t newTauIds = 4;
        size_t nTauIds = tau.tauIDs().size();
        std::vector<pat::Tau::IdPair> tauIds(nTauIds+newTauIds);
        
        for(size_t q = 0; q < nTauIds; ++q){
            tauIds[q] = tau.tauIDs().at(q);
        }
        size_t q = nTauIds;
        tauIds[q].first="chargedIsoPtSumNoOverLap";
        tauIds[q].second= chargedPtIsoSum;
        q=q+1;
        
        tauIds[q].first="chargedIsoPtSumdR03NoOverLap";
        tauIds[q].second= chargedPtIsoSum03;
        q=q+1;
        
        tauIds[q].first="neutralIsoPtSumNoOverLap";
        tauIds[q].second= neutralPtIsoSum;
        q=q+1;
        
        tauIds[q].first="neutralIsoPtSumdR03NoOverLap";
        tauIds[q].second= neutralPtIsoSum03;
        q=q+1;
        
        tau.setTauIDs(tauIds);
        
        
    }// iterate over tau is finished
    
    //This is for Adding new tau Id
    evt.put(std::move(out));
}

DEFINE_FWK_MODULE(PATBoostedTauCleaner);
