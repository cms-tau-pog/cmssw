// -*- C++ -*-
//
// Package:    LowPtTauProducer/lowPtTaus
// Class:      lowPtTaus
// 
/**\class lowPtTaus lowPtTaus.cc LowPtTauProducer/lowPtTaus/plugins/lowPtTaus.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yuta Takahashi (UniZ) [ytakahas]
//         Created:  Mon, 24 Jun 2019 08:26:02 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/lowPtTau/interface/lowPtTau.h"






// This is neede for SV
//#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
//#include "SimDataFormats/JetMatching/interface/JetFlavour.h"

//#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"

//#include "JetMETCorrections/Modules/interface/JetResolution.h"
//#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>


// for vertexing
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"                                                                                             
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
     
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"                                                                                                       
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h" 
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/MuonReco/interface/MuonQuality.h"
             
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"                                                                                                 
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include <tuple>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm> 
#include <stdio.h>
#include <stdlib.h>
#include <map>

#include "TLorentzVector.h"








class lowPtTaus : public edm::stream::EDProducer<> {
public:
  explicit lowPtTaus(const edm::ParameterSet&);
  ~lowPtTaus();  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  std::tuple<Float_t, TransientVertex> vertexProb( const std::vector<reco::TransientTrack>& tracks);
  
private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  edm::Handle< reco::VertexCollection >  vertices_;
  edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
  edm::Handle< reco::GenParticleCollection >  genParticles_;
  
  edm::EDGetTokenT<reco::VertexCollection>                  vtxToken_           ;
  edm::EDGetTokenT<pat::PackedCandidateCollection>          packedpfcandidatesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>             genparticleToken_   ;  

  edm::ESHandle<TransientTrackBuilder> builder;  
  // ----------member data ---------------------------
  Int_t m_abstaucharge; 
  Bool_t m_isMC;

};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
lowPtTaus::lowPtTaus(const edm::ParameterSet& iConfig)
{

  vtxToken_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  packedpfcandidatesToken_ = consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedpfcandidates"));
  genparticleToken_  = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"));

  m_abstaucharge   = iConfig.getParameter<int>("abscharge");
  m_isMC   = iConfig.getParameter<bool>("isMC");

  std::cout << "tau charge = " << m_abstaucharge << std::endl;
  std::cout << "isMC = " << m_isMC << std::endl;
  
  //register your products
  produces<lowPtTauCollection>().setBranchAlias("lowPtTauCollection");



/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed
  
}


lowPtTaus::~lowPtTaus()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}




std::tuple<Float_t, TransientVertex> lowPtTaus::vertexProb( const std::vector<reco::TransientTrack>& tracks){

  Float_t vprob = -1;
  
  KalmanVertexFitter kalman_fitter;
  TransientVertex vertex;

  try{
    vertex = kalman_fitter.vertex(tracks);
  }catch(std::exception e){
    //    std::cout << "No vertex found ... return" << std::endl;
    return std::forward_as_tuple(-9, vertex);
  }

  if(vertex.isValid()){

    vprob =  TMath::Prob(vertex.totalChiSquared(), vertex.degreesOfFreedom());

    //    if(vprob==0) vprob = -1;
    //    vx = vertex.position().x();
    //    vy = vertex.position().y();
    //    vz = vertex.position().z();

    return std::forward_as_tuple(vprob, vertex);
  }else{

    //    std::cout << "No vertex found !!" << std::endl;
    return std::forward_as_tuple(-9, vertex);
  }
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
lowPtTaus::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  
  using namespace edm;
  
  ////////////////////////////////////////////////
  
  iEvent.getByToken( vtxToken_   , vertices_     );
  iEvent.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 

  if(m_isMC){
    iEvent.getByToken(genparticleToken_ , genParticles_); 
  }

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);  


  if(!vertices_.isValid()){
    std::cout << "Vertex is not valid !!!!!" << std::endl;
    return;
  }

  reco::VertexCollection::const_iterator firstGoodVertex = vertices_->begin();

  Float_t pvx = firstGoodVertex->position().x();
  Float_t pvy = firstGoodVertex->position().y();
  Float_t pvz = firstGoodVertex->position().z();


  /* Filter out PF candidates */

  std::vector<pat::PackedCandidate> pfcollection; 
  std::vector<Int_t> pfidcollection; 
  pfcollection.clear();
  pfidcollection.clear();


  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
    
    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
    //    float pf_pt = pf.pt();   
    //    std::cout << pf_pt << std::endl;
    
    if(pf.pt() < 0.5) continue;
    
    if(!pf.hasTrackDetails()) continue;
    
    Float_t precut_dz = pvz - pf.vz();
    if(TMath::Abs(precut_dz) > 0.5) continue;
    
    Bool_t hpflag = pf.trackHighPurity();
    if(!hpflag) continue;
    
    if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
    if(pf.pseudoTrack().normalizedChi2() > 100) continue;
    
    if(TMath::Abs(pf.pdgId())==211 && TMath::Abs(pf.eta()) < 2.3){
      pfcollection.push_back(pf);
      pfidcollection.push_back(ii);
    }
  }
    

  std::cout << "\t # of PF candidates = " << pfcollection.size() << std::endl;

  Int_t numOfch = (int)pfcollection.size();


  /* gen collection */

  std::vector<std::vector<TLorentzVector>> gps;
  std::vector<Int_t> ppdgId;
  //  Bool_t isgen3 = false;
  //  Bool_t isgen3matched = true;


  if(m_isMC){

    for( unsigned p=0; p < genParticles_->size(); ++p){
      
      if(TMath::Abs((*genParticles_)[p].pdgId())!=15) continue;
      if(TMath::Abs((*genParticles_)[p].status())!=2) continue;
      
      std::cout << "\t Tau found with # of daughters = " << (*genParticles_)[p].numberOfDaughters() << " with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
      
      std::vector<TLorentzVector> gp;
      //      Bool_t matched = true;

      for( unsigned int d=0; d < (*genParticles_)[p].numberOfDaughters(); ++d ){
	
	Int_t taupdgId = (*genParticles_)[p].daughter(d)->pdgId();
	
	std::cout << "\t\t --> tau decay:" << taupdgId << std::endl;
	
	if(TMath::Abs(taupdgId)!=211) continue;
	  
	TLorentzVector tlv_gen_pion;
	
	tlv_gen_pion.SetPtEtaPhiM((*genParticles_)[p].daughter(d)->pt(),
				  (*genParticles_)[p].daughter(d)->eta(),
				  (*genParticles_)[p].daughter(d)->phi(),
				  (*genParticles_)[p].daughter(d)->mass());
	
	// check matching to reco PF objects
	Float_t min_dr = 999;
	
	for(int kkk = 0; kkk < numOfch; kkk ++){
	  
	  pat::PackedCandidate _pf = pfcollection[kkk];
	  
	  if(_pf.pdgId()!=taupdgId) continue;
	  
	  Float_t _dR = reco::deltaR(
				     tlv_gen_pion.Eta(), tlv_gen_pion.Phi(),
				     _pf.eta(), _pf.phi()
				     );
	  if(_dR < min_dr && _dR < 0.1){
	    min_dr = _dR;
	  }
	}
	
	//////////////////////////////////////
	//	if(min_dr == 999) matched = false;
	
	gp.push_back(tlv_gen_pion);
	
      }

      
      if(gp.size()==3){
	//	std::cout << "\t -----> This has been registered with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
	gps.push_back(gp);
	ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
	if(TMath::Abs((*genParticles_)[p].mother(0)->pdgId())==541){
	  //	  isgen3 = true;
	  //	  isgen3matched = matched;
	}
      }else{
	//	isgen3matched = false;
      }
    }

    std::cout << "\t # of gen. taus with 3prong = " << gps.size() << std::endl;

  }

  std::vector<lowPtTau> cands;

  for(int iii = 0; iii < numOfch; iii ++){

    pat::PackedCandidate pf1 = pfcollection[iii];
    const reco::Track pf1_track = pf1.pseudoTrack();

    for(int jjj = iii+1; jjj < numOfch; jjj ++){
      
      pat::PackedCandidate pf2 = pfcollection[jjj];
      const reco::Track pf2_track = pf2.pseudoTrack();
//      std::vector<reco::TransientTrack> transient_pairs; 
//      transient_pairs.push_back((*builder).build(pf1_track));
//      transient_pairs.push_back((*builder).build(pf2_track));
//
//      Float_t vprob_2;
//      TransientVertex vertex_2;
//      
//      std::tie(vprob_2, vertex_2) = vertexProb(transient_pairs);
//      if(vprob_2 < 0.01) continue;

//      if(jjj==numOfch-1) break;

      for(int kkk = jjj+1; kkk < numOfch; kkk ++){

	pat::PackedCandidate pf3 = pfcollection[kkk];
	const reco::Track pf3_track = pf3.pseudoTrack();
	Int_t tau_charge = pf1.charge() + pf2.charge() + pf3.charge(); 

	if(TMath::Abs(tau_charge) != m_abstaucharge) continue; 

	std::vector<reco::TransientTrack> transient_tracks; 
	transient_tracks.push_back((*builder).build(pf1_track));
	transient_tracks.push_back((*builder).build(pf2_track));
	transient_tracks.push_back((*builder).build(pf3_track));

	Float_t vprob_3 = -9;
	TransientVertex vertex_3;
	
	std::tie(vprob_3, vertex_3) = vertexProb(transient_tracks);
	if(vprob_3 == -9) continue; 	
	TLorentzVector tlv_pion1; 
	TLorentzVector tlv_pion2;
	TLorentzVector tlv_pion3;

	tlv_pion1.SetPtEtaPhiM(pf1.pt(), pf1.eta(), pf1.phi(), pf1.mass());
	tlv_pion2.SetPtEtaPhiM(pf2.pt(), pf2.eta(), pf2.phi(), pf2.mass());
	tlv_pion3.SetPtEtaPhiM(pf3.pt(), pf3.eta(), pf3.phi(), pf3.mass());

	TLorentzVector tlv_tau = tlv_pion1 + tlv_pion2 + tlv_pion3;

	//	if(tlv_tau.Pt() < 3) continue;
	
	//	Float_t taumass = tlv_tau.M();
	//	if(!(taumass > 0.2 && taumass < 1.5)) continue;

	//	std::cout << "iii, jjj, kkk = " << pfidcollection[iii] << " " << pfidcollection[jjj] << " " << pfidcollection[kkk] << std::endl;



	GlobalVector direction(tlv_tau.Px(), tlv_tau.Py(), tlv_tau.Pz()); //To compute sign of IP

	double trivertex_z = -999;

	double flightSig3D = -999; 
	double flightLength3D =  -999; 
	double flightLengthErr3D = -999; 

	double flightSig2D =  -999; 
	double flightLength2D = -999; 
	double flightLengthErr2D = -999; 

	//	if(vprob_3!=-9){
	trivertex_z = vertex_3.position().z();

	flightSig3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex, vertex_3, direction, true).significance();
	flightLength3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex, vertex_3, direction, true).value();
	flightLengthErr3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex, vertex_3, direction, true).error();
	flightSig2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex, vertex_3, direction, true).significance();
	flightLength2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex, vertex_3, direction, true).value();
	flightLengthErr2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex, vertex_3, direction, true).error();
	//	}

	//	double signed_IP2D = IPTools::signedTransverseImpactParameter(transient_tracks, direction, pv).second.value();


	Float_t max_dr_3prong = -1;
	Float_t dR_12 = tlv_pion1.DeltaR(tlv_pion2);
	Float_t dR_13 = tlv_pion1.DeltaR(tlv_pion3);
	Float_t dR_23 = tlv_pion2.DeltaR(tlv_pion3);

	if(max_dr_3prong < dR_12) max_dr_3prong = dR_12;
	if(max_dr_3prong < dR_13) max_dr_3prong = dR_13;
	if(max_dr_3prong < dR_23) max_dr_3prong = dR_23;


	Float_t isopt03 = 0.;
	Float_t isopt04 = 0.;

	for(int lll = 0; lll < numOfch; lll ++){
	  
	  if(lll==iii) continue;
	  if(lll==jjj) continue;
	  if(lll==kkk) continue;

	  pat::PackedCandidate pf_iso = pfcollection[lll];

	  if(TMath::Abs(pf_iso.pdgId())!=211) continue;
	  
	  TLorentzVector tlv_iso;
	  
	  tlv_iso.SetPtEtaPhiM(pf_iso.pt(), pf_iso.eta(), pf_iso.phi(), pf_iso.mass());

	  Float_t dR_iso = tlv_tau.DeltaR(tlv_iso);

	  if(dR_iso > 0.4) continue;

	  if(dR_iso < 0.4){
	    isopt04 += tlv_iso.Pt();
	  }
	  
	  if(dR_iso < 0.3){
	    isopt03 += tlv_iso.Pt();
	  }
	  

	}


	// check if it matched or not ... 

	Bool_t isRight = false; 
	Int_t pid = -999;

	if(m_isMC){

	  for(unsigned int mmm=0; mmm < gps.size(); mmm++){
	    
	    Bool_t isRight1_ = false;
	    Bool_t isRight2_ = false;
	    Bool_t isRight3_ = false;
	    
	    std::vector<TLorentzVector> tlvs = gps[mmm];
	    
	    for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){
	      
	      if(tlv_pion1.DeltaR(tlvs[nnn]) < 0.1) isRight1_ = true;
	      else if(tlv_pion2.DeltaR(tlvs[nnn]) < 0.1) isRight2_ = true;
	      else if(tlv_pion3.DeltaR(tlvs[nnn]) < 0.1) isRight3_ = true;
	      
	    }
	    
	    Bool_t isRight_ = isRight1_ && isRight2_ && isRight3_;
	    
	    if(isRight_){
	      isRight = true;
	      pid = ppdgId[mmm];
	    }
	  }	
	}


	Float_t dz2pv = vertex_3.position().z() - pvz;
	Float_t dxy2pv = TMath::Sqrt(  TMath::Power( (vertex_3.position().x() - pvx), 2) +  TMath::Power( (vertex_3.position().y() - pvy), 2) );

	if(!(vprob_3 > 0.1 && flightSig3D > 3)) continue; // vprob_fsig

	lowPtTau _cand_ = {
	  
	  (Float_t) tlv_tau.Pt(), 
	  (Float_t) tlv_tau.Eta(), 
	  (Float_t) tlv_tau.Phi(), 
	  (Float_t) tlv_tau.M(), 
	  (Int_t) tau_charge,
	  (Float_t) vprob_3,
	  (Float_t) trivertex_z,
	  (Float_t) dz2pv,
	  (Float_t) dxy2pv,
	  (Float_t) flightSig3D,
	  (Float_t) flightLength3D,
	  (Float_t) flightLengthErr3D,
	  (Float_t) flightSig2D,
	  (Float_t) flightLength2D,
	  (Float_t) flightLengthErr2D,
	  iii, 
	  jjj,
	  kkk,
	  (Float_t) tlv_pion1.Pt(),
	  (Float_t) tlv_pion2.Pt(),
	  (Float_t) tlv_pion3.Pt(),	
	  (Float_t) tlv_pion1.Eta(),
	  (Float_t) tlv_pion2.Eta(),
	  (Float_t) tlv_pion3.Eta(),	
	  (Float_t) tlv_pion1.Phi(),
	  (Float_t) tlv_pion2.Phi(),
	  (Float_t) tlv_pion3.Phi(),	
	  (Float_t) tlv_pion1.M(),
	  (Float_t) tlv_pion2.M(),
	  (Float_t) tlv_pion3.M(),	
	  (Float_t) isopt03,
	  (Float_t) isopt04,
	  (Float_t) max_dr_3prong,	
	  (Bool_t) isRight,
	  (Float_t) 1.,
	  (Float_t) 1.,
	  (Float_t) 1.,
	  (Int_t) pid,
	  
	};

	cands.push_back(_cand_);




	//	if(tlv_all.M() < 4.) continue;
	//	if(tlv_all.M() > 5.5) continue;

	//	if(flightSig3D_all < 3) continue;

	//	if(tlv_tau.Pt() < 2.) continue;
	//	if(tlv_all.Pt() < 10.) continue;

	//	if(tlv_all.Pt()  > max_B_pt && tlv_all.Pt() < 10.){

//	if(tlv_tau.Pt()  > max_tau_pt){
//	  
//	  max_tau_pt = tlv_tau.Pt();
//
//	  pft_tau_vprob = vprob_3;
//	  pft_tau_vz = trivertex_z;
//	  pft_tau_flightSig3D = flightSig3D;
//	  pft_tau_flightLength3D = flightLength3D;
//	  pft_tau_flightLengthErr3D = flightLengthErr3D;
//	  pft_tau_flightSig2D = flightSig2D;
//	  pft_tau_flightLength2D = flightLength2D;
//	  pft_tau_flightLengthErr2D = flightLengthErr2D;
//	  pft_tau_pt = tlv_tau.Pt();
//
//	  pft_tau_pt1 = tlv_pion1.Pt();
//	  pft_tau_pt2 = tlv_pion2.Pt();
//	  pft_tau_pt3 = tlv_pion3.Pt();
//
//	  pft_tau_eta = tlv_tau.Eta();
//	  pft_tau_phi = tlv_tau.Phi();
//	  pft_tau_mass = tlv_tau.M();
//
//	  pft_tau_mass12 = (tlv_pion1 + tlv_pion2).M();
//	  pft_tau_mass13 = (tlv_pion1 + tlv_pion3).M();
//	  pft_tau_mass23 = (tlv_pion2 + tlv_pion3).M();
//
//	  pft_tau_iso03 = isopt03;
//	  pft_tau_iso04 = isopt04;
//	  pft_tau_max_dr_3prong = max_dr_3prong;
//	  
//	  pft_b_vprob = vprob_5;
//	  pft_b_vz = pentavertex_z;
//	  
//	  pft_b_pt = tlv_all.Pt();
//	  pft_b_eta = tlv_all.Eta();
//	  pft_b_phi = tlv_all.Phi();
//	  pft_b_mass = tlv_all.M();
//
//	  pft_tau_dphi2met = reco::deltaPhi(tlv_tau.Phi(), metphi);
//	  
//	  //	std::cout << "Inside loop" << _pft_b_pt.size() << " " << _pft_b_eta.size() << " " << _pft_b_phi.size() << " " << _pft_b_mass.size() <<  std::endl;
//	  
//	  pft_b_flightSig3D = flightSig3D_all;
//	  pft_b_flightLength3D = flightLength3D_all;
//	  pft_b_flightLengthErr3D = flightLengthErr3D_all;
//	  pft_b_flightSig2D = flightSig2D_all;
//	  pft_b_flightLength2D = flightLength2D_all;
//	  pft_b_flightLengthErr2D = flightLengthErr2D_all;
//	  
//	}

//	ncomb += 1;
	
      }
    }
  }
  
  sort(cands.begin(), cands.end());


  //  std::vector<Float_t> _pft_pvz;

//  std::vector<Float_t> _pft_tau_vprob;
//  std::vector<Float_t> _pft_tau_vz;
//  std::vector<Float_t> _pft_tau_dz2pv;
//  std::vector<Float_t> _pft_tau_dxy2pv;
//  std::vector<Float_t> _pft_tau_pt;
//  std::vector<Float_t> _pft_tau_pt1;
//  std::vector<Float_t> _pft_tau_pt2;
//  std::vector<Float_t> _pft_tau_pt3;
//  std::vector<Float_t> _pft_tau_eta; 
//  std::vector<Float_t> _pft_tau_phi;
//  std::vector<Float_t> _pft_tau_mass;
//  std::vector<Float_t> _pft_tau_mass12;
//  std::vector<Float_t> _pft_tau_mass13;
//  std::vector<Float_t> _pft_tau_mass23;
//  std::vector<Float_t> _pft_tau_flightSig3D; 
//  std::vector<Float_t> _pft_tau_flightLength3D;
//  std::vector<Float_t> _pft_tau_flightLengthErr3D;
//  std::vector<Float_t> _pft_tau_flightSig2D; 
//  std::vector<Float_t> _pft_tau_flightLength2D;
//  std::vector<Float_t> _pft_tau_flightLengthErr2D;
//  std::vector<Float_t> _pft_tau_iso03; 
//  std::vector<Float_t> _pft_tau_iso04; 
//  std::vector<Float_t> _pft_tau_max_dr_3prong;
//
//  std::vector<Float_t> _pft_b_vprob;
//  std::vector<Float_t> _pft_b_vz;
//  std::vector<Float_t> _pft_b_pt;
//  std::vector<Float_t> _pft_b_eta;
//  std::vector<Float_t> _pft_b_phi;
//  std::vector<Float_t> _pft_b_mass;
//  std::vector<Float_t> _pft_b_flightSig3D; 
//  std::vector<Float_t> _pft_b_flightLength3D;
//  std::vector<Float_t> _pft_b_flightLengthErr3D;
//  std::vector<Float_t> _pft_b_flightSig2D; 
//  std::vector<Float_t> _pft_b_flightLength2D;
//  std::vector<Float_t> _pft_b_flightLengthErr2D;

  std::vector<Int_t> dict_idx;
  //  Int_t ncomb = 0;


  std::unique_ptr<lowPtTauCollection> result(new lowPtTauCollection);
  
  for(int ic=0; ic < (int)cands.size(); ic++){

    Int_t _idx1 = cands[ic].pfid1;
    Int_t _idx2 = cands[ic].pfid2;
    Int_t _idx3 = cands[ic].pfid3;

    bool flag_overlap = false;
    for(int idc=0; idc<(int) dict_idx.size(); idc++){
      
      if(_idx1 == dict_idx[idc] || 
	 _idx2 == dict_idx[idc] || 
	 _idx3 == dict_idx[idc])
	 
	flag_overlap = true;
    }

    if(flag_overlap) continue; 

    //    std::cout << "sorted:" << ic << " " << cands[ic].cand_tau_pt << " " << cands[ic].cand_tau_id1 << " " << cands[ic].cand_tau_id2 << " " << cands[ic].cand_tau_id3 << " " << cands[ic].cand_tau_vprob << " " << cands[ic].cand_tau_isRight << std::endl;

    //    std::cout << "\t ----> kept ! " << std::endl;
    dict_idx.push_back(_idx1);
    dict_idx.push_back(_idx2);
    dict_idx.push_back(_idx3);

    //    ncomb += 1;

    result->push_back(cands[ic]);


  }


    ///////// filling ... 
    
//    nBranches_->pft_tau_vprob.push_back(cands[ic].cand_tau_vprob);
//    nBranches_->pft_tau_vz.push_back(cands[ic].cand_tau_vz);
//    nBranches_->pft_tau_dz2pv.push_back(cands[ic].cand_tau_dz2pv);
//    nBranches_->pft_tau_dxy2pv.push_back(cands[ic].cand_tau_dxy2pv);
//    nBranches_->pft_tau_flightSig3D.push_back(cands[ic].cand_tau_flightSig3D);
//    nBranches_->pft_tau_flightLength3D.push_back(cands[ic].cand_tau_flightLength3D);
//    nBranches_->pft_tau_flightLengthErr3D.push_back(cands[ic].cand_tau_flightLengthErr3D);
//    nBranches_->pft_tau_flightSig2D.push_back(cands[ic].cand_tau_flightSig2D);
//    nBranches_->pft_tau_flightLength2D.push_back(cands[ic].cand_tau_flightLength2D);
//    nBranches_->pft_tau_flightLengthErr2D.push_back(cands[ic].cand_tau_flightLengthErr2D);
//    nBranches_->pft_tau_pt.push_back(cands[ic].cand_tau_pt);
//    nBranches_->pft_tau_pt1.push_back(cands[ic].cand_tau_pt1);
//    nBranches_->pft_tau_pt2.push_back(cands[ic].cand_tau_pt2);
//    nBranches_->pft_tau_pt3.push_back(cands[ic].cand_tau_pt3);
//    nBranches_->pft_tau_eta.push_back(cands[ic].cand_tau_eta);
//    nBranches_->pft_tau_phi.push_back(cands[ic].cand_tau_phi);
//    nBranches_->pft_tau_mass.push_back(cands[ic].cand_tau_mass);
//    nBranches_->pft_tau_mass12.push_back(cands[ic].cand_tau_mass12);
//    nBranches_->pft_tau_mass13.push_back(cands[ic].cand_tau_mass13);
//    nBranches_->pft_tau_mass23.push_back(cands[ic].cand_tau_mass23);
//    nBranches_->pft_tau_iso03.push_back(cands[ic].cand_tau_isopt03);
//    nBranches_->pft_tau_iso04.push_back(cands[ic].cand_tau_isopt04);
//    nBranches_->pft_tau_max_dr_3prong.push_back(cands[ic].cand_tau_max_dr_3prong);	
//    nBranches_->pft_tau_charge.push_back(cands[ic].cand_tau_charge);
//    nBranches_->pft_tau_isRight.push_back(cands[ic].cand_tau_isRight);
//    nBranches_->pft_tau_ppdgId.push_back(cands[ic].cand_tau_ppdgId);
//
//    nBranches_->pft_b_vprob.push_back(cands[ic].cand_b_vprob);
//    nBranches_->pft_b_vz.push_back(cands[ic].cand_b_vz);
//    nBranches_->pft_b_pt.push_back(cands[ic].cand_b_pt);
//    nBranches_->pft_b_eta.push_back(cands[ic].cand_b_eta);
//    nBranches_->pft_b_phi.push_back(cands[ic].cand_b_phi);
//    nBranches_->pft_b_mass.push_back(cands[ic].cand_b_mass);
//    nBranches_->pft_b_flightSig3D.push_back(cands[ic].cand_b_flightSig3D_all);
//    nBranches_->pft_b_flightLength3D.push_back(cands[ic].cand_b_flightLength3D_all);
//    nBranches_->pft_b_flightLengthErr3D.push_back(cands[ic].cand_b_flightLengthErr3D_all);
//    nBranches_->pft_b_flightSig2D.push_back(cands[ic].cand_b_flightSig2D_all);
//    nBranches_->pft_b_flightLength2D.push_back(cands[ic].cand_b_flightLength2D_all);
//    nBranches_->pft_b_flightLengthErr2D.push_back(cands[ic].cand_b_flightLengthErr2D_all);

//  lowPtTau cand = {1.};
  //   cands.push_back(cand);






//  if(vprob_jpsi!=-9){
//    nBranches_->pft_mu1_pt.push_back(muoncollection[idx_mu1].pt());
//    nBranches_->pft_mu1_eta.push_back(muoncollection[idx_mu1].eta());
//    nBranches_->pft_mu1_phi.push_back(muoncollection[idx_mu1].phi());
//    nBranches_->pft_mu1_dxy.push_back(TMath::Abs(muoncollection[idx_mu1].muonBestTrack()->dxy(firstGoodVertex->position())));
//    nBranches_->pft_mu1_dz.push_back(TMath::Abs(muoncollection[idx_mu1].muonBestTrack()->dz(firstGoodVertex->position())));
//    
//    
//    nBranches_->pft_mu2_pt.push_back(muoncollection[idx_mu2].pt());
//    nBranches_->pft_mu2_eta.push_back(muoncollection[idx_mu2].eta());
//    nBranches_->pft_mu2_phi.push_back(muoncollection[idx_mu2].phi());
//    nBranches_->pft_mu2_dxy.push_back(TMath::Abs(muoncollection[idx_mu2].muonBestTrack()->dxy(firstGoodVertex->position())));
//    nBranches_->pft_mu2_dz.push_back(TMath::Abs(muoncollection[idx_mu2].muonBestTrack()->dz(firstGoodVertex->position())));
//    
//    nBranches_->pft_jpsi_pt.push_back(tlv_jpsi_highest.Pt());
//    nBranches_->pft_jpsi_eta.push_back(tlv_jpsi_highest.Eta());
//    nBranches_->pft_jpsi_phi.push_back(tlv_jpsi_highest.Phi());
//    nBranches_->pft_jpsi_mass.push_back(tlv_jpsi_highest.M());
//    
//    nBranches_->pft_jpsi_vprob.push_back(vprob_jpsi);
//    nBranches_->pft_jpsi_vz.push_back(vz_jpsi);
//    
//    nBranches_->pft_jpsi_flightSig3D.push_back(flightSig3D_jpsi);
//    nBranches_->pft_jpsi_flightLength3D.push_back(flightLength3D_jpsi);
//  }

  //  Bool_t isgen3 = false; 
  
  //  for(unsigned int mmm=0; mmm < ppdgId.size(); mmm++){
  //    if(TMath::Abs(ppdgId[mmm]) == 541){
  //      isgen3 = true;
  //    }
  //  }
  
  //  std::cout << isgen3 << " " << isgen3matched << std::endl;
  
//  nBranches_->pft_isgen3.push_back(isgen3);
//  nBranches_->pft_isgen3matched.push_back(isgen3matched);
//  nBranches_->pft_gentaudm.push_back(1);
//  nBranches_->pft_nch.push_back(numOfch);
//  nBranches_->pft_ngentau.push_back(gps.size());


//  nBranches_->pft_tau_vprob = _pft_tau_vprob;
//  nBranches_->pft_tau_vz = _pft_tau_vz;
//  nBranches_->pft_tau_dz2pv = _pft_tau_dz2pv;
//  nBranches_->pft_tau_dxy2pv = _pft_tau_dxy2pv;
//  nBranches_->pft_tau_pt = _pft_tau_pt;
//  nBranches_->pft_tau_pt1 = _pft_tau_pt1;
//  nBranches_->pft_tau_pt2 = _pft_tau_pt2;
//  nBranches_->pft_tau_pt3 = _pft_tau_pt3;
//  nBranches_->pft_tau_eta = _pft_tau_eta; 
//  nBranches_->pft_tau_phi = _pft_tau_phi;
//  nBranches_->pft_tau_mass = _pft_tau_mass;
//  nBranches_->pft_tau_mass12 = _pft_tau_mass12;
//  nBranches_->pft_tau_mass13 = _pft_tau_mass13;
//  nBranches_->pft_tau_mass23 = _pft_tau_mass23;
//  nBranches_->pft_tau_flightSig3D = _pft_tau_flightSig3D; 
//  nBranches_->pft_tau_flightLength3D = _pft_tau_flightLength3D;
//  nBranches_->pft_tau_flightLengthErr3D = _pft_tau_flightLengthErr3D;
//  nBranches_->pft_tau_flightSig2D = _pft_tau_flightSig2D; 
//  nBranches_->pft_tau_flightLength2D = _pft_tau_flightLength2D;
//  nBranches_->pft_tau_flightLengthErr2D = _pft_tau_flightLengthErr2D;
//  nBranches_->pft_tau_iso03 = _pft_tau_iso03; 
//  nBranches_->pft_tau_iso04 = _pft_tau_iso04; 
//  nBranches_->pft_tau_max_dr_3prong = _pft_tau_max_dr_3prong;
//
//  nBranches_->pft_b_vprob = _pft_b_vprob;
//  nBranches_->pft_b_vz = _pft_b_vz;
//  nBranches_->pft_b_pt = _pft_b_pt;
//  nBranches_->pft_b_eta = _pft_b_eta;
//  nBranches_->pft_b_phi = _pft_b_phi;
//  nBranches_->pft_b_mass = _pft_b_mass;
//  nBranches_->pft_b_flightSig3D = _pft_b_flightSig3D; 
//  nBranches_->pft_b_flightLength3D = _pft_b_flightLength3D;
//  nBranches_->pft_b_flightLengthErr3D = _pft_b_flightLengthErr3D;
//  nBranches_->pft_b_flightSig2D = _pft_b_flightSig2D; 
//  nBranches_->pft_b_flightLength2D = _pft_b_flightLength2D;
//  nBranches_->pft_b_flightLengthErr2D = _pft_b_flightLengthErr2D;


  //  pft_tree->Fill();




///  nBranches_->Jpsi_pvz = pvz;
///  nBranches_->Jpsi_met = met;
///  nBranches_->Jpsi_metphi = metphi;
///  nBranches_->Jpsi_tau_ncomb = ncomb;
///  nBranches_->Jpsi_tau_vprob = pft_tau_vprob;
///  nBranches_->Jpsi_tau_vz = pft_tau_vz;
///  nBranches_->Jpsi_tau_flightSig3D = pft_tau_flightSig3D;
///  nBranches_->Jpsi_tau_flightLength3D = pft_tau_flightLength3D;
///  nBranches_->Jpsi_tau_flightLengthErr3D = pft_tau_flightLengthErr3D;
///  nBranches_->Jpsi_tau_flightSig2D = pft_tau_flightSig2D;
///  nBranches_->Jpsi_tau_flightLength2D = pft_tau_flightLength2D;
///  nBranches_->Jpsi_tau_flightLengthErr2D = pft_tau_flightLengthErr2D;
///  nBranches_->Jpsi_tau_pt = pft_tau_pt;
///  nBranches_->Jpsi_tau_pt1 = pft_tau_pt1;
///  nBranches_->Jpsi_tau_pt2 = pft_tau_pt2;
///  nBranches_->Jpsi_tau_pt3 = pft_tau_pt3;
///  nBranches_->Jpsi_tau_eta = pft_tau_eta;
///  nBranches_->Jpsi_tau_phi = pft_tau_phi;
///  nBranches_->Jpsi_tau_mass = pft_tau_mass;
///  nBranches_->Jpsi_tau_mass12 = pft_tau_mass12;
///  nBranches_->Jpsi_tau_mass13 = pft_tau_mass13;
///  nBranches_->Jpsi_tau_mass23 = pft_tau_mass23;
///
///  nBranches_->Jpsi_tau_iso03 = pft_tau_iso03;
///  nBranches_->Jpsi_tau_iso04 = pft_tau_iso04;
///  nBranches_->Jpsi_tau_max_dr_3prong = pft_tau_max_dr_3prong;
///
///  nBranches_->Jpsi_tau_dphi2jpsi = pft_tau_dphi2met;
///  nBranches_->Jpsi_tau_dz2pv = pft_tau_dz2pv;
///  nBranches_->Jpsi_tau_dxy2pv = pft_tau_dxy2pv;
///  nBranches_->Jpsi_tau_dz2jpsi = pft_tau_dz2jpsi;
///  nBranches_->Jpsi_tau_dxy2jpsi = pft_tau_dxy2jpsi;
///  nBranches_->Jpsi_tau_dz2b = pft_tau_dz2b;
///  nBranches_->Jpsi_tau_dxy2b = pft_tau_dxy2b;
///
///  nBranches_->Jpsi_b_vprob = pft_b_vprob;
///  nBranches_->Jpsi_b_vz = pft_b_vz;
///  nBranches_->Jpsi_b_pt = pft_b_pt;
///  nBranches_->Jpsi_b_eta = pft_b_eta;
///  nBranches_->Jpsi_b_phi = pft_b_phi;
///  nBranches_->Jpsi_b_mass = pft_b_mass;
///  nBranches_->Jpsi_b_flightSig3D = pft_b_flightSig3D;
///  nBranches_->Jpsi_b_flightLength3D = pft_b_flightLength3D;
///  nBranches_->Jpsi_b_flightLengthErr3D = pft_b_flightLengthErr3D;
///  nBranches_->Jpsi_b_flightSig2D = pft_b_flightSig2D;
///  nBranches_->Jpsi_b_flightLength2D = pft_b_flightLength2D;
///  nBranches_->Jpsi_b_flightLengthErr2D = pft_b_flightLengthErr2D;
///
///  nBranches_->Jpsi_mu1_pt = muoncollection[idx_mu1].pt();
///  nBranches_->Jpsi_mu1_eta = muoncollection[idx_mu1].eta();
///  nBranches_->Jpsi_mu1_phi = muoncollection[idx_mu1].phi();
///  nBranches_->Jpsi_mu1_dxy = TMath::Abs(muoncollection[idx_mu1].muonBestTrack()->dxy(firstGoodVertex->position()));
///  nBranches_->Jpsi_mu1_dz = TMath::Abs(muoncollection[idx_mu1].muonBestTrack()->dz(firstGoodVertex->position()));
///
///
///  nBranches_->Jpsi_mu2_pt = muoncollection[idx_mu2].pt();
///  nBranches_->Jpsi_mu2_eta = muoncollection[idx_mu2].eta();
///  nBranches_->Jpsi_mu2_phi = muoncollection[idx_mu2].phi();
///  nBranches_->Jpsi_mu2_dxy = TMath::Abs(muoncollection[idx_mu2].muonBestTrack()->dxy(firstGoodVertex->position()));
///  nBranches_->Jpsi_mu2_dz = TMath::Abs(muoncollection[idx_mu2].muonBestTrack()->dz(firstGoodVertex->position()));
///
///  nBranches_->Jpsi_Jpsi_pt = tlv_jpsi_highest.Pt();
///  nBranches_->Jpsi_Jpsi_eta = tlv_jpsi_highest.Eta();
///  nBranches_->Jpsi_Jpsi_phi = tlv_jpsi_highest.Phi();
///  nBranches_->Jpsi_Jpsi_mass = tlv_jpsi_highest.M();
///
///  nBranches_->Jpsi_Jpsi_vprob = vprob_jpsi;
///  nBranches_->Jpsi_Jpsi_vz = vz_jpsi;
///
///  nBranches_->Jpsi_Jpsi_flightSig3D = flightSig3D_jpsi;
///  nBranches_->Jpsi_Jpsi_flightLength3D = flightLength3D_jpsi;








   ////////////////////////////////////////////////

   
   
  iEvent.put(std::move(result));
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
lowPtTaus::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
lowPtTaus::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
lowPtTaus::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
lowPtTaus::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
lowPtTaus::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
lowPtTaus::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
lowPtTaus::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(lowPtTaus);
