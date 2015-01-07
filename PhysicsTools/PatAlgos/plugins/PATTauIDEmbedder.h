/*
 * PATTauIDEmbedder.h
 *
 *  Created on: Nov 28, 2014
 *      Author: nehrkorn
 */

#ifndef PATTAUIDEMBEDDER_H_
#define PATTAUIDEMBEDDER_H_

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"

class PATTauIDEmbedder : public edm::EDProducer
{
 public:

  explicit PATTauIDEmbedder(const edm::ParameterSet&);
  ~PATTauIDEmbedder(){};

  void produce(edm::Event&, const edm::EventSetup&);

 private:

//--- configuration parameters
  edm::InputTag src_;
  typedef std::pair<std::string, edm::InputTag> NameTag;
  std::vector<NameTag> tauIDSrcs_;
  std::vector<edm::EDGetTokenT<pat::PATTauDiscriminator> > patTauIDTokens_;
};

#endif /* PATTAUIDEMBEDDER_H_ */
