#include "RecoTauTag/RecoTau/interface/AntiElectronDeadECAL.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <TMath.h>

AntiElectronDeadECAL::AntiElectronDeadECAL(const edm::ParameterSet& cfg)
  : isFirstEvent_(true)
{
  minStatus_ = cfg.getParameter<uint32_t>("minStatus");
  dR_ = cfg.getParameter<double>("dR");
}

AntiElectronDeadECAL::~AntiElectronDeadECAL()
{}

void AntiElectronDeadECAL::beginEvent(const edm::EventSetup& es)
{
  updateBadTowers(es);
  positionAtECalEntrance_.beginEvent(es);
}

namespace
{
  template <class Id>
  void loopXtals(std::map<uint32_t, unsigned>& nBadCrystals,
		 std::map<uint32_t, unsigned>& maxStatus,
		 std::map<uint32_t, double>& sumEta,
		 std::map<uint32_t, double>& sumPhi ,
		 const EcalChannelStatus* channelStatus,
		 const CaloGeometry* caloGeometry,
		 const EcalTrigTowerConstituentsMap* ttMap,
                 unsigned minStatus,
                 const uint16_t statusMask)
  {
    // NOTE: modified version of SUSY CAF code
    //         UserCode/SusyCAF/plugins/SusyCAF_EcalDeadChannels.cc
    for ( int i = 0; i < Id::kSizeForDenseIndexing; ++i ) {
      Id id = Id::unhashIndex(i);  
      if ( id == Id(0) ) continue;
      EcalChannelStatusMap::const_iterator it = channelStatus->getMap().find(id.rawId());
      unsigned status = ( it == channelStatus->end() ) ? 
	0 : (it->getStatusCode() & statusMask);
      if ( status >= minStatus ) {
	const GlobalPoint& point = caloGeometry->getPosition(id);
	uint32_t key = ttMap->towerOf(id);
	maxStatus[key] = TMath::Max(status, maxStatus[key]);
	++nBadCrystals[key];
	sumEta[key] += point.eta();
	sumPhi[key] += point.phi();
      }
    }
  }
}

void AntiElectronDeadECAL::updateBadTowers(const edm::EventSetup& es) 
{
  // NOTE: modified version of SUSY CAF code
  //         UserCode/SusyCAF/plugins/SusyCAF_EcalDeadChannels.cc
  const uint32_t channelStatusId = es.get<EcalChannelStatusRcd>().cacheIdentifier();
  const uint32_t caloGeometryId  = es.get<CaloGeometryRecord>().cacheIdentifier();
  const uint32_t idealGeometryId = es.get<IdealGeometryRecord>().cacheIdentifier();

  if ( !isFirstEvent_ && channelStatusId == channelStatusId_cache_ && caloGeometryId == caloGeometryId_cache_ && idealGeometryId == idealGeometryId_cache_  ) return;

  edm::ESHandle<EcalChannelStatus> channelStatus;    
  es.get<EcalChannelStatusRcd>().get(channelStatus);
  channelStatusId_cache_ = channelStatusId;  

  edm::ESHandle<CaloGeometry> caloGeometry;         
  es.get<CaloGeometryRecord>().get(caloGeometry);
  caloGeometryId_cache_ = caloGeometryId;  

  edm::ESHandle<EcalTrigTowerConstituentsMap> ttMap;
  es.get<IdealGeometryRecord>().get(ttMap);
  idealGeometryId_cache_ = idealGeometryId;

  std::map<uint32_t,unsigned> nBadCrystals, maxStatus;
  std::map<uint32_t,double> sumEta, sumPhi;
    
  loopXtals<EBDetId>(nBadCrystals, maxStatus, sumEta, sumPhi, channelStatus.product(), caloGeometry.product(), ttMap.product(), minStatus_, statusMask_);
  loopXtals<EEDetId>(nBadCrystals, maxStatus, sumEta, sumPhi, channelStatus.product(), caloGeometry.product(), ttMap.product(), minStatus_, statusMask_);
    
  badTowers_.clear();
  for ( std::map<uint32_t, unsigned>::const_iterator it = nBadCrystals.begin(); 
	it != nBadCrystals.end(); ++it ) {
    uint32_t key = it->first;
    badTowers_.push_back(towerInfo(key, it->second, maxStatus[key], sumEta[key]/it->second, sumPhi[key]/it->second));
  }

  isFirstEvent_ = false;
}


bool AntiElectronDeadECAL::operator()(const reco::Candidate* leadPFChargedHadron) const
{
  bool isNearBadTower = false;
  bool success = false;
  reco::Candidate::Point positionAtECalEntrance = positionAtECalEntrance_(leadPFChargedHadron, success);
  if ( success ) {
    for ( std::vector<towerInfo>::const_iterator badTower = badTowers_.begin();
	  badTower != badTowers_.end(); ++badTower ) {
      if ( deltaR(badTower->eta_, badTower->phi_, positionAtECalEntrance.eta(), positionAtECalEntrance.phi()) < dR_ ) {
        isNearBadTower = true;
      }
    }
  }
  return isNearBadTower;
}
