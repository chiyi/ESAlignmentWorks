// -*- C++ -*-
//
// Package:    ESHitSkimLoose
// Class:      ESHitSkimLoose
//
/**\class ESHitSkimLoose ESHitSkimLoose.cc SkimTool/ESHitSkimLoose/src/ESHitSkimLoose.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Kai-Yi KAO
//         Created:  Sat Jan 30 13:11:07 CET 2010
// $Id: ESHitSkimLoose.cc,v 1.2 2011/05/10 16:13:41 chiyi Exp $
//
//

#include "SkimTool/ESHitSkimLoose/interface/ESHitSkimLoose.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
//#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/SiStripDetId/interface/TECDetId.h"
//#include "DataFormats/EcalDetId/interface/ESDetId.h"
//#include "DataFormats/EcalDetId/interface/EBDetId.h"
//#include "DataFormats/EcalDetId/interface/EEDetId.h"

//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/Math/interface/LorentzVector.h"

//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"

//#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
//#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
//#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
//#include "MagneticField/Engine/interface/MagneticField.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
//#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
//#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
//#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
//#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
//#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

//ROOT includes
//#include <Math/VectorUtil.h>
//#include <TLorentzVector.h>
//#include <TH1F.h>
//#include <TH1D.h>
//#include <TH2F.h>

#define GEN 0

using namespace std;

//
// constructors and destructor
//
ESHitSkimLoose::ESHitSkimLoose(const edm::ParameterSet& iConfig)
//  :
//  verboseLevel_(0),
//  helper_(iConfig)
{
 std::cout<<"In ESHitSkimLoose Constructor\n";
  _evt_run = 0;

}


ESHitSkimLoose::~ESHitSkimLoose()
{
 std::cout<<"In ESHitSkimLoose destructor\n";
}


//
// member functions
//

// ------------ method called to for each event  ------------
bool ESHitSkimLoose::filter(edm::Event& evt, const edm::EventSetup& iSetup)
//void ESHitSkimLoose::analyze( edm::Event & evt, const edm::EventSetup & iSetup )
{
  using namespace edm;
  using namespace std;

  //if ( verboseLevel_ > 10 )
    std::cout << "ESHitSkimLoose:: in analyze()." << std::endl;
  
  _runNum = evt.id().run();
  _evtNum = evt.id().event();

  _evt_run++;
  Ntrack = Nesrh = 0; 

  for ( int i=0; i<10; i++ ) {
    _TrackPt[i] = _TrackEta[i] = _TrackPhi[i] =0.; _TrackCharge[i] = 0;
    _TrackVx[i] = _TrackVy[i] = _TrackVz[i] =0.;
    _TrackNHit[i]=0; 
    _TrackNChi2[i]=0; 
    _Trackd0[i] = 0.; 
    _TrackPtError[i]=0.; _TrackQuality[i]=-1; 
    _TrackOuterZ[i]=0.; _TrackOuterEta[i]=0.; _TrackOuterPhi[i]=0.;

  }

  //for cluster shape  
  ESHandle<CaloGeometry> caloGeometry;
  iSetup.get<CaloGeometryRecord>().get(caloGeometry);
  const CaloGeometry *caloGeom = caloGeometry.product();
  //ESHandle<MagneticField> theMagField; 
  //iSetup.get<IdealMagneticFieldRecord>().get(theMagField); 
  //ESHandle<GlobalTrackingGeometry> theTrackingGeometry; 
  //iSetup.get<GlobalTrackingGeometryRecord>().get( theTrackingGeometry ); 

 
  // Get ES rechits
  edm::Handle<EcalRecHitCollection> PreshowerRecHits;
  evt.getByLabel(InputTag("ecalPreshowerRecHit","EcalRecHitsES"), PreshowerRecHits);
  if( PreshowerRecHits.isValid() ) 
    EcalRecHitMetaCollection preshowerHits(*PreshowerRecHits);  

  const ESRecHitCollection *ESRH = PreshowerRecHits.product();
  EcalRecHitCollection::const_iterator esrh_it;
  for ( esrh_it = ESRH->begin(); esrh_it != ESRH->end(); esrh_it++)
  { Nesrh++; break; }


  edm::Handle<reco::TrackCollection>   TrackCol;
  evt.getByLabel( "generalTracks",      TrackCol );

  cout << " number of tracks " << TrackCol->size() << endl;

  for(reco::TrackCollection::const_iterator itTrack = TrackCol->begin();
      itTrack != TrackCol->end(); ++itTrack)
  {    
    if ( itTrack->charge()!=0 )
    {
	_TrackPt[Ntrack] = itTrack->pt(); 
	_TrackEta[Ntrack] = itTrack->eta(); 
	_TrackPhi[Ntrack] = itTrack->phi(); 
	_TrackVx[Ntrack] = itTrack->vx(); 
	_TrackVy[Ntrack] = itTrack->vy(); 
	_TrackVz[Ntrack] = itTrack->vz(); 
	_TrackCharge[Ntrack] = itTrack->charge(); 
        _Trackd0[Ntrack]=itTrack->d0(); 
        _TrackNHit[Ntrack]=itTrack->found(); 
        _TrackNChi2[Ntrack]=itTrack->normalizedChi2(); 
        _TrackPtError[Ntrack]=itTrack->ptError();
        _TrackQuality[Ntrack]=itTrack->qualityMask();
        _TrackOuterZ[Ntrack]=itTrack->outerZ();
        _TrackOuterEta[Ntrack]=itTrack->outerEta();
        _TrackOuterPhi[Ntrack]=itTrack->outerPhi();

     if( _TrackPt[Ntrack]>1. //&& _TrackPtError[Ntrack]<2.
         //&& itTrack->d0Error()<0.003
        &&fabs(itTrack->outerZ())>260&&fabs(itTrack->outerZ())<280
        &&fabs(_TrackEta[Ntrack])<2.3&&fabs(_TrackEta[Ntrack])>1.7
        &&_TrackNHit[Ntrack]>=10
        &&((_TrackQuality[Ntrack])%8)>=4
        //&&_TrackNHit[Ntrack]>=5
       )
     {
      Ntrack++; break;
     }//end if TrackPt>1
    }//charge!=0
  }

  if( Nesrh>0 && Ntrack>0 ) return true;  

  return false;

}





// ------------ method called once each job just before starting event loop  ------------
void
ESHitSkimLoose::beginJob()
{

 std::cout<<"In ESHitSkimLoose.beginJob\n";

}



// ------------ method called once each job just after ending the event loop  ------------
void
ESHitSkimLoose::endJob() {
  std::cout<<"In ESHitSkimLoose.endJob\n";

  cout << endl;
  cout << " --------------------------------------------- " << endl;
  cout << " number of events processed  " << _evt_run << endl; 
  cout << " Last Event number # " << _evtNum << endl; 
  cout << " --------------------------------------------- " << endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ESHitSkimLoose);
