// -*- C++ -*-
//
// Package:    ESAlignTool
// Class:      ESAlignTool
// 
/**\class ESAlignTool ESAlignTool.cc AlignmentTool/ESAlignTool/src/ESAlignTool.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kai-Yi Kao,27 2-020,+41227674870,
//         Created:  Wed Jun 16 09:36:53 CEST 2010
// $Id: ESAlignTool.cc,v 1.5 2011/05/26 12:27:33 chiyi Exp $
//
//


// system include files
#include <memory>

// user include files
#include "AlignmentTool/ESAlignTool/interface/ESAlignTool.h"
#include "AlignmentTool/ESAlignTool/plugins/ESAlignSelections.C"
#include "AlignmentTool/ESAlignTool/plugins/ESAlignDeadZone.C"
#include "AlignmentTool/ESAlignTool/func/BadSensorLogic.C"
#include "AlignmentTool/ESAlignTool/func/AlignmentCalculation.C"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
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

//#include "DataFormats/EgammaReco/interface/SuperCluster.h"
//#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
//#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/DetId/interface/DetId.h"
//#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
//#include "DataFormats/EcalDetId/interface/EBDetId.h"
//#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/Math/interface/LorentzVector.h"

//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
//#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"

//#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
//#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
//#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
//#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
//#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
//ROOT includes
#include "TMath.h"
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>

//
// constructors and destructor
//
ESAlignTool::ESAlignTool(const edm::ParameterSet& iConfig)
{
 std::cout << "In ESAlignTool Constructor\n" ;
 initAllPara(iConfig);
}

ESAlignTool::~ESAlignTool()
{
  
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
 std::cout << "In ESAlignTool destructor\n" ;
} 


void ESAlignTool::initAllPara(const edm::ParameterSet& iConfig)
{
 ESpF_Printed=0; ESpR_Printed=0;
 ESmF_Printed=0; ESmR_Printed=0;

 for(int i=0;i<2;i++)
 {for(int j=0;j<2;j++)
  {
   ES_M11[i][j]=0.; ES_M12[i][j]=0.; ES_M13[i][j]=0.;
   ES_M14[i][j]=0.; ES_M15[i][j]=0.; ES_M16[i][j]=0.;
   ES_M22[i][j]=0.; ES_M23[i][j]=0.; ES_M24[i][j]=0.;
   ES_M25[i][j]=0.; ES_M26[i][j]=0.;
   ES_M33[i][j]=0.; ES_M34[i][j]=0.; ES_M35[i][j]=0.;
   ES_M36[i][j]=0.;
   ES_M44[i][j]=0.; ES_M45[i][j]=0.; ES_M46[i][j]=0.;
   ES_M55[i][j]=0.; ES_M56[i][j]=0.;
   ES_M66[i][j]=0.;

   ES_P1[i][j]=0.; ES_P2[i][j]=0.; ES_P3[i][j]=0.;
   ES_P4[i][j]=0.; ES_P5[i][j]=0.; ES_P6[i][j]=0.;
   ES_CHI2[i][j]=0.;
   ES_NTracks[i][j]=0;
   ES_dX[i][j]=0.; ES_dY[i][j]=0.; ES_dZ[i][j]=0.;
   ES_dXerr[i][j]=0.; ES_dYerr[i][j]=0.; ES_dZerr[i][j]=0.;
   ES_dAlpha[i][j]=0.; ES_dBeta[i][j]=0.; ES_dGamma[i][j]=0.;
   ES_dAlphaerr[i][j]=0.; ES_dBetaerr[i][j]=0.; ES_dGammaerr[i][j]=0.;

   ES_O_Alpha[i][j]=0.; ES_O_Beta[i][j]=0.; ES_O_Gamma[i][j]=0.;

   ES_M31Err2[i][j]=0.; ES_M32Err2[i][j]=0.; ES_M33Err2[i][j]=0.;
   ES_P1Err2[i][j]=0.; ES_P2Err2[i][j]=0.; ES_P3Err2[i][j]=0.;

   ES_InvM11Err2[i][j]=0.; ES_InvM12Err2[i][j]=0.; ES_InvM13Err2[i][j]=0.;
   ES_InvM22Err2[i][j]=0.; ES_InvM23Err2[i][j]=0.; ES_InvM33Err2[i][j]=0.;

   ES_R11[i][j]=0.; ES_R12[i][j]=0.; ES_R13[i][j]=0.;
   ES_R21[i][j]=0.; ES_R22[i][j]=0.; ES_R23[i][j]=0.;
   ES_R31[i][j]=0.; ES_R32[i][j]=0.; ES_R33[i][j]=0.;

   ES_O_R11[i][j]=1.; ES_O_R12[i][j]=0.; ES_O_R13[i][j]=0.;
   ES_O_R21[i][j]=0.; ES_O_R22[i][j]=1.; ES_O_R23[i][j]=0.;
   ES_O_R31[i][j]=0.; ES_O_R32[i][j]=0.; ES_O_R33[i][j]=1.;
 }}

 e_xxlimit = iConfig.getParameter<double>("e_xxlimit");
 e_yylimit = iConfig.getParameter<double>("e_yylimit");
 e_yxlimit = iConfig.getParameter<double>("e_yxlimit");
 winlimit = iConfig.getParameter<double>("winlimit");

 b_wRotation = iConfig.getParameter<bool>("withRotation");
 b_DrawMagField = iConfig.getParameter<bool>("DrawMagField");
 b_PrintPosition = iConfig.getParameter<bool>("PrintPosition");
 b_ReSetRfromOutside = iConfig.getParameter<bool>("ReSetRfromOutside");
 b_fromRefitter = iConfig.getParameter<bool>("fromRefitter");

 Cal_ESorigin_from_Geometry = iConfig.getParameter<bool>("Cal_ESorigin_from_Geometry");
 Cal_ESaxes_from_Geometry = iConfig.getParameter<bool>("Cal_ESaxes_from_Geometry");
 b_Overwrite_RotationMatrix_fromGeometry = iConfig.getParameter<bool>("Overwrite_RotationMatrix_fromGeometry");

 b_PrintMatrix = iConfig.getParameter<bool>("PrintMatrix");
 Selected_idee = iConfig.getParameter<unsigned int>("Selected_idee");
 Selected_RUNmin = iConfig.getParameter<int>("Selected_RUNmin");
 Selected_RUNmax = iConfig.getParameter<int>("Selected_RUNmax");

 iterN = iConfig.getParameter<unsigned int>("IterN");
 if(iterN<0||iterN>11) std::cout<<"Error : Out of Range for iterN!!!!\n";
 
 char buf[20];
 for(int iterN_idx=1;iterN_idx<iterN;iterN_idx++)
 {
  sprintf(buf,"Iter%i_ESpFdX",iterN_idx);  iter_ESpFdX[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpFdY",iterN_idx);  iter_ESpFdY[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpFdZ",iterN_idx);  iter_ESpFdZ[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpRdX",iterN_idx);  iter_ESpRdX[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpRdY",iterN_idx);  iter_ESpRdY[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpRdZ",iterN_idx);  iter_ESpRdZ[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmFdX",iterN_idx);  iter_ESmFdX[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmFdY",iterN_idx);  iter_ESmFdY[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmFdZ",iterN_idx);  iter_ESmFdZ[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmRdX",iterN_idx);  iter_ESmRdX[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmRdY",iterN_idx);  iter_ESmRdY[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmRdZ",iterN_idx);  iter_ESmRdZ[iterN_idx-1] = iConfig.getParameter<double>(buf);

  sprintf(buf,"Iter%i_ESpFdAlpha",iterN_idx);  iter_ESpFdAlpha[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpFdBeta",iterN_idx);   iter_ESpFdBeta[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpFdGamma",iterN_idx);  iter_ESpFdGamma[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpRdAlpha",iterN_idx);  iter_ESpRdAlpha[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpRdBeta",iterN_idx);   iter_ESpRdBeta[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESpRdGamma",iterN_idx);  iter_ESpRdGamma[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmFdAlpha",iterN_idx);  iter_ESmFdAlpha[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmFdBeta",iterN_idx);   iter_ESmFdBeta[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmFdGamma",iterN_idx);  iter_ESmFdGamma[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmRdAlpha",iterN_idx);  iter_ESmRdAlpha[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmRdBeta",iterN_idx);   iter_ESmRdBeta[iterN_idx-1] = iConfig.getParameter<double>(buf);
  sprintf(buf,"Iter%i_ESmRdGamma",iterN_idx);  iter_ESmRdGamma[iterN_idx-1] = iConfig.getParameter<double>(buf);

 }

 for(int iterN_idx=1;iterN_idx<=iterN;iterN_idx++)
 {
  ES_dX[0][0] -= iter_ESmFdX[iterN_idx-1]; ES_dY[0][0] -= iter_ESmFdY[iterN_idx-1]; ES_dZ[0][0] -= iter_ESmFdZ[iterN_idx-1];
  ES_dX[0][1] -= iter_ESmRdX[iterN_idx-1]; ES_dY[0][1] -= iter_ESmRdY[iterN_idx-1]; ES_dZ[0][1] -= iter_ESmRdZ[iterN_idx-1];
  ES_dX[1][0] -= iter_ESpFdX[iterN_idx-1]; ES_dY[1][0] -= iter_ESpFdY[iterN_idx-1]; ES_dZ[1][0] -= iter_ESpFdZ[iterN_idx-1];
  ES_dX[1][1] -= iter_ESpRdX[iterN_idx-1]; ES_dY[1][1] -= iter_ESpRdY[iterN_idx-1]; ES_dZ[1][1] -= iter_ESpRdZ[iterN_idx-1];

  ES_dAlpha[0][0] -= iter_ESmFdAlpha[iterN_idx-1]; ES_dBeta[0][0] -= iter_ESmFdBeta[iterN_idx-1]; ES_dGamma[0][0] -= iter_ESmFdGamma[iterN_idx-1];
  ES_dAlpha[0][1] -= iter_ESmRdAlpha[iterN_idx-1]; ES_dBeta[0][1] -= iter_ESmRdBeta[iterN_idx-1]; ES_dGamma[0][1] -= iter_ESmRdGamma[iterN_idx-1];
  ES_dAlpha[1][0] -= iter_ESpFdAlpha[iterN_idx-1]; ES_dBeta[1][0] -= iter_ESpFdBeta[iterN_idx-1]; ES_dGamma[1][0] -= iter_ESpFdGamma[iterN_idx-1];
  ES_dAlpha[1][1] -= iter_ESpRdAlpha[iterN_idx-1]; ES_dBeta[1][1] -= iter_ESpRdBeta[iterN_idx-1]; ES_dGamma[1][1] -= iter_ESpRdGamma[iterN_idx-1];
 }

  _evt_run = 0;

  woRotate=new RotationType();
  //ESpF_O=new PositionType(0.,0.,303.376+ESpFdZ);
  //ESpR_O=new PositionType(0.,0.,307.836+ESpRdZ);
  //ESmF_O=new PositionType(0.,0.,-303.376+ESmFdZ);
  //ESmR_O=new PositionType(0.,0.,-307.836+ESmRdZ);
  ESpF_O=new PositionType(0.,0.,304.186); //default
  ESpR_O=new PositionType(0.,0.,308.646); //default
  ESmF_O=new PositionType(0.,0.,-304.316);//default 
  ESmR_O=new PositionType(0.,0.,-308.776);//default
  ESpF_Oap=new PositionType(ES_dX[1][0],ES_dY[1][0],304.186+ES_dZ[1][0]); //default
  ESpR_Oap=new PositionType(ES_dX[1][1],ES_dY[1][1],308.646+ES_dZ[1][1]); //default
  ESmF_Oap=new PositionType(ES_dX[0][0],ES_dY[0][0],-304.316+ES_dZ[0][0]);//default 
  ESmR_Oap=new PositionType(ES_dX[0][1],ES_dY[0][1],-308.776+ES_dZ[0][1]);//default

  init_RotationMatrices();
  //init_RotationMatrices(0);

  if(b_ReSetRfromOutside)
  {
   ES_R21[0][0]=-0.00225996; ES_R22[0][0]=0.999996;  ES_R23[0][0]=0.00159956;
   ES_R31[0][0]=0.000904; ES_R32[0][0]=-0.001598; ES_R33[0][0]=0.999998;

   ES_R11[0][1]=0.999998; ES_R12[0][1]=0.00173077;  ES_R13[0][1]=-0.000499602;
   ES_R21[0][1]=-0.00173008; ES_R22[0][1]=0.999997;  ES_R23[0][1]=0.00169866;
   ES_R31[0][1]=0.000503; ES_R32[0][1]=-0.001698; ES_R33[0][1]=0.999998;

   ES_R11[1][0]=0.999999; ES_R12[1][0]=0.00164022;  ES_R13[1][0]=0.;
   ES_R21[1][0]=-0.00164007; ES_R22[1][0]=0.999997;  ES_R23[1][0]=-0.00169612;
   ES_R31[1][0]=-0.000003; ES_R32[1][0]=0.001700; ES_R33[1][0]=0.999999;

   ES_R11[1][1]=0.999999; ES_R12[1][1]=0.00137978;  ES_R13[1][1]=0.0001500904;
   ES_R21[1][1]=-0.00138015; ES_R22[1][1]=0.999998;  ES_R23[1][1]=-0.00169865;
   ES_R31[1][1]=-0.000152; ES_R32[1][1]=0.001698; ES_R33[1][1]=0.999999;
  } 
  ESpF_wRotateap=new RotationType(ES_R11[1][0],ES_R12[1][0],ES_R13[1][0],ES_R21[1][0],ES_R22[1][0],ES_R23[1][0],ES_R31[1][0],ES_R32[1][0],ES_R33[1][0]);
  ESpR_wRotateap=new RotationType(ES_R11[1][1],ES_R12[1][1],ES_R13[1][1],ES_R21[1][1],ES_R22[1][1],ES_R23[1][1],ES_R31[1][1],ES_R32[1][1],ES_R33[1][1]);
  ESmF_wRotateap=new RotationType(ES_R11[0][0],ES_R12[0][0],ES_R13[0][0],ES_R21[0][0],ES_R22[0][0],ES_R23[0][0],ES_R31[0][0],ES_R32[0][0],ES_R33[0][0]);
  ESmR_wRotateap=new RotationType(ES_R11[0][1],ES_R12[0][1],ES_R13[0][1],ES_R21[0][1],ES_R22[0][1],ES_R23[0][1],ES_R31[0][1],ES_R32[0][1],ES_R33[0][1]);


  ESpF_residualX=new TH1D("ESpF_residualX","ES+Front residualX",300,-15,15); 
  ESpF_residualY=new TH1D("ESpF_residualY","ES+Front residualY",300,-15,15);
  ESpR_residualX=new TH1D("ESpR_residualX","ES+Rear residualX",300,-15,15);
  ESpR_residualY=new TH1D("ESpR_residualY","ES+Rear residualY",300,-15,15);
  ESmF_residualX=new TH1D("ESmF_residualX","ES-Front residualX",300,-15,15);
  ESmF_residualY=new TH1D("ESmF_residualY","ES-Front residualY",300,-15,15);
  ESmR_residualX=new TH1D("ESmR_residualX","ES-Rear residualX",300,-15,15);
  ESmR_residualY=new TH1D("ESmR_residualY","ES-Rear residualY",300,-15,15);

}

void ESAlignTool::init_RotationMatrices()
{
  for(int i=0;i<2;i++)
  {for(int j=0;j<2;j++)
   {
    ES_R11[i][j] = cos(ES_dBeta[i][j])*cos(ES_dGamma[i][j]) - sin(ES_dAlpha[i][j])*sin(ES_dBeta[i][j])*sin(ES_dGamma[i][j]);
    ES_R12[i][j] = cos(ES_dBeta[i][j])*sin(ES_dGamma[i][j]) + sin(ES_dAlpha[i][j])*sin(ES_dBeta[i][j])*cos(ES_dGamma[i][j]);
    ES_R13[i][j] = -cos(ES_dAlpha[i][j])*sin(ES_dBeta[i][j]);
    ES_R21[i][j] = -cos(ES_dAlpha[i][j])*sin(ES_dGamma[i][j]);
    ES_R22[i][j] = cos(ES_dAlpha[i][j])*cos(ES_dGamma[i][j]);
    ES_R23[i][j] = sin(ES_dAlpha[i][j]);
    ES_R31[i][j] = sin(ES_dBeta[i][j])*cos(ES_dGamma[i][j]) + sin(ES_dAlpha[i][j])*cos(ES_dBeta[i][j])*sin(ES_dGamma[i][j]);
    ES_R32[i][j] = sin(ES_dBeta[i][j])*sin(ES_dGamma[i][j]) - sin(ES_dAlpha[i][j])*cos(ES_dBeta[i][j])*cos(ES_dGamma[i][j]);
    ES_R33[i][j] = cos(ES_dAlpha[i][j])*cos(ES_dBeta[i][j]);
  }}
 
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
ESAlignTool::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

 //1.initialize for each events
 std::cout << "In ESAlignTool.analyze\n";
 std::cout << "ESAlignTool:: in analyze().\n";
 
 _runNum = iEvent.id().run();
 _evtNum = iEvent.id().event();
 init_perEvent();

 int Run = iEvent.id().run();
 int RunMin=Selected_RUNmin;
 int RunMax=Selected_RUNmax;
 //std::cout<<"RunMin="<<RunMin<<", RunMax="<<RunMax<<", Run="<<Run<<"\n";
 if(  Selected_RUNmin==0 || Selected_RUNmax==0
    || ( RunMin!=0 && RunMax!=0 && Run >= RunMin && Run <= RunMax )
   )   
 {
  std::cout<<"Pass EvtRun Selection.\n";
 
  ESHandle<CaloGeometry> caloGeometry;
  iSetup.get<CaloGeometryRecord>().get(caloGeometry);
  const CaloGeometry *caloGeom = caloGeometry.product();
  ESHandle<MagneticField> theMagField; 
  iSetup.get<IdealMagneticFieldRecord>().get(theMagField); 
  ESHandle<GlobalTrackingGeometry> theTrackingGeometry; 
  iSetup.get<GlobalTrackingGeometryRecord>().get( theTrackingGeometry ); 
  //for MagneticField on ES
  if(b_DrawMagField)  LoadMagFieldonES(_evt_run,caloGeom,theMagField);
  //for Log some position
  if(b_PrintPosition) PrintPosition(_evt_run,caloGeom);
 
  
  //2.fill es hits
  edm::Handle<EcalRecHitCollection> PreshowerRecHits;
  iEvent.getByLabel(InputTag("ecalPreshowerRecHit","EcalRecHitsES"), PreshowerRecHits);
  fill_esRecHit(caloGeom,PreshowerRecHits);
 
 
  //3.fill Track
  edm::Handle<reco::TrackCollection>   TrackCol;
  if(b_fromRefitter)
   iEvent.getByLabel( "TrackRefitter",      TrackCol );
  else
   iEvent.getByLabel( "generalTracks",      TrackCol );
 
  std::cout << " number of tracks " << TrackCol->size() << std::endl;
 
  fill_tracks(TrackCol);
  //3.1 Association for ESRecHits and tracks
  fill_esRecHitAssociation();
 
 
  //4.locate local ES origin
  set_ESorigin(_evt_run,caloGeom);
  set_ESaxes(_evt_run,caloGeom);
  if(b_Overwrite_RotationMatrix_fromGeometry) Overwrite_RotationMatrix(_evt_run);
 
 
  //5.fill PredictionState on ES
  int iz=-1; int ip=1; int idee=0; //dee=1 -X/Y ; dee=2 +X/Y
  if(!b_wRotation)
  {
   iz=1;ip=1;idee=0; fill_PredictionState(iz,ip,idee,TrackCol,theMagField,theTrackingGeometry);
   iz=1;ip=2;idee=0; fill_PredictionState(iz,ip,idee,TrackCol,theMagField,theTrackingGeometry);
   iz=-1;ip=1;idee=0; fill_PredictionState(iz,ip,idee,TrackCol,theMagField,theTrackingGeometry);
   iz=-1;ip=2;idee=0; fill_PredictionState(iz,ip,idee,TrackCol,theMagField,theTrackingGeometry);
  }
  if(b_wRotation)
  {
   iz=1;ip=1;idee=0; fill_PredictionState_wRotation(iz,ip,idee,TrackCol,theMagField,theTrackingGeometry);
   iz=1;ip=2;idee=0; fill_PredictionState_wRotation(iz,ip,idee,TrackCol,theMagField,theTrackingGeometry);
   iz=-1;ip=1;idee=0; fill_PredictionState_wRotation(iz,ip,idee,TrackCol,theMagField,theTrackingGeometry);
   iz=-1;ip=2;idee=0; fill_PredictionState_wRotation(iz,ip,idee,TrackCol,theMagField,theTrackingGeometry);
  }
 
  //6.Residual and Calculation
  if(!b_wRotation)
  {
   iz=1; fill_residual(iz);
   iz=-1; fill_residual(iz);
  }
 
  if(b_wRotation)
  {
   //iz=1; fill_residual(iz);
   //iz=-1; fill_residual(iz);
   iz=1; fill_residual_wRotation_v2(iz);
   iz=-1; fill_residual_wRotation_v2(iz);
  }
 
  t_ESAlign->Fill();

 }//end if seleted RUN
}


// ------------ method called once each job just before starting event loop  ------------
void ESAlignTool::init_perEvent()
{
 _evt_run++; Ntrack = Nesrh = 0;

 for ( int i=0; i<10000; i++)
 {
   _esRecHit_E[i]=-1; _esRecHit_X[i]=0; _esRecHit_Y[i]=0; _esRecHit_Z[i]=0;
   _esRecHit_siZ[i]=0;   _esRecHit_siP[i]=0; 
   _esRecHit_siX[i]=0;   _esRecHit_siY[i]=0;  _esRecHit_Strip[i]=0; 
   _esRecHit_Noisy[i]=0; _esRecHit_MatchedTrk_fromOuter[i]=-1;
 }//end init esRecHit

 for ( int i=0; i<2000; i++ )
 {
  _TrackPt[i] = _TrackEta[i] = _TrackPhi[i] =0.; _TrackCharge[i] = 0;
  _TrackVx[i] = _TrackVy[i] = _TrackVz[i] =0.;
  _TrackNHit[i]=0; 
  _TrackNChi2[i]=0; 
  _Trackd0[i] = 0.; 
  _TrackPtError[i]=0.; _TrackQuality[i]=-1; 
  _TrackOuterZ[i]=0.; _TrackOuterEta[i]=0.; _TrackOuterPhi[i]=0.;
  _TrackInnerX[i]=0.; _TrackInnerY[i]=0.; _TrackInnerZ[i]=0.;

  for(int j=0;j<2;j++)
  {for(int k=0;k<2;k++)
   {
    PredictionState_iz[i][j][k]=0;
    PredictionState_ip[i][j][k]=0;
    PredictionState_valid[i][j][k]=0;
    PredictionState_X[i][j][k]=0.;
    PredictionState_Y[i][j][k]=0.;
    PredictionState_Z[i][j][k]=0.;
    PredictionState_Px[i][j][k]=0.;
    PredictionState_Py[i][j][k]=0.;
    PredictionState_Pz[i][j][k]=0.;
    PredictionState_Bx[i][j][k]=0.;
    PredictionState_By[i][j][k]=0.;
    PredictionState_Bz[i][j][k]=0.;
    PredictionState_Exx[i][j][k]=0.;
    PredictionState_Eyx[i][j][k]=0.;
    PredictionState_Eyy[i][j][k]=0.;
    PredictionState_E44[i][j][k]=0.;
    PredictionState_E55[i][j][k]=0.;
    PredictionState_E66[i][j][k]=0.;
    PredictionState_E45[i][j][k]=0.;
    PredictionState_E46[i][j][k]=0.;
    PredictionState_E56[i][j][k]=0.;
    PredictionState_E14[i][j][k]=0.;
    PredictionState_E15[i][j][k]=0.;
    PredictionState_E16[i][j][k]=0.;
    PredictionState_E24[i][j][k]=0.;
    PredictionState_E25[i][j][k]=0.;
    PredictionState_E26[i][j][k]=0.;
    PredictionState_delX[i][j][k]=0.;
    PredictionState_delY[i][j][k]=0.;
    PredictionState_MatchedRec[i][j][k]=-1;
    PredictionState_resiX[i][j][k]=0.;
    PredictionState_resiY[i][j][k]=0.;
  }}
 }//end init Track n PredictionState

}

void ESAlignTool::LoadMagFieldonES(Long64_t _evt_run, const CaloGeometry *caloGeom, edm::ESHandle<MagneticField> theMagField )
{
 if(_evt_run==1)
 {
  for(int ix=1;ix<=40;ix++)
  {for(int iy=1;iy<=40;iy++)
   {for(int is=1;is<=32;is++)
    {
     ESpF_B_x[ix-1][iy-1][is-1]=0.;
     ESpF_B_y[ix-1][iy-1][is-1]=0.;
     ESpF_B_z[ix-1][iy-1][is-1]=0.;
     ESpR_B_x[ix-1][iy-1][is-1]=0.;
     ESpR_B_y[ix-1][iy-1][is-1]=0.;
     ESpR_B_z[ix-1][iy-1][is-1]=0.;
     ESmF_B_x[ix-1][iy-1][is-1]=0.;
     ESmF_B_y[ix-1][iy-1][is-1]=0.;
     ESmF_B_z[ix-1][iy-1][is-1]=0.;
     ESmR_B_x[ix-1][iy-1][is-1]=0.;
     ESmR_B_y[ix-1][iy-1][is-1]=0.;
     ESmR_B_z[ix-1][iy-1][is-1]=0.;
  }}}
 }

 if(_evt_run==1)
 {
  for(int iz=-1;iz<2;iz=iz+2)
  {for(int ip=1;ip<=2;ip++)
   {for(int ix=1;ix<=40;ix++)
    {for(int iy=1;iy<=40;iy++)
     {for(int is=1;is<=32;is++)
      {
       if(ESDetId::validDetId(is,ix,iy,ip,iz))
       {
        ESDetId SaveMag(is,ix,iy,ip,iz);

        if(SaveMag.zside()==1&&SaveMag.plane()==1)
        {
         ESpF_B_x[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).x();
         ESpF_B_y[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).y();
         ESpF_B_z[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).z();
        }
        if(SaveMag.zside()==1&&SaveMag.plane()==2)
        {
         ESpR_B_x[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).x();
         ESpR_B_y[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).y();
         ESpR_B_z[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).z();
        }
        if(SaveMag.zside()==-1&&SaveMag.plane()==1)
        {
         ESmF_B_x[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).x();
         ESmF_B_y[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).y();
         ESmF_B_z[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).z();
        }
        if(SaveMag.zside()==-1&&SaveMag.plane()==2)
        {
         ESmR_B_x[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).x();
         ESmR_B_y[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).y();
         ESmR_B_z[SaveMag.six()-1][SaveMag.siy()-1][SaveMag.strip()-1]=theMagField->inTesla(caloGeom->getPosition(SaveMag)).z();
        }
  }}}}}}

  t_ESField->Fill();
 }
}

void ESAlignTool::PrintPosition(Long64_t _evt_run, const CaloGeometry *caloGeom)
{
 if(_evt_run==1)
 {
  for(int iz=-1;iz<2;iz=iz+2)
  {for(int ip=1;ip<=2;ip++)
   {for(int ix=1;ix<=40;ix++)
    {for(int iy=1;iy<=40;iy++)
     {for(int is=1;is<=32;is++)
      {
       if(ESDetId::validDetId(is,ix,iy,ip,iz))
       {
        ESDetId SaveMag(is,ix,iy,ip,iz);
        //For logout position //START
        if(  (  SaveMag.zside()==1&&SaveMag.plane()==1&&SaveMag.six()==5
              &&SaveMag.siy()==20&&SaveMag.strip()==1 
             )
           ||(  SaveMag.zside()==1&&SaveMag.plane()==1&&SaveMag.six()==35
              &&SaveMag.siy()==20&&SaveMag.strip()==1
             )
           ||(  SaveMag.zside()==1&&SaveMag.plane()==1&&SaveMag.six()==5
              &&SaveMag.siy()==20&&SaveMag.strip()==32
             )
           ||(  SaveMag.zside()==1&&SaveMag.plane()==1&&SaveMag.six()==5
              &&SaveMag.siy()==21&&SaveMag.strip()==1
             )
           ||(  SaveMag.zside()==1&&SaveMag.plane()==1&&SaveMag.six()==5
              &&SaveMag.siy()==21&&SaveMag.strip()==32
             )
           ||(  SaveMag.zside()==-1&&SaveMag.plane()==1&&SaveMag.six()==5
              &&SaveMag.siy()==20&&SaveMag.strip()==1
             )
           ||(  SaveMag.zside()==-1&&SaveMag.plane()==1&&SaveMag.six()==5
              &&SaveMag.siy()==20&&SaveMag.strip()==32
             )
           ||(  SaveMag.zside()==-1&&SaveMag.plane()==1&&SaveMag.six()==5
              &&SaveMag.siy()==21&&SaveMag.strip()==1
             )
           ||(  SaveMag.zside()==-1&&SaveMag.plane()==1&&SaveMag.six()==5
              &&SaveMag.siy()==21&&SaveMag.strip()==32
             )
           ||(  SaveMag.zside()==1&&SaveMag.plane()==2&&SaveMag.six()==20
              &&SaveMag.siy()==5&&SaveMag.strip()==1
             )
           ||(  SaveMag.zside()==1&&SaveMag.plane()==2&&SaveMag.six()==20
              &&SaveMag.siy()==5&&SaveMag.strip()==32
             )
           ||(  SaveMag.zside()==1&&SaveMag.plane()==2&&SaveMag.six()==21
              &&SaveMag.siy()==5&&SaveMag.strip()==1
             )
           ||(  SaveMag.zside()==1&&SaveMag.plane()==2&&SaveMag.six()==21
              &&SaveMag.siy()==5&&SaveMag.strip()==32
             )
           ||(  SaveMag.zside()==-1&&SaveMag.plane()==2&&SaveMag.six()==20
              &&SaveMag.siy()==5&&SaveMag.strip()==1
             )
           ||(  SaveMag.zside()==-1&&SaveMag.plane()==2&&SaveMag.six()==20
              &&SaveMag.siy()==5&&SaveMag.strip()==32
             )
           ||(  SaveMag.zside()==-1&&SaveMag.plane()==2&&SaveMag.six()==21
              &&SaveMag.siy()==5&&SaveMag.strip()==1
             )
           ||(  SaveMag.zside()==-1&&SaveMag.plane()==2&&SaveMag.six()==21
              &&SaveMag.siy()==5&&SaveMag.strip()==32
             )
            )
        {
         std::cout<<"iz="<<SaveMag.zside()<<",ip="<<SaveMag.plane();
         std::cout<<",ix="<<SaveMag.six()<<",iy="<<SaveMag.siy();
         std::cout<<",is="<<SaveMag.strip()<<"\n";
         std::cout<<"(X,Y,Z)=("<<caloGeom->getPosition(SaveMag).x();
         std::cout<<","<<caloGeom->getPosition(SaveMag).y();
         std::cout<<","<<caloGeom->getPosition(SaveMag).z();
         std::cout<<")\n";
        }//For log out position //END
 }}}}}}}
}

void ESAlignTool::fill_esRecHit(const CaloGeometry *caloGeom, edm::Handle<EcalRecHitCollection> PreshowerRecHits)
{
 const ESRecHitCollection *ESRH = PreshowerRecHits.product();
 EcalRecHitCollection::const_iterator esrh_it;

 for ( esrh_it = ESRH->begin(); esrh_it != ESRH->end(); esrh_it++)
 {
  Double_t esrh_x = caloGeom->getPosition(esrh_it->id()).x();
  Double_t esrh_y = caloGeom->getPosition(esrh_it->id()).y();
  Double_t esrh_z = caloGeom->getPosition(esrh_it->id()).z();
  Double_t esrh_eta = caloGeom->getPosition(esrh_it->id()).eta();
  Double_t esrh_phi = caloGeom->getPosition(esrh_it->id()).phi();
  ESDetId esdetid = ESDetId(esrh_it->id());
  if(Nesrh>=10000)
  {
   edm::LogWarning("fill_esRecHit")<<"Too many ES RecHits.\n";
   continue;
  }
  _esRecHit_E[Nesrh] = esrh_it->energy(); 
  _esRecHit_X[Nesrh] = esrh_x; 
  _esRecHit_Y[Nesrh] = esrh_y; 
  _esRecHit_Z[Nesrh] = esrh_z; 
  _esRecHit_Eta[Nesrh] = esrh_eta; 
  _esRecHit_Phi[Nesrh] = esrh_phi; 
  _esRecHit_siZ[Nesrh] = esdetid.zside(); 
  _esRecHit_siP[Nesrh] = esdetid.plane(); 
  _esRecHit_siX[Nesrh] = esdetid.six(); 
  _esRecHit_siY[Nesrh] = esdetid.siy(); 
  _esRecHit_Strip[Nesrh] = esdetid.strip(); 
  Nesrh++;
 }

 for(int iesrh=0;iesrh<Nesrh;iesrh++)
 {for(int jesrh=iesrh+1;jesrh<Nesrh;jesrh++)
  {
   if(
        (_esRecHit_siZ[iesrh]==_esRecHit_siZ[jesrh])
      &&(_esRecHit_siP[iesrh]==_esRecHit_siP[jesrh])
      &&(_esRecHit_siX[iesrh]==_esRecHit_siX[jesrh])
      &&(_esRecHit_siY[iesrh]==_esRecHit_siY[jesrh]) 
     )
   { _esRecHit_Noisy[iesrh]=1; _esRecHit_Noisy[jesrh]=1; }
 }}//end find noisy sensor / ambiguous for matching

}

void ESAlignTool::fill_tracks(edm::Handle<reco::TrackCollection> TrackCol)
{
 for(reco::TrackCollection::const_iterator itTrack = TrackCol->begin();
     itTrack != TrackCol->end(); ++itTrack)
 {    
  if(Ntrack>=2000)
  {
   edm::LogWarning("fill_tracks")<<"Too many selected tracks.\n";
   continue;
  }
  if ( itTrack->charge()!=0 )
  {
   if( pass_TrackSelection(itTrack) )
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

    if(itTrack->innerOk())
    {
     _TrackInnerX[Ntrack] = itTrack->innerPosition().X();
     _TrackInnerY[Ntrack] = itTrack->innerPosition().Y();
     _TrackInnerZ[Ntrack] = itTrack->innerPosition().Z();
    }
    Ntrack++;
   }//end pass_TrackSelection
  }//charge!=0
 }//TrackCollection
}

void ESAlignTool::fill_esRecHitAssociation()
{
 for(int itrk=0;itrk<Ntrack;itrk++)
 {
  for(int i=0;i<Nesrh;i++)
  {
   double r = deltaR(_TrackOuterEta[itrk],_TrackOuterPhi[itrk],_esRecHit_Eta[i],_esRecHit_Phi[i]);
   if(r<0.5) _esRecHit_MatchedTrk_fromOuter[i]=itrk;
  }
 }
}

void ESAlignTool::set_ESorigin(int _evt_run, const CaloGeometry *caloGeom)
{
 int iz=-1; int ip=1; int ix=0; int iy=0; int is=0;
 Double_t X=0.; Double_t Y=0.; Double_t Z=0.;
 ESDetId esid;
 if( Cal_ESorigin_from_Geometry && _evt_run==1)
 {
  std::cout<<"\nSetting ES local origins... \n";
  ip=1; iz=-1; //ES-Front
  ix=2;iy=15;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=2;iy=15;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=2;iy=26;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=2;iy=26;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=39;iy=15;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=39;iy=15;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=39;iy=26;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=39;iy=26;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=16;iy=1;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=16;iy=1;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=16;iy=40;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=16;iy=40;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=25;iy=1;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=25;iy=1;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=25;iy=40;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=25;iy=40;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  X/=16.;Y/=16.;Z/=16.;
  delete ESmF_O;  delete ESmF_Oap;
  ESmF_O=new PositionType(X,Y,Z);
  ESmF_Oap=new PositionType(X+ES_dX[0][0],Y+ES_dY[0][0],Z+ES_dZ[0][0]);
  std::cout<<"ESmF_O("<<ESmF_O->x()<<","<<ESmF_O->y()<<","<<ESmF_O->z()<<")\n";

//  ESmF_O_X=ESmF_O->x(); ESmF_O_Y=ESmF_O->y(); ESmF_O_Z=ESmF_O->z();
//  ESmF_Oap_X=ESmF_Oap->x(); ESmF_Oap_Y=ESmF_Oap->y(); ESmF_Oap_Z=ESmF_Oap->z();

  X=0.0; Y=0.0; Z=0.0;
  ip=1; iz=1; //ES+Front
  ix=2;iy=15;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=2;iy=15;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=2;iy=26;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=2;iy=26;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=39;iy=15;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=39;iy=15;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=39;iy=26;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=39;iy=26;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=16;iy=1;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=16;iy=1;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=16;iy=40;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=16;iy=40;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=25;iy=1;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=25;iy=1;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=25;iy=40;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=25;iy=40;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  X/=16.;Y/=16.;Z/=16.;
  delete ESpF_O;  delete ESpF_Oap;
  ESpF_O=new PositionType(X,Y,Z);
  ESpF_Oap=new PositionType(X+ES_dX[1][0],Y+ES_dY[1][0],Z+ES_dZ[1][0]);
  std::cout<<"ESpF_O("<<ESpF_O->x()<<","<<ESpF_O->y()<<","<<ESpF_O->z()<<")\n";
//  ESpF_O_X=ESpF_O->x(); ESpF_O_Y=ESpF_O->y(); ESpF_O_Z=ESpF_O->z();
//  ESpF_Oap_X=ESpF_Oap->x(); ESpF_Oap_Y=ESpF_Oap->y(); ESpF_Oap_Z=ESpF_Oap->z();

  X=0.0; Y=0.0; Z=0.0;
  ip=2; iz=-1; //ES-Rear
  ix=1;iy=25;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=1;iy=25;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=1;iy=16;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=1;iy=16;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=40;iy=25;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=40;iy=25;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=40;iy=16;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=40;iy=16;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=15;iy=2;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=15;iy=2;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=26;iy=2;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=26;iy=2;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=15;iy=39;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=15;iy=39;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=26;iy=39;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=26;iy=39;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  X/=16.;Y/=16.;Z/=16.;
  delete ESmR_O;  delete ESmR_Oap;
  ESmR_O=new PositionType(X,Y,Z);
  ESmR_Oap=new PositionType(X+ES_dX[0][1],Y+ES_dY[0][1],Z+ES_dZ[0][1]);
  std::cout<<"ESmR_O("<<ESmR_O->x()<<","<<ESmR_O->y()<<","<<ESmR_O->z()<<")\n";
//  ESmR_O_X=ESmR_O->x(); ESmR_O_Y=ESmR_O->y(); ESmR_O_Z=ESmR_O->z();
//  ESmR_Oap_X=ESmR_Oap->x(); ESmR_Oap_Y=ESmR_Oap->y(); ESmR_Oap_Z=ESmR_Oap->z();

  X=0.0; Y=0.0; Z=0.0;
  ip=2; iz=1; //ES+Rear
  ix=1;iy=25;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=1;iy=25;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=1;iy=16;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=1;iy=16;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=40;iy=25;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=40;iy=25;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=40;iy=16;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=40;iy=16;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=15;iy=2;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=15;iy=2;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=26;iy=2;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=26;iy=2;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=15;iy=39;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=15;iy=39;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=26;iy=39;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  ix=26;iy=39;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X, Y, Z);
  X/=16.;Y/=16.;Z/=16.;
  delete ESpR_O;  delete ESpR_Oap;
  ESpR_O=new PositionType(X,Y,Z);
  ESpR_Oap=new PositionType(X+ES_dX[1][1],Y+ES_dY[1][1],Z+ES_dZ[1][1]);
  std::cout<<"ESpR_O("<<ESpR_O->x()<<","<<ESpR_O->y()<<","<<ESpR_O->z()<<")\n";
//  ESpR_O_X=ESpR_O->x(); ESpR_O_Y=ESpR_O->y(); ESpR_O_Z=ESpR_O->z();
//  ESpR_Oap_X=ESpR_Oap->x(); ESpR_Oap_Y=ESpR_Oap->y(); ESpR_Oap_Z=ESpR_Oap->z();
  std::cout<<"\n\n";
 }//if( Cal_ESorigin_from_Geometry && _evt_run==1)

 if( _evt_run==1)
 {
  ES_O_X[0][0]=ESmF_O->x(); ES_O_Y[0][0]=ESmF_O->y(); ES_O_Z[0][0]=ESmF_O->z();
  ES_Oap_X[0][0]=ESmF_Oap->x(); ES_Oap_Y[0][0]=ESmF_Oap->y(); ES_Oap_Z[0][0]=ESmF_Oap->z();
  ES_O_X[0][1]=ESmR_O->x(); ES_O_Y[0][1]=ESmR_O->y(); ES_O_Z[0][1]=ESmR_O->z();
  ES_Oap_X[0][1]=ESmR_Oap->x(); ES_Oap_Y[0][1]=ESmR_Oap->y(); ES_Oap_Z[0][1]=ESmR_Oap->z();
  ES_O_X[1][0]=ESpF_O->x(); ES_O_Y[1][0]=ESpF_O->y(); ES_O_Z[1][0]=ESpF_O->z();
  ES_Oap_X[1][0]=ESpF_Oap->x(); ES_Oap_Y[1][0]=ESpF_Oap->y(); ES_Oap_Z[1][0]=ESpF_Oap->z();
  ES_O_X[1][1]=ESpR_O->x(); ES_O_Y[1][1]=ESpR_O->y(); ES_O_Z[1][1]=ESpR_O->z();
  ES_Oap_X[1][1]=ESpR_Oap->x(); ES_Oap_Y[1][1]=ESpR_Oap->y(); ES_Oap_Z[1][1]=ESpR_Oap->z();
 }
}

void ESAlignTool::set_ESaxes(int _evt_run, const CaloGeometry *caloGeom)
{
 int iz=-1; int ip=1; int ix=0; int iy=0; int is=0;
 Double_t X=0.; Double_t Y=0.; Double_t Z=0.;
 Double_t X1=0.; Double_t Y1=0.; Double_t Z1=0.;
 Double_t X2=0.; Double_t Y2=0.; Double_t Z2=0.;
 ESDetId esid;
 if( Cal_ESaxes_from_Geometry && _evt_run==1)
 {
  std::cout<<"\nSetting ES local axes... \n";
  ip=1; iz=-1; //ES-Front
  X1=0.;  Y1=0.;  Z1=0.;  X2=0.;  Y2=0.;  Z2=0.;
  ix=2;iy=15;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=2;iy=15;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=2;iy=26;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=2;iy=26;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=39;iy=15;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=39;iy=15;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=39;iy=26;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=39;iy=26;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  X=(X2-X1)/8.;  Y=(Y2-Y1)/8.;  Z=(Z2-Z1)/8.;   normalize(X,Y,Z);
  ES_O_R11[0][0]=X;  ES_O_R12[0][0]=Y;  ES_O_R13[0][0]=Z;

  X1=0.;  Y1=0.;  Z1=0.;  X2=0.;  Y2=0.;  Z2=0.;
  ix=16;iy=1;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=16;iy=1;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=25;iy=1;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=25;iy=1;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=16;iy=40;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=16;iy=40;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=25;iy=40;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=25;iy=40;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  X=(X2-X1)/8.;  Y=(Y2-Y1)/8.;  Z=(Z2-Z1)/8.;   normalize(X,Y,Z);
  ES_O_R21[0][0]=X;  ES_O_R22[0][0]=Y;  ES_O_R23[0][0]=Z;

  X=ES_O_R12[0][0]*ES_O_R23[0][0]-ES_O_R13[0][0]*ES_O_R22[0][0];
  Y=ES_O_R13[0][0]*ES_O_R21[0][0]-ES_O_R11[0][0]*ES_O_R23[0][0];
  Z=ES_O_R11[0][0]*ES_O_R22[0][0]-ES_O_R12[0][0]*ES_O_R21[0][0];
  normalize(X,Y,Z);
  ES_O_R31[0][0]=X;  ES_O_R32[0][0]=Y;  ES_O_R33[0][0]=Z;

  std::cout<<"ESmF_O_R:\n";
  std::cout<<ES_O_R11[0][0]<<", "<<ES_O_R12[0][0]<<", "<<ES_O_R13[0][0]<<")\n";
  std::cout<<ES_O_R21[0][0]<<", "<<ES_O_R22[0][0]<<", "<<ES_O_R23[0][0]<<")\n";
  std::cout<<ES_O_R31[0][0]<<", "<<ES_O_R32[0][0]<<", "<<ES_O_R33[0][0]<<")\n";

  ip=2; iz=-1; //ES-Rear
  X1=0.;  Y1=0.;  Z1=0.;  X2=0.;  Y2=0.;  Z2=0.;
  ix=1;iy=16;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=1;iy=16;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=1;iy=25;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=1;iy=25;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=40;iy=16;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=40;iy=16;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=40;iy=25;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=40;iy=25;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  X=(X2-X1)/8.;  Y=(Y2-Y1)/8.;  Z=(Z2-Z1)/8.;   normalize(X,Y,Z);
  ES_O_R11[0][1]=X;  ES_O_R12[0][1]=Y;  ES_O_R13[0][1]=Z;

  X1=0.;  Y1=0.;  Z1=0.;  X2=0.;  Y2=0.;  Z2=0.;
  ix=15;iy=2;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=15;iy=2;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=26;iy=2;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=26;iy=2;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=15;iy=39;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=15;iy=39;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=26;iy=39;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=26;iy=39;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  X=(X2-X1)/8.;  Y=(Y2-Y1)/8.;  Z=(Z2-Z1)/8.;   normalize(X,Y,Z);
  ES_O_R21[0][1]=X;  ES_O_R22[0][1]=Y;  ES_O_R23[0][1]=Z;

  X=ES_O_R12[0][1]*ES_O_R23[0][1]-ES_O_R13[0][1]*ES_O_R22[0][1];
  Y=ES_O_R13[0][1]*ES_O_R21[0][1]-ES_O_R11[0][1]*ES_O_R23[0][1];
  Z=ES_O_R11[0][1]*ES_O_R22[0][1]-ES_O_R12[0][1]*ES_O_R21[0][1];
  normalize(X,Y,Z);
  ES_O_R31[0][1]=X;  ES_O_R32[0][1]=Y;  ES_O_R33[0][1]=Z;

  std::cout<<"ESmR_O_R:\n";
  std::cout<<ES_O_R11[0][1]<<", "<<ES_O_R12[0][1]<<", "<<ES_O_R13[0][1]<<")\n";
  std::cout<<ES_O_R21[0][1]<<", "<<ES_O_R22[0][1]<<", "<<ES_O_R23[0][1]<<")\n";
  std::cout<<ES_O_R31[0][1]<<", "<<ES_O_R32[0][1]<<", "<<ES_O_R33[0][1]<<")\n";

  ip=1; iz=1; //ES+Front
  X1=0.;  Y1=0.;  Z1=0.;  X2=0.;  Y2=0.;  Z2=0.;
  ix=2;iy=15;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=2;iy=15;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=2;iy=26;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=2;iy=26;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=39;iy=15;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=39;iy=15;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=39;iy=26;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=39;iy=26;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  X=(X2-X1)/8.;  Y=(Y2-Y1)/8.;  Z=(Z2-Z1)/8.;   normalize(X,Y,Z);
  ES_O_R11[1][0]=X;  ES_O_R12[1][0]=Y;  ES_O_R13[1][0]=Z;

  X1=0.;  Y1=0.;  Z1=0.;  X2=0.;  Y2=0.;  Z2=0.;
  ix=16;iy=1;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=16;iy=1;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=25;iy=1;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=25;iy=1;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=16;iy=40;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=16;iy=40;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=25;iy=40;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=25;iy=40;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  X=(X2-X1)/8.;  Y=(Y2-Y1)/8.;  Z=(Z2-Z1)/8.;   normalize(X,Y,Z);
  ES_O_R21[1][0]=X;  ES_O_R22[1][0]=Y;  ES_O_R23[1][0]=Z;

  X=ES_O_R12[1][0]*ES_O_R23[1][0]-ES_O_R13[1][0]*ES_O_R22[1][0];
  Y=ES_O_R13[1][0]*ES_O_R21[1][0]-ES_O_R11[1][0]*ES_O_R23[1][0];
  Z=ES_O_R11[1][0]*ES_O_R22[1][0]-ES_O_R12[1][0]*ES_O_R21[1][0];
  normalize(X,Y,Z);
  ES_O_R31[1][0]=X;  ES_O_R32[1][0]=Y;  ES_O_R33[1][0]=Z;

  std::cout<<"ESpF_O_R:\n";
  std::cout<<ES_O_R11[1][0]<<", "<<ES_O_R12[1][0]<<", "<<ES_O_R13[1][0]<<")\n";
  std::cout<<ES_O_R21[1][0]<<", "<<ES_O_R22[1][0]<<", "<<ES_O_R23[1][0]<<")\n";
  std::cout<<ES_O_R31[1][0]<<", "<<ES_O_R32[1][0]<<", "<<ES_O_R33[1][0]<<")\n";

  ip=2; iz=1; //ES+Rear
  X1=0.;  Y1=0.;  Z1=0.;  X2=0.;  Y2=0.;  Z2=0.;
  ix=1;iy=16;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=1;iy=16;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=1;iy=25;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=1;iy=25;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=40;iy=16;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=40;iy=16;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=40;iy=25;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=40;iy=25;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  X=(X2-X1)/8.;  Y=(Y2-Y1)/8.;  Z=(Z2-Z1)/8.;   normalize(X,Y,Z);
  ES_O_R11[1][1]=X;  ES_O_R12[1][1]=Y;  ES_O_R13[1][1]=Z;

  X1=0.;  Y1=0.;  Z1=0.;  X2=0.;  Y2=0.;  Z2=0.;
  ix=15;iy=2;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=15;iy=2;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=26;iy=2;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=26;iy=2;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X1, Y1, Z1);
  ix=15;iy=39;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=15;iy=39;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=26;iy=39;is=1;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  ix=26;iy=39;is=32;  esid=ESDetId(is,ix,iy,ip,iz);
  Sum_ESPosition(caloGeom, esid, X2, Y2, Z2);
  X=(X2-X1)/8.;  Y=(Y2-Y1)/8.;  Z=(Z2-Z1)/8.;   normalize(X,Y,Z);
  ES_O_R21[1][1]=X;  ES_O_R22[1][1]=Y;  ES_O_R23[1][1]=Z;

  X=ES_O_R12[1][1]*ES_O_R23[1][1]-ES_O_R13[1][1]*ES_O_R22[1][1];
  Y=ES_O_R13[1][1]*ES_O_R21[1][1]-ES_O_R11[1][1]*ES_O_R23[1][1];
  Z=ES_O_R11[1][1]*ES_O_R22[1][1]-ES_O_R12[1][1]*ES_O_R21[1][1];
  normalize(X,Y,Z);
  ES_O_R31[1][1]=X;  ES_O_R32[1][1]=Y;  ES_O_R33[1][1]=Z;

  std::cout<<"ESpR_O_R:\n";
  std::cout<<ES_O_R11[1][1]<<", "<<ES_O_R12[1][1]<<", "<<ES_O_R13[1][1]<<")\n";
  std::cout<<ES_O_R21[1][1]<<", "<<ES_O_R22[1][1]<<", "<<ES_O_R23[1][1]<<")\n";
  std::cout<<ES_O_R31[1][1]<<", "<<ES_O_R32[1][1]<<", "<<ES_O_R33[1][1]<<")\n";

 }
}

void ESAlignTool::Sum_ESPosition(const CaloGeometry *caloGeom, ESDetId esid, Double_t &X, Double_t &Y, Double_t &Z)
{
  X += ( caloGeom->getPosition(esid).x() );
  Y += ( caloGeom->getPosition(esid).y() );
  Z += ( caloGeom->getPosition(esid).z() );
}

void ESAlignTool::normalize(Double_t &X, Double_t &Y, Double_t &Z)
{
 Double_t buf=X*X+Y*Y+Z*Z;
 buf=TMath::Sqrt(buf);
 X/=buf; Y/=buf; Z/=buf;
}

void ESAlignTool::Overwrite_RotationMatrix(int _evt_run)
{
 if( _evt_run==1 )
 {
  for(int i=0;i<2;i++)
  {
   for(int j=0;j<2;j++)
   {
    ES_R11[i][j]=ES_O_R11[i][j]; ES_R12[i][j]=ES_O_R12[i][j]; ES_R13[i][j]=ES_O_R13[i][j];
    ES_R21[i][j]=ES_O_R21[i][j]; ES_R22[i][j]=ES_O_R22[i][j]; ES_R23[i][j]=ES_O_R23[i][j];
    ES_R31[i][j]=ES_O_R31[i][j]; ES_R32[i][j]=ES_O_R32[i][j]; ES_R33[i][j]=ES_O_R33[i][j];
  }}
 }
}

void ESAlignTool::fill_PredictionState(int iz, int ip, int idee, edm::Handle<reco::TrackCollection> TrackCol, edm::ESHandle<MagneticField> theMagField, edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry)
{
 PositionType *ES_Oap;

 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
   int a,b;
  iz==-1 ? a=0 : a=1;
  b=ip-1;
  ES_Oap=new PositionType(ES_Oap_X[a][b],ES_Oap_Y[a][b],ES_Oap_Z[a][b]);
  SteppingHelixPropagator ESpropa(theMagField.product());
  Plane::PlanePointer Plane_ES = Plane::build(*ES_Oap,*woRotate);


  int iTrk=0;
  for(reco::TrackCollection::const_iterator itTrack = TrackCol->begin();
      itTrack != TrackCol->end(); ++itTrack)
  {    
   if(iTrk>=2000)
   {
    edm::LogWarning("fill_tracks")<<"Too many selected tracks.\n";
    continue;
   }
   if ( itTrack->charge()!=0 )
   {
    if( pass_TrackSelection(itTrack) )
    {
     reco::TransientTrack TTrack(*itTrack,theMagField.product(),theTrackingGeometry);
     TrajectoryStateOnSurface ipTSOS=TTrack.impactPointState();
     if(ipTSOS.isValid()&&TTrack.impactPointStateAvailable())
     {
      TrajectoryStateOnSurface ES_prediction=ESpropa.propagate(*(ipTSOS.freeState()),*Plane_ES);
      if( (ES_prediction.isValid())&&(ES_prediction.hasError()) )
      {
       PredictionState_iz[iTrk][a][b]=iz;
       PredictionState_ip[iTrk][a][b]=ip;
       PredictionState_valid[iTrk][a][b]=1;
       PredictionState_X[iTrk][a][b]=ES_prediction.globalPosition().x();
       PredictionState_Y[iTrk][a][b]=ES_prediction.globalPosition().y();
       PredictionState_Z[iTrk][a][b]=ES_prediction.globalPosition().z();
       PredictionState_Px[iTrk][a][b]=ES_prediction.globalMomentum().x();
       PredictionState_Py[iTrk][a][b]=ES_prediction.globalMomentum().y();
       PredictionState_Pz[iTrk][a][b]=ES_prediction.globalMomentum().z();
       PredictionState_Bx[iTrk][a][b]=theMagField->inTesla(ES_prediction.globalPosition()).x();
       PredictionState_By[iTrk][a][b]=theMagField->inTesla(ES_prediction.globalPosition()).y();
       PredictionState_Bz[iTrk][a][b]=theMagField->inTesla(ES_prediction.globalPosition()).z();
       PredictionState_Exx[iTrk][a][b]=ES_prediction.cartesianError().position().cxx();
       PredictionState_Eyx[iTrk][a][b]=ES_prediction.cartesianError().position().cyx();
       PredictionState_Eyy[iTrk][a][b]=ES_prediction.cartesianError().position().cyy();
       PredictionState_E44[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(3,3));
       PredictionState_E55[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(4,4));
       PredictionState_E66[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(5,5));
       PredictionState_E45[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(3,4));
       PredictionState_E46[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(3,5));
       PredictionState_E56[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(4,5));
       PredictionState_E14[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(0,3));
       PredictionState_E15[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(0,4));
       PredictionState_E16[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(0,5));
       PredictionState_E24[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(1,3));
       PredictionState_E25[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(1,4));
       PredictionState_E26[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(1,5));
       PredictionState_delX[iTrk][a][b]=PredictionState_Px[iTrk][a][b]/PredictionState_Pz[iTrk][a][b];
       PredictionState_delY[iTrk][a][b]=PredictionState_Py[iTrk][a][b]/PredictionState_Pz[iTrk][a][b];
      }//end if ES_prediction Valid and hasError
     }//end if innerTSOS.isValid()
     iTrk++;
    }//pass trk selections
   }//trk charge !=0
  }//loop trks
 }//if one of 4 planes
}

void ESAlignTool::fill_PredictionState_wRotation(int iz, int ip, int idee, edm::Handle<reco::TrackCollection> TrackCol, edm::ESHandle<MagneticField> theMagField, edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry)
{
 PositionType *ES_Oap;
 RotationType *ES_Rotation;

 Double_t R11,R12,R13,R21,R22,R23,R31,R32,R33;
 Double_t cxx,cyy,czz,cyx,czx,czy;

 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  int a,b;
  iz==-1 ? a=0 : a=1;
  b=ip-1;
  ES_Oap=new PositionType(ES_Oap_X[a][b],ES_Oap_Y[a][b],ES_Oap_Z[a][b]);
  ES_Rotation=new RotationType(ES_R11[a][b],ES_R12[a][b],ES_R13[a][b], ES_R21[a][b],ES_R22[a][b],ES_R23[a][b],ES_R31[a][b],ES_R32[a][b],ES_R33[a][b]);

  R11=ES_R11[a][b]; R12=ES_R12[a][b]; R13=ES_R13[a][b];
  R21=ES_R21[a][b]; R22=ES_R22[a][b]; R23=ES_R23[a][b];
  R31=ES_R31[a][b]; R32=ES_R32[a][b]; R33=ES_R33[a][b];

  if(b_PrintMatrix)
  {
   if(ESpF_Printed==0)
   {
    if(iz==1&&ip==1)
    {
     std::cout << "Printing Rotation Matrix :\n";
     std::cout << "ES+F :\n";
     std::cout << R11 << ", " << R12 << ", " << R13 << ",\n";
     std::cout << R21 << ", " << R22 << ", " << R23 << ",\n";
     std::cout << R31 << ", " << R32 << ", " << R33 << ",\n";
     ESpF_Printed=1;
    }
   }
   if(ESpR_Printed==0)
   {
    if(iz==1&&ip==2)
    {
     std::cout << "Printing Rotation Matrix :\n";
     std::cout << "ES+R :\n";
     std::cout << R11 << ", " << R12 << ", " << R13 << ",\n";
     std::cout << R21 << ", " << R22 << ", " << R23 << ",\n";
     std::cout << R31 << ", " << R32 << ", " << R33 << ",\n";
     ESpR_Printed=1;
    }
   }
   if(ESmF_Printed==0)
   {
    if(iz==-1&&ip==1)
    {
     std::cout << "Printing Rotation Matrix :\n";
     std::cout << "ES-F :\n";
     std::cout << R11 << ", " << R12 << ", " << R13 << ",\n";
     std::cout << R21 << ", " << R22 << ", " << R23 << ",\n";
     std::cout << R31 << ", " << R32 << ", " << R33 << ",\n";
     ESmF_Printed=1;
    }
   }
   if(ESmR_Printed==0)
   {
    if(iz==-1&&ip==2)
    {
     std::cout << "Printing Rotation Matrix :\n";
     std::cout << "ES-R :\n";
     std::cout << R11 << ", " << R12 << ", " << R13 << ",\n";
     std::cout << R21 << ", " << R22 << ", " << R23 << ",\n";
     std::cout << R31 << ", " << R32 << ", " << R33 << ",\n";
     ESmR_Printed=1;
    }
   }
  }

  SteppingHelixPropagator ESpropa(theMagField.product());
  Plane::PlanePointer Plane_ES = Plane::build(*ES_Oap,*ES_Rotation);

  int iTrk=0;
  for(reco::TrackCollection::const_iterator itTrack = TrackCol->begin();
      itTrack != TrackCol->end(); ++itTrack)
  {    
   if(iTrk>=2000)
   {
    edm::LogWarning("fill_tracks")<<"Too many selected tracks.\n";
    continue;
   }
   if ( itTrack->charge()!=0 )
   {
    if( pass_TrackSelection(itTrack) )
    {
     reco::TransientTrack TTrack(*itTrack,theMagField.product(),theTrackingGeometry);
     TrajectoryStateOnSurface ipTSOS=TTrack.impactPointState();
     if(ipTSOS.isValid()&&TTrack.impactPointStateAvailable())
     {
      TrajectoryStateOnSurface ES_prediction=ESpropa.propagate(*(ipTSOS.freeState()),*Plane_ES);
      if( (ES_prediction.isValid())&&(ES_prediction.hasError()) )
      {
       PredictionState_iz[iTrk][a][b]=iz;
       PredictionState_ip[iTrk][a][b]=ip;
       PredictionState_valid[iTrk][a][b]=1;
       PredictionState_X[iTrk][a][b]=ES_prediction.globalPosition().x();
       PredictionState_Y[iTrk][a][b]=ES_prediction.globalPosition().y();
       PredictionState_Z[iTrk][a][b]=ES_prediction.globalPosition().z();
       PredictionState_Px[iTrk][a][b]=ES_prediction.globalMomentum().x();
       PredictionState_Py[iTrk][a][b]=ES_prediction.globalMomentum().y();
       PredictionState_Pz[iTrk][a][b]=ES_prediction.globalMomentum().z();
       PredictionState_Bx[iTrk][a][b]=theMagField->inTesla(ES_prediction.globalPosition()).x();
       PredictionState_By[iTrk][a][b]=theMagField->inTesla(ES_prediction.globalPosition()).y();
       PredictionState_Bz[iTrk][a][b]=theMagField->inTesla(ES_prediction.globalPosition()).z();
 
       //Original position error matrix
       cxx=ES_prediction.cartesianError().position().cxx();
       cyy=ES_prediction.cartesianError().position().cyy();
       czz=ES_prediction.cartesianError().position().czz();
       cyx=ES_prediction.cartesianError().position().cyx();
       czx=ES_prediction.cartesianError().position().czx();
       czy=ES_prediction.cartesianError().position().czy();
       //Transform this error matrix on local coordinate
       PredictionState_Exx[iTrk][a][b]= R11*R11*cxx+R12*R12*cyy+R13*R13*czz+2.*R11*R12*cyx+2.*R11*R13*czx+2.*R12*R13*czy;
       PredictionState_Eyy[iTrk][a][b]= R21*R21*cxx+R22*R22*cyy+R23*R23*czz+2.*R21*R22*cyx+2.*R21*R23*czx+2.*R22*R23*czy;
       PredictionState_Eyx[iTrk][a][b]= R11*R21*cxx+R12*R22*cyy+R13*R23*czz+(R11*R22+R21*R12)*cyx+(R11*R23+R21*R13)*czx+(R12*R23+R22*R13)*czy;

       PredictionState_E44[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(3,3));
       PredictionState_E55[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(4,4));
       PredictionState_E66[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(5,5));
       PredictionState_E45[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(3,4));
       PredictionState_E46[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(3,5));
       PredictionState_E56[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(4,5));
       PredictionState_E14[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(0,3));
       PredictionState_E15[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(0,4));
       PredictionState_E16[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(0,5));
       PredictionState_E24[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(1,3));
       PredictionState_E25[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(1,4));
       PredictionState_E26[iTrk][a][b]=((ES_prediction.cartesianError().matrix())(1,5));
       PredictionState_delX[iTrk][a][b]=PredictionState_Px[iTrk][a][b]/PredictionState_Pz[iTrk][a][b];
       PredictionState_delY[iTrk][a][b]=PredictionState_Py[iTrk][a][b]/PredictionState_Pz[iTrk][a][b];
      }//end if ES_prediction Valid and hasError
     }//end if innerTSOS.isValid()
     iTrk++;
    }//pass trk selections
   }//trk charge !=0
  }//loop trks
 }//if one of 4 planes
}

void ESAlignTool::fill_residual(int iz)
{
 if(iz==1||iz==-1)
 {
  int a; iz==-1 ? a=0 : a=1;
  for(int iTrk=0;iTrk<Ntrack;iTrk++)
  {
   if(PredictionState_valid[iTrk][a][0]!=1) continue; //ES+/- F
   if(PredictionState_valid[iTrk][a][1]!=1) continue; //ES+/- R
   Double_t eF_xx= (PredictionState_Exx[iTrk][a][0]);
   Double_t eF_yy= (PredictionState_Eyy[iTrk][a][0]);
   Double_t eF_yx= (PredictionState_Eyx[iTrk][a][0]);
   //e_yx=0.;
   if( eF_xx*eF_yy-eF_yx*eF_yx <= 0 ) continue;
   if( eF_xx > e_xxlimit ) continue;
   if( eF_yy > e_yylimit ) continue;
   if( fabs(eF_yx) > e_yxlimit ) continue;
   eF_xx+=pow(0.055029,2.); eF_yy+=pow(1.760918,2.);
   Double_t eR_yy= (PredictionState_Eyy[iTrk][a][1]);
   Double_t eR_xx= (PredictionState_Exx[iTrk][a][1]);
   Double_t eR_yx= (PredictionState_Eyx[iTrk][a][1]);
   //e_yx=0.;
   if( eR_xx*eR_yy-eR_yx*eR_yx <= 0 ) continue;
   if( eR_xx > e_xxlimit ) continue;
   if( eR_yy > e_yylimit ) continue;
   if( fabs(eR_yx) > e_yxlimit ) continue;
   eR_yy+=pow(0.055029,2.); eR_xx+=pow(1.760918,2.);
   Double_t X0F=(PredictionState_X[iTrk][a][0]-ES_Oap_X[a][0]);
   Double_t Y0F=(PredictionState_Y[iTrk][a][0]-ES_Oap_Y[a][0]);
   if(check_DeadZone(iz,1,X0F,Y0F)==0) continue; 
   if(check_Radius(X0F,Y0F)==0) continue;
   Double_t X0R=(PredictionState_X[iTrk][a][1]-ES_Oap_X[a][1]);
   Double_t Y0R=(PredictionState_Y[iTrk][a][1]-ES_Oap_Y[a][1]);
   if(check_DeadZone(iz,2,X0R,Y0R)==0) continue;
   if(check_Radius(X0R,Y0R)==0) continue;

   Double_t disF=(winlimit*winlimit); int indF=-1;
   for(int irec=0;irec<Nesrh;irec++)
   {
    if(_esRecHit_siZ[irec]!=iz||_esRecHit_siP[irec]!=1) continue;
    if(_esRecHit_X[irec]==0.&&_esRecHit_Y[irec]==0.) continue;
    if(_esRecHit_MatchedTrk_fromOuter[irec]!=iTrk) continue;
    Double_t X=_esRecHit_X[irec]-ES_O_X[a][0];  
    Double_t Y=_esRecHit_Y[irec]-ES_O_Y[a][0];
    if( fabs(X-X0F)>winlimit ) continue;
    if( fabs(Y-Y0F)>winlimit ) continue;
    Double_t buf = pow(X-X0F,2.)+pow(Y-Y0F,2.);
    if(buf<disF)
    { indF=irec; disF=buf; }
   }//end for-loop ESrechit
   Double_t disR=(winlimit*winlimit); int indR=-1;
   for(int irec=0;irec<Nesrh;irec++)
   {
    if(_esRecHit_siZ[irec]!=iz||_esRecHit_siP[irec]!=2) continue;
    if(_esRecHit_X[irec]==0.&&_esRecHit_Y[irec]==0.) continue;
    if(_esRecHit_MatchedTrk_fromOuter[irec]!=iTrk) continue;
    Double_t X=_esRecHit_X[irec]-ES_O_X[a][1];  
    Double_t Y=_esRecHit_Y[irec]-ES_O_Y[a][1];
    if( fabs(X-X0R)>winlimit) continue;
    if( fabs(Y-Y0R)>winlimit) continue;
    Double_t buf = pow(X-X0R,2.)+pow(Y-Y0R,2.);
    if(buf<disR)
    { indR=irec; disR=buf; }
   }//end for-loop ESrechita

   if(indF>-1&&_esRecHit_Noisy[indF]==0&&indR>-1&&_esRecHit_Noisy[indR]==0
      && BadSensor(_esRecHit_siZ[indF],_esRecHit_siP[indF],_esRecHit_siX[indF],_esRecHit_siY[indF])==0
      && BadSensor(_esRecHit_siZ[indR],_esRecHit_siP[indR],_esRecHit_siX[indR],_esRecHit_siY[indR])==0
     )
   {
    PredictionState_MatchedRec[iTrk][a][0]=indF;
    PredictionState_resiX[iTrk][a][0]=PredictionState_X[iTrk][a][0]-ES_dX[a][0]-_esRecHit_X[indF];
    PredictionState_resiY[iTrk][a][0]=PredictionState_Y[iTrk][a][0]-ES_dY[a][0]-_esRecHit_Y[indF];
    if(a==1)
    {
     ESpF_residualX->Fill(PredictionState_resiX[iTrk][a][0]);
     ESpF_residualY->Fill(PredictionState_resiY[iTrk][a][0]);
    }
    if(a==0)
    {
     ESmF_residualX->Fill(PredictionState_resiX[iTrk][a][0]);
     ESmF_residualY->Fill(PredictionState_resiY[iTrk][a][0]);
    }
    Double_t XDF=PredictionState_resiX[iTrk][a][0];
    Double_t YDF=PredictionState_resiY[iTrk][a][0];
    Double_t determinant= eF_xx*eF_yy-eF_yx*eF_yx;
    Double_t derX=(PredictionState_delX[iTrk][a][0]);
    Double_t derY=(PredictionState_delY[iTrk][a][0]);

    Cal_MatrixM(iz,1,eF_xx,eF_yx,eF_yy,derX,derY,determinant);
    Cal_VectorP(iz,1,eF_xx,eF_yx,eF_yy,derX,derY,determinant,XDF,YDF);
    Cal_CHI2(iz,1,eF_xx,eF_yx,eF_yy,XDF,YDF);
    ES_NTracks[a][0] += 1;

    Double_t Exx=eF_xx;Double_t Eyx=eF_yx;Double_t Eyy=eF_yy;
    Double_t Px=PredictionState_Px[iTrk][a][0];
    Double_t Py=PredictionState_Py[iTrk][a][0];
    Double_t Pz=PredictionState_Pz[iTrk][a][0];
    Double_t sPxx=PredictionState_E44[iTrk][a][0];
    Double_t sPxy=PredictionState_E45[iTrk][a][0];
    Double_t sPxz=PredictionState_E46[iTrk][a][0];
    Double_t sPyy=PredictionState_E55[iTrk][a][0];
    Double_t sPyz=PredictionState_E56[iTrk][a][0];
    Double_t sPzz=PredictionState_E66[iTrk][a][0];
    Double_t sE14=PredictionState_E14[iTrk][a][0];
    Double_t sE15=PredictionState_E15[iTrk][a][0];
    Double_t sE16=PredictionState_E16[iTrk][a][0];
    Double_t sE24=PredictionState_E24[iTrk][a][0];
    Double_t sE25=PredictionState_E25[iTrk][a][0];
    Double_t sE26=PredictionState_E26[iTrk][a][0];
    Double_t Rxx=pow(0.055029,2.); Double_t Ryy=pow(1.760918,2.);
    Double_t Xpre=X0F;  Double_t Ypre=Y0F;
    Double_t Xrec=X0F-XDF;  Double_t Yrec=Y0F-YDF;

    Cal_MatrixMErr2(iz,1,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz);
    Cal_VectorPErr2(iz,1,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz,Rxx,Ryy,Xpre,Ypre,Xrec,Yrec,sE14,sE15,sE16,sE24,sE25,sE26,determinant);

    PredictionState_MatchedRec[iTrk][a][1]=indR;
    PredictionState_resiX[iTrk][a][1]=PredictionState_X[iTrk][a][1]-ES_dX[a][1]-_esRecHit_X[indR];
    PredictionState_resiY[iTrk][a][1]=PredictionState_Y[iTrk][a][1]-ES_dY[a][1]-_esRecHit_Y[indR];
    if(a==1)
    {
     ESpR_residualX->Fill(PredictionState_resiX[iTrk][a][1]);
     ESpR_residualY->Fill(PredictionState_resiY[iTrk][a][1]);
    }
    if(a==0)
    {
     ESmR_residualX->Fill(PredictionState_resiX[iTrk][a][1]);
     ESmR_residualY->Fill(PredictionState_resiY[iTrk][a][1]);
    }
    Double_t XDR=PredictionState_resiX[iTrk][a][1];
    Double_t YDR=PredictionState_resiY[iTrk][a][1];
    determinant= eR_xx*eR_yy-eR_yx*eR_yx;
    derX=(PredictionState_delX[iTrk][a][1]);
    derY=(PredictionState_delY[iTrk][a][1]);

    Cal_MatrixM(iz,2,eR_xx,eR_yx,eR_yy,derX,derY,determinant);
    Cal_VectorP(iz,2,eR_xx,eR_yx,eR_yy,derX,derY,determinant,XDR,YDR);
    Cal_CHI2(iz,2,eR_xx,eR_yx,eR_yy,XDR,YDR);
    ES_NTracks[a][1] += 1;

    Exx=eF_xx;Eyx=eF_yx;Eyy=eF_yy;
    Px=PredictionState_Px[iTrk][a][1];
    Py=PredictionState_Py[iTrk][a][1];
    Pz=PredictionState_Pz[iTrk][a][1];
    sPxx=PredictionState_E44[iTrk][a][1];
    sPxy=PredictionState_E45[iTrk][a][1];
    sPxz=PredictionState_E46[iTrk][a][1];
    sPyy=PredictionState_E55[iTrk][a][1];
    sPyz=PredictionState_E56[iTrk][a][1];
    sPzz=PredictionState_E66[iTrk][a][1];
    sE14=PredictionState_E14[iTrk][a][1];
    sE15=PredictionState_E15[iTrk][a][1];
    sE16=PredictionState_E16[iTrk][a][1];
    sE24=PredictionState_E24[iTrk][a][1];
    sE25=PredictionState_E25[iTrk][a][1];
    sE26=PredictionState_E26[iTrk][a][1];
    Ryy=pow(0.055029,2.); Rxx=pow(1.760918,2.);
    Xpre=X0F;  Ypre=Y0F;
    Xrec=X0F-XDF;  Yrec=Y0F-YDF;

    Cal_MatrixMErr2(iz,2,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz);
    Cal_VectorPErr2(iz,2,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz,Rxx,Ryy,Xpre,Ypre,Xrec,Yrec,sE14,sE15,sE16,sE24,sE25,sE26,determinant);

   }//end if good matching
  }//end for-loop Trk

 }//iz==1 or -1
}

void ESAlignTool::fill_residual_wRotation_v2(int iz)
{
 if(iz==1||iz==-1)
 {
  int a; iz==-1 ? a=0 : a=1;
  for(int iTrk=0;iTrk<Ntrack;iTrk++)
  {
   if(PredictionState_valid[iTrk][a][0]!=1) continue; //ES+/- F
   if(PredictionState_valid[iTrk][a][1]!=1) continue; //ES+/- R
   Double_t eF_xx= (PredictionState_Exx[iTrk][a][0]);
   Double_t eF_yy= (PredictionState_Eyy[iTrk][a][0]);
   Double_t eF_yx= (PredictionState_Eyx[iTrk][a][0]);
   //e_yx=0.;
   if( eF_xx*eF_yy-eF_yx*eF_yx <= 0 ) continue;
   if( eF_xx > e_xxlimit ) continue;
   if( eF_yy > e_yylimit ) continue;
   if( fabs(eF_yx) > e_yxlimit ) continue;
   eF_xx+=pow(0.055029,2.); eF_yy+=pow(1.760918,2.);
   Double_t eR_yy= (PredictionState_Eyy[iTrk][a][1]);
   Double_t eR_xx= (PredictionState_Exx[iTrk][a][1]);
   Double_t eR_yx= (PredictionState_Eyx[iTrk][a][1]);
   //e_yx=0.;
   if( eR_xx*eR_yy-eR_yx*eR_yx <= 0 ) continue;
   if( eR_xx > e_xxlimit ) continue;
   if( eR_yy > e_yylimit ) continue;
   if( fabs(eR_yx) > e_yxlimit ) continue;
   eR_yy+=pow(0.055029,2.); eR_xx+=pow(1.760918,2.);
   //Prediction Point on Local coordinate
   Double_t X0F,Y0F;
   X0F
   = (PredictionState_X[iTrk][a][0]-ES_Oap_X[a][0])*ES_R11[a][0]
    +(PredictionState_Y[iTrk][a][0]-ES_Oap_Y[a][0])*ES_R12[a][0]
    +(PredictionState_Z[iTrk][a][0]-ES_Oap_Z[a][0])*ES_R13[a][0];
   Y0F
   = (PredictionState_X[iTrk][a][0]-ES_Oap_X[a][0])*ES_R21[a][0]
    +(PredictionState_Y[iTrk][a][0]-ES_Oap_Y[a][0])*ES_R22[a][0]
    +(PredictionState_Z[iTrk][a][0]-ES_Oap_Z[a][0])*ES_R23[a][0];
   if(check_DeadZone(iz,1,X0F,Y0F)==0) continue; 
   if(check_Radius(X0F,Y0F)==0) continue;
   Double_t X0R,Y0R;
   X0R
   = (PredictionState_X[iTrk][a][1]-ES_Oap_X[a][1])*ES_R11[a][1]
    +(PredictionState_Y[iTrk][a][1]-ES_Oap_Y[a][1])*ES_R12[a][1]
    +(PredictionState_Z[iTrk][a][1]-ES_Oap_Z[a][1])*ES_R13[a][1];
   Y0R
   = (PredictionState_X[iTrk][a][1]-ES_Oap_X[a][1])*ES_R21[a][1]
    +(PredictionState_Y[iTrk][a][1]-ES_Oap_Y[a][1])*ES_R22[a][1]
    +(PredictionState_Z[iTrk][a][1]-ES_Oap_Z[a][1])*ES_R23[a][1];
   if(check_DeadZone(iz,2,X0R,Y0R)==0) continue;
   if(check_Radius(X0R,Y0R)==0) continue;

   Double_t disF=(winlimit*winlimit); int indF=-1;
   Double_t resiXF,resiYF,resiXR,resiYR;
   resiXF=30.; resiYF=30.; resiXR=30.; resiYR=30.;
   for(int irec=0;irec<Nesrh;irec++)
   {
    if(_esRecHit_siZ[irec]!=iz||_esRecHit_siP[irec]!=1) continue;
    if(_esRecHit_X[irec]==0.&&_esRecHit_Y[irec]==0.) continue;
    if(_esRecHit_MatchedTrk_fromOuter[irec]!=iTrk) continue;
    Double_t X,Y;
    X= (_esRecHit_X[irec]-ES_O_X[a][0])*ES_O_R11[a][0]
      +(_esRecHit_Y[irec]-ES_O_Y[a][0])*ES_O_R12[a][0]
      +(_esRecHit_Z[irec]-ES_O_Z[a][0])*ES_O_R13[a][0];
    Y= (_esRecHit_X[irec]-ES_O_X[a][0])*ES_O_R21[a][0]
      +(_esRecHit_Y[irec]-ES_O_Y[a][0])*ES_O_R22[a][0]
      +(_esRecHit_Z[irec]-ES_O_Z[a][0])*ES_O_R23[a][0];
    if( fabs(X-X0F)>winlimit ) continue;
    if( fabs(Y-Y0F)>winlimit ) continue;
    Double_t buf = pow(X-X0F,2.)+pow(Y-Y0F,2.);
    if(buf<disF)
    { indF=irec; disF=buf; resiXF=X0F-X; resiYF=Y0F-Y; }
   }//end for-loop ESrechit
   Double_t disR=(winlimit*winlimit); int indR=-1;
   for(int irec=0;irec<Nesrh;irec++)
   {
    if(_esRecHit_siZ[irec]!=iz||_esRecHit_siP[irec]!=2) continue;
    if(_esRecHit_X[irec]==0.&&_esRecHit_Y[irec]==0.) continue;
    if(_esRecHit_MatchedTrk_fromOuter[irec]!=iTrk) continue;
    Double_t X,Y;
    X= (_esRecHit_X[irec]-ES_O_X[a][1])*ES_O_R11[a][1]
      +(_esRecHit_Y[irec]-ES_O_Y[a][1])*ES_O_R12[a][1]
      +(_esRecHit_Z[irec]-ES_O_Z[a][1])*ES_O_R13[a][1];
    Y= (_esRecHit_X[irec]-ES_O_X[a][1])*ES_O_R21[a][1]
      +(_esRecHit_Y[irec]-ES_O_Y[a][1])*ES_O_R22[a][1]
      +(_esRecHit_Z[irec]-ES_O_Z[a][1])*ES_O_R23[a][1];
    if( fabs(X-X0R)>winlimit) continue;
    if( fabs(Y-Y0R)>winlimit) continue;
    Double_t buf = pow(X-X0R,2.)+pow(Y-Y0R,2.);
    if(buf<disR)
    { indR=irec; disR=buf; resiXR=X0R-X; resiYR=Y0R-Y; }
   }//end for-loop ESrechita

   if(indF>-1&&_esRecHit_Noisy[indF]==0&&indR>-1&&_esRecHit_Noisy[indR]==0
      && BadSensor(_esRecHit_siZ[indF],_esRecHit_siP[indF],_esRecHit_siX[indF],_esRecHit_siY[indF])==0
      && BadSensor(_esRecHit_siZ[indR],_esRecHit_siP[indR],_esRecHit_siX[indR],_esRecHit_siY[indR])==0
     )
   {
    PredictionState_MatchedRec[iTrk][a][0]=indF;
    PredictionState_resiX[iTrk][a][0]=resiXF;
    PredictionState_resiY[iTrk][a][0]=resiYF;
    if(a==1)
    {
     ESpF_residualX->Fill(PredictionState_resiX[iTrk][a][0]);
     ESpF_residualY->Fill(PredictionState_resiY[iTrk][a][0]);
    }
    if(a==0)
    {
     ESmF_residualX->Fill(PredictionState_resiX[iTrk][a][0]);
     ESmF_residualY->Fill(PredictionState_resiY[iTrk][a][0]);
    }
    Double_t XDF=PredictionState_resiX[iTrk][a][0];
    Double_t YDF=PredictionState_resiY[iTrk][a][0];
    Double_t determinant= eF_xx*eF_yy-eF_yx*eF_yx;
    Double_t derX=(PredictionState_delX[iTrk][a][0]);
    Double_t derY=(PredictionState_delY[iTrk][a][0]);

    if(  Selected_idee==0
       || (Selected_idee==1&&_esRecHit_X[indF]-ES_O_X[a][0]>6.1 )
       || (Selected_idee==2&&_esRecHit_X[indF]-ES_O_X[a][0]<-6.1 )
      )
    {
     Cal_MatrixM_wRotation(iTrk,iz,1,eF_xx,eF_yx,eF_yy);
     Cal_VectorP_wRotation(iTrk,iz,1,eF_xx,eF_yx,eF_yy);
     Cal_CHI2(iz,1,eF_xx,eF_yx,eF_yy,XDF,YDF);
     ES_NTracks[a][0] += 1;
    }

    /*
    Double_t Exx=eF_xx;Double_t Eyx=eF_yx;Double_t Eyy=eF_yy;
    Double_t Px=PredictionState_Px[iTrk][a][0];
    Double_t Py=PredictionState_Py[iTrk][a][0];
    Double_t Pz=PredictionState_Pz[iTrk][a][0];
    Double_t sPxx=PredictionState_E44[iTrk][a][0];
    Double_t sPxy=PredictionState_E45[iTrk][a][0];
    Double_t sPxz=PredictionState_E46[iTrk][a][0];
    Double_t sPyy=PredictionState_E55[iTrk][a][0];
    Double_t sPyz=PredictionState_E56[iTrk][a][0];
    Double_t sPzz=PredictionState_E66[iTrk][a][0];
    Double_t sE14=PredictionState_E14[iTrk][a][0];
    Double_t sE15=PredictionState_E15[iTrk][a][0];
    Double_t sE16=PredictionState_E16[iTrk][a][0];
    Double_t sE24=PredictionState_E24[iTrk][a][0];
    Double_t sE25=PredictionState_E25[iTrk][a][0];
    Double_t sE26=PredictionState_E26[iTrk][a][0];
    Double_t Rxx=pow(0.055029,2.); Double_t Ryy=pow(1.760918,2.);
    Double_t Xpre=X0F;  Double_t Ypre=Y0F;
    Double_t Xrec=X0F-XDF;  Double_t Yrec=Y0F-YDF;

    Cal_MatrixMErr2(iz,1,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz);
    Cal_VectorPErr2(iz,1,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz,Rxx,Ryy,Xpre,Ypre,Xrec,Yrec,sE14,sE15,sE16,sE24,sE25,sE26,determinant);
    */

    PredictionState_MatchedRec[iTrk][a][1]=indR;
    PredictionState_resiX[iTrk][a][1]=resiXR;
    PredictionState_resiY[iTrk][a][1]=resiYR;

    if(a==1)
    {
     ESpR_residualX->Fill(PredictionState_resiX[iTrk][a][1]);
     ESpR_residualY->Fill(PredictionState_resiY[iTrk][a][1]);
    }
    if(a==0)
    {
     ESmR_residualX->Fill(PredictionState_resiX[iTrk][a][1]);
     ESmR_residualY->Fill(PredictionState_resiY[iTrk][a][1]);
    }
    Double_t XDR=PredictionState_resiX[iTrk][a][1];
    Double_t YDR=PredictionState_resiY[iTrk][a][1];
    determinant= eR_xx*eR_yy-eR_yx*eR_yx;
    derX=(PredictionState_delX[iTrk][a][1]);
    derY=(PredictionState_delY[iTrk][a][1]);

    if(  Selected_idee==0
       || (Selected_idee==1&&_esRecHit_Y[indR]-ES_O_Y[a][0]>6.1 )
       || (Selected_idee==2&&_esRecHit_Y[indR]-ES_O_Y[a][0]<-6.1 )
      )
    {
     Cal_MatrixM_wRotation(iTrk,iz,2,eR_xx,eR_yx,eR_yy);
     Cal_VectorP_wRotation(iTrk,iz,2,eR_xx,eR_yx,eR_yy);
     Cal_CHI2(iz,2,eR_xx,eR_yx,eR_yy,XDR,YDR);
     ES_NTracks[a][1] += 1;
    }

    /*
    Exx=eF_xx;Eyx=eF_yx;Eyy=eF_yy;
    Px=PredictionState_Px[iTrk][a][1];
    Py=PredictionState_Py[iTrk][a][1];
    Pz=PredictionState_Pz[iTrk][a][1];
    sPxx=PredictionState_E44[iTrk][a][1];
    sPxy=PredictionState_E45[iTrk][a][1];
    sPxz=PredictionState_E46[iTrk][a][1];
    sPyy=PredictionState_E55[iTrk][a][1];
    sPyz=PredictionState_E56[iTrk][a][1];
    sPzz=PredictionState_E66[iTrk][a][1];
    sE14=PredictionState_E14[iTrk][a][1];
    sE15=PredictionState_E15[iTrk][a][1];
    sE16=PredictionState_E16[iTrk][a][1];
    sE24=PredictionState_E24[iTrk][a][1];
    sE25=PredictionState_E25[iTrk][a][1];
    sE26=PredictionState_E26[iTrk][a][1];
    Ryy=pow(0.055029,2.); Rxx=pow(1.760918,2.);
    Xpre=X0F;  Ypre=Y0F;
    Xrec=X0F-XDF;  Yrec=Y0F-YDF;

    Cal_MatrixMErr2(iz,2,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz);
    Cal_VectorPErr2(iz,2,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz,Rxx,Ryy,Xpre,Ypre,Xrec,Yrec,sE14,sE15,sE16,sE24,sE25,sE26,determinant);
    */

   }//end if good matching
  }//end for-loop Trk

 }//iz==1 or -1
}

void ESAlignTool::fill_residual_wRotation(int iz)
{
 if(iz==1||iz==-1)
 {
  int a; iz==-1 ? a=0 : a=1;
  for(int iTrk=0;iTrk<Ntrack;iTrk++)
  {
   if(PredictionState_valid[iTrk][a][0]!=1) continue; //ES+/- F
   if(PredictionState_valid[iTrk][a][1]!=1) continue; //ES+/- R
   Double_t eF_xx= (PredictionState_Exx[iTrk][a][0]);
   Double_t eF_yy= (PredictionState_Eyy[iTrk][a][0]);
   Double_t eF_yx= (PredictionState_Eyx[iTrk][a][0]);
   //e_yx=0.;
   if( eF_xx*eF_yy-eF_yx*eF_yx <= 0 ) continue;
   if( eF_xx > e_xxlimit ) continue;
   if( eF_yy > e_yylimit ) continue;
   if( fabs(eF_yx) > e_yxlimit ) continue;
   eF_xx+=pow(0.055029,2.); eF_yy+=pow(1.760918,2.);
   Double_t eR_yy= (PredictionState_Eyy[iTrk][a][1]);
   Double_t eR_xx= (PredictionState_Exx[iTrk][a][1]);
   Double_t eR_yx= (PredictionState_Eyx[iTrk][a][1]);
   //e_yx=0.;
   if( eR_xx*eR_yy-eR_yx*eR_yx <= 0 ) continue;
   if( eR_xx > e_xxlimit ) continue;
   if( eR_yy > e_yylimit ) continue;
   if( fabs(eR_yx) > e_yxlimit ) continue;
   eR_yy+=pow(0.055029,2.); eR_xx+=pow(1.760918,2.);
   //Prediction Point on Local coordinate
   Double_t X0F,Y0F;
   X0F
   = (PredictionState_X[iTrk][a][0]-ES_Oap_X[a][0])*ES_R11[a][0]
    +(PredictionState_Y[iTrk][a][0]-ES_Oap_Y[a][0])*ES_R12[a][0]
    +(PredictionState_Z[iTrk][a][0]-ES_Oap_Z[a][0])*ES_R13[a][0];
   Y0F
   = (PredictionState_X[iTrk][a][0]-ES_Oap_X[a][0])*ES_R21[a][0]
    +(PredictionState_Y[iTrk][a][0]-ES_Oap_Y[a][0])*ES_R22[a][0]
    +(PredictionState_Z[iTrk][a][0]-ES_Oap_Z[a][0])*ES_R23[a][0];
   if(check_DeadZone(iz,1,X0F,Y0F)==0) continue; 
   if(check_Radius(X0F,Y0F)==0) continue;
   Double_t X0R,Y0R;
   X0R
   = (PredictionState_X[iTrk][a][1]-ES_Oap_X[a][1])*ES_R11[a][1]
    +(PredictionState_Y[iTrk][a][1]-ES_Oap_Y[a][1])*ES_R12[a][1]
    +(PredictionState_Z[iTrk][a][1]-ES_Oap_Z[a][1])*ES_R13[a][1];
   Y0R
   = (PredictionState_X[iTrk][a][1]-ES_Oap_X[a][1])*ES_R21[a][1]
    +(PredictionState_Y[iTrk][a][1]-ES_Oap_Y[a][1])*ES_R22[a][1]
    +(PredictionState_Z[iTrk][a][1]-ES_Oap_Z[a][1])*ES_R23[a][1];
   if(check_DeadZone(iz,2,X0R,Y0R)==0) continue;
   if(check_Radius(X0R,Y0R)==0) continue;

   Double_t disF=(winlimit*winlimit); int indF=-1;
   for(int irec=0;irec<Nesrh;irec++)
   {
    if(_esRecHit_siZ[irec]!=iz||_esRecHit_siP[irec]!=1) continue;
    if(_esRecHit_X[irec]==0.&&_esRecHit_Y[irec]==0.) continue;
    if(_esRecHit_MatchedTrk_fromOuter[irec]!=iTrk) continue;
    Double_t X=_esRecHit_X[irec]-ES_O_X[a][0];  
    Double_t Y=_esRecHit_Y[irec]-ES_O_Y[a][0];
    if( fabs(X-X0F)>winlimit ) continue;
    if( fabs(Y-Y0F)>winlimit ) continue;
    Double_t buf = pow(X-X0F,2.)+pow(Y-Y0F,2.);
    if(buf<disF)
    { indF=irec; disF=buf; }
   }//end for-loop ESrechit
   Double_t disR=(winlimit*winlimit); int indR=-1;
   for(int irec=0;irec<Nesrh;irec++)
   {
    if(_esRecHit_siZ[irec]!=iz||_esRecHit_siP[irec]!=2) continue;
    if(_esRecHit_X[irec]==0.&&_esRecHit_Y[irec]==0.) continue;
    if(_esRecHit_MatchedTrk_fromOuter[irec]!=iTrk) continue;
    Double_t X=_esRecHit_X[irec]-ES_O_X[a][1];  
    Double_t Y=_esRecHit_Y[irec]-ES_O_Y[a][1];
    if( fabs(X-X0R)>winlimit) continue;
    if( fabs(Y-Y0R)>winlimit) continue;
    Double_t buf = pow(X-X0R,2.)+pow(Y-Y0R,2.);
    if(buf<disR)
    { indR=irec; disR=buf; }
   }//end for-loop ESrechita

   if(indF>-1&&_esRecHit_Noisy[indF]==0&&indR>-1&&_esRecHit_Noisy[indR]==0
      && BadSensor(_esRecHit_siZ[indF],_esRecHit_siP[indF],_esRecHit_siX[indF],_esRecHit_siY[indF])==0
      && BadSensor(_esRecHit_siZ[indR],_esRecHit_siP[indR],_esRecHit_siX[indR],_esRecHit_siY[indR])==0
     )
   {
    PredictionState_MatchedRec[iTrk][a][0]=indF;
    PredictionState_resiX[iTrk][a][0]=X0F-(_esRecHit_X[indF]-ES_O_X[a][0]);
    PredictionState_resiY[iTrk][a][0]=Y0F-(_esRecHit_Y[indF]-ES_O_Y[a][0]);
    if(a==1)
    {
     ESpF_residualX->Fill(PredictionState_resiX[iTrk][a][0]);
     ESpF_residualY->Fill(PredictionState_resiY[iTrk][a][0]);
    }
    if(a==0)
    {
     ESmF_residualX->Fill(PredictionState_resiX[iTrk][a][0]);
     ESmF_residualY->Fill(PredictionState_resiY[iTrk][a][0]);
    }
    Double_t XDF=PredictionState_resiX[iTrk][a][0];
    Double_t YDF=PredictionState_resiY[iTrk][a][0];
    Double_t determinant= eF_xx*eF_yy-eF_yx*eF_yx;
    Double_t derX=(PredictionState_delX[iTrk][a][0]);
    Double_t derY=(PredictionState_delY[iTrk][a][0]);

    if(  Selected_idee==0
       || (Selected_idee==1&&_esRecHit_X[indF]-ES_O_X[a][0]>6.1 )
       || (Selected_idee==2&&_esRecHit_X[indF]-ES_O_X[a][0]<-6.1 )
      )
    {
     Cal_MatrixM_wRotation(iTrk,iz,1,eF_xx,eF_yx,eF_yy);
     Cal_VectorP_wRotation(iTrk,iz,1,eF_xx,eF_yx,eF_yy);
     Cal_CHI2(iz,1,eF_xx,eF_yx,eF_yy,XDF,YDF);
     ES_NTracks[a][0] += 1;
    }

    /*
    Double_t Exx=eF_xx;Double_t Eyx=eF_yx;Double_t Eyy=eF_yy;
    Double_t Px=PredictionState_Px[iTrk][a][0];
    Double_t Py=PredictionState_Py[iTrk][a][0];
    Double_t Pz=PredictionState_Pz[iTrk][a][0];
    Double_t sPxx=PredictionState_E44[iTrk][a][0];
    Double_t sPxy=PredictionState_E45[iTrk][a][0];
    Double_t sPxz=PredictionState_E46[iTrk][a][0];
    Double_t sPyy=PredictionState_E55[iTrk][a][0];
    Double_t sPyz=PredictionState_E56[iTrk][a][0];
    Double_t sPzz=PredictionState_E66[iTrk][a][0];
    Double_t sE14=PredictionState_E14[iTrk][a][0];
    Double_t sE15=PredictionState_E15[iTrk][a][0];
    Double_t sE16=PredictionState_E16[iTrk][a][0];
    Double_t sE24=PredictionState_E24[iTrk][a][0];
    Double_t sE25=PredictionState_E25[iTrk][a][0];
    Double_t sE26=PredictionState_E26[iTrk][a][0];
    Double_t Rxx=pow(0.055029,2.); Double_t Ryy=pow(1.760918,2.);
    Double_t Xpre=X0F;  Double_t Ypre=Y0F;
    Double_t Xrec=X0F-XDF;  Double_t Yrec=Y0F-YDF;

    Cal_MatrixMErr2(iz,1,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz);
    Cal_VectorPErr2(iz,1,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz,Rxx,Ryy,Xpre,Ypre,Xrec,Yrec,sE14,sE15,sE16,sE24,sE25,sE26,determinant);
    */

    PredictionState_MatchedRec[iTrk][a][1]=indR;
    PredictionState_resiX[iTrk][a][1]=X0R-(_esRecHit_X[indR]-ES_O_X[a][1]);
    PredictionState_resiY[iTrk][a][1]=Y0R-(_esRecHit_Y[indR]-ES_O_Y[a][1]);

    if(a==1)
    {
     ESpR_residualX->Fill(PredictionState_resiX[iTrk][a][1]);
     ESpR_residualY->Fill(PredictionState_resiY[iTrk][a][1]);
    }
    if(a==0)
    {
     ESmR_residualX->Fill(PredictionState_resiX[iTrk][a][1]);
     ESmR_residualY->Fill(PredictionState_resiY[iTrk][a][1]);
    }
    Double_t XDR=PredictionState_resiX[iTrk][a][1];
    Double_t YDR=PredictionState_resiY[iTrk][a][1];
    determinant= eR_xx*eR_yy-eR_yx*eR_yx;
    derX=(PredictionState_delX[iTrk][a][1]);
    derY=(PredictionState_delY[iTrk][a][1]);

    if(  Selected_idee==0
       || (Selected_idee==1&&_esRecHit_Y[indR]-ES_O_Y[a][0]>6.1 )
       || (Selected_idee==2&&_esRecHit_Y[indR]-ES_O_Y[a][0]<-6.1 )
      )
    {
     Cal_MatrixM_wRotation(iTrk,iz,2,eR_xx,eR_yx,eR_yy);
     Cal_VectorP_wRotation(iTrk,iz,2,eR_xx,eR_yx,eR_yy);
     Cal_CHI2(iz,2,eR_xx,eR_yx,eR_yy,XDR,YDR);
     ES_NTracks[a][1] += 1;
    }

    /*
    Exx=eF_xx;Eyx=eF_yx;Eyy=eF_yy;
    Px=PredictionState_Px[iTrk][a][1];
    Py=PredictionState_Py[iTrk][a][1];
    Pz=PredictionState_Pz[iTrk][a][1];
    sPxx=PredictionState_E44[iTrk][a][1];
    sPxy=PredictionState_E45[iTrk][a][1];
    sPxz=PredictionState_E46[iTrk][a][1];
    sPyy=PredictionState_E55[iTrk][a][1];
    sPyz=PredictionState_E56[iTrk][a][1];
    sPzz=PredictionState_E66[iTrk][a][1];
    sE14=PredictionState_E14[iTrk][a][1];
    sE15=PredictionState_E15[iTrk][a][1];
    sE16=PredictionState_E16[iTrk][a][1];
    sE24=PredictionState_E24[iTrk][a][1];
    sE25=PredictionState_E25[iTrk][a][1];
    sE26=PredictionState_E26[iTrk][a][1];
    Ryy=pow(0.055029,2.); Rxx=pow(1.760918,2.);
    Xpre=X0F;  Ypre=Y0F;
    Xrec=X0F-XDF;  Yrec=Y0F-YDF;

    Cal_MatrixMErr2(iz,2,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz);
    Cal_VectorPErr2(iz,2,Exx,Eyx,Eyy,Px,Py,Pz,sPxx,sPyy,sPzz,sPxy,sPxz,sPyz,Rxx,Ryy,Xpre,Ypre,Xrec,Yrec,sE14,sE15,sE16,sE24,sE25,sE26,determinant);
    */

   }//end if good matching
  }//end for-loop Trk

 }//iz==1 or -1
}

int ESAlignTool::check_Radius(Double_t X, Double_t Y)
{
 int res=0;
 if(  (X*X+Y*Y)>60.*60. && (X*X+Y*Y)<110.*110.) res=1;
 return res;
}


void 
ESAlignTool::beginJob()
{
  std::cout << "In ESAlignTool.beginJob\n";

  f=new TFile("AlignmentFile.root","RECREATE");
  t_ESAlign = new TTree("ESAlign","tree",0);

  t_ESAlign->Branch("runNumber", &_runNum, "runNumber/L");
  t_ESAlign->Branch("evtNumber", &_evtNum, "evtNumber/L");
  
  t_ESAlign->Branch("Nesrh",    &Nesrh,       "Nesrh/I"); 
  t_ESAlign->Branch("_esRecHit_E",   &_esRecHit_E[0],   "_esRecHit_E[Nesrh]/D"); 
  t_ESAlign->Branch("_esRecHit_X",  &_esRecHit_X[0],  "_esRecHit_X[Nesrh]/D");
  t_ESAlign->Branch("_esRecHit_Y", &_esRecHit_Y[0], "_esRecHit_Y[Nesrh]/D");
  t_ESAlign->Branch("_esRecHit_Z", &_esRecHit_Z[0], "_esRecHit_Z[Nesrh]/D");
  t_ESAlign->Branch("_esRecHit_Eta",  &_esRecHit_Eta[0],  "_esRecHit_Eta[Nesrh]/D");
  t_ESAlign->Branch("_esRecHit_Phi",  &_esRecHit_Phi[0],  "_esRecHit_Phi[Nesrh]/D");
  t_ESAlign->Branch("_esRecHit_siZ", &_esRecHit_siZ[0], "_esRecHit_siZ[Nesrh]/I");
  t_ESAlign->Branch("_esRecHit_siP", &_esRecHit_siP[0], "_esRecHit_siP[Nesrh]/I");
  t_ESAlign->Branch("_esRecHit_siX", &_esRecHit_siX[0], "_esRecHit_siX[Nesrh]/IL");
  t_ESAlign->Branch("_esRecHit_siY", &_esRecHit_siY[0], "_esRecHit_siY[Nesrh]/IL");
  t_ESAlign->Branch("_esRecHit_Strip", &_esRecHit_Strip[0], "_esRecHit_Strip[Nesrh]/I");
  t_ESAlign->Branch("_esRecHit_Noisy",_esRecHit_Noisy,"_esRecHit_Noisy[Nesrh]/S");
  t_ESAlign->Branch("_esRecHit_MatchedTrk_fromOuter",_esRecHit_MatchedTrk_fromOuter,"_esRecHit_MatchedTrk_fromOuter[Nesrh]/I");
  t_ESAlign->Branch("Ntrack"          , &Ntrack          , "Ntrack/I");
  t_ESAlign->Branch("TrackPt"         , &_TrackPt[0]      , "TrackPt[Ntrack]/D");
  t_ESAlign->Branch("TrackEta"        , &_TrackEta[0]     , "TrackEta[Ntrack]/D");
  t_ESAlign->Branch("TrackPhi"        , &_TrackPhi[0]     , "TrackPhi[Ntrack]/D");
  t_ESAlign->Branch("TrackVx"         , &_TrackVx[0]      , "TrackVx[Ntrack]/D");
  t_ESAlign->Branch("TrackVy"         , &_TrackVy[0]      , "TrackVy[Ntrack]/D");
  t_ESAlign->Branch("TrackVz"         , &_TrackVz[0]      , "TrackVz[Ntrack]/D");
  t_ESAlign->Branch("Trackd0"     , &_Trackd0[0]  , "Trackd0[Ntrack]/D");
  t_ESAlign->Branch("TrackInnerX"     , &_TrackInnerX[0]  , "TrackInnerX[Ntrack]/D");
  t_ESAlign->Branch("TrackInnerY"     , &_TrackInnerY[0]  , "TrackInnerY[Ntrack]/D");
  t_ESAlign->Branch("TrackInnerZ"     , &_TrackInnerZ[0]  , "TrackInnerZ[Ntrack]/D");
  t_ESAlign->Branch("TrackOuterZ"     , &_TrackOuterZ[0]  , "TrackOuterZ[Ntrack]/D");
  t_ESAlign->Branch("TrackOuterEta"     , &_TrackOuterEta[0]  , "TrackOuterEta[Ntrack]/D");
  t_ESAlign->Branch("TrackOuterPhi"     , &_TrackOuterPhi[0]  , "TrackOuterPhi[Ntrack]/D");
  t_ESAlign->Branch("TrackCharge"     , &_TrackCharge[0]  , "TrackCharge[Ntrack]/S");
  t_ESAlign->Branch("TrackNHit"     , &_TrackNHit[0]  , "TrackNHit[Ntrack]/I");
  t_ESAlign->Branch("TrackNChi2"     , &_TrackNChi2[0]  , "TrackNChi2[Ntrack]/D");
  t_ESAlign->Branch("TrackPtError",_TrackPtError,"TrackPtError[Ntrack]/D");
  t_ESAlign->Branch("TrackQuality",_TrackQuality,"TrackQuality[Ntrack]/I");

  t_ESAlign->Branch("PredictionState_iz",PredictionState_iz,"PredictionState_iz[Ntrack][2][2]/I");
  t_ESAlign->Branch("PredictionState_ip",PredictionState_ip,"PredictionState_ip[Ntrack][2][2]/I");
  t_ESAlign->Branch("PredictionState_valid",PredictionState_valid,"PredictionState_valid[Ntrack][2][2]/S");
  t_ESAlign->Branch("PredictionState_X",PredictionState_X,"PredictionState_X[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Y",PredictionState_Y,"PredictionState_Y[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Z",PredictionState_Z,"PredictionState_Z[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Px",PredictionState_Px,"PredictionState_Px[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Py",PredictionState_Py,"PredictionState_Py[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Pz",PredictionState_Pz,"PredictionState_Pz[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Bx",PredictionState_Bx,"PredictionState_Bx[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_By",PredictionState_By,"PredictionState_By[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Bz",PredictionState_Bz,"PredictionState_Bz[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Exx",PredictionState_Exx,"PredictionState_Exx[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Eyx",PredictionState_Eyx,"PredictionState_Eyx[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_Eyy",PredictionState_Eyy,"PredictionState_Eyy[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_delX",PredictionState_delX,"PredictionState_delX[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_delY",PredictionState_delY,"PredictionState_delY[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_MatchedRec",PredictionState_MatchedRec,"PredictionState_MatchedRec[Ntrack][2][2]/I");
  t_ESAlign->Branch("PredictionState_resiX",PredictionState_resiX,"PredictionState_resiX[Ntrack][2][2]/D");
  t_ESAlign->Branch("PredictionState_resiY",PredictionState_resiY,"PredictionState_resiY[Ntrack][2][2]/D");

  t_ESAlign->Branch("ES_CHI2",ES_CHI2,"ES_CHI2[2][2]/D");
  t_ESAlign->Branch("ES_NTracks",ES_NTracks,"ES_NTracks[2][2]/L");
  t_ESAlign->Branch("ES_M11",ES_M11,"ES_M11[2][2]/D");
  t_ESAlign->Branch("ES_M12",ES_M12,"ES_M12[2][2]/D");
  t_ESAlign->Branch("ES_M13",ES_M13,"ES_M13[2][2]/D");
  t_ESAlign->Branch("ES_M14",ES_M14,"ES_M14[2][2]/D");
  t_ESAlign->Branch("ES_M15",ES_M15,"ES_M15[2][2]/D");
  t_ESAlign->Branch("ES_M16",ES_M16,"ES_M16[2][2]/D");
  t_ESAlign->Branch("ES_M22",ES_M22,"ES_M22[2][2]/D");
  t_ESAlign->Branch("ES_M23",ES_M23,"ES_M23[2][2]/D");
  t_ESAlign->Branch("ES_M24",ES_M24,"ES_M24[2][2]/D");
  t_ESAlign->Branch("ES_M25",ES_M25,"ES_M25[2][2]/D");
  t_ESAlign->Branch("ES_M26",ES_M26,"ES_M26[2][2]/D");
  t_ESAlign->Branch("ES_M33",ES_M33,"ES_M33[2][2]/D");
  t_ESAlign->Branch("ES_M34",ES_M34,"ES_M34[2][2]/D");
  t_ESAlign->Branch("ES_M35",ES_M35,"ES_M35[2][2]/D");
  t_ESAlign->Branch("ES_M36",ES_M36,"ES_M36[2][2]/D");
  t_ESAlign->Branch("ES_M44",ES_M44,"ES_M44[2][2]/D");
  t_ESAlign->Branch("ES_M45",ES_M45,"ES_M45[2][2]/D");
  t_ESAlign->Branch("ES_M46",ES_M46,"ES_M46[2][2]/D");
  t_ESAlign->Branch("ES_M55",ES_M55,"ES_M55[2][2]/D");
  t_ESAlign->Branch("ES_M56",ES_M56,"ES_M56[2][2]/D");
  t_ESAlign->Branch("ES_M66",ES_M66,"ES_M66[2][2]/D");
  t_ESAlign->Branch("ES_P1",ES_P1,"ES_P1[2][2]/D");
  t_ESAlign->Branch("ES_P2",ES_P2,"ES_P2[2][2]/D");
  t_ESAlign->Branch("ES_P3",ES_P3,"ES_P3[2][2]/D");
  t_ESAlign->Branch("ES_P4",ES_P4,"ES_P4[2][2]/D");
  t_ESAlign->Branch("ES_P5",ES_P5,"ES_P5[2][2]/D");
  t_ESAlign->Branch("ES_P6",ES_P6,"ES_P6[2][2]/D");
  t_ESAlign->Branch("ES_M31Err2",ES_M31Err2,"ES_M31Err2[2][2]/D");
  t_ESAlign->Branch("ES_M32Err2",ES_M32Err2,"ES_M32Err2[2][2]/D");
  t_ESAlign->Branch("ES_M33Err2",ES_M33Err2,"ES_M33Err2[2][2]/D");
  t_ESAlign->Branch("ES_P1Err2",ES_P1Err2,"ES_P1Err2[2][2]/D");
  t_ESAlign->Branch("ES_P2Err2",ES_P2Err2,"ES_P2Err2[2][2]/D");
  t_ESAlign->Branch("ES_P3Err2",ES_P3Err2,"ES_P3Err2[2][2]/D");

  if(b_DrawMagField)
  {
   t_ESField = new TTree("ESField","tree",0);
   t_ESField->Branch("ESpF_B_x",ESpF_B_x,"ESpF_B_x[40][40][32]/D");
   t_ESField->Branch("ESpF_B_y",ESpF_B_y,"ESpF_B_y[40][40][32]/D");
   t_ESField->Branch("ESpF_B_z",ESpF_B_z,"ESpF_B_z[40][40][32]/D");
   t_ESField->Branch("ESpR_B_x",ESpR_B_x,"ESpR_B_x[40][40][32]/D");
   t_ESField->Branch("ESpR_B_y",ESpR_B_y,"ESpR_B_y[40][40][32]/D");
   t_ESField->Branch("ESpR_B_z",ESpR_B_z,"ESpR_B_z[40][40][32]/D");
   t_ESField->Branch("ESmF_B_x",ESmF_B_x,"ESmF_B_x[40][40][32]/D");
   t_ESField->Branch("ESmF_B_y",ESmF_B_y,"ESmF_B_y[40][40][32]/D");
   t_ESField->Branch("ESmF_B_z",ESmF_B_z,"ESmF_B_z[40][40][32]/D");
   t_ESField->Branch("ESmR_B_x",ESmR_B_x,"ESmR_B_x[40][40][32]/D");
   t_ESField->Branch("ESmR_B_y",ESmR_B_y,"ESmR_B_y[40][40][32]/D");
   t_ESField->Branch("ESmR_B_z",ESmR_B_z,"ESmR_B_z[40][40][32]/D");

  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ESAlignTool::endJob() {
  std::cout << "In ESAlignTool.endJob\n";

  f->cd();
  t_ESAlign->Write(); 
  if(b_DrawMagField)  t_ESField->Write();

  ESpF_residualX->Write();  ESpF_residualY->Write();
  ESpR_residualX->Write();  ESpR_residualY->Write();
  ESmF_residualX->Write();  ESmF_residualY->Write();
  ESmR_residualX->Write();  ESmR_residualY->Write();


  std::cout << std::endl;
  std::cout << " --------------------------------------------- " << std::endl;
  std::cout << " number of events processed  " << _evt_run << std::endl; 
  std::cout << " Last Event number # " << _evtNum << std::endl; 
  std::cout << " --------------------------------------------- " << std::endl;


  delete woRotate;
  delete ESpF_wRotateap; delete ESpR_wRotateap;
  delete ESmF_wRotateap; delete ESmR_wRotateap;
  delete ESpF_O;  delete ESpR_O;  delete ESmF_O;  delete ESmR_O;
  delete ESpF_Oap;  delete ESpR_Oap;  delete ESmF_Oap;  delete ESmR_Oap;
 
  delete ESpF_residualX;  delete ESpF_residualY;
  delete ESpR_residualX;  delete ESpR_residualY;
  delete ESmF_residualX;  delete ESmF_residualY;
  delete ESmR_residualX;  delete ESmR_residualY;

/*
  cout<<"Alignment on ES+Front :\n";
  if(GetIteration(1,1,ESpFdX,ESpFdY,ESpFdZ)==1)
  {
   cout<<"dX="<<ESpFdX<<", dY="<<ESpFdY<<", dZ="<<ESpFdZ<<"\n";
   GetIterationError(1,1,ESpFdXerr,ESpFdYerr,ESpFdZerr);
   cout<<"dXerr="<<ESpFdXerr<<", dYerr="<<ESpFdYerr<<", dZerr="<<ESpFdZerr<<"\n";
   cout<<"normalized CHI2 before iteration="<<ESpFCHI2/ESpFNTracks<<" ";
   cout<<"# of Tracks="<<ESpFNTracks<<"\n";
  }
  cout<<"Alignment on ES+Rear :\n";
  if(GetIteration(1,2,ESpRdX,ESpRdY,ESpRdZ)==1)
  {
   cout<<"dX="<<ESpRdX<<", dY="<<ESpRdY<<", dZ="<<ESpRdZ<<"\n";
   GetIterationError(1,2,ESpRdXerr,ESpRdYerr,ESpRdZerr);
   cout<<"dXerr="<<ESpRdXerr<<", dYerr="<<ESpRdYerr<<", dZerr="<<ESpRdZerr<<"\n";
   cout<<"normalized CHI2 before iteration="<<ESpRCHI2/ESpRNTracks<<" ";
   cout<<"# of Tracks="<<ESpRNTracks<<"\n";
  }
  cout<<"Alignment on ES-Front :\n";
  if(GetIteration(-1,1,ESmFdX,ESmFdY,ESmFdZ)==1)
  {
   cout<<"dX="<<ESmFdX<<", dY="<<ESmFdY<<", dZ="<<ESmFdZ<<"\n";
   GetIterationError(-1,1,ESmFdXerr,ESmFdYerr,ESmFdZerr);
   cout<<"dXerr="<<ESmFdXerr<<", dYerr="<<ESmFdYerr<<", dZerr="<<ESmFdZerr<<"\n";
   cout<<"normalized CHI2 before iteration="<<ESmFCHI2/ESmFNTracks<<" ";
   cout<<"# of Tracks="<<ESmFNTracks<<"\n";
  }
  cout<<"Alignment on ES-Rear :\n";
  if(GetIteration(-1,2,ESmRdX,ESmRdY,ESmRdZ)==1)
  {
   cout<<"dX="<<ESmRdX<<", dY="<<ESmRdY<<", dZ="<<ESmRdZ<<"\n";
   GetIterationError(-1,2,ESmRdXerr,ESmRdYerr,ESmRdZerr);
   cout<<"dXerr="<<ESmRdXerr<<", dYerr="<<ESmRdYerr<<", dZerr="<<ESmRdZerr<<"\n";
   cout<<"normalized CHI2 before iteration="<<ESmRCHI2/ESmRNTracks<<" ";
   cout<<"# of Tracks="<<ESmRNTracks<<"\n";
  }
*/


/*
  ofstream file;
  file.open("file.txt",ios::out|ios::app);
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESpFdX = cms.double("<<ESpFdX<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESpFdY = cms.double("<<ESpFdY<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESpFdZ = cms.double("<<ESpFdZ<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESpRdX = cms.double("<<ESpRdX<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESpRdY = cms.double("<<ESpRdY<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESpRdZ = cms.double("<<ESpRdZ<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESmFdX = cms.double("<<ESmFdX<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESmFdY = cms.double("<<ESmFdY<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESmFdZ = cms.double("<<ESmFdZ<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESmRdX = cms.double("<<ESmRdX<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESmRdY = cms.double("<<ESmRdY<<")\n";
  file<<"process.patAlignTool3.Iter"<<iterN<<"_ESmRdZ = cms.double("<<ESmRdZ<<")\n";
  file<<"\n";
  file.close();

  ofstream file2;
  file2.open("file2.txt",ios::out|ios::app);
  file2<<"idx="<<iterN-1<<";\n";
  file2<<"ESpFdX-="<<ESpFdX<<"; ESpFdY-="<<ESpFdY<<"; ESpFdZ-="<<ESpFdZ<<";\n";
  file2<<"ESpFdXArray[idx]=ESpFdX; ESpFdYArray[idx]=ESpFdY; ESpFdZArray[idx]=ESpFdZ;\n";
  file2<<"ESpFchisqArray[idx]="<<ESpFCHI2/ESpFNTracks<<"; //";
  file2<<ESpFNTracks<<"\n";

  file2<<"ESpRdX-="<<ESpRdX<<"; ESpRdY-="<<ESpRdY<<"; ESpRdZ-="<<ESpRdZ<<";\n";
  file2<<"ESpRdXArray[idx]=ESpRdX; ESpRdYArray[idx]=ESpRdY; ESpRdZArray[idx]=ESpRdZ;\n";
  file2<<"ESpRchisqArray[idx]="<<ESpRCHI2/ESpRNTracks<<"; //";
  file2<<ESpRNTracks<<"\n";

  file2<<"ESmFdX-="<<ESmFdX<<"; ESmFdY-="<<ESmFdY<<"; ESmFdZ-="<<ESmFdZ<<";\n";
  file2<<"ESmFdXArray[idx]=ESmFdX; ESmFdYArray[idx]=ESmFdY; ESmFdZArray[idx]=ESmFdZ;\n";
  file2<<"ESmFchisqArray[idx]="<<ESmFCHI2/ESmFNTracks<<"; //";
  file2<<ESmFNTracks<<"\n";

  file2<<"ESmRdX-="<<ESmRdX<<"; ESmRdY-="<<ESmRdY<<"; ESmRdZ-="<<ESmRdZ<<";\n";
  file2<<"ESmRdXArray[idx]=ESmRdX; ESmRdYArray[idx]=ESmRdY; ESmRdZArray[idx]=ESmRdZ;\n";
  file2<<"ESmRchisqArray[idx]="<<ESmRCHI2/ESmRNTracks<<"; //";
  file2<<ESmRNTracks<<"\n";


  file2<<"\n";
  file2.close();
*/
 t_ESAlign->Delete();
 if(b_DrawMagField)  t_ESField->Delete();

 f->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(ESAlignTool);
