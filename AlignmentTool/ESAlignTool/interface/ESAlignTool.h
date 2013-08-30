#ifndef AlignmentTool_ESAlignTool_h
#define AlignmentTool_ESAlignTool_h



// -*- C++ -*-
//// -*- C++ -*-
//
// Package:    ESAlignTool
// Class:      ESAlignTool
//
#include <memory>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"


#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
//#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
//#include "DataFormats/GeometrySurface/interface/Surface.h"
//#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"

// ROOT include files
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>

//
// class declaration
//

class ESAlignTool : public edm::EDAnalyzer  // public edm::EDAnalyzer
{
public:
  explicit ESAlignTool(const edm::ParameterSet&);
  virtual ~ESAlignTool();
    // beginJob
  virtual void beginJob() ;
  // produce is where the ntuples are made
  virtual void analyze(const edm::Event &, const edm::EventSetup & );
  // endJob
  virtual void endJob() ;


protected:
  void initAllPara(const edm::ParameterSet & );
  void init_perEvent();
  void init_RotationMatrices();
  void LoadMagFieldonES(Long64_t ,const CaloGeometry *, edm::ESHandle<MagneticField>);
  void PrintPosition(Long64_t , const CaloGeometry *);
  void fill_esRecHit(const CaloGeometry *, edm::Handle<EcalRecHitCollection>);
  void fill_tracks(edm::Handle<reco::TrackCollection>);
  bool pass_TrackSelection(reco::TrackCollection::const_iterator);
  void fill_esRecHitAssociation();
  void set_ESorigin(int, const CaloGeometry *);
  void set_ESaxes(int, const CaloGeometry *);
  void Sum_ESPosition(const CaloGeometry *, ESDetId, Double_t &, Double_t &, Double_t &);
  void Overwrite_RotationMatrix(int);
  void normalize(Double_t &, Double_t &, Double_t &);
  void fill_PredictionState(int, int, int, edm::Handle<reco::TrackCollection>, edm::ESHandle<MagneticField>, edm::ESHandle<GlobalTrackingGeometry>);
  void fill_PredictionState_wRotation(int, int, int, edm::Handle<reco::TrackCollection>, edm::ESHandle<MagneticField>, edm::ESHandle<GlobalTrackingGeometry>);
  void fill_residual(int);
  void fill_residual_wRotation(int);
  void fill_residual_wRotation_v2(int);
 
  // The main sub-object which does the real work
  // Verbosity
  Long64_t _evt_run; 
  Long64_t _runNum, _evtNum;
  Int_t Nesrh; 
  Double_t _esRecHit_E[10000];
  Double_t _esRecHit_X[10000];  Double_t _esRecHit_Y[10000]; Double_t _esRecHit_Z[10000];
  Int_t _esRecHit_siZ[10000];  Int_t _esRecHit_siP[10000]; 
  Int_t _esRecHit_siX[10000];  Int_t _esRecHit_siY[10000]; 
  Double_t _esRecHit_Eta[10000]; Double_t _esRecHit_Phi[10000];
  Int_t _esRecHit_Strip[10000]; 
  Short_t _esRecHit_Noisy[10000]; 
  Int_t _esRecHit_MatchedTrk_fromOuter[10000];

  Int_t Ntrack; 
  Double_t _TrackPt[2000], _TrackEta[2000], _TrackPhi[2000]; 
  Double_t _TrackVx[2000], _TrackVy[2000], _TrackVz[2000]; 
  Double_t _TrackInnerX[2000], _TrackInnerY[2000], _TrackInnerZ[2000];
  Double_t _TrackOuterZ[2000]; 
  Double_t _TrackOuterEta[2000]; 
  Double_t _TrackOuterPhi[2000]; 
  Int_t _TrackNHit[2000]; 
  Double_t _TrackNChi2[2000];  
  Short_t _TrackCharge[2000]; 
  Double_t _Trackd0[2000];
  Double_t _TrackPtError[2000]; 
  Int_t _TrackQuality[2000]; 

  TFile *f;
  TTree *t_ESAlign;
  TTree *t_ESField;

  Double_t ESpF_B_x[40][40][32];
  Double_t ESpF_B_y[40][40][32];
  Double_t ESpF_B_z[40][40][32];
  Double_t ESpR_B_x[40][40][32];
  Double_t ESpR_B_y[40][40][32];
  Double_t ESpR_B_z[40][40][32];
  Double_t ESmF_B_x[40][40][32];
  Double_t ESmF_B_y[40][40][32];
  Double_t ESmF_B_z[40][40][32];
  Double_t ESmR_B_x[40][40][32];
  Double_t ESmR_B_y[40][40][32];
  Double_t ESmR_B_z[40][40][32];

  Int_t PredictionState_iz[2000][2][2];
  Int_t PredictionState_ip[2000][2][2];
  Short_t PredictionState_valid[2000][2][2];
  Double_t PredictionState_X[2000][2][2];
  Double_t PredictionState_Y[2000][2][2];
  Double_t PredictionState_Z[2000][2][2];
  Double_t PredictionState_Px[2000][2][2];
  Double_t PredictionState_Py[2000][2][2];
  Double_t PredictionState_Pz[2000][2][2];
  Double_t PredictionState_Bx[2000][2][2];
  Double_t PredictionState_By[2000][2][2];
  Double_t PredictionState_Bz[2000][2][2];
  Double_t PredictionState_Exx[2000][2][2];
  Double_t PredictionState_Eyx[2000][2][2];
  Double_t PredictionState_Eyy[2000][2][2];
  Double_t PredictionState_delX[2000][2][2];
  Double_t PredictionState_delY[2000][2][2];
  Double_t PredictionState_E44[2000][2][2];
  Double_t PredictionState_E45[2000][2][2];
  Double_t PredictionState_E46[2000][2][2];
  Double_t PredictionState_E55[2000][2][2];
  Double_t PredictionState_E56[2000][2][2];
  Double_t PredictionState_E66[2000][2][2];
  Double_t PredictionState_E14[2000][2][2];
  Double_t PredictionState_E15[2000][2][2];
  Double_t PredictionState_E16[2000][2][2];
  Double_t PredictionState_E24[2000][2][2];
  Double_t PredictionState_E25[2000][2][2];
  Double_t PredictionState_E26[2000][2][2];
  Int_t PredictionState_MatchedRec[2000][2][2];
  Double_t PredictionState_resiX[2000][2][2];
  Double_t PredictionState_resiY[2000][2][2];


  Double_t ES_M11[2][2],ES_M12[2][2],ES_M13[2][2],ES_M14[2][2],ES_M15[2][2],ES_M16[2][2],ES_M22[2][2],ES_M23[2][2],ES_M24[2][2],ES_M25[2][2],ES_M26[2][2],ES_M33[2][2],ES_M34[2][2],ES_M35[2][2],ES_M36[2][2],ES_M44[2][2],ES_M45[2][2],ES_M46[2][2],ES_M55[2][2],ES_M56[2][2],ES_M66[2][2];
  Double_t ES_P1[2][2],ES_P2[2][2],ES_P3[2][2],ES_P4[2][2],ES_P5[2][2],ES_P6[2][2];
  Double_t ES_CHI2[2][2];
  Long64_t ES_NTracks[2][2];
  Double_t ES_dX[2][2],ES_dY[2][2],ES_dZ[2][2];
  Double_t ES_dXerr[2][2],ES_dYerr[2][2],ES_dZerr[2][2];
  Double_t ES_dAlpha[2][2],ES_dBeta[2][2],ES_dGamma[2][2];
  Double_t ES_dAlphaerr[2][2],ES_dBetaerr[2][2],ES_dGammaerr[2][2];
  Double_t ES_O_Alpha[2][2],ES_O_Beta[2][2],ES_O_Gamma[2][2];
  Double_t ES_M31Err2[2][2],ES_M32Err2[2][2],ES_M33Err2[2][2];
  Double_t ES_P1Err2[2][2],ES_P2Err2[2][2],ES_P3Err2[2][2];
  Double_t ES_InvM11Err2[2][2],ES_InvM12Err2[2][2],ES_InvM13Err2[2][2],ES_InvM22Err2[2][2],ES_InvM23Err2[2][2],ES_InvM33Err2[2][2];


  typedef Surface::PositionType                  PositionType;
  typedef Surface::RotationType                  RotationType;
  RotationType *woRotate;
  PositionType *ESpF_O,*ESpR_O,*ESmF_O,*ESmR_O;
  PositionType *ESpF_Oap,*ESpR_Oap,*ESmF_Oap,*ESmR_Oap;
  RotationType *ESpF_wRotateO,*ESpR_wRotateO,*ESmF_wRotateO,*ESmR_wRotateO;
  RotationType *ESpF_wRotateap,*ESpR_wRotateap,*ESmF_wRotateap,*ESmR_wRotateap;
  //Rotation Matrices
  Double_t ES_R11[2][2],ES_R12[2][2],ES_R13[2][2];
  Double_t ES_R21[2][2],ES_R22[2][2],ES_R23[2][2];
  Double_t ES_R31[2][2],ES_R32[2][2],ES_R33[2][2];
  Double_t ES_O_R11[2][2],ES_O_R12[2][2],ES_O_R13[2][2];
  Double_t ES_O_R21[2][2],ES_O_R22[2][2],ES_O_R23[2][2];
  Double_t ES_O_R31[2][2],ES_O_R32[2][2],ES_O_R33[2][2];


  Double_t ES_O_X[2][2];  Double_t ES_O_Y[2][2];  Double_t ES_O_Z[2][2];
  Double_t ES_Oap_X[2][2];  Double_t ES_Oap_Y[2][2];  Double_t ES_Oap_Z[2][2];

  TH1D *ESpF_residualX;   TH1D *ESpF_residualY;
  TH1D *ESpR_residualX;   TH1D *ESpR_residualY;
  TH1D *ESmF_residualX;   TH1D *ESmF_residualY;
  TH1D *ESmR_residualX;   TH1D *ESmR_residualY;

private:
  /*double detSym(double,double,double,double,double,double);
  double GetIterX(double,double,double,double,double,double,double,double,double);
  double GetIterY(double,double,double,double,double,double,double,double,double);
  double GetIterZ(double,double,double,double,double,double,double,double,double);

  int GetIteration(int,int,double &,double &,double &);
  void setInvMError(int,int,double,double,double,double,double,double,double,double,double);
  int GetIterationError(int,int,double &,double &,double &);*/

  int check_Radius(Double_t, Double_t);
  int check_DeadZone(int, int, Double_t, Double_t);

  int BadSensor(int,int,int,int);

  void Cal_MatrixM(int,int,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t);
  void Cal_MatrixM_wRotation(int,int,int,Double_t,Double_t,Double_t);
  void Cal_VectorP(int,int,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t);
  void Cal_VectorP_wRotation(int,int,int,Double_t,Double_t,Double_t);
  void Cal_CHI2(int,int,Double_t,Double_t,Double_t,Double_t,Double_t);
  void Cal_MatrixMErr2(int,int,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t);
  void Cal_VectorPErr2(int,int,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t);

  Double_t e_xxlimit;
  Double_t e_yylimit;
  Double_t e_yxlimit;
  Double_t winlimit;

  Double_t iter_ESpFdX[11], iter_ESpFdY[11], iter_ESpFdZ[11];
  Double_t iter_ESpRdX[11], iter_ESpRdY[11], iter_ESpRdZ[11];
  Double_t iter_ESmFdX[11], iter_ESmFdY[11], iter_ESmFdZ[11];
  Double_t iter_ESmRdX[11], iter_ESmRdY[11], iter_ESmRdZ[11];

  Double_t iter_ESpFdAlpha[11], iter_ESpFdBeta[11], iter_ESpFdGamma[11];
  Double_t iter_ESpRdAlpha[11], iter_ESpRdBeta[11], iter_ESpRdGamma[11];
  Double_t iter_ESmFdAlpha[11], iter_ESmFdBeta[11], iter_ESmFdGamma[11];
  Double_t iter_ESmRdAlpha[11], iter_ESmRdBeta[11], iter_ESmRdGamma[11];


  int iterN;

  bool b_wRotation;
  bool b_DrawMagField;   bool b_PrintPosition;
  bool b_PrintMatrix; bool b_ReSetRfromOutside;
  bool b_fromRefitter;
  int ESpF_Printed;  int ESpR_Printed;
  int ESmF_Printed;  int ESmR_Printed;
  bool Cal_ESorigin_from_Geometry;
  bool Cal_ESaxes_from_Geometry;
  bool b_Overwrite_RotationMatrix_fromGeometry;
  int Selected_idee;
  int Selected_RUNmin;  int Selected_RUNmax;
};

#endif

