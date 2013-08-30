//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 27 11:07:26 2010 by ROOT version 5.22/00d
// from TTree ESAlign/tree
// found on file: AlignmentFile.root
//////////////////////////////////////////////////////////

#ifndef ESAlign_h
#define ESAlign_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ESAlign {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Long64_t        runNumber;
   Long64_t        evtNumber;
   Int_t           Nesrh;
   Double_t        _esRecHit_E[10000];   //[Nesrh]
   Double_t        _esRecHit_X[10000];   //[Nesrh]
   Double_t        _esRecHit_Y[10000];   //[Nesrh]
   Double_t        _esRecHit_Z[10000];   //[Nesrh]
   Double_t        _esRecHit_Eta[10000];   //[Nesrh]
   Double_t        _esRecHit_Phi[10000];   //[Nesrh]
   Int_t           _esRecHit_siZ[10000];   //[Nesrh]
   Int_t           _esRecHit_siP[10000];   //[Nesrh]
   Int_t           _esRecHit_siX[10000];   //[Nesrh]
   Int_t           _esRecHit_siY[10000];   //[Nesrh]
   Int_t           _esRecHit_Strip[10000];   //[Nesrh]
   Short_t         _esRecHit_Noisy[10000];   //[Nesrh]
   Int_t           _esRecHit_MatchedTrk_fromOuter[10000];   //[Nesrh]
   Int_t           Ntrack;
   Double_t        TrackPt[2000];   //[Ntrack]
   Double_t        TrackEta[2000];   //[Ntrack]
   Double_t        TrackPhi[2000];   //[Ntrack]
   Double_t        TrackVx[2000];   //[Ntrack]
   Double_t        TrackVy[2000];   //[Ntrack]
   Double_t        TrackVz[2000];   //[Ntrack]
   Double_t        Trackd0[2000];   //[Ntrack]
   Double_t        TrackInnerX[2000];   //[Ntrack]
   Double_t        TrackInnerY[2000];   //[Ntrack]
   Double_t        TrackInnerZ[2000];   //[Ntrack]
   Double_t        TrackOuterZ[2000];   //[Ntrack]
   Double_t        TrackOuterEta[2000];   //[Ntrack]
   Double_t        TrackOuterPhi[2000];   //[Ntrack]
   Short_t         TrackCharge[2000];   //[Ntrack]
   Int_t           TrackNHit[2000];   //[Ntrack]
   Double_t        TrackNChi2[2000];   //[Ntrack]
   Double_t        TrackPtError[2000];   //[Ntrack]
   Int_t           TrackQuality[2000];   //[Ntrack]
   Int_t           PredictionState_iz[2000][2][2];   //[Ntrack]
   Int_t           PredictionState_ip[2000][2][2];   //[Ntrack]
   Short_t         PredictionState_valid[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_X[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Y[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Z[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Px[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Py[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Pz[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Bx[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_By[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Bz[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Exx[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Eyx[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_Eyy[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_delX[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_delY[2000][2][2];   //[Ntrack]
   Int_t           PredictionState_MatchedRec[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_resiX[2000][2][2];   //[Ntrack]
   Double_t        PredictionState_resiY[2000][2][2];   //[Ntrack]
   Double_t        ES_CHI2[2][2];
   Long64_t        ES_NTracks[2][2];
   Double_t        ES_M11[2][2];
   Double_t        ES_M12[2][2];
   Double_t        ES_M13[2][2];
   Double_t        ES_M14[2][2];
   Double_t        ES_M15[2][2];
   Double_t        ES_M16[2][2];
   Double_t        ES_M22[2][2];
   Double_t        ES_M23[2][2];
   Double_t        ES_M24[2][2];
   Double_t        ES_M25[2][2];
   Double_t        ES_M26[2][2];
   Double_t        ES_M33[2][2];
   Double_t        ES_M34[2][2];
   Double_t        ES_M35[2][2];
   Double_t        ES_M36[2][2];
   Double_t        ES_M44[2][2];
   Double_t        ES_M45[2][2];
   Double_t        ES_M46[2][2];
   Double_t        ES_M55[2][2];
   Double_t        ES_M56[2][2];
   Double_t        ES_M66[2][2];
   Double_t        ES_P1[2][2];
   Double_t        ES_P2[2][2];
   Double_t        ES_P3[2][2];
   Double_t        ES_P4[2][2];
   Double_t        ES_P5[2][2];
   Double_t        ES_P6[2][2];
   Double_t        ES_M31Err2[2][2];
   Double_t        ES_M32Err2[2][2];
   Double_t        ES_M33Err2[2][2];
   Double_t        ES_P1Err2[2][2];
   Double_t        ES_P2Err2[2][2];
   Double_t        ES_P3Err2[2][2];

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_evtNumber;   //!
   TBranch        *b_Nesrh;   //!
   TBranch        *b__esRecHit_E;   //!
   TBranch        *b__esRecHit_X;   //!
   TBranch        *b__esRecHit_Y;   //!
   TBranch        *b__esRecHit_Z;   //!
   TBranch        *b__esRecHit_Eta;   //!
   TBranch        *b__esRecHit_Phi;   //!
   TBranch        *b__esRecHit_siZ;   //!
   TBranch        *b__esRecHit_siP;   //!
   TBranch        *b__esRecHit_siX;   //!
   TBranch        *b__esRecHit_siY;   //!
   TBranch        *b__esRecHit_Strip;   //!
   TBranch        *b__esRecHit_Noisy;   //!
   TBranch        *b__esRecHit_MatchedTrk_fromOuter;   //!
   TBranch        *b_Ntrack;   //!
   TBranch        *b_TrackPt;   //!
   TBranch        *b_TrackEta;   //!
   TBranch        *b_TrackPhi;   //!
   TBranch        *b_TrackVx;   //!
   TBranch        *b_TrackVy;   //!
   TBranch        *b_TrackVz;   //!
   TBranch        *b_Trackd0;   //!
   TBranch        *b_TrackInnerX;   //!
   TBranch        *b_TrackInnerY;   //!
   TBranch        *b_TrackInnerZ;   //!
   TBranch        *b_TrackOuterZ;   //!
   TBranch        *b_TrackOuterEta;   //!
   TBranch        *b_TrackOuterPhi;   //!
   TBranch        *b_TrackCharge;   //!
   TBranch        *b_TrackNHit;   //!
   TBranch        *b_TrackNChi2;   //!
   TBranch        *b_TrackPtError;   //!
   TBranch        *b_TrackQuality;   //!
   TBranch        *b_PredictionState_iz;   //!
   TBranch        *b_PredictionState_ip;   //!
   TBranch        *b_PredictionState_valid;   //!
   TBranch        *b_PredictionState_X;   //!
   TBranch        *b_PredictionState_Y;   //!
   TBranch        *b_PredictionState_Z;   //!
   TBranch        *b_PredictionState_Px;   //!
   TBranch        *b_PredictionState_Py;   //!
   TBranch        *b_PredictionState_Pz;   //!
   TBranch        *b_PredictionState_Bx;   //!
   TBranch        *b_PredictionState_By;   //!
   TBranch        *b_PredictionState_Bz;   //!
   TBranch        *b_PredictionState_Exx;   //!
   TBranch        *b_PredictionState_Eyx;   //!
   TBranch        *b_PredictionState_Eyy;   //!
   TBranch        *b_PredictionState_delX;   //!
   TBranch        *b_PredictionState_delY;   //!
   TBranch        *b_PredictionState_MatchedRec;   //!
   TBranch        *b_PredictionState_resiX;   //!
   TBranch        *b_PredictionState_resiY;   //!
   TBranch        *b_ES_CHI2;   //!
   TBranch        *b_ES_NTracks;   //!
   TBranch        *b_ES_M11;   //!
   TBranch        *b_ES_M12;   //!
   TBranch        *b_ES_M13;   //!
   TBranch        *b_ES_M14;   //!
   TBranch        *b_ES_M15;   //!
   TBranch        *b_ES_M16;   //!
   TBranch        *b_ES_M22;   //!
   TBranch        *b_ES_M23;   //!
   TBranch        *b_ES_M24;   //!
   TBranch        *b_ES_M25;   //!
   TBranch        *b_ES_M26;   //!
   TBranch        *b_ES_M33;   //!
   TBranch        *b_ES_M34;   //!
   TBranch        *b_ES_M35;   //!
   TBranch        *b_ES_M36;   //!
   TBranch        *b_ES_M44;   //!
   TBranch        *b_ES_M45;   //!
   TBranch        *b_ES_M46;   //!
   TBranch        *b_ES_M55;   //!
   TBranch        *b_ES_M56;   //!
   TBranch        *b_ES_M66;   //!
   TBranch        *b_ES_P1;   //!
   TBranch        *b_ES_P2;   //!
   TBranch        *b_ES_P3;   //!
   TBranch        *b_ES_P4;   //!
   TBranch        *b_ES_P5;   //!
   TBranch        *b_ES_P6;   //!
   TBranch        *b_ES_M31Err2;   //!
   TBranch        *b_ES_M32Err2;   //!
   TBranch        *b_ES_M33Err2;   //!
   TBranch        *b_ES_P1Err2;   //!
   TBranch        *b_ES_P2Err2;   //!
   TBranch        *b_ES_P3Err2;   //!

   ESAlign(TTree *tree=0);
   virtual ~ESAlign();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   void             ReadElements();
   void             ReadElements_wRotation();
};

#endif

#ifdef ESAlign_cxx
ESAlign::ESAlign(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("AlignmentFile.root");
      if (!f) {
         f = new TFile("AlignmentFile.root");
      }
      tree = (TTree*)gDirectory->Get("ESAlign");

   }
   Init(tree);
}

ESAlign::~ESAlign()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ESAlign::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ESAlign::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ESAlign::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("evtNumber", &evtNumber, &b_evtNumber);
   fChain->SetBranchAddress("Nesrh", &Nesrh, &b_Nesrh);
   fChain->SetBranchAddress("_esRecHit_E", _esRecHit_E, &b__esRecHit_E);
   fChain->SetBranchAddress("_esRecHit_X", _esRecHit_X, &b__esRecHit_X);
   fChain->SetBranchAddress("_esRecHit_Y", _esRecHit_Y, &b__esRecHit_Y);
   fChain->SetBranchAddress("_esRecHit_Z", _esRecHit_Z, &b__esRecHit_Z);
   fChain->SetBranchAddress("_esRecHit_Eta", _esRecHit_Eta, &b__esRecHit_Eta);
   fChain->SetBranchAddress("_esRecHit_Phi", _esRecHit_Phi, &b__esRecHit_Phi);
   fChain->SetBranchAddress("_esRecHit_siZ", _esRecHit_siZ, &b__esRecHit_siZ);
   fChain->SetBranchAddress("_esRecHit_siP", _esRecHit_siP, &b__esRecHit_siP);
   fChain->SetBranchAddress("_esRecHit_siX", _esRecHit_siX, &b__esRecHit_siX);
   fChain->SetBranchAddress("_esRecHit_siY", _esRecHit_siY, &b__esRecHit_siY);
   fChain->SetBranchAddress("_esRecHit_Strip", _esRecHit_Strip, &b__esRecHit_Strip);
   fChain->SetBranchAddress("_esRecHit_Noisy", _esRecHit_Noisy, &b__esRecHit_Noisy);
   fChain->SetBranchAddress("_esRecHit_MatchedTrk_fromOuter", _esRecHit_MatchedTrk_fromOuter, &b__esRecHit_MatchedTrk_fromOuter);
   fChain->SetBranchAddress("Ntrack", &Ntrack, &b_Ntrack);
   fChain->SetBranchAddress("TrackPt", TrackPt, &b_TrackPt);
   fChain->SetBranchAddress("TrackEta", TrackEta, &b_TrackEta);
   fChain->SetBranchAddress("TrackPhi", TrackPhi, &b_TrackPhi);
   fChain->SetBranchAddress("TrackVx", TrackVx, &b_TrackVx);
   fChain->SetBranchAddress("TrackVy", TrackVy, &b_TrackVy);
   fChain->SetBranchAddress("TrackVz", TrackVz, &b_TrackVz);
   fChain->SetBranchAddress("Trackd0", Trackd0, &b_Trackd0);
   fChain->SetBranchAddress("TrackInnerX", TrackInnerX, &b_TrackInnerX);
   fChain->SetBranchAddress("TrackInnerY", TrackInnerY, &b_TrackInnerY);
   fChain->SetBranchAddress("TrackInnerZ", TrackInnerZ, &b_TrackInnerZ);
   fChain->SetBranchAddress("TrackOuterZ", TrackOuterZ, &b_TrackOuterZ);
   fChain->SetBranchAddress("TrackOuterEta", TrackOuterEta, &b_TrackOuterEta);
   fChain->SetBranchAddress("TrackOuterPhi", TrackOuterPhi, &b_TrackOuterPhi);
   fChain->SetBranchAddress("TrackCharge", TrackCharge, &b_TrackCharge);
   fChain->SetBranchAddress("TrackNHit", TrackNHit, &b_TrackNHit);
   fChain->SetBranchAddress("TrackNChi2", TrackNChi2, &b_TrackNChi2);
   fChain->SetBranchAddress("TrackPtError", TrackPtError, &b_TrackPtError);
   fChain->SetBranchAddress("TrackQuality", TrackQuality, &b_TrackQuality);
   fChain->SetBranchAddress("PredictionState_iz", PredictionState_iz, &b_PredictionState_iz);
   fChain->SetBranchAddress("PredictionState_ip", PredictionState_ip, &b_PredictionState_ip);
   fChain->SetBranchAddress("PredictionState_valid", PredictionState_valid, &b_PredictionState_valid);
   fChain->SetBranchAddress("PredictionState_X", PredictionState_X, &b_PredictionState_X);
   fChain->SetBranchAddress("PredictionState_Y", PredictionState_Y, &b_PredictionState_Y);
   fChain->SetBranchAddress("PredictionState_Z", PredictionState_Z, &b_PredictionState_Z);
   fChain->SetBranchAddress("PredictionState_Px", PredictionState_Px, &b_PredictionState_Px);
   fChain->SetBranchAddress("PredictionState_Py", PredictionState_Py, &b_PredictionState_Py);
   fChain->SetBranchAddress("PredictionState_Pz", PredictionState_Pz, &b_PredictionState_Pz);
   fChain->SetBranchAddress("PredictionState_Bx", PredictionState_Bx, &b_PredictionState_Bx);
   fChain->SetBranchAddress("PredictionState_By", PredictionState_By, &b_PredictionState_By);
   fChain->SetBranchAddress("PredictionState_Bz", PredictionState_Bz, &b_PredictionState_Bz);
   fChain->SetBranchAddress("PredictionState_Exx", PredictionState_Exx, &b_PredictionState_Exx);
   fChain->SetBranchAddress("PredictionState_Eyx", PredictionState_Eyx, &b_PredictionState_Eyx);
   fChain->SetBranchAddress("PredictionState_Eyy", PredictionState_Eyy, &b_PredictionState_Eyy);
   fChain->SetBranchAddress("PredictionState_delX", PredictionState_delX, &b_PredictionState_delX);
   fChain->SetBranchAddress("PredictionState_delY", PredictionState_delY, &b_PredictionState_delY);
   fChain->SetBranchAddress("PredictionState_MatchedRec", PredictionState_MatchedRec, &b_PredictionState_MatchedRec);
   fChain->SetBranchAddress("PredictionState_resiX", PredictionState_resiX, &b_PredictionState_resiX);
   fChain->SetBranchAddress("PredictionState_resiY", PredictionState_resiY, &b_PredictionState_resiY);
   fChain->SetBranchAddress("ES_CHI2", ES_CHI2, &b_ES_CHI2);
   fChain->SetBranchAddress("ES_NTracks", ES_NTracks, &b_ES_NTracks);
   fChain->SetBranchAddress("ES_M11", ES_M11, &b_ES_M11);
   fChain->SetBranchAddress("ES_M12", ES_M12, &b_ES_M12);
   fChain->SetBranchAddress("ES_M13", ES_M13, &b_ES_M13);
   fChain->SetBranchAddress("ES_M14", ES_M14, &b_ES_M14);
   fChain->SetBranchAddress("ES_M15", ES_M15, &b_ES_M15);
   fChain->SetBranchAddress("ES_M16", ES_M16, &b_ES_M16);
   fChain->SetBranchAddress("ES_M22", ES_M22, &b_ES_M22);
   fChain->SetBranchAddress("ES_M23", ES_M23, &b_ES_M23);
   fChain->SetBranchAddress("ES_M24", ES_M24, &b_ES_M24);
   fChain->SetBranchAddress("ES_M25", ES_M25, &b_ES_M25);
   fChain->SetBranchAddress("ES_M26", ES_M26, &b_ES_M26);
   fChain->SetBranchAddress("ES_M33", ES_M33, &b_ES_M33);
   fChain->SetBranchAddress("ES_M34", ES_M34, &b_ES_M34);
   fChain->SetBranchAddress("ES_M35", ES_M35, &b_ES_M35);
   fChain->SetBranchAddress("ES_M36", ES_M36, &b_ES_M36);
   fChain->SetBranchAddress("ES_M44", ES_M44, &b_ES_M44);
   fChain->SetBranchAddress("ES_M45", ES_M45, &b_ES_M45);
   fChain->SetBranchAddress("ES_M46", ES_M46, &b_ES_M46);
   fChain->SetBranchAddress("ES_M55", ES_M55, &b_ES_M55);
   fChain->SetBranchAddress("ES_M56", ES_M56, &b_ES_M56);
   fChain->SetBranchAddress("ES_M66", ES_M66, &b_ES_M66);
   fChain->SetBranchAddress("ES_P1", ES_P1, &b_ES_P1);
   fChain->SetBranchAddress("ES_P2", ES_P2, &b_ES_P2);
   fChain->SetBranchAddress("ES_P3", ES_P3, &b_ES_P3);
   fChain->SetBranchAddress("ES_P4", ES_P4, &b_ES_P4);
   fChain->SetBranchAddress("ES_P5", ES_P5, &b_ES_P5);
   fChain->SetBranchAddress("ES_P6", ES_P6, &b_ES_P6);
   fChain->SetBranchAddress("ES_M31Err2", ES_M31Err2, &b_ES_M31Err2);
   fChain->SetBranchAddress("ES_M32Err2", ES_M32Err2, &b_ES_M32Err2);
   fChain->SetBranchAddress("ES_M33Err2", ES_M33Err2, &b_ES_M33Err2);
   fChain->SetBranchAddress("ES_P1Err2", ES_P1Err2, &b_ES_P1Err2);
   fChain->SetBranchAddress("ES_P2Err2", ES_P2Err2, &b_ES_P2Err2);
   fChain->SetBranchAddress("ES_P3Err2", ES_P3Err2, &b_ES_P3Err2);
   Notify();
}

Bool_t ESAlign::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ESAlign::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ESAlign::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ESAlign_cxx
