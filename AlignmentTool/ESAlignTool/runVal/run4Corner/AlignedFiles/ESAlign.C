#define ESAlign_cxx

#include <iostream>
#include <fstream>

#include "ESAlign.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include "tdrstyle.C"


void ESAlign::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L ESAlign.C
//      Root > ESAlign t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

void ESAlign::ReadElements()
{

 fChain->GetEntry(fChain->GetEntries()-1);
 ofstream oText;
 oText.open("outputText.txt",ios::out|ios::app);

 char buf[20]; Long64_t valI=0; double valD=0.;
 for(int i=0;i<2;i++)
 {for(int j=0;j<2;j++)
  {
   valD=ES_M11[i][j]; sprintf(buf,"ES_M11[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M12[i][j]; sprintf(buf,"ES_M12[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M13[i][j]; sprintf(buf,"ES_M13[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M22[i][j]; sprintf(buf,"ES_M22[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M23[i][j]; sprintf(buf,"ES_M23[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M33[i][j]; sprintf(buf,"ES_M33[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_CHI2[i][j]; sprintf(buf,"ES_CHI2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M31Err2[i][j]; sprintf(buf,"ES_M31Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M32Err2[i][j]; sprintf(buf,"ES_M32Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M33Err2[i][j]; sprintf(buf,"ES_M33Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P1Err2[i][j]; sprintf(buf,"ES_P1Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P2Err2[i][j]; sprintf(buf,"ES_P2Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P3Err2[i][j]; sprintf(buf,"ES_P3Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valI=ES_NTracks[i][j]; sprintf(buf,"ES_NTracks[%i][%i]+= (",i,j);
   oText<<buf<<valI<<");";
   valD=ES_P1[i][j]; sprintf(buf,"ES_P1[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P2[i][j]; sprintf(buf,"ES_P2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P3[i][j]; sprintf(buf,"ES_P3[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P1Err2[i][j]; sprintf(buf,"ES_P1Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P2Err2[i][j]; sprintf(buf,"ES_P2Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P3Err2[i][j]; sprintf(buf,"ES_P3Err2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";

 }}

oText.close();

}

void ESAlign::ReadElements_wRotation()
{

 fChain->GetEntry(fChain->GetEntries()-1);
 ofstream oText;
 oText.open("outputText.txt",ios::out|ios::app);

 char buf[20]; Long64_t valI=0; double valD=0.;
 for(int i=0;i<2;i++)
 {for(int j=0;j<2;j++)
  {
   valD=ES_M11[i][j]; sprintf(buf,"ES_M11[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M12[i][j]; sprintf(buf,"ES_M12[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M13[i][j]; sprintf(buf,"ES_M13[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M14[i][j]; sprintf(buf,"ES_M14[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M15[i][j]; sprintf(buf,"ES_M15[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M16[i][j]; sprintf(buf,"ES_M16[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M22[i][j]; sprintf(buf,"ES_M22[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M23[i][j]; sprintf(buf,"ES_M23[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M24[i][j]; sprintf(buf,"ES_M24[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M25[i][j]; sprintf(buf,"ES_M25[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M26[i][j]; sprintf(buf,"ES_M26[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M33[i][j]; sprintf(buf,"ES_M33[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M34[i][j]; sprintf(buf,"ES_M34[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M35[i][j]; sprintf(buf,"ES_M35[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M36[i][j]; sprintf(buf,"ES_M36[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M44[i][j]; sprintf(buf,"ES_M44[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M45[i][j]; sprintf(buf,"ES_M45[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M46[i][j]; sprintf(buf,"ES_M46[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M55[i][j]; sprintf(buf,"ES_M55[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M56[i][j]; sprintf(buf,"ES_M56[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_M66[i][j]; sprintf(buf,"ES_M66[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";

   valI=ES_NTracks[i][j]; sprintf(buf,"ES_NTracks[%i][%i]+= (",i,j);
   oText<<buf<<valI<<");";

   valD=ES_P1[i][j]; sprintf(buf,"ES_P1[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P2[i][j]; sprintf(buf,"ES_P2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P3[i][j]; sprintf(buf,"ES_P3[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P4[i][j]; sprintf(buf,"ES_P4[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P5[i][j]; sprintf(buf,"ES_P5[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";
   valD=ES_P6[i][j]; sprintf(buf,"ES_P6[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");\n";

 }}

oText.close();

}

void ESAlign::Draw_vsXY()
{
 //gROOT->ProcessLine(".L ./format/tdrstyle.C");
 setTDRStyle();
 gStyle->SetPalette(1);
   
// TChain *EcalRecHit = new TChain("EcalRecHit");
// EcalRecHit->Add("AlignmentFile_*.root");
 double X=0.0;double Y=0.0;

 TH2D *h_ESpF_resiX_vs_recoX=new TH2D("h_ESpF_resiX_vs_recoX","",100,-110,110,100,-3,3);
 TH2D *h_ESpR_resiY_vs_recoY=new TH2D("h_ESpR_resiY_vs_recoY","",100,-110,110,100,-3,3);
 TH2D *h_ESmF_resiX_vs_recoX=new TH2D("h_ESmF_resiX_vs_recoX","",100,-110,110,100,-3,3);
 TH2D *h_ESmR_resiY_vs_recoY=new TH2D("h_ESmR_resiY_vs_recoY","",100,-110,110,100,-3,3);
 h_ESpF_resiX_vs_recoX->SetXTitle("RecHit X(cm)"); h_ESpF_resiX_vs_recoX->SetYTitle("ESpF Residual X(cm)");
 h_ESpR_resiY_vs_recoY->SetXTitle("RecHit Y(cm)"); h_ESpR_resiY_vs_recoY->SetYTitle("ESpR Residual Y(cm)");
 h_ESmF_resiX_vs_recoX->SetXTitle("RecHit X(cm)"); h_ESmF_resiX_vs_recoX->SetYTitle("ESmF Residual X(cm)");
 h_ESmR_resiY_vs_recoY->SetXTitle("RecHit Y(cm)"); h_ESmR_resiY_vs_recoY->SetYTitle("ESmR Residual Y(cm)");

 TH2D *h_ESpF_resiY_vs_recoY=new TH2D("h_ESpF_resiY_vs_recoY","",100,-110,110,100,-3,3);
 TH2D *h_ESpR_resiX_vs_recoX=new TH2D("h_ESpR_resiX_vs_recoX","",100,-110,110,100,-3,3);
 TH2D *h_ESmF_resiY_vs_recoY=new TH2D("h_ESmF_resiY_vs_recoY","",100,-110,110,100,-3,3);
 TH2D *h_ESmR_resiX_vs_recoX=new TH2D("h_ESmR_resiX_vs_recoX","",100,-110,110,100,-3,3);
 h_ESpF_resiY_vs_recoY->SetXTitle("RecHit Y(cm)"); h_ESpF_resiY_vs_recoY->SetYTitle("ESpF Residual Y(cm)");
 h_ESpR_resiX_vs_recoX->SetXTitle("RecHit X(cm)"); h_ESpR_resiX_vs_recoX->SetYTitle("ESpR Residual X(cm)");
 h_ESmF_resiY_vs_recoY->SetXTitle("RecHit Y(cm)"); h_ESmF_resiY_vs_recoY->SetYTitle("ESmF Residual Y(cm)");
 h_ESmR_resiX_vs_recoX->SetXTitle("RecHit X(cm)"); h_ESmR_resiX_vs_recoX->SetYTitle("ESmR Residual X(cm)");

 TH2D *h_ESpF_resiX_vs_recoY=new TH2D("h_ESpF_resiX_vs_recoY","",100,-110,110,100,-3,3);
 TH2D *h_ESpR_resiY_vs_recoX=new TH2D("h_ESpR_resiY_vs_recoX","",100,-110,110,100,-3,3);
 TH2D *h_ESmF_resiX_vs_recoY=new TH2D("h_ESmF_resiX_vs_recoY","",100,-110,110,100,-3,3);
 TH2D *h_ESmR_resiY_vs_recoX=new TH2D("h_ESmR_resiY_vs_recoX","",100,-110,110,100,-3,3);
 h_ESpF_resiX_vs_recoY->SetXTitle("RecHit Y(cm)"); h_ESpF_resiX_vs_recoY->SetYTitle("ESpF Residual X(cm)");
 h_ESpR_resiY_vs_recoX->SetXTitle("RecHit X(cm)"); h_ESpR_resiY_vs_recoX->SetYTitle("ESpR Residual Y(cm)");
 h_ESmF_resiX_vs_recoY->SetXTitle("RecHit Y(cm)"); h_ESmF_resiX_vs_recoY->SetYTitle("ESmF Residual X(cm)");
 h_ESmR_resiY_vs_recoX->SetXTitle("RecHit X(cm)"); h_ESmR_resiY_vs_recoX->SetYTitle("ESmR Residual Y(cm)");



 TH2D *h_ESpF_resiX_vs_recoXY=new TH2D("h_ESpF_resiX_vs_recoXY","ESpF residual on reco XY-plane",36,-114.5005,114.5005,36,-114.5005,114.5005);
 TH2I *h_ESpF_N_vs_recoXY=new TH2I("h_ESpF_N_vs_recoXY","N on reco XY-plane",36,-114.5005,114.5005,36,-114.5005,114.5005);
 TH2D *h_ESmF_resiX_vs_recoXY=new TH2D("h_ESmF_resiX_vs_recoXY","ESmF residual on reco XY-plane",36,-114.5005,114.5005,36,-114.5005,114.5005);
 TH2I *h_ESmF_N_vs_recoXY=new TH2I("h_ESmF_N_vs_recoXY","N on reco XY-plane",36,-114.5005,114.5005,36,-114.5005,114.5005);
 TH2D *h_ESpR_resiY_vs_recoXY=new TH2D("h_ESpR_resiY_vs_recoXY","ESpR residual on reco XY-plane",36,-114.5005,114.5005,36,-114.5005,114.5005);
 TH2I *h_ESpR_N_vs_recoXY=new TH2I("h_ESpR_N_vs_recoXY","N on reco XY-plane",36,-114.5005,114.5005,36,-114.5005,114.5005);
 TH2D *h_ESmR_resiY_vs_recoXY=new TH2D("h_ESmR_resiY_vs_recoXY","ESmR residual on reco XY-plane",36,-114.5005,114.5005,36,-114.5005,114.5005);
 TH2I *h_ESmR_N_vs_recoXY=new TH2I("h_ESmR_N_vs_recoXY","N on reco XY-plane",36,-114.5005,114.5005,36,-114.5005,114.5005);


 TCanvas *v1=new TCanvas("v1","",1600,800);
 v1->Divide(2,2);
 
 int NEntries=fChain->GetEntries();
 for(int i=0;i< NEntries ;i++)
 {
  X=0.0; Y=0.0;
  //LoadTree(i);
  fChain->GetEntry(i);
  if( (i<10) || ( i<100 && i%10==0) || ( i<1000 && i%100==0) || i%1000==0 )
   std::cout<< i << "/" << NEntries << "\n";
  

  for(int a=0;a<2;a++)
  {
   for(int b=0;b<2;b++)
   {
    for(int iTrk=0;iTrk<Ntrack;iTrk++)
    {
     if(a==0&&PredictionState_iz[iTrk][a][b]!=-1) continue;
     if(a==1&&PredictionState_iz[iTrk][a][b]!=1) continue;
     if(b==0&&PredictionState_ip[iTrk][a][b]!=1) continue;
     if(b==1&&PredictionState_ip[iTrk][a][b]!=2) continue;
     if(PredictionState_MatchedRec[iTrk][a][b]>-1)
     {
      if(    _esRecHit_siZ[PredictionState_MatchedRec[iTrk][a][b]]
           != PredictionState_iz[iTrk][a][b]
         ||  _esRecHit_siP[PredictionState_MatchedRec[iTrk][a][b]]
           !=PredictionState_ip[iTrk][a][b]
        )
      {
       std::cout<<"Error Found : MatchedRec !!\n";
       std::cout<<"entry="<<i<<"\n";
       std::cout<<"siZ="<<_esRecHit_siZ[PredictionState_MatchedRec[iTrk][a][b]]<<"\n";
       std::cout<<"siP="<<_esRecHit_siP[PredictionState_MatchedRec[iTrk][a][b]]<<"\n";
       std::cout<<"index="<<PredictionState_MatchedRec[iTrk][a][b]<<"\n";
       
       break;
      }

      //X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
      //  +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
      //  +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
      //Y= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R21[a][b]
      //  +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R22[a][b]
      //  +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R23[a][b];


      //X=_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b];
      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
      Y=PredictionState_resiX[iTrk][a][b];
      
      if(a==0&&b==0) h_ESmF_resiX_vs_recoX->Fill(X,Y);
      if(a==0&&b==1) h_ESmR_resiX_vs_recoX->Fill(X,Y);
      if(a==1&&b==0) h_ESpF_resiX_vs_recoX->Fill(X,Y);
      if(a==1&&b==1) h_ESpR_resiX_vs_recoX->Fill(X,Y);
  
      //X=_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b];
      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
      Y=PredictionState_resiY[iTrk][a][b];
      if(a==0&&b==0) h_ESmF_resiY_vs_recoY->Fill(X,Y);
      if(a==0&&b==1) h_ESmR_resiY_vs_recoY->Fill(X,Y);
      if(a==1&&b==0) h_ESpF_resiY_vs_recoY->Fill(X,Y);
      if(a==1&&b==1) h_ESpR_resiY_vs_recoY->Fill(X,Y);
  
      //X=_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b];
      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R21[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R22[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R23[a][b];
      Y=PredictionState_resiX[iTrk][a][b];
      if(a==0&&b==0) h_ESmF_resiX_vs_recoY->Fill(X,Y);
      if(a==1&&b==0) h_ESpF_resiX_vs_recoY->Fill(X,Y);

      //X=_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b];
      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
      Y=PredictionState_resiY[iTrk][a][b];
      if(a==0&&b==1) h_ESmR_resiY_vs_recoX->Fill(X,Y);
      if(a==1&&b==1) h_ESpR_resiY_vs_recoX->Fill(X,Y);

      //X=_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b];
      //Y=_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b];
      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
      Y= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R21[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R22[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R23[a][b];
      if(a==0&&b==0)
      {
       h_ESmF_resiX_vs_recoXY->Fill(X,Y,PredictionState_resiX[iTrk][a][b]);
       h_ESmF_N_vs_recoXY->Fill(X,Y);
      }
      if(a==0&&b==1)
      {
       h_ESmR_resiY_vs_recoXY->Fill(X,Y,PredictionState_resiY[iTrk][a][b]);
       h_ESmR_N_vs_recoXY->Fill(X,Y);
      }
      if(a==1&&b==0)
      {
       h_ESpF_resiX_vs_recoXY->Fill(X,Y,PredictionState_resiX[iTrk][a][b]);
       h_ESpF_N_vs_recoXY->Fill(X,Y);
      }
      if(a==1&&b==1)
      {
       h_ESpR_resiY_vs_recoXY->Fill(X,Y,PredictionState_resiY[iTrk][a][b]);
       h_ESpR_N_vs_recoXY->Fill(X,Y);
      }

     }//end if MatchedRec>-1
    }//end all trks  
  }}

 }//loop all events

 v1->cd(1)->SetRightMargin(0.3); h_ESpF_resiX_vs_recoX->Draw("colz");
 v1->cd(2)->SetRightMargin(0.3); h_ESpR_resiY_vs_recoY->Draw("colz");
 v1->cd(3)->SetRightMargin(0.3); h_ESmF_resiX_vs_recoX->Draw("colz");
 v1->cd(4)->SetRightMargin(0.3); h_ESmR_resiY_vs_recoY->Draw("colz");
 v1->SaveAs("Residual_vs_Strip.png");

 v1->cd(1)->SetRightMargin(0.3); h_ESpF_resiY_vs_recoY->Draw("colz");
 v1->cd(2)->SetRightMargin(0.3); h_ESpR_resiX_vs_recoX->Draw("colz");
 v1->cd(3)->SetRightMargin(0.3); h_ESmF_resiY_vs_recoY->Draw("colz");
 v1->cd(4)->SetRightMargin(0.3); h_ESmR_resiX_vs_recoX->Draw("colz");
 v1->SaveAs("Residual_vs_Sensor.png");

 v1->cd(1)->SetRightMargin(0.3); h_ESpF_resiX_vs_recoY->Draw("colz");
 v1->cd(2)->SetRightMargin(0.3); h_ESpR_resiY_vs_recoX->Draw("colz");
 v1->cd(3)->SetRightMargin(0.3); h_ESmF_resiX_vs_recoY->Draw("colz");
 v1->cd(4)->SetRightMargin(0.3); h_ESmR_resiY_vs_recoX->Draw("colz");
 v1->SaveAs("Residual_vs_reco90d.png");

 int N=0; double val=0.;
 int Nthre=50;
 for(int ix=1;ix<=36;ix++)
 {
  for(int iy=1;iy<=36;iy++)
  {
   N=0; val=0.;
   N=h_ESpF_N_vs_recoXY->GetBinContent(ix,iy);
   val=h_ESpF_resiX_vs_recoXY->GetBinContent(ix,iy);
   if(N>Nthre) h_ESpF_resiX_vs_recoXY->SetBinContent(ix,iy,10.+val/(double)N);
   if(N<=Nthre) h_ESpF_resiX_vs_recoXY->SetBinContent(ix,iy,0);

   N=0; val=0.;
   N=h_ESpR_N_vs_recoXY->GetBinContent(ix,iy);
   val=h_ESpR_resiY_vs_recoXY->GetBinContent(ix,iy);
   if(N>Nthre) h_ESpR_resiY_vs_recoXY->SetBinContent(ix,iy,10.+val/(double)N);
   if(N<=Nthre) h_ESpR_resiY_vs_recoXY->SetBinContent(ix,iy,0);

   N=0; val=0.;
   N=h_ESmF_N_vs_recoXY->GetBinContent(ix,iy);
   val=h_ESmF_resiX_vs_recoXY->GetBinContent(ix,iy);
   if(N>Nthre) h_ESmF_resiX_vs_recoXY->SetBinContent(ix,iy,10.+val/(double)N);
   if(N<=Nthre) h_ESmF_resiX_vs_recoXY->SetBinContent(ix,iy,0);

   N=0; val=0.;
   N=h_ESmR_N_vs_recoXY->GetBinContent(ix,iy);
   val=h_ESmR_resiY_vs_recoXY->GetBinContent(ix,iy);
   if(N>Nthre) h_ESmR_resiY_vs_recoXY->SetBinContent(ix,iy,10.+val/(double)N);
   if(N<=Nthre) h_ESmR_resiY_vs_recoXY->SetBinContent(ix,iy,0);
 }}
 h_ESpF_resiX_vs_recoXY->SetMaximum(10.+0.3); h_ESpF_resiX_vs_recoXY->SetMinimum(10.-0.3);
 h_ESpR_resiY_vs_recoXY->SetMaximum(10.+0.3); h_ESpR_resiY_vs_recoXY->SetMinimum(10.-0.3);
 h_ESmF_resiX_vs_recoXY->SetMaximum(10.+0.3); h_ESmF_resiX_vs_recoXY->SetMinimum(10.-0.3);
 h_ESmR_resiY_vs_recoXY->SetMaximum(10.+0.3); h_ESmR_resiY_vs_recoXY->SetMinimum(10.-0.3);

 TH2D *h2d_empty=new TH2D("h2d_empty","",36,-114.5005,114.5005,36,-114.5005,114.5005);
 h2d_empty->SetMaximum(0.3); h2d_empty->SetMinimum(-0.3);
 h2d_empty->Fill(-110,0);

 v1->cd(1)->SetRightMargin(0.3);
  h2d_empty->Draw("colz"); gPad->Update();
  h_ESpF_resiX_vs_recoXY->Draw("col");
  h2d_empty->Draw("axiszsame");
 v1->cd(2)->SetRightMargin(0.3);
  h2d_empty->Draw("colz"); gPad->Update();
  h_ESpR_resiY_vs_recoXY->Draw("col");
  h2d_empty->Draw("axiszsame");
 v1->cd(3)->SetRightMargin(0.3);
  h2d_empty->Draw("colz"); gPad->Update();
  h_ESmF_resiX_vs_recoXY->Draw("col");
  h2d_empty->Draw("axiszsame");
 v1->cd(4)->SetRightMargin(0.3);
  h2d_empty->Draw("colz"); gPad->Update();
  h_ESmR_resiY_vs_recoXY->Draw("col");
  h2d_empty->Draw("axiszsame");
 v1->SaveAs("residual_onPlane.png");


 v1->cd(1)->SetRightMargin(0.3); h_ESpF_N_vs_recoXY->Draw("colz");
 v1->cd(2)->SetRightMargin(0.3); h_ESpR_N_vs_recoXY->Draw("colz");
 v1->cd(3)->SetRightMargin(0.3); h_ESmF_N_vs_recoXY->Draw("colz");
 v1->cd(4)->SetRightMargin(0.3); h_ESmR_N_vs_recoXY->Draw("colz");
 v1->SaveAs("N_onPlane.png");


}

void ESAlign::Draw_4quad()
{
 //gROOT->ProcessLine(".L ./format/tdrstyle.C");
 setTDRStyle();
 gStyle->SetPalette(1);
   
// TChain *EcalRecHit = new TChain("EcalRecHit");
// EcalRecHit->Add("AlignmentFile_*.root");
 double X=0.0;double Y=0.0;

 TH1D *h_ESpF_resiX_quad1=new TH1D("h_ESpF_resiX_quad1","quad1",60,-3,3);
 TH1D *h_ESpF_resiX_quad2=new TH1D("h_ESpF_resiX_quad2","quad2",60,-3,3);
 TH1D *h_ESpF_resiX_quad3=new TH1D("h_ESpF_resiX_quad3","quad3",60,-3,3);
 TH1D *h_ESpF_resiX_quad4=new TH1D("h_ESpF_resiX_quad4","quad4",60,-3,3);
 TH1D *h_ESpR_resiY_quad1=new TH1D("h_ESpR_resiY_quad1","quad1",60,-3,3);
 TH1D *h_ESpR_resiY_quad2=new TH1D("h_ESpR_resiY_quad2","quad2",60,-3,3);
 TH1D *h_ESpR_resiY_quad3=new TH1D("h_ESpR_resiY_quad3","quad3",60,-3,3);
 TH1D *h_ESpR_resiY_quad4=new TH1D("h_ESpR_resiY_quad4","quad4",60,-3,3);
 TH1D *h_ESmF_resiX_quad1=new TH1D("h_ESmF_resiX_quad1","quad1",60,-3,3);
 TH1D *h_ESmF_resiX_quad2=new TH1D("h_ESmF_resiX_quad2","quad2",60,-3,3);
 TH1D *h_ESmF_resiX_quad3=new TH1D("h_ESmF_resiX_quad3","quad3",60,-3,3);
 TH1D *h_ESmF_resiX_quad4=new TH1D("h_ESmF_resiX_quad4","quad4",60,-3,3);
 TH1D *h_ESmR_resiY_quad1=new TH1D("h_ESmR_resiY_quad1","quad1",60,-3,3);
 TH1D *h_ESmR_resiY_quad2=new TH1D("h_ESmR_resiY_quad2","quad2",60,-3,3);
 TH1D *h_ESmR_resiY_quad3=new TH1D("h_ESmR_resiY_quad3","quad3",60,-3,3);
 TH1D *h_ESmR_resiY_quad4=new TH1D("h_ESmR_resiY_quad4","quad4",60,-3,3);

 int NEntries=fChain->GetEntries();
 for(int i=0;i< NEntries ;i++)
 {
  fChain->GetEntry(i);

  for(int a=0;a<2;a++)
  {
   for(int b=0;b<2;b++)
   {
    for(int iTrk=0;iTrk<Ntrack;iTrk++)
    {
     if(a==0&&PredictionState_iz[iTrk][a][b]!=-1) continue;
     if(a==1&&PredictionState_iz[iTrk][a][b]!=1) continue;
     if(b==0&&PredictionState_ip[iTrk][a][b]!=1) continue;
     if(b==1&&PredictionState_ip[iTrk][a][b]!=2) continue;
     if(PredictionState_MatchedRec[iTrk][a][b]>-1)
     {
      //X=_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b];
      //Y=_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b];
      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
      Y= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R21[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R22[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R23[a][b];
      if(a==0&&b==0)
      {
       if(X>6.39&&Y>0.)
        h_ESmF_resiX_quad1->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X<6.39&&Y>0.)
        h_ESmF_resiX_quad2->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X<-6.39&&Y<0.)
        h_ESmF_resiX_quad3->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X>-6.39&&Y<0.)
        h_ESmF_resiX_quad4->Fill(PredictionState_resiX[iTrk][a][b]);
      }//end if ES-F
      if(a==0&&b==1)
      {
       if(X>0.&&Y>0.)
        h_ESmR_resiY_quad1->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X<0.&&Y>0.)
        h_ESmR_resiY_quad2->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X<0.&&Y<0.)
        h_ESmR_resiY_quad3->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X>0.&&Y<0.)
        h_ESmR_resiY_quad4->Fill(PredictionState_resiY[iTrk][a][b]);
      }//end if ES-R
      if(a==1&&b==0)
      {
       if(X>-6.39&&Y>0.)
        h_ESpF_resiX_quad1->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X<-6.39&&Y>0.)
        h_ESpF_resiX_quad2->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X<6.39&&Y<0.)
        h_ESpF_resiX_quad3->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X>6.39&&Y<0.)
        h_ESpF_resiX_quad4->Fill(PredictionState_resiX[iTrk][a][b]);
      }//end if ES+F
      if(a==1&&b==1)
      {
       if(X>0.&&Y>0.)
        h_ESpR_resiY_quad1->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X<0.&&Y>0.)
        h_ESpR_resiY_quad2->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X<0.&&Y<0.)
        h_ESpR_resiY_quad3->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X>0.&&Y<0.)
        h_ESpR_resiY_quad4->Fill(PredictionState_resiY[iTrk][a][b]);
      }//end if ES+R
     }//end if matched
    }//end for-loop iTrk
   }//end for-loop ip
  }//end for-loop iz

 }//end loop all entries


 double mean=0.; double sig=0.;

 TFile *f1,*f2,*f3,*f4;

 if(iterN==1)
 {
  f1=new TFile("ESmF_4quad_iter1.root","RECREATE");
  f2=new TFile("ESmR_4quad_iter1.root","RECREATE");
  f3=new TFile("ESpF_4quad_iter1.root","RECREATE");
  f4=new TFile("ESpR_4quad_iter1.root","RECREATE");
 }
 if(iterN==11)
 {
  f1=new TFile("ESmF_4quad_iter11.root","RECREATE");
  f2=new TFile("ESmR_4quad_iter11.root","RECREATE");
  f3=new TFile("ESpF_4quad_iter11.root","RECREATE");
  f4=new TFile("ESpR_4quad_iter11.root","RECREATE");
 }
 if(iterN!=1&&iterN!=11)
 {
  f1=new TFile("ESmF_4quad.root","RECREATE");
  f2=new TFile("ESmR_4quad.root","RECREATE");
  f3=new TFile("ESpF_4quad.root","RECREATE");
  f4=new TFile("ESpR_4quad.root","RECREATE");
 }

 f1->cd();
 h_ESmF_resiX_quad1->Write();
 h_ESmF_resiX_quad2->Write();
 h_ESmF_resiX_quad3->Write();
 h_ESmF_resiX_quad4->Write();
 f2->cd();
 h_ESmR_resiY_quad1->Write();
 h_ESmR_resiY_quad2->Write();
 h_ESmR_resiY_quad3->Write();
 h_ESmR_resiY_quad4->Write();
 f3->cd();
 h_ESpF_resiX_quad1->Write();
 h_ESpF_resiX_quad2->Write();
 h_ESpF_resiX_quad3->Write();
 h_ESpF_resiX_quad4->Write();
 f4->cd();
 h_ESpR_resiY_quad1->Write();
 h_ESpR_resiY_quad2->Write();
 h_ESpR_resiY_quad3->Write();
 h_ESpR_resiY_quad4->Write();

}

void ESAlign::Draw_6Corner()
{
 //gROOT->ProcessLine(".L ./format/tdrstyle.C");
 setTDRStyle();
 gStyle->SetPalette(1);
   
// TChain *EcalRecHit = new TChain("EcalRecHit");
// EcalRecHit->Add("AlignmentFile_*.root");
 double X=0.0;double Y=0.0;

 TH1D *h_ESpF_resiX_Right=new TH1D("h_ESpF_resiX_Right","Right",60,-3,3);
 TH1D *h_ESpF_resiX_Left=new TH1D("h_ESpF_resiX_Left","Left",60,-3,3);
 TH1D *h_ESpF_resiX_TopL=new TH1D("h_ESpF_resiX_TopL","TopL",60,-3,3);
 TH1D *h_ESpF_resiX_BottomL=new TH1D("h_ESpF_resiX_BottomL","BottomL",60,-3,3);
 TH1D *h_ESpR_resiY_Right=new TH1D("h_ESpR_resiY_Right","Right",60,-3,3);
 TH1D *h_ESpR_resiY_Left=new TH1D("h_ESpR_resiY_Left","Left",60,-3,3);
 TH1D *h_ESpR_resiY_TopL=new TH1D("h_ESpR_resiY_TopL","TopL",60,-3,3);
 TH1D *h_ESpR_resiY_BottomL=new TH1D("h_ESpR_resiY_BottomL","BottomL",60,-3,3);
 TH1D *h_ESmF_resiX_Right=new TH1D("h_ESmF_resiX_Right","Right",60,-3,3);
 TH1D *h_ESmF_resiX_Left=new TH1D("h_ESmF_resiX_Left","Left",60,-3,3);
 TH1D *h_ESmF_resiX_TopL=new TH1D("h_ESmF_resiX_TopL","TopL",60,-3,3);
 TH1D *h_ESmF_resiX_BottomL=new TH1D("h_ESmF_resiX_BottomL","BottomL",60,-3,3);
 TH1D *h_ESmR_resiY_Right=new TH1D("h_ESmR_resiY_Right","Right",60,-3,3);
 TH1D *h_ESmR_resiY_Left=new TH1D("h_ESmR_resiY_Left","Left",60,-3,3);
 TH1D *h_ESmR_resiY_TopL=new TH1D("h_ESmR_resiY_TopL","TopL",60,-3,3);
 TH1D *h_ESmR_resiY_BottomL=new TH1D("h_ESmR_resiY_BottomL","BottomL",60,-3,3);
 TH1D *h_ESpF_resiX_TopR=new TH1D("h_ESpF_resiX_TopR","TopR",60,-3,3);
 TH1D *h_ESpF_resiX_BottomR=new TH1D("h_ESpF_resiX_BottomR","BottomR",60,-3,3);
 TH1D *h_ESpR_resiY_TopR=new TH1D("h_ESpR_resiY_TopR","TopR",60,-3,3);
 TH1D *h_ESpR_resiY_BottomR=new TH1D("h_ESpR_resiY_BottomR","BottomR",60,-3,3);
 TH1D *h_ESmF_resiX_TopR=new TH1D("h_ESmF_resiX_TopR","TopR",60,-3,3);
 TH1D *h_ESmF_resiX_BottomR=new TH1D("h_ESmF_resiX_BottomR","BottomR",60,-3,3);
 TH1D *h_ESmR_resiY_TopR=new TH1D("h_ESmR_resiY_TopR","TopR",60,-3,3);
 TH1D *h_ESmR_resiY_BottomR=new TH1D("h_ESmR_resiY_BottomR","BottomR",60,-3,3);

 int NEntries=fChain->GetEntries();
 for(int i=0;i< NEntries ;i++)
 {
  fChain->GetEntry(i);

  for(int a=0;a<2;a++)
  {
   for(int b=0;b<2;b++)
   {
    for(int iTrk=0;iTrk<Ntrack;iTrk++)
    {
     if(a==0&&PredictionState_iz[iTrk][a][b]!=-1) continue;
     if(a==1&&PredictionState_iz[iTrk][a][b]!=1) continue;
     if(b==0&&PredictionState_ip[iTrk][a][b]!=1) continue;
     if(b==1&&PredictionState_ip[iTrk][a][b]!=2) continue;
     if(PredictionState_MatchedRec[iTrk][a][b]>-1)
     {
      //X=_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b];
      //Y=_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b];
      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
      Y= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R21[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R22[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R23[a][b];
      if(a==0&&b==0)
      {
       if(X>85.&&fabs(Y)<20.)
        h_ESmF_resiX_Right->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X<-85.&&fabs(Y)<20.)
        h_ESmF_resiX_Left->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.&&X>6.39)
        h_ESmF_resiX_TopR->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.&&X>-6.39)
        h_ESmF_resiX_BottomR->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.&&X<6.39)
        h_ESmF_resiX_TopL->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.&&X<-6.39)
        h_ESmF_resiX_BottomL->Fill(PredictionState_resiX[iTrk][a][b]);
      }//end if ES-F
      if(a==0&&b==1)
      {
       if(X>85.&&fabs(Y)<20.)
        h_ESmR_resiY_Right->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X<-85.&&fabs(Y)<20.)
        h_ESmR_resiY_Left->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.&&X>0.)
        h_ESmR_resiY_TopR->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.&&X>0.)
        h_ESmR_resiY_BottomR->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.&&X<0.)
        h_ESmR_resiY_TopL->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.&&X<0.)
        h_ESmR_resiY_BottomL->Fill(PredictionState_resiY[iTrk][a][b]);
      }//end if ES-R
      if(a==1&&b==0)
      {
       if(X>85.&&fabs(Y)<20.)
        h_ESpF_resiX_Right->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X<-85.&&fabs(Y)<20.)
        h_ESpF_resiX_Left->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.&&X>-6.39)
        h_ESpF_resiX_TopR->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.&&X>6.39)
        h_ESpF_resiX_BottomR->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.&&X<-6.39)
        h_ESpF_resiX_TopL->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.&&X<6.39)
        h_ESpF_resiX_BottomL->Fill(PredictionState_resiX[iTrk][a][b]);
      }//end if ES+F
      if(a==1&&b==1)
      {
       if(X>85.&&fabs(Y)<20.)
        h_ESpR_resiY_Right->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X<-85.&&fabs(Y)<20.)
        h_ESpR_resiY_Left->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.&&X>0.)
        h_ESpR_resiY_TopR->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.&&X>0.)
        h_ESpR_resiY_BottomR->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.&&X<0.)
        h_ESpR_resiY_TopL->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.&&X<0.)
        h_ESpR_resiY_BottomL->Fill(PredictionState_resiY[iTrk][a][b]);
      }//end if ES+R
     }//end if matched
    }//end for-loop iTrk
   }//end for-loop ip
  }//end for-loop iz

 }//end loop all entries


 double mean=0.; double sig=0.;

 TFile *f1,*f2,*f3,*f4;

 if(iterN==1)
 {
  f1=new TFile("ESmF_6Corner_iter1.root","RECREATE");
  f2=new TFile("ESmR_6Corner_iter1.root","RECREATE");
  f3=new TFile("ESpF_6Corner_iter1.root","RECREATE");
  f4=new TFile("ESpR_6Corner_iter1.root","RECREATE");
 }
 if(iterN==11)
 {
  f1=new TFile("ESmF_6Corner_iter11.root","RECREATE");
  f2=new TFile("ESmR_6Corner_iter11.root","RECREATE");
  f3=new TFile("ESpF_6Corner_iter11.root","RECREATE");
  f4=new TFile("ESpR_6Corner_iter11.root","RECREATE");
 }
 if(iterN!=1&&iterN!=11)
 {
  f1=new TFile("ESmF_6Corner.root","RECREATE");
  f2=new TFile("ESmR_6Corner.root","RECREATE");
  f3=new TFile("ESpF_6Corner.root","RECREATE");
  f4=new TFile("ESpR_6Corner.root","RECREATE");
 }

 f1->cd();
 h_ESmF_resiX_Right->Write();
 h_ESmF_resiX_Left->Write();
 h_ESmF_resiX_TopL->Write();
 h_ESmF_resiX_BottomL->Write();
 h_ESmF_resiX_TopR->Write();
 h_ESmF_resiX_BottomR->Write();

 f2->cd();
 h_ESmR_resiY_Right->Write();
 h_ESmR_resiY_Left->Write();
 h_ESmR_resiY_TopL->Write();
 h_ESmR_resiY_BottomL->Write();
 h_ESmR_resiY_TopR->Write();
 h_ESmR_resiY_BottomR->Write();

 f3->cd();
 h_ESpF_resiX_Right->Write();
 h_ESpF_resiX_Left->Write();
 h_ESpF_resiX_TopL->Write();
 h_ESpF_resiX_BottomL->Write();
 h_ESpF_resiX_TopR->Write();
 h_ESpF_resiX_BottomR->Write();

 f4->cd();
 h_ESpR_resiY_Right->Write();
 h_ESpR_resiY_Left->Write();
 h_ESpR_resiY_TopL->Write();
 h_ESpR_resiY_BottomL->Write();
 h_ESpR_resiY_TopR->Write();
 h_ESpR_resiY_BottomR->Write();

}

void ESAlign::Draw_4Corner()
{
 //gROOT->ProcessLine(".L ./format/tdrstyle.C");
 //setTDRStyle();
 //gStyle->SetPalette(1);
   
// TChain *EcalRecHit = new TChain("EcalRecHit");
// EcalRecHit->Add("AlignmentFile_*.root");
 double X=0.0;double Y=0.0;

 TH1D *h_ESpF_resiX_Right=new TH1D("h_ESpF_resiX_Right","Right",60,-3,3);
 TH1D *h_ESpF_resiX_Left=new TH1D("h_ESpF_resiX_Left","Left",60,-3,3);
 TH1D *h_ESpF_resiX_Top=new TH1D("h_ESpF_resiX_Top","Top",60,-3,3);
 TH1D *h_ESpF_resiX_Bottom=new TH1D("h_ESpF_resiX_Bottom","Bottom",60,-3,3);
 TH1D *h_ESpR_resiY_Right=new TH1D("h_ESpR_resiY_Right","Right",60,-3,3);
 TH1D *h_ESpR_resiY_Left=new TH1D("h_ESpR_resiY_Left","Left",60,-3,3);
 TH1D *h_ESpR_resiY_Top=new TH1D("h_ESpR_resiY_Top","Top",60,-3,3);
 TH1D *h_ESpR_resiY_Bottom=new TH1D("h_ESpR_resiY_Bottom","Bottom",60,-3,3);
 TH1D *h_ESmF_resiX_Right=new TH1D("h_ESmF_resiX_Right","Right",60,-3,3);
 TH1D *h_ESmF_resiX_Left=new TH1D("h_ESmF_resiX_Left","Left",60,-3,3);
 TH1D *h_ESmF_resiX_Top=new TH1D("h_ESmF_resiX_Top","Top",60,-3,3);
 TH1D *h_ESmF_resiX_Bottom=new TH1D("h_ESmF_resiX_Bottom","Bottom",60,-3,3);
 TH1D *h_ESmR_resiY_Right=new TH1D("h_ESmR_resiY_Right","Right",60,-3,3);
 TH1D *h_ESmR_resiY_Left=new TH1D("h_ESmR_resiY_Left","Left",60,-3,3);
 TH1D *h_ESmR_resiY_Top=new TH1D("h_ESmR_resiY_Top","Top",60,-3,3);
 TH1D *h_ESmR_resiY_Bottom=new TH1D("h_ESmR_resiY_Bottom","Bottom",60,-3,3);

 int NEntries=fChain->GetEntries();
 for(int i=0;i< NEntries ;i++)
 {
  fChain->GetEntry(i);

  for(int a=0;a<2;a++)
  {
   for(int b=0;b<2;b++)
   {
    for(int iTrk=0;iTrk<Ntrack;iTrk++)
    {
     if(a==0&&PredictionState_iz[iTrk][a][b]!=-1) continue;
     if(a==1&&PredictionState_iz[iTrk][a][b]!=1) continue;
     if(b==0&&PredictionState_ip[iTrk][a][b]!=1) continue;
     if(b==1&&PredictionState_ip[iTrk][a][b]!=2) continue;
     if(PredictionState_MatchedRec[iTrk][a][b]>-1)
     {
      //X=_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b];
      //Y=_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b];
      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
      Y= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R21[a][b]
        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R22[a][b]
        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R23[a][b];
      if(a==0&&b==0)
      {
       if(X>85.&&fabs(Y)<20.)
        h_ESmF_resiX_Right->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X<-85.&&fabs(Y)<20.)
        h_ESmF_resiX_Left->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.)
        h_ESmF_resiX_Top->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.)
        h_ESmF_resiX_Bottom->Fill(PredictionState_resiX[iTrk][a][b]);
      }//end if ES-F
      if(a==0&&b==1)
      {
       if(X>85.&&fabs(Y)<20.)
        h_ESmR_resiY_Right->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X<-85.&&fabs(Y)<20.)
        h_ESmR_resiY_Left->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.)
        h_ESmR_resiY_Top->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.)
        h_ESmR_resiY_Bottom->Fill(PredictionState_resiY[iTrk][a][b]);
      }//end if ES-R
      if(a==1&&b==0)
      {
       if(X>85.&&fabs(Y)<20.)
        h_ESpF_resiX_Right->Fill(PredictionState_resiX[iTrk][a][b]);
       if(X<-85.&&fabs(Y)<20.)
        h_ESpF_resiX_Left->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.)
        h_ESpF_resiX_Top->Fill(PredictionState_resiX[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.)
        h_ESpF_resiX_Bottom->Fill(PredictionState_resiX[iTrk][a][b]);
      }//end if ES+F
      if(a==1&&b==1)
      {
       if(X>85.&&fabs(Y)<20.)
        h_ESpR_resiY_Right->Fill(PredictionState_resiY[iTrk][a][b]);
       if(X<-85.&&fabs(Y)<20.)
        h_ESpR_resiY_Left->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y>85.&&fabs(X)<20.)
        h_ESpR_resiY_Top->Fill(PredictionState_resiY[iTrk][a][b]);
       if(Y<-85.&&fabs(X)<20.)
        h_ESpR_resiY_Bottom->Fill(PredictionState_resiY[iTrk][a][b]);
      }//end if ES+R
     }//end if matched
    }//end for-loop iTrk
   }//end for-loop ip
  }//end for-loop iz

 }//end loop all entries


 double mean=0.; double sig=0.;

 TFile *f1,*f2,*f3,*f4;

 if(iterN==1)
 {
  f1=new TFile("ESmF_4Corner_iter1.root","RECREATE");
  f2=new TFile("ESmR_4Corner_iter1.root","RECREATE");
  f3=new TFile("ESpF_4Corner_iter1.root","RECREATE");
  f4=new TFile("ESpR_4Corner_iter1.root","RECREATE");
 }
 if(iterN==11)
 {
  f1=new TFile("ESmF_4Corner_iter11.root","RECREATE");
  f2=new TFile("ESmR_4Corner_iter11.root","RECREATE");
  f3=new TFile("ESpF_4Corner_iter11.root","RECREATE");
  f4=new TFile("ESpR_4Corner_iter11.root","RECREATE");
 }
 if(iterN!=1&&iterN!=11)
 {
  f1=new TFile("ESmF_4Corner.root","RECREATE");
  f2=new TFile("ESmR_4Corner.root","RECREATE");
  f3=new TFile("ESpF_4Corner.root","RECREATE");
  f4=new TFile("ESpR_4Corner.root","RECREATE");
 }

 /*
 =new TFile("ESmF_4Corner_iter1.root","RECREATE");
 TFile *f2=new TFile("ESmR_4Corner_iter1.root","RECREATE");
 TFile *f3=new TFile("ESpF_4Corner_iter1.root","RECREATE");
 TFile *f4=new TFile("ESpR_4Corner_iter1.root","RECREATE");*/

// h_ESmF_resiX_Right->Draw();
// mean=h_ESmF_resiX_Right->GetMean();
// h_ESmF_resiX_Right->Fit("gaus","","",mean-1.5,mean+1.5);
// h_ESmF_resiX_Right->GetFunction("gaus")->SetLineStyle(2);
// mean=h_ESmF_resiX_Right->GetFunction("gaus")->GetParameter(1);
// sig=h_ESmF_resiX_Right->GetFunction("gaus")->GetParameter(2);
// h_ESmF_resiX_Right->Fit("gaus","","pE1",mean-2*sig,mean+2*sig);

 f1->cd();
 h_ESmF_resiX_Right->Write();
 h_ESmF_resiX_Left->Write();
 h_ESmF_resiX_Top->Write();
 h_ESmF_resiX_Bottom->Write();

 f2->cd();
 h_ESmR_resiY_Right->Write();
 h_ESmR_resiY_Left->Write();
 h_ESmR_resiY_Top->Write();
 h_ESmR_resiY_Bottom->Write();

 f3->cd();
 h_ESpF_resiX_Right->Write();
 h_ESpF_resiX_Left->Write();
 h_ESpF_resiX_Top->Write();
 h_ESpF_resiX_Bottom->Write();

 f4->cd();
 h_ESpR_resiY_Right->Write();
 h_ESpR_resiY_Left->Write();
 h_ESpR_resiY_Top->Write();
 h_ESpR_resiY_Bottom->Write();

}

//void ESAlign::Draw_HistoOnPlanes()
//{
// //gROOT->ProcessLine(".L ./format/tdrstyle.C");
// setTDRStyle();
// gStyle->SetPalette(1);
//   
//// TChain *EcalRecHit = new TChain("EcalRecHit");
//// EcalRecHit->Add("AlignmentFile_*.root");
//
// // XP1YP5 | XP2YP5 | XP3YP5          |                  110
// //        |        |                 |
// //----------------------------------------------------   85
// // XP1YP3 | XP2YP3 | XP3YP4 | XP4YP4 | XP5YP3
// //        |        |------------------                  22.5
// //        |        | XP3YP3 | XP4YP3 |
// //----------------------------------------------------   40
// //                 | XP3YP2          | XP5YP2
// //----------------------------------------------------   20
// //                 | XP3YP1          | XP5YP1
// //----------------------------------------------------
// //     20     40     22.5   85            110
//  
//
// double X=0.0;double Y=0.0;
//
// TH1D *h_ESpF_resiX_XP1YP3=new TH1D("h_ESpF_resiX_XP1YP3","XP1YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XP1YP5=new TH1D("h_ESpF_resiX_XP1YP5","XP1YP5",60,-3,3);
// TH1D *h_ESpF_resiX_XP2YP3=new TH1D("h_ESpF_resiX_XP2YP3","XP2YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XP2YP5=new TH1D("h_ESpF_resiX_XP2YP5","XP2YP5",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YP1=new TH1D("h_ESpF_resiX_XP3YP1","XP3YP1",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YP2=new TH1D("h_ESpF_resiX_XP3YP2","XP3YP2",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YP3=new TH1D("h_ESpF_resiX_XP3YP3","XP3YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YP4=new TH1D("h_ESpF_resiX_XP3YP4","XP3YP4",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YP5=new TH1D("h_ESpF_resiX_XP3YP5","XP3YP5",60,-3,3);
// TH1D *h_ESpF_resiX_XP4YP3=new TH1D("h_ESpF_resiX_XP4YP3","XP4YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XP4YP4=new TH1D("h_ESpF_resiX_XP4YP4","XP4YP4",60,-3,3);
// TH1D *h_ESpF_resiX_XP5YP1=new TH1D("h_ESpF_resiX_XP5YP1","XP5YP1",60,-3,3);
// TH1D *h_ESpF_resiX_XP5YP2=new TH1D("h_ESpF_resiX_XP5YP2","XP5YP2",60,-3,3);
// TH1D *h_ESpF_resiX_XP5YP3=new TH1D("h_ESpF_resiX_XP5YP3","XP5YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XP1YN3=new TH1D("h_ESpF_resiX_XP1YN3","XP1YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XP1YN5=new TH1D("h_ESpF_resiX_XP1YN5","XP1YN5",60,-3,3);
// TH1D *h_ESpF_resiX_XP2YN3=new TH1D("h_ESpF_resiX_XP2YN3","XP2YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XP2YN5=new TH1D("h_ESpF_resiX_XP2YN5","XP2YN5",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YN1=new TH1D("h_ESpF_resiX_XP3YN1","XP3YN1",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YN2=new TH1D("h_ESpF_resiX_XP3YN2","XP3YN2",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YN3=new TH1D("h_ESpF_resiX_XP3YN3","XP3YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YN4=new TH1D("h_ESpF_resiX_XP3YN4","XP3YN4",60,-3,3);
// TH1D *h_ESpF_resiX_XP3YN5=new TH1D("h_ESpF_resiX_XP3YN5","XP3YN5",60,-3,3);
// TH1D *h_ESpF_resiX_XP4YN3=new TH1D("h_ESpF_resiX_XP4YN3","XP4YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XP4YN4=new TH1D("h_ESpF_resiX_XP4YN4","XP4YN4",60,-3,3);
// TH1D *h_ESpF_resiX_XP5YN1=new TH1D("h_ESpF_resiX_XP5YN1","XP5YN1",60,-3,3);
// TH1D *h_ESpF_resiX_XP5YN2=new TH1D("h_ESpF_resiX_XP5YN2","XP5YN2",60,-3,3);
// TH1D *h_ESpF_resiX_XP5YN3=new TH1D("h_ESpF_resiX_XP5YN3","XP5YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XN1YP3=new TH1D("h_ESpF_resiX_XN1YP3","XN1YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XN1YP5=new TH1D("h_ESpF_resiX_XN1YP5","XN1YP5",60,-3,3);
// TH1D *h_ESpF_resiX_XN2YP3=new TH1D("h_ESpF_resiX_XN2YP3","XN2YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XN2YP5=new TH1D("h_ESpF_resiX_XN2YP5","XN2YP5",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YP1=new TH1D("h_ESpF_resiX_XN3YP1","XN3YP1",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YP2=new TH1D("h_ESpF_resiX_XN3YP2","XN3YP2",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YP3=new TH1D("h_ESpF_resiX_XN3YP3","XN3YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YP4=new TH1D("h_ESpF_resiX_XN3YP4","XN3YP4",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YP5=new TH1D("h_ESpF_resiX_XN3YP5","XN3YP5",60,-3,3);
// TH1D *h_ESpF_resiX_XN4YP3=new TH1D("h_ESpF_resiX_XN4YP3","XN4YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XN4YP4=new TH1D("h_ESpF_resiX_XN4YP4","XN4YP4",60,-3,3);
// TH1D *h_ESpF_resiX_XN5YP1=new TH1D("h_ESpF_resiX_XN5YP1","XN5YP1",60,-3,3);
// TH1D *h_ESpF_resiX_XN5YP2=new TH1D("h_ESpF_resiX_XN5YP2","XN5YP2",60,-3,3);
// TH1D *h_ESpF_resiX_XN5YP3=new TH1D("h_ESpF_resiX_XN5YP3","XN5YP3",60,-3,3);
// TH1D *h_ESpF_resiX_XN1YN3=new TH1D("h_ESpF_resiX_XN1YN3","XN1YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XN1YN5=new TH1D("h_ESpF_resiX_XN1YN5","XN1YN5",60,-3,3);
// TH1D *h_ESpF_resiX_XN2YN3=new TH1D("h_ESpF_resiX_XN2YN3","XN2YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XN2YN5=new TH1D("h_ESpF_resiX_XN2YN5","XN2YN5",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YN1=new TH1D("h_ESpF_resiX_XN3YN1","XN3YN1",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YN2=new TH1D("h_ESpF_resiX_XN3YN2","XN3YN2",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YN3=new TH1D("h_ESpF_resiX_XN3YN3","XN3YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YN4=new TH1D("h_ESpF_resiX_XN3YN4","XN3YN4",60,-3,3);
// TH1D *h_ESpF_resiX_XN3YN5=new TH1D("h_ESpF_resiX_XN3YN5","XN3YN5",60,-3,3);
// TH1D *h_ESpF_resiX_XN4YN3=new TH1D("h_ESpF_resiX_XN4YN3","XN4YN3",60,-3,3);
// TH1D *h_ESpF_resiX_XN4YN4=new TH1D("h_ESpF_resiX_XN4YN4","XN4YN4",60,-3,3);
// TH1D *h_ESpF_resiX_XN5YN1=new TH1D("h_ESpF_resiX_XN5YN1","XN5YN1",60,-3,3);
// TH1D *h_ESpF_resiX_XN5YN2=new TH1D("h_ESpF_resiX_XN5YN2","XN5YN2",60,-3,3);
// TH1D *h_ESpF_resiX_XN5YN3=new TH1D("h_ESpF_resiX_XN5YN3","XN5YN3",60,-3,3);
//
//
//
// int NEntries=fChain->GetEntries();
// for(int i=0;i< NEntries ;i++)
// {
//  fChain->GetEntry(i);
//
//  for(int a=0;a<2;a++)
//  {
//   for(int b=0;b<2;b++)
//   {
//    for(int iTrk=0;iTrk<Ntrack;iTrk++)
//    {
//     if(a==0&&PredictionState_iz[iTrk][a][b]!=-1) continue;
//     if(a==1&&PredictionState_iz[iTrk][a][b]!=1) continue;
//     if(b==0&&PredictionState_ip[iTrk][a][b]!=1) continue;
//     if(b==1&&PredictionState_ip[iTrk][a][b]!=2) continue;
//     if(PredictionState_MatchedRec[iTrk][a][b]>-1)
//     {
//      //X=_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b];
//      //Y=_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b];
//      X= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R11[a][b]
//        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R12[a][b]
//        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R13[a][b];
//      Y= (_esRecHit_X[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_X[a][b])*ES_O_R21[a][b]
//        +(_esRecHit_Y[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Y[a][b])*ES_O_R22[a][b]
//        +(_esRecHit_Z[PredictionState_MatchedRec[iTrk][a][b]]-ES_O_Z[a][b])*ES_O_R23[a][b];
//      if(a==0&&b==0)
//      {
//       if(X>85.&&fabs(Y)<20.)
//        h_ESmF_resiX_Right->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(X<-85.&&fabs(Y)<20.)
//        h_ESmF_resiX_Left->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(Y>85.&&fabs(X)<20.&&X>6.39)
//        h_ESmF_resiX_TopR->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(Y<-85.&&fabs(X)<20.&&X>-6.39)
//        h_ESmF_resiX_BottomR->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(Y>85.&&fabs(X)<20.&&X<6.39)
//        h_ESmF_resiX_TopL->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(Y<-85.&&fabs(X)<20.&&X<-6.39)
//        h_ESmF_resiX_BottomL->Fill(PredictionState_resiX[iTrk][a][b]);
//      }//end if ES-F
//      if(a==0&&b==1)
//      {
//       if(X>85.&&fabs(Y)<20.)
//        h_ESmR_resiY_Right->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(X<-85.&&fabs(Y)<20.)
//        h_ESmR_resiY_Left->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(Y>85.&&fabs(X)<20.&&X>0.)
//        h_ESmR_resiY_TopR->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(Y<-85.&&fabs(X)<20.&&X>0.)
//        h_ESmR_resiY_BottomR->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(Y>85.&&fabs(X)<20.&&X<0.)
//        h_ESmR_resiY_TopL->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(Y<-85.&&fabs(X)<20.&&X<0.)
//        h_ESmR_resiY_BottomL->Fill(PredictionState_resiY[iTrk][a][b]);
//      }//end if ES-R
//      if(a==1&&b==0)
//      {
//       if(X>85.&&fabs(Y)<20.)
//        h_ESpF_resiX_Right->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(X<-85.&&fabs(Y)<20.)
//        h_ESpF_resiX_Left->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(Y>85.&&fabs(X)<20.&&X>-6.39)
//        h_ESpF_resiX_TopR->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(Y<-85.&&fabs(X)<20.&&X>6.39)
//        h_ESpF_resiX_BottomR->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(Y>85.&&fabs(X)<20.&&X<-6.39)
//        h_ESpF_resiX_TopL->Fill(PredictionState_resiX[iTrk][a][b]);
//       if(Y<-85.&&fabs(X)<20.&&X<6.39)
//        h_ESpF_resiX_BottomL->Fill(PredictionState_resiX[iTrk][a][b]);
//      }//end if ES+F
//      if(a==1&&b==1)
//      {
//       if(X>85.&&fabs(Y)<20.)
//        h_ESpR_resiY_Right->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(X<-85.&&fabs(Y)<20.)
//        h_ESpR_resiY_Left->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(Y>85.&&fabs(X)<20.&&X>0.)
//        h_ESpR_resiY_TopR->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(Y<-85.&&fabs(X)<20.&&X>0.)
//        h_ESpR_resiY_BottomR->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(Y>85.&&fabs(X)<20.&&X<0.)
//        h_ESpR_resiY_TopL->Fill(PredictionState_resiY[iTrk][a][b]);
//       if(Y<-85.&&fabs(X)<20.&&X<0.)
//        h_ESpR_resiY_BottomL->Fill(PredictionState_resiY[iTrk][a][b]);
//      }//end if ES+R
//     }//end if matched
//    }//end for-loop iTrk
//   }//end for-loop ip
//  }//end for-loop iz
//
// }//end loop all entries
//
//
// double mean=0.; double sig=0.;
//
// TFile *f1,*f2,*f3,*f4;
//
// if(iterN==1)
// {
//  f1=new TFile("ESmF_6Corner_iter1.root","RECREATE");
//  f2=new TFile("ESmR_6Corner_iter1.root","RECREATE");
//  f3=new TFile("ESpF_6Corner_iter1.root","RECREATE");
//  f4=new TFile("ESpR_6Corner_iter1.root","RECREATE");
// }
// if(iterN==11)
// {
//  f1=new TFile("ESmF_6Corner_iter11.root","RECREATE");
//  f2=new TFile("ESmR_6Corner_iter11.root","RECREATE");
//  f3=new TFile("ESpF_6Corner_iter11.root","RECREATE");
//  f4=new TFile("ESpR_6Corner_iter11.root","RECREATE");
// }
// if(iterN!=1&&iterN!=11)
// {
//  f1=new TFile("ESmF_6Corner.root","RECREATE");
//  f2=new TFile("ESmR_6Corner.root","RECREATE");
//  f3=new TFile("ESpF_6Corner.root","RECREATE");
//  f4=new TFile("ESpR_6Corner.root","RECREATE");
// }
//
// f1->cd();
// h_ESmF_resiX_Right->Write();
// h_ESmF_resiX_Left->Write();
// h_ESmF_resiX_TopL->Write();
// h_ESmF_resiX_BottomL->Write();
// h_ESmF_resiX_TopR->Write();
// h_ESmF_resiX_BottomR->Write();
//
// f2->cd();
// h_ESmR_resiY_Right->Write();
// h_ESmR_resiY_Left->Write();
// h_ESmR_resiY_TopL->Write();
// h_ESmR_resiY_BottomL->Write();
// h_ESmR_resiY_TopR->Write();
// h_ESmR_resiY_BottomR->Write();
//
// f3->cd();
// h_ESpF_resiX_Right->Write();
// h_ESpF_resiX_Left->Write();
// h_ESpF_resiX_TopL->Write();
// h_ESpF_resiX_BottomL->Write();
// h_ESpF_resiX_TopR->Write();
// h_ESpF_resiX_BottomR->Write();
//
// f4->cd();
// h_ESpR_resiY_Right->Write();
// h_ESpR_resiY_Left->Write();
// h_ESpR_resiY_TopL->Write();
// h_ESpR_resiY_BottomL->Write();
// h_ESpR_resiY_TopR->Write();
// h_ESpR_resiY_BottomR->Write();
//
//}
