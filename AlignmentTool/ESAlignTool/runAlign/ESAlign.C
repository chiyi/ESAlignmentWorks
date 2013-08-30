#define ESAlign_cxx

#include <iostream>
#include <fstream>

#include "ESAlign.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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
   valD=ES_CHI2[i][j]; sprintf(buf,"ES_CHI2[%i][%i]+= (",i,j);
   oText<<buf<<valD<<");";

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
