#include "ESAlign.C"

void runPlot()
{
 TChain * chn1 = new TChain("ESAlign","");

 chn1->Add("AlignmentFile_*.root");


 ESAlign a(chn1);
 a.iterN=11;
 //a.Loop();
 //a.Draw_4quad();
 //a.Draw_vsXY();
 a.Draw_4Corner();

}


