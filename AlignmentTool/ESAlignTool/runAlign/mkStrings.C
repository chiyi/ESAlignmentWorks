{
 Double_t ES_Alpha0[2][2]; Double_t ES_Beta0[2][2]; Double_t ES_Gamma0[2][2];
 Double_t ES_dAlpha[2][2]; Double_t ES_dBeta[2][2]; Double_t ES_dGamma[2][2];
 //ES+Front Plane
 ES_Alpha0[1][0]=-0.0017; ES_Beta0[1][0]=0.0000; ES_Gamma0[1][0]=0.00164;
 //ES+Rear Plane
 ES_Alpha0[1][1]=-0.0017; ES_Beta0[1][1]=0.0001; ES_Gamma0[1][1]=0.00138;
 //ES-Front Plane
 ES_Alpha0[0][0]=0.0016; ES_Beta0[0][0]=0.0009; ES_Gamma0[0][0]=0.00226;
 //ES-Rear Plane
 ES_Alpha0[0][1]=0.0017; ES_Beta0[0][1]=0.0005; ES_Gamma0[0][1]=0.00173;

 for(int i=0;i<2;i++)
 {
  for(int j=0;j<2;j++)
  {
   ES_dAlpha[i][j]=-ES_Alpha0[i][j];
   ES_dBeta[i][j]=-ES_Beta0[i][j];
   ES_dGamma[i][j]=-ES_Gamma0[i][j];
 }}

  int iterN=1;
  ofstream file;
  //file.open("tmp.txt",ios::out|ios::app);
  file.open("file.txt_Sol1",ios::out);
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESpFdAlpha = cms.double("<<ES_dAlpha[1][0]<<")\n";
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESpFdBeta = cms.double("<<ES_dBeta[1][0]<<")\n";
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESpFdGamma = cms.double("<<ES_dGamma[1][0]<<")\n";

  file<<"process.esAlignTool.Iter"<<iterN<<"_ESpRdAlpha = cms.double("<<ES_dAlpha[1][1]<<")\n";
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESpRdBeta = cms.double("<<ES_dBeta[1][1]<<")\n";
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESpRdGamma = cms.double("<<ES_dGamma[1][1]<<")\n";

  file<<"process.esAlignTool.Iter"<<iterN<<"_ESmFdAlpha = cms.double("<<ES_dAlpha[0][0]<<")\n";
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESmFdBeta = cms.double("<<ES_dBeta[0][0]<<")\n";
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESmFdGamma = cms.double("<<ES_dGamma[0][0]<<")\n";

  file<<"process.esAlignTool.Iter"<<iterN<<"_ESmRdAlpha = cms.double("<<ES_dAlpha[0][1]<<")\n";
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESmRdBeta = cms.double("<<ES_dBeta[0][1]<<")\n";
  file<<"process.esAlignTool.Iter"<<iterN<<"_ESmRdGamma = cms.double("<<ES_dGamma[0][1]<<")\n";

 file<<"\n";
 file.close();

  ofstream file2;
  //file2.open("file2.txt",ios::out|ios::app);
  file2.open("file2.txt_Sol1",ios::out);
  file2<<"ESpFdAlpha-="<<ES_dAlpha[1][0]<<"; ESpFdBeta-="<<ES_dBeta[1][0]<<"; ESpFdGamma-="<<ES_dGamma[1][0]<<";\n";
  file2<<"ESpRdAlpha-="<<ES_dAlpha[1][1]<<"; ESpRdBeta-="<<ES_dBeta[1][1]<<"; ESpRdGamma-="<<ES_dGamma[1][1]<<";\n";
  file2<<"ESmFdAlpha-="<<ES_dAlpha[0][0]<<"; ESmFdBeta-="<<ES_dBeta[0][0]<<"; ESmFdGamma-="<<ES_dGamma[0][0]<<";\n";
  file2<<"ESmRdAlpha-="<<ES_dAlpha[0][1]<<"; ESmRdBeta-="<<ES_dBeta[0][1]<<"; ESmRdGamma-="<<ES_dGamma[0][1]<<";\n";

  file2<<"\n";
  file2.close();

}
