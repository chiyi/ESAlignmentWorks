
void ESAlignTool::Cal_MatrixM(int iz,int ip,Double_t eF_xx,Double_t eF_yx,Double_t eF_yy,Double_t derX,Double_t derY,Double_t determinant)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;
 
  ES_M11[a][b] += (eF_yy/determinant);
  ES_M12[a][b] += (-eF_yx/determinant);
  ES_M13[a][b] += ( (-eF_yy*derX+eF_yx*derY)/determinant );
  ES_M22[a][b] += (eF_xx/determinant);
  ES_M23[a][b] += ( (-eF_xx*derY+eF_yx*derX)/determinant );
  ES_M33[a][b] += ( ( eF_yy*pow(derX,2.)+eF_xx*pow(derY,2.)-2.*eF_yx*derX*derY
                    )/determinant);
 }
}

void ESAlignTool::Cal_VectorP(int iz,int ip,Double_t eF_xx,Double_t eF_yx,Double_t eF_yy,Double_t derX,Double_t derY,Double_t determinant,Double_t XDF,Double_t YDF)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  ES_P1[a][b] += ( (-eF_yy*XDF+eF_yx*YDF)/determinant );
  ES_P2[a][b] += ( (eF_yx*XDF-eF_xx*YDF)/determinant );
  ES_P3[a][b] += ( ( (eF_yy*derX-eF_yx*derY)*XDF + (eF_xx*derY-eF_yx*derX)*YDF )
                   /determinant );
 }
}

void ESAlignTool::Cal_CHI2(int iz,int ip,Double_t eF_xx,Double_t eF_yx,Double_t eF_yy,Double_t XDF,Double_t YDF)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  ES_CHI2[a][b] += (( eF_yy*XDF*XDF-2.*eF_yx*XDF*YDF+eF_xx*YDF*YDF )/(eF_xx*eF_yy-eF_yx*eF_yx));
 }
}

void ESAlignTool::Cal_MatrixMErr2(int iz,int ip,Double_t Exx,Double_t Eyx,Double_t Eyy,Double_t Px,Double_t Py,Double_t Pz,Double_t sPxx,Double_t sPyy,Double_t sPzz,Double_t sPxy,Double_t sPxz,Double_t sPyz)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  ES_M31Err2[a][b] +=
  ( (Eyy*Eyy*(Pz*Pz*sPxx - 2.*Px*Pz*sPxz + Px*Px*sPzz) 
   - 2.*Eyx*Eyy*(Pz*Pz*sPxy - Pz*(Py*sPxz + Px*sPyz) 
   + Px*Py*sPzz) + Eyx*Eyx*(Pz*Pz*sPyy - 2.*Py*Pz*sPyz 
   + Py*Py*sPzz))/(Pz*Pz*Pz*Pz) );
  ES_M32Err2[a][b] +=
  ( (Eyx*Eyx*(Pz*Pz*sPxx - 2.*Px*Pz*sPxz + Px*Px*sPzz)
   - 2.*Exx*Eyx*(Pz*Pz*sPxy - Pz*(Py*sPxz + Px*sPyz)
   + Px*Py*sPzz) + Exx*Exx*(Pz*Pz*sPyy - 2.*Py*Pz*sPyz 
   + Py*Py*sPzz))/(Pz*Pz*Pz*Pz) );
  ES_M33Err2[a][b] +=
  ( (4.*(Eyy*Eyy*Px*Px*(Pz*Pz*sPxx - 2.*Px*Pz*sPxz + Px*Px*sPzz) 
    + Exx*Exx*Py*Py*(Pz*Pz*sPyy - 2.*Py*Pz*sPyz + Py*Py*sPzz) 
    + 2.*Exx*Eyx*Py*(-(Px*Pz*Pz*sPyy) + Py*Pz*(-(Pz*sPxy) + 3.*Px*sPyz) 
    + Py*Py*(Pz*sPxz - 2.*Px*sPzz)) 
    + Eyx*Eyx*(Px*Px*Pz*Pz*sPyy + 2*Px*Py*Pz*(Pz*sPxy - 2*Px*sPyz) 
    + Py*Py*(Pz*Pz*sPxx - 4.*Px*Pz*sPxz + 4.*Px*Px*sPzz)) 
    - 2.*Eyy*Px*(Exx*Py*(-(Pz*Pz*sPxy) + Py*Pz*sPxz + Px*Pz*sPyz
    - Px*Py*sPzz) + Eyx*(Px*Pz*(Pz*sPxy - Px*sPyz) 
    + Py*(Pz*Pz*sPxx - 3.*Px*Pz*sPxz + 2.*Px*Px*sPzz)))))/(Pz*Pz*Pz*Pz*Pz*Pz)
  );
 } 
}

void ESAlignTool::Cal_VectorPErr2(int iz,int ip,Double_t Exx,Double_t Eyx,Double_t Eyy,Double_t Px,Double_t Py,Double_t Pz,Double_t sPxx,Double_t sPyy,Double_t sPzz,Double_t sPxy,Double_t sPxz,Double_t sPyz,Double_t Rxx,Double_t Ryy,Double_t Xpre,Double_t Ypre,Double_t Xrec,Double_t Yrec,Double_t sE14,Double_t sE15,Double_t sE16,Double_t sE24,Double_t sE25,Double_t sE26,Double_t determinant)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  ES_P1Err2[a][b] +=
  ( - (Eyx*Eyx*Eyy)/pow(determinant,2.) 
    + (Exx*Eyy*Eyy)/pow(determinant,2.) 
    + (Eyy*Eyy*Rxx)/pow(determinant,2.) 
    + (Eyx*Eyx*Ryy)/pow(determinant,2.)
  );
  ES_P2Err2[a][b] +=
  ( - (Exx*Eyx*Eyx)/pow(determinant,2.) 
    + (Exx*Exx*Eyy)/pow(determinant,2.) 
    + (Eyx*Eyx*Rxx)/pow(determinant,2.) 
    + (Exx*Exx*Ryy)/pow(determinant,2.)
  );
  ES_P3Err2[a][b] +=
  ( (Eyy*pow(Eyx*Px - Exx*Py,2.)*Pz*Pz + Exx*pow(Eyy*Px - Eyx*Py,2.)*Pz*Pz 
   - 2.*Eyx*(Eyx*Px - Exx*Py)*(-(Eyy*Px) + Eyx*Py)*Pz*Pz 
   + pow(Eyy*Px - Eyx*Py,2.)*Pz*Pz*Rxx + pow(Eyx*Px - Exx*Py,2.)*Pz*Pz*Ryy
   - 2.*(Eyy*Px - Eyx*Py)*Pz*Pz*sE15*(Eyx*(Xpre - Xrec) + Exx*(Ypre - Yrec))
   - 2.*(Eyx*Px - Exx*Py)*Pz*Pz*sE25*(Eyx*(Xpre - Xrec) + Exx*(Ypre - Yrec))
   + Pz*Pz*sPyy*pow(Eyx*(Xpre - Xrec) + Exx*(Ypre - Yrec),2.) 
   + 2.*(Eyy*Px - Eyx*Py)*Pz*Pz*sE14*(Eyy*(Xpre - Xrec) 
   + Eyx*(Ypre - Yrec)) + 2.*(Eyx*Px - Exx*Py)*Pz*Pz*sE24*(Eyy*(Xpre - Xrec)
   + Eyx*(Ypre - Yrec))
   + Pz*Pz*sPxx*pow(Eyy*(Xpre - Xrec) + Eyx*(Ypre - Yrec),2.) 
   - 2.*(Eyy*Px - Eyx*Py)*Pz*sE16*(Eyy*Px*(Xpre - Xrec) 
   + Exx*Py*(-Ypre + Yrec) + Eyx*(-(Py*Xpre) + Py*Xrec + Px*Ypre - Px*Yrec))
   - 2.*(Eyx*Px - Exx*Py)*Pz*sE26*(Eyy*Px*(Xpre - Xrec)
   + Exx*Py*(-Ypre + Yrec) + Eyx*(-(Py*Xpre) + Py*Xrec + Px*Ypre - Px*Yrec))
   + sPzz*pow( Eyy*Px*(Xpre - Xrec) + Exx*Py*(-Ypre + Yrec)
              + Eyx*( -(Py*Xpre) + Py*Xrec + Px*Ypre - Px*Yrec)
          ,2.)
     )/(pow(Eyx*Eyx - Exx*Eyy,2.)*Pz*Pz*Pz*Pz)
  );

 }
}

/*
void ESAlignTool::Cal_JacobianMatrix_forHuman(int iz,int ip,Double_t e_xx,Double_t e_yx,Double_t e_yy,Double_t determinant,int iTrk,Double_t &J11,Double_t &J12,Double_t &J21,Double_t &J22,Double_t &J31,Double_t &J32,Double_t &J41,Double_t &J42,Double_t &J51,Double_t &J52,Double_t &J61,Double_t &J62)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  Double_t R11=ES_R11[a][b];
  Double_t R12=ES_R12[a][b];
  Double_t R13=ES_R13[a][b];
  Double_t R21=ES_R21[a][b];
  Double_t R22=ES_R22[a][b];
  Double_t R23=ES_R23[a][b];
  Double_t R31=ES_R31[a][b];
  Double_t R32=ES_R32[a][b];
  Double_t R33=ES_R33[a][b];
  
  Double_t dAlpha=ES_dAlpha[a][b];
  Double_t dBeta=ES_dBeta[a][b];
  Double_t dGamma=ES_dGamma[a][b];
  
  Double_t DR11_Alpha = -(cos(dAlpha)*sin(dBeta)*sin(dGamma));
  Double_t DR11_Beta = -(cos(dGamma)*sin(dBeta)) - cos(dBeta)*sin(dAlpha)*sin(dGamma);
  Double_t DR11_Gamma = -(cos(dGamma)*sin(dAlpha)*sin(dBeta)) - cos(dBeta)*sin(dGamma);
  Double_t DR12_Alpha = cos(dAlpha)*cos(dGamma)*sin(dBeta);
  Double_t DR12_Beta = cos(dBeta)*cos(dGamma)*sin(dAlpha) - sin(dBeta)*sin(dGamma);
  Double_t DR12_Gamma = cos(dBeta)*cos(dGamma) - sin(dAlpha)*sin(dBeta)*sin(dGamma);
  Double_t DR13_Alpha = sin(dAlpha)*sin(dBeta);
  Double_t DR13_Beta = -(cos(dAlpha)*cos(dBeta));
  Double_t DR13_Gamma = 0;
  Double_t DR21_Alpha = sin(dAlpha)*sin(dGamma);
  Double_t DR21_Beta = 0;
  Double_t DR21_Gamma = -(cos(dAlpha)*cos(dGamma));
  Double_t DR22_Alpha = -(cos(dGamma)*sin(dAlpha));
  Double_t DR22_Beta = 0;
  Double_t DR22_Gamma = -(cos(dAlpha)*sin(dGamma));
  Double_t DR23_Alpha = cos(dAlpha);
  Double_t DR23_Beta = 0;
  Double_t DR23_Gamma = 0;

  Double_t DR31_Alpha= sin(dAlpha)*sin(dGamma);
  Double_t DR31_Beta=0;
  Double_t DR31_Gamma=-( cos(dAlpha)*cos(dGamma));
  Double_t DR32_Alpha=-( cos(dGamma)*sin(dAlpha));
  Double_t DR32_Beta=0;
  Double_t DR32_Gamma=-( cos(dAlpha)*sin(dGamma));
  Double_t DR33_Alpha= cos(dAlpha);
  Double_t DR33_Beta=0;
  Double_t DR33_Gamma=0;

  Double_t Xpre = PredictionState_X[iTrk][a][b];
  Double_t Ypre = PredictionState_Y[iTrk][a][b];
  Double_t Zpre = PredictionState_Z[iTrk][a][b];
  Double_t Px = PredictionState_Px[iTrk][a][b];
  Double_t Py = PredictionState_Py[iTrk][a][b];
  Double_t Pz = PredictionState_Pz[iTrk][a][b];
  Double_t Xap = ES_Oap_X[a][b];
  Double_t Yap = ES_Oap_Y[a][b];
  Double_t Zap = ES_Oap_Z[a][b];

  Double_t PzLocal = R31*Px + R32*Py + R33*Pz;

  J11
  = -R11 + ((Px*R11 + Py*R12 + Pz*R13)*R31)/PzLocal;

  J12 
  = -R21 + ((Px*R21 + Py*R22 + Pz*R23)*R31)/PzLocal;

  J21 
  = -R12 + ((Px*R11 + Py*R12 + Pz*R13)*R32)/PzLocal;

  J22 
  = -R22 + ((Px*R21 + Py*R22 + Pz*R23)*R32)/PzLocal;

  J31 
  = -R13 + ((Px*R11 + Py*R12 + Pz*R13)*R33)/PzLocal;

  J32 
  = -R23 + ((Px*R21 + Py*R22 + Pz*R23)*R33)/PzLocal;

  J41 
  = (DR11_Alpha*PzLocal*(Xpre - Xap) - (Px*R11 + Py*R12 + Pz*R13)*
       (DR31_Alpha*(Xap - Xpre) + DR32_Alpha*Yap - DR32_Alpha*Ypre + 
        DR33_Alpha*Zap - DR33_Alpha*Zpre) + PzLocal*(DR12_Alpha*(Ypre - Yap) + 
        DR13_Alpha*(-Zap + Zpre)))/PzLocal;

  J42 
  = (DR21_Alpha*PzLocal*(Xpre - Xap) - (Px*R21 + Py*R22 + Pz*R23)*
       (DR31_Alpha*(Xap - Xpre) + DR32_Alpha*Yap - DR32_Alpha*Ypre + 
        DR33_Alpha*Zap - DR33_Alpha*Zpre) + PzLocal*(DR22_Alpha*(Ypre - Yap) + 
        DR23_Alpha*(-Zap + Zpre)))/PzLocal;

  J51 
  = (DR11_Beta*PzLocal*(Xpre - Xap) - (Px*R11 + Py*R12 + Pz*R13)*
       (DR31_Beta*(Xap - Xpre) + DR32_Beta*Yap - DR32_Beta*Ypre + 
        DR33_Beta*Zap - DR33_Beta*Zpre) + PzLocal*(DR12_Beta*(Ypre - Yap) + 
        DR13_Beta*(-Zap + Zpre)))/PzLocal;

  J52 
  = (DR21_Beta*PzLocal*(Xpre - Xap) - (Px*R21 + Py*R22 + Pz*R23)*
       (DR31_Beta*(Xap - Xpre) + DR32_Beta*Yap - DR32_Beta*Ypre + 
        DR33_Beta*Zap - DR33_Beta*Zpre) + PzLocal*(DR22_Beta*(Ypre - Yap) + 
        DR23_Beta*(-Zap + Zpre)))/PzLocal;

  J61 
  = (DR11_Gamma*PzLocal*(Xpre - Xap) - (Px*R11 + Py*R12 + Pz*R13)*
       (DR31_Gamma*(Xap - Xpre) + DR32_Gamma*Yap - DR32_Gamma*Ypre + 
        DR33_Gamma*Zap - DR33_Gamma*Zpre) + PzLocal*(DR12_Gamma*(Ypre - Yap) + 
        DR13_Gamma*(-Zap + Zpre)))/PzLocal;

  J62 
  = (DR21_Gamma*PzLocal*(Xpre - Xap) - (Px*R21 + Py*R22 + Pz*R23)*
       (DR31_Gamma*(Xap - Xpre) + DR32_Gamma*Yap - DR32_Gamma*Ypre + 
        DR33_Gamma*Zap - DR33_Gamma*Zpre) + PzLocal*(DR22_Gamma*(Ypre - Yap) + 
        DR23_Gamma*(-Zap + Zpre)))/PzLocal;

 }
}
*/

/*
void ESAlignTool::Cal_MatrixM_wRotation_forHuman(int iz,int ip,Double_t e_xx,Double_t e_yx,Double_t e_yy,Double_t derX,Double_t derY,Double_t determinant,Double_t J11,Double_t J12,Double_t J21,Double_t J22,Double_t J31,Double_t J32,Double_t J41,Double_t J42,Double_t J51,Double_t J52,Double_t J61,Double_t J62)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  //Adding MatrixM
  //M_ij*determinant
  //  = J_i1*e_yy*J_j1 - e_yx*(J_i1*J_j2+J_i2*J_j1) + J_i2*e_xx*J_j2;
  Double_t buf=0.;
  buf = J11*e_yy*J11 - (e_yx)*( J11*J12 + J12*J11 ) + J12*e_xx*J12 ;
  buf/=determinant; ES_M11[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J21 - (e_yx)*( J11*J22 + J12*J21 ) + J12*e_xx*J22 ;
  buf/=determinant; ES_M12[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J31 - (e_yx)*( J11*J32 + J12*J31 ) + J12*e_xx*J32 ;
  buf/=determinant; ES_M13[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J41 - (e_yx)*( J11*J42 + J12*J41 ) + J12*e_xx*J42 ;
  buf/=determinant; ES_M14[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J51 - (e_yx)*( J11*J52 + J12*J51 ) + J12*e_xx*J52 ;
  buf/=determinant; ES_M15[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J61 - (e_yx)*( J11*J62 + J12*J61 ) + J12*e_xx*J62 ;
  buf/=determinant; ES_M16[a][b]+=buf;

  buf=0.;
  buf = J21*e_yy*J21 - (e_yx)*( J21*J22 + J22*J21 ) + J22*e_xx*J22 ;
  buf/=determinant; ES_M22[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*J31 - (e_yx)*( J21*J32 + J22*J31 ) + J22*e_xx*J32 ;
  buf/=determinant; ES_M23[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*J41 - (e_yx)*( J21*J42 + J22*J41 ) + J22*e_xx*J42 ;
  buf/=determinant; ES_M24[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*J51 - (e_yx)*( J21*J52 + J22*J51 ) + J22*e_xx*J52 ;
  buf/=determinant; ES_M25[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*J61 - (e_yx)*( J21*J62 + J22*J61 ) + J22*e_xx*J62 ;
  buf/=determinant; ES_M26[a][b]+=buf;

  buf=0.;
  buf = J31*e_yy*J31 - (e_yx)*( J31*J32 + J32*J31 ) + J32*e_xx*J32 ;
  buf/=determinant; ES_M33[a][b]+=buf;
  buf=0.;
  buf = J31*e_yy*J41 - (e_yx)*( J31*J42 + J32*J41 ) + J32*e_xx*J42 ;
  buf/=determinant; ES_M34[a][b]+=buf;
  buf=0.;
  buf = J31*e_yy*J51 - (e_yx)*( J31*J52 + J32*J51 ) + J32*e_xx*J52 ;
  buf/=determinant; ES_M35[a][b]+=buf;
  buf=0.;
  buf = J31*e_yy*J61 - (e_yx)*( J31*J62 + J32*J61 ) + J32*e_xx*J62 ;
  buf/=determinant; ES_M36[a][b]+=buf;

  buf=0.;
  buf = J41*e_yy*J41 - (e_yx)*( J41*J42 + J42*J41 ) + J42*e_xx*J42 ;
  buf/=determinant; ES_M44[a][b]+=buf;
  buf=0.;
  buf = J41*e_yy*J51 - (e_yx)*( J41*J52 + J42*J51 ) + J42*e_xx*J52 ;
  buf/=determinant; ES_M45[a][b]+=buf;
  buf=0.;
  buf = J41*e_yy*J61 - (e_yx)*( J41*J62 + J42*J61 ) + J42*e_xx*J62 ;
  buf/=determinant; ES_M46[a][b]+=buf;

  buf=0.;
  buf = J51*e_yy*J51 - (e_yx)*( J51*J52 + J52*J51 ) + J52*e_xx*J52 ;
  buf/=determinant; ES_M55[a][b]+=buf;
  buf=0.;
  buf = J51*e_yy*J61 - (e_yx)*( J51*J62 + J52*J61 ) + J52*e_xx*J62 ;
  buf/=determinant; ES_M56[a][b]+=buf;

  buf=0.;
  buf = J61*e_yy*J61 - (e_yx)*( J61*J62 + J62*J61 ) + J62*e_xx*J62 ;
  buf/=determinant; ES_M66[a][b]+=buf;

 }
}
*/

/*
void ESAlignTool::Cal_VectorP_wRotation_forHuman(int iz,int ip,Double_t e_xx,Double_t e_yx,Double_t e_yy,Double_t derX,Double_t derY,Double_t determinant,Double_t XD,Double_t YD,Double_t J11,Double_t J12,Double_t J21,Double_t J22,Double_t J31,Double_t J32,Double_t J41,Double_t J42,Double_t J51,Double_t J52,Double_t J61,Double_t J62)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  //Adding MatrixM
  //vp_i*determinant
  //  = J_i1*e_yy*XD - e_yx*(J_i1*YD+J_i2*XD) + J_i2*e_xx*YD;
  Double_t buf=0.;
  buf = J11*e_yy*XD - e_yx*(J11*YD+J12*XD) + J12*e_xx*YD;
  buf/=determinant;  ES_P1[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*XD - e_yx*(J21*YD+J22*XD) + J22*e_xx*YD;
  buf/=determinant;  ES_P2[a][b]+=buf;
  buf=0.;
  buf = J31*e_yy*XD - e_yx*(J31*YD+J32*XD) + J32*e_xx*YD;
  buf/=determinant;  ES_P3[a][b]+=buf;
  buf=0.;
  buf = J41*e_yy*XD - e_yx*(J41*YD+J42*XD) + J42*e_xx*YD;
  buf/=determinant;  ES_P4[a][b]+=buf;
  buf=0.;
  buf = J51*e_yy*XD - e_yx*(J51*YD+J52*XD) + J52*e_xx*YD;
  buf/=determinant;  ES_P5[a][b]+=buf;
  buf=0.;
  buf = J61*e_yy*XD - e_yx*(J61*YD+J62*XD) + J62*e_xx*YD;
  buf/=determinant;  ES_P6[a][b]+=buf;
 }
}
*/

/*
void ESAlignTool::Cal_MatrixM_wRotation(int iTrk,int iz,int ip,Double_t e_xx,Double_t e_yx,Double_t e_yy)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  Double_t Px = PredictionState_Px[iTrk][a][b];
  Double_t Py = PredictionState_Py[iTrk][a][b];
  Double_t Pz = PredictionState_Pz[iTrk][a][b];
  Double_t Xap = ES_Oap_X[a][b];
  Double_t Yap = ES_Oap_Y[a][b];
  Double_t Zap = ES_Oap_Z[a][b];
  Double_t Xpre = PredictionState_X[iTrk][a][b];
  Double_t Ypre = PredictionState_Y[iTrk][a][b];
  Double_t Zpre = PredictionState_Z[iTrk][a][b];
  Double_t dAlpha = ES_dAlpha[a][b];
  Double_t dBeta = ES_dBeta[a][b];
  Double_t dGamma = ES_dGamma[a][b];

  Double_t R11=ES_R11[a][b];
  Double_t R12=ES_R12[a][b];
  Double_t R13=ES_R13[a][b];
  Double_t R21=ES_R21[a][b];
  Double_t R22=ES_R22[a][b];
  Double_t R23=ES_R23[a][b];
  Double_t R31=ES_R31[a][b];
  Double_t R32=ES_R32[a][b];
  Double_t R33=ES_R33[a][b];
  
  Double_t DR11_Alpha = -(cos(dAlpha)*sin(dBeta)*sin(dGamma));
  Double_t DR11_Beta = -(cos(dGamma)*sin(dBeta)) - cos(dBeta)*sin(dAlpha)*sin(dGamma);
  Double_t DR11_Gamma = -(cos(dGamma)*sin(dAlpha)*sin(dBeta)) - cos(dBeta)*sin(dGamma);
  Double_t DR12_Alpha = cos(dAlpha)*cos(dGamma)*sin(dBeta);
  Double_t DR12_Beta = cos(dBeta)*cos(dGamma)*sin(dAlpha) - sin(dBeta)*sin(dGamma);
  Double_t DR12_Gamma = cos(dBeta)*cos(dGamma) - sin(dAlpha)*sin(dBeta)*sin(dGamma);
  Double_t DR13_Alpha = sin(dAlpha)*sin(dBeta);
  Double_t DR13_Beta = -(cos(dAlpha)*cos(dBeta));
  Double_t DR13_Gamma = 0.;
  Double_t DR21_Alpha = sin(dAlpha)*sin(dGamma);
  Double_t DR21_Beta = 0.;
  Double_t DR21_Gamma = -(cos(dAlpha)*cos(dGamma));
  Double_t DR22_Alpha = -(cos(dGamma)*sin(dAlpha));
  Double_t DR22_Beta = 0.;
  Double_t DR22_Gamma = -(cos(dAlpha)*sin(dGamma));
  Double_t DR23_Alpha = cos(dAlpha);
  Double_t DR23_Beta = 0.;
  Double_t DR23_Gamma = 0.;

  Double_t DR31_Alpha= sin(dAlpha)*sin(dGamma);
  Double_t DR31_Beta=0.;
  Double_t DR31_Gamma=-( cos(dAlpha)*cos(dGamma));
  Double_t DR32_Alpha=-( cos(dGamma)*sin(dAlpha));
  Double_t DR32_Beta=0.;
  Double_t DR32_Gamma=-( cos(dAlpha)*sin(dGamma));
  Double_t DR33_Alpha= cos(dAlpha);
  Double_t DR33_Beta=0.;
  Double_t DR33_Gamma=0.;

  Double_t PzLocal = R31*Px + R32*Py + R33*Pz;
  Double_t J11,J12,J21,J22,J31,J32,J41,J42,J51,J52,J61,J62;
    
  J11
  = -R11 - ((Px*R11 + Py*R12 + Pz*R13)*R31)/PzLocal;

  J12 
  = -R21 - ((Px*R21 + Py*R22 + Pz*R23)*R31)/PzLocal;

  J21 
  = -R12 - ((Px*R11 + Py*R12 + Pz*R13)*R32)/PzLocal;

  J22 
  = -R22 - ((Px*R21 + Py*R22 + Pz*R23)*R32)/PzLocal;

  J31 
  = -R13 - ((Px*R11 + Py*R12 + Pz*R13)*R33)/PzLocal;

  J32 
  = -R23 - ((Px*R21 + Py*R22 + Pz*R23)*R33)/PzLocal;

  J41 = DR11_Alpha*(Xpre-Xap) + DR12_Alpha*(Ypre-Yap) + DR13_Alpha*(Zpre-Zap);
  J42 = DR21_Alpha*(Xpre-Xap) + DR22_Alpha*(Ypre-Yap) + DR23_Alpha*(Zpre-Zap);
  J51 = DR11_Beta*(Xpre-Xap) + DR12_Beta*(Ypre-Yap) + DR13_Beta*(Zpre-Zap);
  J52 = DR21_Beta*(Xpre-Xap) + DR22_Beta*(Ypre-Yap) + DR23_Beta*(Zpre-Zap);
  J61 = DR11_Gamma*(Xpre-Xap) + DR12_Gamma*(Ypre-Yap) + DR13_Gamma*(Zpre-Zap);
  J62 = DR21_Gamma*(Xpre-Xap) + DR22_Gamma*(Ypre-Yap) + DR23_Gamma*(Zpre-Zap);

  Double_t determinant = e_xx*e_yy - e_yx*e_yx;
  //  = J_i1*e_yy*J_j1 - e_yx*(J_i1*J_j2+J_i2*J_j1) + J_i2*e_xx*J_j2;
  Double_t buf=0.;
  buf = J11*e_yy*J11 - (e_yx)*( J11*J12 + J12*J11 ) + J12*e_xx*J12 ;
  buf/=determinant; ES_M11[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J21 - (e_yx)*( J11*J22 + J12*J21 ) + J12*e_xx*J22 ;
  buf/=determinant; ES_M12[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J31 - (e_yx)*( J11*J32 + J12*J31 ) + J12*e_xx*J32 ;
  buf/=determinant; ES_M13[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J41 - (e_yx)*( J11*J42 + J12*J41 ) + J12*e_xx*J42 ;
  buf/=determinant; ES_M14[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J51 - (e_yx)*( J11*J52 + J12*J51 ) + J12*e_xx*J52 ;
  buf/=determinant; ES_M15[a][b]+=buf;
  buf=0.;
  buf = J11*e_yy*J61 - (e_yx)*( J11*J62 + J12*J61 ) + J12*e_xx*J62 ;
  buf/=determinant; ES_M16[a][b]+=buf;

  buf=0.;
  buf = J21*e_yy*J21 - (e_yx)*( J21*J22 + J22*J21 ) + J22*e_xx*J22 ;
  buf/=determinant; ES_M22[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*J31 - (e_yx)*( J21*J32 + J22*J31 ) + J22*e_xx*J32 ;
  buf/=determinant; ES_M23[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*J41 - (e_yx)*( J21*J42 + J22*J41 ) + J22*e_xx*J42 ;
  buf/=determinant; ES_M24[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*J51 - (e_yx)*( J21*J52 + J22*J51 ) + J22*e_xx*J52 ;
  buf/=determinant; ES_M25[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*J61 - (e_yx)*( J21*J62 + J22*J61 ) + J22*e_xx*J62 ;
  buf/=determinant; ES_M26[a][b]+=buf;

  buf=0.;
  buf = J31*e_yy*J31 - (e_yx)*( J31*J32 + J32*J31 ) + J32*e_xx*J32 ;
  buf/=determinant; ES_M33[a][b]+=buf;
  buf=0.;
  buf = J31*e_yy*J41 - (e_yx)*( J31*J42 + J32*J41 ) + J32*e_xx*J42 ;
  buf/=determinant; ES_M34[a][b]+=buf;
  buf=0.;
  buf = J31*e_yy*J51 - (e_yx)*( J31*J52 + J32*J51 ) + J32*e_xx*J52 ;
  buf/=determinant; ES_M35[a][b]+=buf;
  buf=0.;
  buf = J31*e_yy*J61 - (e_yx)*( J31*J62 + J32*J61 ) + J32*e_xx*J62 ;
  buf/=determinant; ES_M36[a][b]+=buf;  
  
  buf=0.;
  buf = J41*e_yy*J41 - (e_yx)*( J41*J42 + J42*J41 ) + J42*e_xx*J42 ;
  buf/=determinant; ES_M44[a][b]+=buf;
  buf=0.;
  buf = J41*e_yy*J51 - (e_yx)*( J41*J52 + J42*J51 ) + J42*e_xx*J52 ;
  buf/=determinant; ES_M45[a][b]+=buf;
  buf=0.;
  buf = J41*e_yy*J61 - (e_yx)*( J41*J62 + J42*J61 ) + J42*e_xx*J62 ;
  buf/=determinant; ES_M46[a][b]+=buf;  
  
  buf=0.;
  buf = J51*e_yy*J51 - (e_yx)*( J51*J52 + J52*J51 ) + J52*e_xx*J52 ;
  buf/=determinant; ES_M55[a][b]+=buf;
  buf=0.;
  buf = J51*e_yy*J61 - (e_yx)*( J51*J62 + J52*J61 ) + J52*e_xx*J62 ;
  buf/=determinant; ES_M56[a][b]+=buf;  
  
  buf=0.;
  buf = J61*e_yy*J61 - (e_yx)*( J61*J62 + J62*J61 ) + J62*e_xx*J62 ;
  buf/=determinant; ES_M66[a][b]+=buf;  
  
 }//end if (iz ip initialized )
}

void ESAlignTool::Cal_VectorP_wRotation(int iTrk,int iz,int ip,Double_t e_xx,Double_t e_yx,Double_t e_yy)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  Double_t XD=PredictionState_resiX[iTrk][a][b];
  Double_t YD=PredictionState_resiY[iTrk][a][b];
  
  Double_t Px = PredictionState_Px[iTrk][a][b];
  Double_t Py = PredictionState_Py[iTrk][a][b];
  Double_t Pz = PredictionState_Pz[iTrk][a][b];
  Double_t Xap = ES_Oap_X[a][b];
  Double_t Yap = ES_Oap_Y[a][b];
  Double_t Zap = ES_Oap_Z[a][b];
  Double_t Xpre = PredictionState_X[iTrk][a][b];
  Double_t Ypre = PredictionState_Y[iTrk][a][b];
  Double_t Zpre = PredictionState_Z[iTrk][a][b];
  Double_t dAlpha = ES_dAlpha[a][b];
  Double_t dBeta = ES_dBeta[a][b];
  Double_t dGamma = ES_dGamma[a][b];

  Double_t R11=ES_R11[a][b];
  Double_t R12=ES_R12[a][b];
  Double_t R13=ES_R13[a][b];
  Double_t R21=ES_R21[a][b];
  Double_t R22=ES_R22[a][b];
  Double_t R23=ES_R23[a][b];
  Double_t R31=ES_R31[a][b];
  Double_t R32=ES_R32[a][b];
  Double_t R33=ES_R33[a][b];
  
  Double_t DR11_Alpha = -(cos(dAlpha)*sin(dBeta)*sin(dGamma));
  Double_t DR11_Beta = -(cos(dGamma)*sin(dBeta)) - cos(dBeta)*sin(dAlpha)*sin(dGamma);
  Double_t DR11_Gamma = -(cos(dGamma)*sin(dAlpha)*sin(dBeta)) - cos(dBeta)*sin(dGamma);
  Double_t DR12_Alpha = cos(dAlpha)*cos(dGamma)*sin(dBeta);
  Double_t DR12_Beta = cos(dBeta)*cos(dGamma)*sin(dAlpha) - sin(dBeta)*sin(dGamma);
  Double_t DR12_Gamma = cos(dBeta)*cos(dGamma) - sin(dAlpha)*sin(dBeta)*sin(dGamma);
  Double_t DR13_Alpha = sin(dAlpha)*sin(dBeta);
  Double_t DR13_Beta = -(cos(dAlpha)*cos(dBeta));
  Double_t DR13_Gamma = 0.;
  Double_t DR21_Alpha = sin(dAlpha)*sin(dGamma);
  Double_t DR21_Beta = 0.;
  Double_t DR21_Gamma = -(cos(dAlpha)*cos(dGamma));
  Double_t DR22_Alpha = -(cos(dGamma)*sin(dAlpha));
  Double_t DR22_Beta = 0.;
  Double_t DR22_Gamma = -(cos(dAlpha)*sin(dGamma));
  Double_t DR23_Alpha = cos(dAlpha);
  Double_t DR23_Beta = 0.;
  Double_t DR23_Gamma = 0.;

  Double_t DR31_Alpha= sin(dAlpha)*sin(dGamma);
  Double_t DR31_Beta=0.;
  Double_t DR31_Gamma=-( cos(dAlpha)*cos(dGamma));
  Double_t DR32_Alpha=-( cos(dGamma)*sin(dAlpha));
  Double_t DR32_Beta=0.;
  Double_t DR32_Gamma=-( cos(dAlpha)*sin(dGamma));
  Double_t DR33_Alpha= cos(dAlpha);
  Double_t DR33_Beta=0.;
  Double_t DR33_Gamma=0.;

  Double_t PzLocal = R31*Px + R32*Py + R33*Pz;
  Double_t J11,J12,J21,J22,J31,J32,J41,J42,J51,J52,J61,J62;
    
  J11
  = -R11 - ((Px*R11 + Py*R12 + Pz*R13)*R31)/PzLocal;

  J12 
  = -R21 - ((Px*R21 + Py*R22 + Pz*R23)*R31)/PzLocal;

  J21 
  = -R12 - ((Px*R11 + Py*R12 + Pz*R13)*R32)/PzLocal;

  J22 
  = -R22 - ((Px*R21 + Py*R22 + Pz*R23)*R32)/PzLocal;

  J31 
  = -R13 - ((Px*R11 + Py*R12 + Pz*R13)*R33)/PzLocal;

  J32 
  = -R23 - ((Px*R21 + Py*R22 + Pz*R23)*R33)/PzLocal;

  J41 = DR11_Alpha*(Xpre-Xap) + DR12_Alpha*(Ypre-Yap) + DR13_Alpha*(Zpre-Zap);
  J42 = DR21_Alpha*(Xpre-Xap) + DR22_Alpha*(Ypre-Yap) + DR23_Alpha*(Zpre-Zap);
  J51 = DR11_Beta*(Xpre-Xap) + DR12_Beta*(Ypre-Yap) + DR13_Beta*(Zpre-Zap);
  J52 = DR21_Beta*(Xpre-Xap) + DR22_Beta*(Ypre-Yap) + DR23_Beta*(Zpre-Zap);
  J61 = DR11_Gamma*(Xpre-Xap) + DR12_Gamma*(Ypre-Yap) + DR13_Gamma*(Zpre-Zap);
  J62 = DR21_Gamma*(Xpre-Xap) + DR22_Gamma*(Ypre-Yap) + DR23_Gamma*(Zpre-Zap);

  Double_t determinant = e_xx*e_yy - e_yx*e_yx;
  //vp_i*determinant
  //  = J_i1*e_yy*XD - e_yx*(J_i1*YD+J_i2*XD) + J_i2*e_xx*YD;
  Double_t buf=0.; 
  buf = J11*e_yy*XD - e_yx*(J11*YD+J12*XD) + J12*e_xx*YD;
  buf/=determinant;  ES_P1[a][b]+=buf;
  buf=0.;
  buf = J21*e_yy*XD - e_yx*(J21*YD+J22*XD) + J22*e_xx*YD;
  buf/=determinant;  ES_P2[a][b]+=buf;
  buf=0.;
  buf = J31*e_yy*XD - e_yx*(J31*YD+J32*XD) + J32*e_xx*YD;
  buf/=determinant;  ES_P3[a][b]+=buf;
  buf=0.;
  buf = J41*e_yy*XD - e_yx*(J41*YD+J42*XD) + J42*e_xx*YD;
  buf/=determinant;  ES_P4[a][b]+=buf;
  buf=0.;
  buf = J51*e_yy*XD - e_yx*(J51*YD+J52*XD) + J52*e_xx*YD;
  buf/=determinant;  ES_P5[a][b]+=buf;
  buf=0.;
  buf = J61*e_yy*XD - e_yx*(J61*YD+J62*XD) + J62*e_xx*YD;
  buf/=determinant;  ES_P6[a][b]+=buf;
  
 }//end if (iz ip initialized )
}
*/

void ESAlignTool::Cal_MatrixM_wRotation(int iTrk,int iz,int ip,Double_t e_xx,Double_t e_yx,Double_t e_yy)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  Double_t Px = PredictionState_Px[iTrk][a][b];
  Double_t Py = PredictionState_Py[iTrk][a][b];
  Double_t Pz = PredictionState_Pz[iTrk][a][b];
  Double_t Xap = ES_Oap_X[a][b];
  Double_t Yap = ES_Oap_Y[a][b];
  Double_t Zap = ES_Oap_Z[a][b];
  Double_t Xpre = PredictionState_X[iTrk][a][b];
  Double_t Ypre = PredictionState_Y[iTrk][a][b];
  Double_t Zpre = PredictionState_Z[iTrk][a][b];
  Double_t dAlpha = ES_dAlpha[a][b];
  Double_t dBeta = ES_dBeta[a][b];
  Double_t dGamma = ES_dGamma[a][b];

  //Adding MatrixM
  Double_t buf=0.;
  buf=(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*(-((e_yx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));
  ES_M11[a][b] += (buf);

  buf=(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*(-((e_yx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));
  ES_M12[a][b] += (buf);

  buf=(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*(-((e_yx*(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));

  ES_M22[a][b] += (buf);

  buf=(-sin(dAlpha)+(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(cos(dAlpha)*sin(dBeta)+(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*(-((e_yx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));
  ES_M13[a][b] += (buf);

  buf=(-sin(dAlpha)+(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(cos(dAlpha)*sin(dBeta)+(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*(-((e_yx*(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));
  ES_M23[a][b] += (buf);

  buf=(-sin(dAlpha)-(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(-sin(dAlpha)-(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(cos(dAlpha)*sin(dBeta)-(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(cos(dAlpha)*sin(dBeta)-(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*(-((e_yx*(-sin(dAlpha)-(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(cos(dAlpha)*sin(dBeta)-(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));
  ES_M33[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+((((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M14[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(-(cos(dAlpha)*cos(dGamma))-((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)-((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+((((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(-(cos(dAlpha)*cos(dGamma))-((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)-((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M24[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(-sin(dAlpha)+(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(cos(dAlpha)*sin(dBeta)+(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+((((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(-sin(dAlpha)+(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(cos(dAlpha)*sin(dBeta)+(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M34[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*(((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+((((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*(((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M44[a][b] += (buf);

  buf=((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*((e_xx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M15[a][b] += (buf);

  buf=((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*((e_xx*(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M25[a][b] += (buf);

  buf=((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*((e_xx*(-sin(dAlpha)-(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(cos(dAlpha)*sin(dBeta)-(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(-sin(dAlpha)-(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(cos(dAlpha)*sin(dBeta)-(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M35[a][b] += (buf);

  buf=((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*((e_xx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*(((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*(((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M45[a][b] += (buf);

  buf=((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*((e_xx*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M55[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M16[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(-(cos(dAlpha)*cos(dGamma))-((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)-((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(-(cos(dAlpha)*cos(dGamma))-((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)-((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M26[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(-sin(dAlpha)+(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(cos(dAlpha)*sin(dBeta)+(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(-sin(dAlpha)+(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(cos(dAlpha)*sin(dBeta)+(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M36[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*(((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*(((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M46[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M56[a][b] += (buf);

  buf=(((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))*((e_xx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))*(-((e_yx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)));
  ES_M66[a][b] += (buf);
 }
}

void ESAlignTool::Cal_VectorP_wRotation(int iTrk,int iz,int ip,Double_t e_xx,Double_t e_yx,Double_t e_yy)
{
 int a,b;
 if( (iz==1||iz==-1) && (ip==1||ip==2) )
 {
  if(iz==-1) a=0;
  if(iz==1) a=1;
  b=ip-1;

  Double_t Px = PredictionState_Px[iTrk][a][b];
  Double_t Py = PredictionState_Py[iTrk][a][b];
  Double_t Pz = PredictionState_Pz[iTrk][a][b];
  Double_t Xap = ES_Oap_X[a][b];
  Double_t Yap = ES_Oap_Y[a][b];
  Double_t Zap = ES_Oap_Z[a][b];
  Double_t Xpre = PredictionState_X[iTrk][a][b];
  Double_t Ypre = PredictionState_Y[iTrk][a][b];
  Double_t Zpre = PredictionState_Z[iTrk][a][b];
  Double_t dAlpha = ES_dAlpha[a][b];
  Double_t dBeta = ES_dBeta[a][b];
  Double_t dGamma = ES_dGamma[a][b];
  Double_t XD
=PredictionState_resiX[iTrk][a][b];
  Double_t YD
=PredictionState_resiY[iTrk][a][b];

  Double_t buf=0.;
  buf=YD*((e_xx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+XD*(-((e_yx*(cos(dAlpha)*sin(dGamma)+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dBeta)*cos(dGamma))+sin(dAlpha)*sin(dBeta)*sin(dGamma)+((cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));
  ES_P1[a][b] += (buf) ;

  buf=YD*((e_xx*(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+XD*(-((e_yx*(-(cos(dAlpha)*cos(dGamma))+((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma)+((-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));
  ES_P2[a][b] += (buf) ;

  buf=YD*((e_xx*(-sin(dAlpha)+(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy)-(e_yx*(cos(dAlpha)*sin(dBeta)+(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+XD*(-((e_yx*(-sin(dAlpha)+(cos(dAlpha)*cos(dBeta)*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma)))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy))+(e_yy*(cos(dAlpha)*sin(dBeta)+(cos(dAlpha)*cos(dBeta)*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/(-e_yx*e_yx+e_xx*e_yy));
  ES_P3[a][b] += (buf) ;

  buf=YD*((e_xx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*(((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+XD*(-((e_yx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))+((-Zap+Zpre)*cos(dAlpha)-(-Yap+Ypre)*cos(dGamma)*sin(dAlpha))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Xap+Xpre)*sin(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*(((-Yap+Ypre)*cos(dAlpha)*cos(dGamma)*sin(dBeta)+(-Zap+Zpre)*sin(dAlpha)*sin(dBeta))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*sin(dBeta)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Yap*cos(dAlpha)*cos(dBeta)*cos(dGamma))+Ypre*cos(dAlpha)*cos(dBeta)*cos(dGamma)-Zap*cos(dBeta)*sin(dAlpha)+Zpre*cos(dBeta)*sin(dAlpha)+(Xap-Xpre)*cos(dAlpha)*cos(dBeta)*sin(dGamma))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))));
  ES_P4[a][b] += (buf) ;

  buf=YD*((e_xx*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+XD*(-((e_yx*(Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*((-Xap+Xpre)*(-(cos(dGamma)*sin(dBeta))-cos(dBeta)*sin(dAlpha)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-((-Zap+Zpre)*cos(dAlpha)*cos(dBeta))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))*(-(Zap*cos(dAlpha)*sin(dBeta))+Zpre*cos(dAlpha)*sin(dBeta)+Yap*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))));
  ES_P5[a][b] += (buf) ;

  buf=YD*((e_xx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma))))-(e_yx*((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+XD*(-((e_yx*((Py*cos(dAlpha)*cos(dGamma)+Pz*sin(dAlpha)-Px*cos(dAlpha)*sin(dGamma))*(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))-(-Xap+Xpre)*cos(dAlpha)*cos(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))-(-Yap+Ypre)*cos(dAlpha)*sin(dGamma)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))))+(e_yy*((-Xap+Xpre)*(-(cos(dGamma)*sin(dAlpha)*sin(dBeta))-cos(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(-Yap+Ypre)*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma))*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))+(Yap*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))-Ypre*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+(Xap-Xpre)*(cos(dBeta)*cos(dGamma)*sin(dAlpha)-sin(dBeta)*sin(dGamma)))*(-(Pz*cos(dAlpha)*sin(dBeta))+Py*(cos(dGamma)*sin(dAlpha)*sin(dBeta)+cos(dBeta)*sin(dGamma))+Px*(cos(dBeta)*cos(dGamma)-sin(dAlpha)*sin(dBeta)*sin(dGamma)))))/((-e_yx*e_yx+e_xx*e_yy)*(Pz*cos(dAlpha)*cos(dBeta)+Px*(cos(dGamma)*sin(dBeta)+cos(dBeta)*sin(dAlpha)*sin(dGamma))+Py*(-(cos(dBeta)*cos(dGamma)*sin(dAlpha))+sin(dBeta)*sin(dGamma)))));
  ES_P6[a][b] += (buf) ;
 }
}

