#include "TCanvas.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLine.h"
#include "TNtupleD.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TFile.h"

void readMap(Char_t* filename);

const Double_t c   = 0.299792458; // speed of light in vacuum in m/ns

/*
 *    Simulated beam parameters
 */
  Double_t mu_x   = 0;
  Double_t mu_xp  = 0;
  Double_t mu_y   = 0;
  Double_t mu_yp  = 0;
  Double_t mu_l   = 0;
  Double_t mu_dkk = 0;
  Double_t sig_x   = 0*1e-2;
  Double_t sig_xp  = 0*1e-4;
  Double_t sig_y   = 0*1e-2;
  Double_t sig_yp  = 0*1e-4;
  Double_t sig_l   = 0*10;
  Double_t sig_dkk = 0*1e-2;


static const UInt_t maxElements = 2000;
static Double_t mapx  [maxElements][7]; UInt_t nx  = 0;
static Double_t mapxp [maxElements][7]; UInt_t nxp = 0;
static Double_t mapy  [maxElements][7]; UInt_t ny  = 0;
static Double_t mapyp [maxElements][7]; UInt_t nyp = 0;
static Double_t mapl  [maxElements][7]; UInt_t nll = 0;
static Double_t mapd  [maxElements][7]; UInt_t nd  = 0;
static UInt_t ns[3][3]={{0,0,0},{0,0,0},{0,0,0}};
static Double_t sxsx[maxElements][7];
static Double_t sxsy[maxElements][7];
static Double_t sxsz[maxElements][7];
static Double_t sysx[maxElements][7];
static Double_t sysy[maxElements][7];
static Double_t sysz[maxElements][7];
static Double_t szsx[maxElements][7];
static Double_t szsy[maxElements][7];
static Double_t szsz[maxElements][7];
static UInt_t maxOrder;
static const UInt_t maxMaxOrder = 10;
static Double_t M0,p0,K0,gamma0,beta0,G,agamma,L0;

Double_t calcTime(Double_t* vs);

void readMap(Char_t* filename)
{
  FILE* file = fopen(filename,"r");
  if (!file)
  {
    fprintf(stderr,"ERROR: cannot open file '%s'\n",filename);
    return;
  }
  printf("SUCCESS: opened file '%s'\n",filename);

  Char_t line[200];

  // extract particle information
  fgets(line,200,file);
  printf(">>%s<<\n",line);
  if (sscanf(line,"Muon Mass =   %lf MeV/c^2",&M0) != 1) return;
  printf("Particle mass      = %lf MeV/c^2\n",M0);
  fgets(line,200,file);
  if (sscanf(line,"Muon Momentum =   %lf MeV/c",&p0) != 1) return;
  printf("Reference momentum = %lf MeV/c\n",p0);
  fgets(line,200,file);
  if (sscanf(line,"Muon Kinetic Energy =   %lf MeV",&K0) != 1) return;
  printf("Reference Kin. E.  = %lf MeV\n",K0);
  fgets(line,200,file);
  if (sscanf(line,"Muon gamma =   %lf",&gamma0) != 1) return;
  printf("Reference gamma    = %lf \n",gamma0);
  fgets(line,200,file);
  if (sscanf(line,"Muon beta =  %lf",&beta0) != 1) return;
  printf("Reference velocity = %lf * c\n",beta0);
  fgets(line,200,file);
  if (sscanf(line,"Muon Anomaly G = %lf",&G) != 1) return;
  printf("Spin anomaly       = %lf\n",G);
  fgets(line,200,file);
  if (sscanf(line,"Muon Spin Tune G.gamma = %lf",&agamma) != 1) return;
  printf("Spin tune          = %lf\n",agamma);

  // extract length of reference orbit
  fgets(line,200,file);

  if (!sscanf(line," L %lf",&L0))
  {
    fprintf(stderr,"ERROR: wrong format\n");
    return;
  }
  printf("Reference orbit    = %f m\n",L0);

  // read `P', whatever it is ....
  fgets(line,200,file);
  if (line[1] != 'P')
  {
    fprintf(stderr,"ERROR: wrong format\n");
    return;
  }

  // read `A', whatever it is ....
  fgets(line,200,file);
  if (line[1] != 'A')
  {
    fprintf(stderr,"ERROR: wrong format\n");
    return;
  }

  // now read coordinate map data
  UInt_t dum,order,px,pa,py,pb,pt,pd;
  Double_t coef;

  // reset size counters
  nx  = 0;
  nxp = 0;
  ny  = 0;
  nyp = 0;
  nll = 0;
  nd  = 0;

  ns[0][0] = 0;
  ns[0][1] = 0;
  ns[0][2] = 0;
  ns[1][0] = 0;
  ns[1][1] = 0;
  ns[1][2] = 0;
  ns[2][0] = 0;
  ns[2][1] = 0;
  ns[2][2] = 0;
  maxOrder = 0;

  // dependences of x
  do
  {
    fgets(line,200,file);
    UInt_t nm = sscanf(line, "%d %lf %d %d %d %d %d %d %d",
      &dum,&coef,&order,&px,&pa,&py,&pb,&pt,&pd);
    if (nm!=9) continue;
    if (order>maxOrder) maxOrder = order;
    mapx[nx][0] = coef;
    mapx[nx][1] = px;
    mapx[nx][2] = pa;
    mapx[nx][3] = py;
    mapx[nx][4] = pb;
    mapx[nx][5] = pt;
    mapx[nx][6] = pd;
    nx++;
  } while (!strstr(line,"-------------"));


  // dependences of x'
  do
  {
    fgets(line,200,file);
    UInt_t nm = sscanf(line, "%d %lf %d %d %d %d %d %d %d",
      &dum,&coef,&order,&px,&pa,&py,&pb,&pt,&pd);
    if (nm!=9) continue;
    if (order>maxOrder) maxOrder = order;
    mapxp[nxp][0] = coef;
    mapxp[nxp][1] = px;
    mapxp[nxp][2] = pa;
    mapxp[nxp][3] = py;
    mapxp[nxp][4] = pb;
    mapxp[nxp][5] = pt;
    mapxp[nxp][6] = pd;
    nxp++;
  } while (!strstr(line,"-------------"));

  // dependences of y
  do
  {
    fgets(line,200,file);
    UInt_t nm = sscanf(line, "%d %lf %d %d %d %d %d %d %d",
      &dum,&coef,&order,&px,&pa,&py,&pb,&pt,&pd);
    if (nm!=9) continue;
    if (order>maxOrder) maxOrder = order;
    mapy[ny][0] = coef;
    mapy[ny][1] = px;
    mapy[ny][2] = pa;
    mapy[ny][3] = py;
    mapy[ny][4] = pb;
    mapy[ny][5] = pt;
    mapy[ny][6] = pd;
    ny++;
  } while (!strstr(line,"-------------"));

  // dependences of y'
  do
  {
    fgets(line,200,file);
    UInt_t nm = sscanf(line, "%d %lf %d %d %d %d %d %d %d",
      &dum,&coef,&order,&px,&pa,&py,&pb,&pt,&pd);
    if (nm!=9) continue;
    if (order>maxOrder) maxOrder = order;
    mapyp[nyp][0] = coef;
    mapyp[nyp][1] = px;
    mapyp[nyp][2] = pa;
    mapyp[nyp][3] = py;
    mapyp[nyp][4] = pb;
    mapyp[nyp][5] = pt;
    mapyp[nyp][6] = pd;
    nyp++;
  } while (!strstr(line,"-------------"));

  // dependences of dl
  do
  {
    fgets(line,200,file);
    UInt_t nm = sscanf(line, "%d %lf %d %d %d %d %d %d %d",
      &dum,&coef,&order,&px,&pa,&py,&pb,&pt,&pd);
    if (nm!=9) continue;
    if (order>maxOrder) maxOrder = order;
    mapl[nll][0] = coef;
    mapl[nll][1] = px;
    mapl[nll][2] = pa;
    mapl[nll][3] = py;
    mapl[nll][4] = pb;
    mapl[nll][5] = pt;
    mapl[nll][6] = pd;
    nll++;
  } while (!strstr(line,"-------------"));

  // dependences of dK/K
  do
  {
    fgets(line,200,file);
    UInt_t nm = sscanf(line, "%d %lf %d %d %d %d %d %d %d",
      &dum,&coef,&order,&px,&pa,&py,&pb,&pt,&pd);
    if (nm!=9) continue;
    if (order>maxOrder) maxOrder = order;
    mapd[nd][0] = coef;
    mapd[nd][1] = px;
    mapd[nd][2] = pa;
    mapd[nd][3] = py;
    mapd[nd][4] = pb;
    mapd[nd][5] = pt;
    mapd[nd][6] = pd;
    nd++;
  } while (!strstr(line,"-------------"));

  // skip the lines for the mass and charge dependence, assuming they are all zero!
  fgets(line,200,file);
  fgets(line,200,file);
  fgets(line,200,file);
  fgets(line,200,file);

  // spin <sx|sx> <sx|sy> <sx|sz>
  Double_t c1,c2,c3;
  UInt_t pow;
  do
  {
    fgets(line,200,file);
    if (!sscanf(line,"%lf %lf %lf %d", &c1, &c2, &c3, &pow)) continue;
    UInt_t ppow = pow;
    if (c1!=0)
    {
      sxsx[ns[0][0]][0] = c1;
      for (UInt_t l = 6; l>0; l--) 
      {
        sxsx[ns[0][0]][l] = pow%10;
        pow/=10;
      }
      (ns[0][0])++;
    }
    pow=ppow;
    if (c2!=0)
    {
      sxsy[ns[0][1]][0] = c2;
      for (UInt_t l = 6; l>0; l--)
      {
        sxsy[ns[0][1]][l] = pow%10;
        pow/=10;
      }
      (ns[0][1])++;
    }
    pow=ppow;
    if (c3!=0)
    { 
      sxsz[ns[0][2]][0] = c3;
      for (UInt_t l = 6; l>0; l--)
      {
        sxsz[ns[0][2]][l] = pow%10;
        pow/=10;
      }
      (ns[0][2])++;
    }
  } while (!strstr(line,"-------------"));
  // spin <sy|sx> <sy|sy> <sy|sz>
  do
  {
    fgets(line,200,file);
    if (!sscanf(line,"%lf %lf %lf %d", &c1, &c2, &c3, &pow)) continue;
    UInt_t ppow=pow;
    if (c1!=0)
    {
      sysx[ns[1][0]][0] = c1;
      for (UInt_t l = 6; l>0; l--)
      {
        sysx[ns[1][0]][l] = pow%10;
        pow/=10;
      }
      (ns[1][0])++;
    }
    pow=ppow;
    if (c2!=0)
    {
      sysy[ns[1][1]][0] = c2;
      for (UInt_t l = 6; l>0; l--)
      {
        sysy[ns[1][1]][l] = pow%10;
        pow/=10;
      }
      (ns[1][1])++;
    }
    pow=ppow;
    if (c3!=0)
    { 
      sysz[ns[1][2]][0] = c3;
      for (UInt_t l = 6; l>0; l--)
      {
        sysz[ns[1][2]][l] = pow%10;
        pow/=10;
      }
      (ns[1][2])++;
    }
  } while (!strstr(line,"-------------"));
  // spin <sz|sx> <sz|sy> <sz|sz>
  do
  {
    fgets(line,200,file);
    if (!sscanf(line,"%lf %lf %lf %d", &c1, &c2, &c3, &pow)) continue;
    UInt_t ppow=pow;
    if (c1!=0)
    {
      szsx[ns[2][0]][0] = c1;
      for (UInt_t l = 6; l>0; l--)
      {
        szsx[ns[2][0]][l] = pow%10;
        pow/=10;
      }
      (ns[2][0])++;
    }
    pow=ppow;
    if (c2!=0)
    {
      szsy[ns[2][1]][0] = c2;
      for (UInt_t l = 6; l>0; l--)
      {
        szsy[ns[2][1]][l] = pow%10;
        pow/=10;
      }
      (ns[2][1])++;
    }
    pow=ppow;
    if (c3!=0)
    { 
      szsz[ns[2][2]][0] = c3;
      for (UInt_t l = 6; l>0; l--)
      {
        szsz[ns[2][2]][l] = pow%10;
        pow/=10;
      }
      (ns[2][2])++;
    }
  } while (!strstr(line,"-------------"));

  printf("------------------------------------------\n");
  printf("Nr. of Beam Dynamics Elements read\n");
  printf("%d %d %d %d %d %d\n",nx,nxp,ny,nyp,nll,nd);
  printf("%d %d %d\n",ns[0][0],ns[0][1],ns[0][2]);
  printf("%d %d %d\n",ns[1][0],ns[1][1],ns[1][2]);
  printf("%d %d %d\n",ns[2][0],ns[2][1],ns[2][2]);
  printf("------------------------------------------\n");
  printf("Maximum Order = %d\n",maxOrder);

  fclose(file);
}

void printMap()
{
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<nx; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  mapx[i][0],
  (UInt_t)mapx[i][1],
  (UInt_t)mapx[i][2],
  (UInt_t)mapx[i][3],
  (UInt_t)mapx[i][4],
  (UInt_t)mapx[i][5],
  (UInt_t)mapx[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<nxp; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  mapxp[i][0],
  (UInt_t)mapxp[i][1],
  (UInt_t)mapxp[i][2],
  (UInt_t)mapxp[i][3],
  (UInt_t)mapxp[i][4],
  (UInt_t)mapxp[i][5],
  (UInt_t)mapxp[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ny; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  mapy[i][0],
  (UInt_t)mapy[i][1],
  (UInt_t)mapy[i][2],
  (UInt_t)mapy[i][3],
  (UInt_t)mapy[i][4],
  (UInt_t)mapy[i][5],
  (UInt_t)mapy[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<nyp; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  mapyp[i][0],
  (UInt_t)mapyp[i][1],
  (UInt_t)mapyp[i][2],
  (UInt_t)mapyp[i][3],
  (UInt_t)mapyp[i][4],
  (UInt_t)mapyp[i][5],
  (UInt_t)mapyp[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<nll; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  mapl[i][0],
  (UInt_t)mapl[i][1],
  (UInt_t)mapl[i][2],
  (UInt_t)mapl[i][3],
  (UInt_t)mapl[i][4],
  (UInt_t)mapl[i][5],
  (UInt_t)mapl[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<nd; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  mapd[i][0],
  (UInt_t)mapd[i][1],
  (UInt_t)mapd[i][2],
  (UInt_t)mapd[i][3],
  (UInt_t)mapd[i][4],
  (UInt_t)mapd[i][5],
  (UInt_t)mapd[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[0][0]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  sxsx[i][0],
  (UInt_t)sxsx[i][1],
  (UInt_t)sxsx[i][2],
  (UInt_t)sxsx[i][3],
  (UInt_t)sxsx[i][4],
  (UInt_t)sxsx[i][5],
  (UInt_t)sxsx[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[0][1]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  sxsy[i][0],
  (UInt_t)sxsy[i][1],
  (UInt_t)sxsy[i][2],
  (UInt_t)sxsy[i][3],
  (UInt_t)sxsy[i][4],
  (UInt_t)sxsy[i][5],
  (UInt_t)sxsy[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[0][2]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  sxsz[i][0],
  (UInt_t)sxsz[i][1],
  (UInt_t)sxsz[i][2],
  (UInt_t)sxsz[i][3],
  (UInt_t)sxsz[i][4],
  (UInt_t)sxsz[i][5],
  (UInt_t)sxsz[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[1][0]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  sysx[i][0],
  (UInt_t)sysx[i][1],
  (UInt_t)sysx[i][2],
  (UInt_t)sysx[i][3],
  (UInt_t)sysx[i][4],
  (UInt_t)sysx[i][5],
  (UInt_t)sysx[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[1][1]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  sysy[i][0],
  (UInt_t)sysy[i][1],
  (UInt_t)sysy[i][2],
  (UInt_t)sysy[i][3],
  (UInt_t)sysy[i][4],
  (UInt_t)sysy[i][5],
  (UInt_t)sysy[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[1][2]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  sysz[i][0],
  (UInt_t)sysz[i][1],
  (UInt_t)sysz[i][2],
  (UInt_t)sysz[i][3],
  (UInt_t)sysz[i][4],
  (UInt_t)sysz[i][5],
  (UInt_t)sysz[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[2][0]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  szsx[i][0],
  (UInt_t)szsx[i][1],
  (UInt_t)szsx[i][2],
  (UInt_t)szsx[i][3],
  (UInt_t)szsx[i][4],
  (UInt_t)szsx[i][5],
  (UInt_t)szsx[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[2][1]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  szsy[i][0],
  (UInt_t)szsy[i][1],
  (UInt_t)szsy[i][2],
  (UInt_t)szsy[i][3],
  (UInt_t)szsy[i][4],
  (UInt_t)szsy[i][5],
  (UInt_t)szsy[i][6]);
  printf("---------------------------------------------\n");
  for (UInt_t i = 0; i<ns[2][2]; i++)
  printf("%25.15e %d %d  %d %d  %d %d\n",
  szsz[i][0],
  (UInt_t)szsz[i][1],
  (UInt_t)szsz[i][2],
  (UInt_t)szsz[i][3],
  (UInt_t)szsz[i][4],
  (UInt_t)szsz[i][5],
  (UInt_t)szsz[i][6]);
  printf("---------------------------------------------\n");
}

Double_t xn[6][maxMaxOrder];
void applyMap(Double_t* v, Double_t* s, Bool_t withSpin = kTRUE)
{
  Double_t w[6] = {0,0,0,0,0,0}; // new coordinates
  Double_t ss[3] = {0,0,0}; // new spin
  Double_t thisTerm;

  for (UInt_t i = 0 ; i<6; i++) // loop over six coordinates
  {
    xn[i][0] = 1;  // v[i]^0
    for (UInt_t j = 1; j<=maxOrder; j++) xn[i][j] = xn[i][j-1]*v[i];
  }

  //  x
  for (UInt_t i = 0; i<nx; i++)
  {
    thisTerm = mapx[i][0];
    for (UInt_t c = 0; c<6; c++) // loop over all six coordinates
    {
      thisTerm*= xn[c][(UInt_t)(mapx[i][c+1])];
    }
    w[0] += thisTerm;
  }

  // xp
  for (UInt_t i = 0; i<nxp; i++)
  {
    thisTerm = mapxp[i][0];
    for (UInt_t c = 0; c<6; c++)
    {
      thisTerm*= xn[c][(UInt_t)(mapxp[i][c+1])];
    }
    w[1] += thisTerm;
  }

  // y
  for (UInt_t i = 0; i<ny; i++)
  {
    thisTerm = mapy[i][0];
    for (UInt_t c = 0; c<6; c++)
    {
      thisTerm*= xn[c][(UInt_t)(mapy[i][c+1])];
    }
    w[2] += thisTerm;
  }

  // yp
  for (UInt_t i = 0; i<nyp; i++)
  {
    thisTerm = mapyp[i][0];
    for (UInt_t c = 0; c<6; c++)
    {
      thisTerm*= xn[c][(UInt_t)(mapyp[i][c+1])];
    }
    w[3] += thisTerm;
  }

  // dl
  for (UInt_t i = 0; i<nll; i++)
  {
    thisTerm = mapl[i][0];
    for (UInt_t c = 0; c<6; c++)
    {
      thisTerm*= xn[c][(UInt_t)(mapl[i][c+1])];
    }
    w[4] += thisTerm;
  }

  // dkk
  for (UInt_t i = 0; i<nd; i++)
  {
    thisTerm = mapd[i][0];
    for (UInt_t c = 0; c<6; c++)
    {
      thisTerm*= xn[c][(UInt_t)(mapd[i][c+1])];
    }
    w[5] += thisTerm;
  }

  Double_t msxsx, msxsy, msxsz, msysx, msysy, msysz, mszsx, mszsy, mszsz; // spin rotation matrix
  msxsx=msxsy=msxsz=msysx=msysy=msysz=mszsx=mszsy=mszsz=0;  // init. to zero

  if (withSpin)
  {
    // spin matrix
    msxsx = 0;
    for (UInt_t i = 0; i<ns[0][0]; i++)
    {
      thisTerm = sxsx[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(sxsx[i][c+1])];
      }
      msxsx += thisTerm;
    }

    msxsy = 0;
    for (UInt_t i = 0; i<ns[0][1]; i++)
    {
      thisTerm = sxsy[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(sxsy[i][c+1])];
      }
      msxsy += thisTerm;
    }

    msxsz = 0;
    for (UInt_t i = 0; i<ns[0][2]; i++)
    {
      thisTerm = sxsz[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(sxsz[i][c+1])];
      }
      msxsz += thisTerm;
    }

    msysx = 0;
    for (UInt_t i = 0; i<ns[1][0]; i++)
    {
      thisTerm = sysx[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(sysx[i][c+1])];
      }
      msysx += thisTerm;
    }

    msysy = 0;
    for (UInt_t i = 0; i<ns[1][1]; i++)
    {
      thisTerm = sysy[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(sysy[i][c+1])];
      }
      msysy += thisTerm;
    }

    msysz = 0;
    for (UInt_t i = 0; i<ns[1][2]; i++)
    {
      thisTerm = sysz[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(sysz[i][c+1])];
      }
      msysz += thisTerm;
    }

    mszsx = 0;
    for (UInt_t i = 0; i<ns[2][0]; i++)
    {
      thisTerm = szsx[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(szsx[i][c+1])];
      }
      mszsx += thisTerm;
    }

    mszsy = 0;
    for (UInt_t i = 0; i<ns[2][1]; i++)
    {
      thisTerm = szsy[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(szsy[i][c+1])];
      }
      mszsy += thisTerm;
    }

    mszsz = 0;
    for (UInt_t i = 0; i<ns[2][2]; i++)
    {
      thisTerm = szsz[i][0];
      for (UInt_t c = 0; c<6; c++)
      {
        thisTerm*= xn[c][(UInt_t)(szsz[i][c+1])];
      }
      mszsz += thisTerm;
    }
  }

  // copy new coordinates
  for (UInt_t i = 0; i<6; i++) v[i]=w[i];

  if (withSpin)
  {
    ss[0] = msxsx*s[0] + msxsy*s[1] + msxsz*s[2] ;
    ss[1] = msysx*s[0] + msysy*s[1] + msysz*s[2] ;
    ss[2] = mszsx*s[0] + mszsy*s[1] + mszsz*s[2] ;
    Double_t sl = 1 ;// sqrt(ss[0]*ss[0]+ss[1]*ss[1]+ss[2]*ss[2]);

    s[0] = ss[0]/sl;
    s[1] = ss[1]/sl;
    s[2] = ss[2]/sl;
  }
}

Bool_t tooFarOut(Double_t* vs)
{
  Double_t r = sqrt(vs[0]*vs[0]+vs[2]*vs[2]);
  if (r>0.045) return kTRUE; // circular aperture
  return kFALSE;
}

void track(UInt_t maxturn = 7000, UInt_t npart = 1, Char_t* name = "thirdOrderMap.dat", Bool_t withSpin=kTRUE, UInt_t
saveOneIn=1, UInt_t marker=0)
{
  if (gRandom) gRandom->Delete();
  gRandom = new TRandom3(0);

  if (name) readMap(name);
  // components 1-6: beam     7-9: spin    10: turn  11: absolute time  12: particle id
  // +12 : initial values
  Double_t vs[25] = {0,0,0,0,0,0,
                     0,0,0,
                     0,0,0};
  vs[24] = marker;

  TNtupleD* nt = (TNtupleD*)gDirectory->Get("nt");
  if (!nt)
  {
    nt = new TNtupleD("nt","tracking","x:xp:y:yp:dl:dkk:sx:sy:sz:i:t:id:x0:xp0:y0:yp0:dl0:dkk0:sx0:sy0:sz0:i0:t0:id0:marker");
    nt->SetMarkerStyle(20);
    nt->SetMarkerSize(0.5);
    nt->SetMarkerColor(2);
  }

  UInt_t turn = 0;
  UInt_t particle = 0;
  ULong_t particleTurns = 0;
  TStopwatch* st = new TStopwatch();
  st->Start();
  for (particle = 0; particle<npart; particle++)
  {
    if (!(particle%100)) printf("%d of %d\n",particle,npart);
    vs[0] = gRandom->Gaus(mu_x,sig_x);
    vs[1] = gRandom->Gaus(mu_xp,sig_xp);
    vs[2] = gRandom->Gaus(mu_y,sig_y);
    vs[3] = gRandom->Gaus(mu_yp,sig_yp);
    vs[4] = gRandom->Gaus(mu_l,sig_l);
    vs[5] = gRandom->Gaus(mu_dkk,sig_dkk);
    vs[6] = 0;
    vs[7] = 0;
    vs[8] = 1;
    vs[9] = 0; // turn number
    vs[10] = 0; // time
    vs[11] = particle;

    for (UInt_t i = 0; i<12; i++)
    {
      vs[12+i] = vs[i];
    }

    for (turn = 0; turn<maxturn; turn++)
    {
      vs[9] = turn;
      vs[10] = calcTime(vs);
      if ( (turn%saveOneIn)==0) nt->Fill(vs);
      applyMap(vs,vs+6,withSpin);
      if (tooFarOut(vs)) break;
      particleTurns++;
    }
    vs[9] = turn;
    vs[10] = calcTime(vs);
    if (turn>0) nt->Fill(vs);
  }
  st->Stop();
  st->Print();
  printf("Time per particle.turn = %f mus\n",st->RealTime()/particleTurns*1e6);
  printf("Requested # particles = %d\n",npart);
  printf("Requested # turns     = %d\n",maxturn);
  printf("Realized # part.turns = %ld (%.0f percent)\n",particleTurns,(particleTurns*100.)/(npart*maxturn));
}

Double_t calcTime(Double_t* vs)
{
  Double_t gamma = (1+vs[5])*K0/M0+1;
  return -(gamma+1)/gamma/(beta0*c)*vs[4] + vs[9]*L0/(beta0*c);
  // return vs[9]*L0/(beta0*c);
}

void doit(UInt_t oneIn = 1, UInt_t nturn = 100001, UInt_t npart = 100)
{
  TFile* outfile = new TFile("fullring.root","recreate");
  Char_t inname[100];
  sprintf(inname,"ring.map");

  mu_x   = 0;
  mu_xp  = 0;
  mu_y   = 0;
  mu_yp  = 0;
  mu_l   = 0;
  mu_dkk = 0;

  sig_x   = 0*0.02;
  sig_xp  = 0*0.002;
  sig_y   = 0*0.02;
  sig_yp  = 0.02;
  sig_l   = 0*7.2;
  sig_dkk = 0*1e-3;

  track(nturn,npart,inname,kTRUE,oneIn,1);

  outfile->Write();
  outfile->Close();
  outfile = new TFile("fullring.root");
}
