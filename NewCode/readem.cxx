#include "TNtupleD.h"
#include "TRandom.h"
#include "TTimeStamp.h"
#include "TFile.h"


Double_t calcTime(Double_t* vs, UInt_t turn);
static Double_t M0,p0,K0,gamma0,beta0,G,agamma,L0;
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

struct davar
{
  Double_t amp; // amplitude
  Int_t    p[6]; // powers
  davar*   nextElement; // pointer to next element

  // constructors
  davar() // default
  {
    amp = 0;
    memset(p,0,sizeof(p));
    nextElement = NULL;
  }

  davar(Double_t v, UInt_t n0, UInt_t n1, UInt_t n2, UInt_t n3, UInt_t n4, UInt_t n5)
  {
    amp = v;
    p[0] = n0;
    p[1] = n1;
    p[2] = n2;
    p[3] = n3;
    p[4] = n4;
    p[5] = n5;
    nextElement = NULL;
  }

  // destructor
  ~davar()
  {
    /*
    if (nextElement)
    {
      delete nextElement;
      nextElement = NULL;
    }
    */
  }

  void add(davar* u)
  {
    if (nextElement == NULL)
    {
      nextElement = u;
    }
    else
    {
      (*nextElement).add(u);
    }
  }

  Double_t eval(Double_t* x)
  {
    Double_t res = amp;
    for (UInt_t i = 0; i<6; i++)
    {
      for (UInt_t u = 0; u<p[i]; u++) res *= x[i];
    }
    if (nextElement) res += (*nextElement).eval(x);
    return res;
  }

  void print()
  {
    printf("DAV : %15.7e",amp);
    for (UInt_t i = 0 ; i<6; i++)
    {
      printf(" : %2d",p[i]);
    }
    printf("\n");
    if (nextElement) (*nextElement).print();
  }
};

struct damap
{
  davar map[6];

  void apply(Double_t* xin)
  {
    Double_t xout[6];
    for (UInt_t i = 0; i<6; i++)
    {
      xout[i] = map[i].eval(xin);
    }
    for (UInt_t i = 0; i<6; i++)
    {
      xin[i] = xout[i];
    }
  }

  void print()
  {
    printf("--- X ---\n");
    map[0].print();
    printf("--- XP --\n");
    map[1].print();
    printf("--- Y ---\n");
    map[2].print();
    printf("--- YP --\n");
    map[3].print();
    printf("--- T ---\n");
    map[4].print();
    printf("--- dKK -\n");
    map[5].print();
  }
};

struct smap
{
  davar map[3][3];

  void print()
  {
    for (UInt_t i = 0; i<3; i++)
    {
      for (UInt_t j = 0; j<3; j++)
      {
        map[i][j].print();
      }
    }
  }

  void apply(Double_t* xin, Double_t* spinin)
  {
    Double_t matrix[3][3];
    for (UInt_t i = 0; i<3; i++)
    {
      for (UInt_t j = 0; j<3; j++)
      {
        matrix[i][j] = map[i][j].eval(xin);
      }
    }
    Double_t spinout[3];
    for (UInt_t i = 0; i<3; i++)
    {
      spinout[i] = 0;
      for (UInt_t j = 0; j<3; j++)
      {
        spinout[i] += matrix[i][j]*spinin[j];
      }
    }
    Double_t sum = 0;
    for (UInt_t i = 0; i<3; i++)
    {
      sum += spinout[i]*spinout[i];
    }
    for (UInt_t i = 0; i<3; i++)
    {
      spinin[i] = spinout[i]/sqrt(sum);
    }
  }
};

void readMap(Char_t* filename, damap* map, smap* smap)
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

  printf("Revolution Period  = %lf ns\n", L0/(beta0*c));
  printf("Precession Period  = %lf ns\n", L0/(beta0*c)/agamma);

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

  // map data
  UInt_t dum,order,px,pa,py,pb,pt,pd;
  Double_t coef;

  // dependences of 6D spatial coordinates
  for (UInt_t dim = 0; dim<6; dim++)
  {
    printf("doing %d ....",dim);
    UInt_t nl = 0;
    do
    {
      fgets(line,200,file);
      UInt_t nm = sscanf(line, "%d %lf %d %d %d %d %d %d %d",
        &dum,&coef,&order,&px,&pa,&py,&pb,&pt,&pd);
      if (nm!=9) continue;
      davar* cc = new davar(coef,px,pa,py,pb,pt,pd);
      (*map).map[dim].add(cc);
      nl++;
    } while (!strstr(line,"-------------"));
    printf(" done! (%d lines)\n",nl);
  }

  // skip the lines for the mass and charge dependence, assuming they are all zero!
  fgets(line,200,file);
  fgets(line,200,file);
  fgets(line,200,file);
  fgets(line,200,file);

  // dependence on spin matrix
  Double_t c1,c2,c3;
  UInt_t pow;
  UInt_t pww[6];
  for (UInt_t i = 0; i<3; i++)
  {
    do
    {
      fgets(line,200,file);
      if (!sscanf(line,"%lf %lf %lf %d", &c1, &c2, &c3, &pow)) continue;

      UInt_t ppow = pow;
      if (c1!=0)
      {
        for (UInt_t l = 6; l>0; l--) 
        {
          pww[l-1] = pow%10;
          pow/=10;
        }
        davar* cc = new davar(c1,pww[0],pww[1],pww[2],pww[3],pww[4],pww[5]);
        (*smap).map[i][0].add(cc);
      }

      pow=ppow;
      if (c2!=0)
      {
        for (UInt_t l = 6; l>0; l--)
        {
          pww[l-1] = pow%10;
          pow/=10;
        }
        davar* cc = new davar(c2,pww[0],pww[1],pww[2],pww[3],pww[4],pww[5]);
        (*smap).map[i][1].add(cc);
      }

      pow=ppow;
      if (c3!=0)
      { 
        for (UInt_t l = 6; l>0; l--)
        {
          pww[l-1] = pow%10;
          pow/=10;
        }
        davar* cc = new davar(c3,pww[0],pww[1],pww[2],pww[3],pww[4],pww[5]);
        (*smap).map[i][2].add(cc);
      }
    } while (!strstr(line,"-------------"));
  }
}

void testit(UInt_t nParticles = 10)
{
  TFile* outfile = new TFile("outfile.root","recreate");
  TNtupleD* nt = new TNtupleD("nt","nt","i:t0:x0:xp0:y0:yp0:phi0:dkk0:sx0:sy0:sz0:t:x:xp:y:yp:phi:dkk:sx:sy:sz");
  Double_t vs[21];
  nt->SetMarkerStyle(20);
  damap oneturn;
  smap  spin1turn;
  readMap("mapsy3.dat",&oneturn, &spin1turn);

  mu_x   = 0;
  mu_xp  = 0;
  mu_y   = 0;
  mu_yp  = 0;
  mu_l   = 0;
  mu_dkk = 0;

  sig_x   = 0* 0.02;
  sig_xp  = 0* 0.002;
  sig_y   = 0* 0.02;
  sig_yp  = 0.02;
  sig_l   = 0* 7.2;
  sig_dkk = 0* 1e-3;

  TTimeStamp* ts = new TTimeStamp();
  ts->Set();
  Double_t tbegin = ts->AsDouble();
  for (UInt_t partNr = 0; partNr<nParticles; partNr++)
  {
    if (!(partNr%10)) printf("%d\n",partNr);
    memset(vs,0,sizeof(vs));
    vs[0] = 0; // turn Nr
    vs[2] = gRandom->Gaus(mu_x,sig_x);
    vs[3] = gRandom->Gaus(mu_xp,sig_xp);
    vs[4] = gRandom->Gaus(mu_y,sig_y);
    vs[5] = gRandom->Gaus(mu_yp,sig_yp);
    vs[6] = gRandom->Gaus(mu_l,sig_l);
    vs[7] = gRandom->Gaus(mu_dkk,sig_dkk);
    vs[8] = 0; // sx
    vs[9] = 0; // sy
    vs[10] = 1; // sz
    vs[1] = calcTime(vs+2,vs[0]); // time in ns
    memcpy(vs+11,vs+1,11*sizeof(Double_t));
    // vs[11] = vs[1];
    // for (UInt_t k = 0; k<6; k++) vs[12+k] = init[k];
    // for (UInt_t k = 0; k<3; k++) vs[18+k] = spin[k];
    nt->Fill(vs);
    vs[0]++;

    for (UInt_t turn = 0; turn<4001; turn++)
    {
      spin1turn.apply(vs+12,vs+18);
      oneturn.apply(vs+12);
      vs[11] = calcTime(vs+12,vs[0]);
      if (!(turn%100)) nt->Fill(vs);
      vs[0]++;
    }
  }
  ts->Set();
  Double_t tend = ts->AsDouble();
  printf("Total time = %lf s\n",tend-tbegin);
  printf("Time per Turn = %lf mus/particle.turn\n",(tend-tbegin)*1e6/nParticles/4000);

  Double_t T0 = 4231.778309; // ns
  // nt->Draw("sx*sin(2*TMath::Pi()/4231.778309*t)-sz*cos(2*TMath::Pi()/4231.778309*t):t");
  nt->Draw("sx*cos(2*TMath::Pi()/4231.778309*t)+sz*sin(2*TMath::Pi()/4231.778309*t):t:yp0:yp0","","");

  outfile->Write();
}

Double_t calcTime(Double_t* vs, UInt_t turn)
{
  Double_t gamma = (1+vs[5])*K0/M0+1;
  Double_t time =  -(gamma+1)/gamma/(beta0*c)*vs[4] + turn*L0/(beta0*c);
  return time;
  // return vs[9]*L0/(beta0*c);
}
