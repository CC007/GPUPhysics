struct vv {
  Double_t amp;
  Int_t    p[6];

  vv() {
    amp = 1;
    memset(p,0,sizeof(p));
  }

  Double_t eval(Double_t* x)
  {
    Double_t res = 1;
    for (UInt_t i = 0; i<6; i++)
    {
      res *= pow(x[i],p[i]);
    }
    return res*amp;
  }

  void print()
  {
    printf("Amp : %e\n",amp);
    for (UInt_t i = 0 ; i<6; i++)
    {
      printf("N%d  : %d\n",i,p[i]);
    }
  }
};

void test()
{
  vv das[6];

  das[0].amp = 2;
  das[1].amp = 4;
  das[2].amp = 8;
  das[3].amp = 16;
  das[4].amp = 32;
  das[5].amp = 64;
  for (UInt_t i = 0; i<6; i++)
  {
    das[i].p[0] = 1;
  }

  Double_t x[6];

  for (UInt_t i = 0 ; i<6; i++)
  {
    memset(x,0,sizeof(x));
    x[i] = i+1;
    for (UInt_t j = 0; j<6; j++)
    {
      printf("%d: %f\n",j, v(x,das[j]));
    }
  }

}
