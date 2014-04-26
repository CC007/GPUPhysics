#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <string.h>

#define ITER 16

typedef struct Map{
  int length;
  double *A;
  int *x;
  int *dx;
  int *y;
  int *dy;
  int *delta;
  int *phi;
}Map;

typedef struct Coefs{
  int length;
  double *x;
  double *dx;
  double *y;
  double *dy;
  double *delta;
  double *phi;
}Coofs;

typedef struct Vars{
  double mass;
  double momentum;
  double kinEn;
  double gamma;
  double beta;
  double mAnomalyG;
  double spinTuneGgamma;
  double lRefOrbit;
}Vars;

void mallocMap(Map *m, int p){
  (*m).length = p;
  if(p>0){
    (*m).A = (double*)calloc(p,sizeof(double));
    (*m).x = (int*)calloc(p,sizeof(int));
    (*m).dx = (int*)calloc(p,sizeof(int));
    (*m).y = (int*)calloc(p,sizeof(int));
    (*m).dy = (int*)calloc(p,sizeof(int));
    (*m).delta = (int*)calloc(p,sizeof(int));
    (*m).phi = (int*)calloc(p,sizeof(int));
  }
}

void freeMap(Map *m){
    if((*m).length > 0){
      free((*m).A);
      free((*m).x);
      free((*m).dx);
      free((*m).y);
      free((*m).dy);
      free((*m).delta);
      free((*m).phi);
    }
}

void mallocCoefs(Coofs *c, int p){
  (*c).length = p;
  if(p>0){
    (*c).x = (double*)calloc(p,sizeof(double));
    (*c).dx = (double*)calloc(p,sizeof(double));
    (*c).y = (double*)calloc(p,sizeof(double));
    (*c).dy = (double*)calloc(p,sizeof(double));
    (*c).delta = (double*)calloc(p,sizeof(double));
    (*c).phi = (double*)calloc(p,sizeof(double));
  }
}

void freeCoefs(Coofs *c){
    if((*c).length > 0){
      free((*c).x);
      free((*c).dx);
      free((*c).y);
      free((*c).dy);
      free((*c).delta);
      free((*c).phi);
    }
}

void readMap(FILE *fp, Map *m){
  char* line = (char*)malloc(200*sizeof(char));
  line = fgets(line, 200, fp);
  int dum1, dum2;
  if(strncmp(line, "     I  COEFFICIENT            ORDER EXPONENTS", 46)!=0){
     if(strncmp(line, "     ALL COMPONENTS ZERO ", 24)!=0) exit(EXIT_FAILURE);
  }
  for(int i=0;!strstr((line = fgets(line, 200, fp)), "------");i++){
    sscanf(line,"%d %lf %d %d %d %d %d %d %d",
           &dum1,
           &((*m).A[i]),
           &dum2,
           &((*m).x[i]),
           &((*m).dx[i]),
           &((*m).y[i]),
           &((*m).dy[i]),
           &((*m).delta[i]),
           &((*m).phi[i])
    );
  }
  free(line);
}

void readVars(FILE *fp, Vars *v){
  char* line = (char*)malloc(200*sizeof(char));
  line = fgets(line, 200, fp);
  printf("check and store muon mass\n");
  if (sscanf(line,"Muon Mass =   %lf MeV/c^2",&((*v).mass)) != 1) exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  printf("check and store muon momentum\n");
  if (sscanf(line,"Muon Momentum =   %lf MeV/c",&((*v).momentum)) != 1) exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  printf("check and store muon kin energy\n");
  if (sscanf(line,"Muon Kinetic Energy =   %lf MeV",&((*v).kinEn)) != 1) exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  printf("check and store Muon gamma\n");
  if (sscanf(line,"Muon gamma =   %lf",&((*v).gamma)) != 1) exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  printf("check and store Muon beta\n");
  if (sscanf(line,"Muon beta =  %lf",&((*v).beta)) != 1) exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  printf("check and store Muon Anomaly G\n");
  if (sscanf(line,"Muon Anomaly G =  %lf",&((*v).mAnomalyG)) != 1) exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  printf("check and store Muon Spin Tune G.gamma\n");
  if (sscanf(line,"Muon Spin Tune G.gamma =  %lf",&((*v).spinTuneGgamma)) != 1) exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  if (sscanf(line," L    %lf",&((*v).lRefOrbit)) != 1) exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  if (line[1] !='P') exit(EXIT_FAILURE);
  line = fgets(line, 200, fp);
  if (line[1] !='A') exit(EXIT_FAILURE);
  free(line);
}

void sumArrayHelper(double *nums, int length, int interval){
  int index = 0;
  int next = interval/2;
  do{
    if(next < length){
      nums[index] += nums[next];
    }
    index += interval;
    next += interval;
  } while (index < length);
}

double sumArray(double *nums, int length){
  if(length <= 0){
    return 0;
  }
  int interval = 2;
  while(interval < length*2){
    sumArrayHelper(nums, length, interval);
    interval *= 2;
  }
  return nums[0];
}

void getCoefs(Coefs *c){
  printf("Geef de 6 coefficienten: ");

  scanf("%lf %lf %lf %lf %lf %lf",
        &((*c).x[0]),
        &((*c).dx[0]),
        &((*c).y[0]),
        &((*c).dy[0]),
        &((*c).delta[0]),
        &((*c).phi[0])
  );
}

void calcCoefs(Coefs *c, int idx, Map *m, double *newValue){
  double *nums = (double*)calloc((*m).length,sizeof(double));

  for(int i = 0; i < (*m).length; i++) {
    nums[i] = (*m).A[i] * pow((*c).x[idx],(*m).x[i])
                        * pow((*c).dx[idx],(*m).dx[i])
                        * pow((*c).y[idx],(*m).y[i])
                        * pow((*c).dy[idx],(*m).dy[i])
                        * pow((*c).delta[idx],(*m).delta[i])
                        * pow((*c).phi[idx],(*m).phi[i]);
  }

  *newValue = sumArray(nums, (*m).length);
  free(nums);
}

int main(int argc, char **argv){
  char *fileName = "Test1.dat";
  Map x;
  Map dx;
  Map y;
  Map dy;
  Map delta;
  Map phi;
  Coefs c;
  Vars v;
  printf("map x\n");
  mallocMap(&x, 17);
  printf("map dx\n");
  mallocMap(&dx, 6);
  printf("map y\n");
  mallocMap(&y, 11);
  printf("map dy\n");
  mallocMap(&dy, 1);
  printf("map delta\n");
  mallocMap(&delta, 0);
  printf("map phi\n");
  mallocMap(&phi, 0);
  printf("map coefs\n");
  mallocCoefs(&c, ITER);

  printf("open file\n");
  FILE *fp = fopen(fileName, "r");
  printf("check if file is NULL\n");
  if( fp == NULL ){
    printf("Error while opening the file.\n");
    exit(EXIT_FAILURE);
  }
  printf("read vars");
  readVars(fp,&v);
  printf("read x\n");
  readMap(fp, &x);
  printf("read dx\n");
  readMap(fp, &dx);
  printf("read y\n");
  readMap(fp, &y);
  printf("read dy\n");
  readMap(fp, &dy);
  printf("read delta\n");
  readMap(fp, &delta);
  printf("read phi\n");
  readMap(fp, &phi);

  printf("read coefs\n");
  getCoefs(&c);
  for(int i = 0;i < ITER;i++){
    printf("%f %f %f %f %f %f\n", c.x[i], c.dx[i], c.y[i], c.dy[i], c.delta[i], c.phi[i]);
    calcCoefs(&c, i, &x, &(c.x[i+1]));
    calcCoefs(&c, i, &dx, &(c.dx[i+1]));
    calcCoefs(&c, i, &y, &(c.y[i+1]));
    calcCoefs(&c, i, &dy, &(c.dy[i+1]));
    calcCoefs(&c, i, &delta, &(c.delta[i+1]));
    calcCoefs(&c, i, &phi, &(c.phi[i+1]));
  }
  fclose(fp);
  freeMap(&x);
  freeMap(&dx);
  freeMap(&y);
  freeMap(&dy);
  freeMap(&delta);
  freeMap(&phi);
  freeCoefs(&c);
  getch();
  return 0;
}
