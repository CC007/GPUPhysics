#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define ITER 4000

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
}Coefs;

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

void scanFile(FILE* fp, int *size){
	char* line = (char*)malloc(200*sizeof(char));
	line = fgets(line, 200, fp);
	if(strncmp(line, "     I  COEFFICIENT            ORDER EXPONENTS", 46)!=0){
		if(strncmp(line, "     ALL COMPONENTS ZERO ", 25)!=0){
			exit(EXIT_FAILURE);
		}else{
			*size=1;
		}
	}else{
		for((*size)=0;!strstr((line = fgets(line, 200, fp)), "------");(*size)++);
	}
	free(line);
}

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

void mallocCoefs(Coefs **c, int iter, int p){
	if(iter>0){
		int i;
		(*c) = (Coefs*) malloc(p*sizeof(Coefs));
		for(i=0;i<p;i++){
			(*c)[i].length = iter;
			(*c)[i].x = (double*)calloc(iter,sizeof(double));
			(*c)[i].dx = (double*)calloc(iter,sizeof(double));
			(*c)[i].y = (double*)calloc(iter,sizeof(double));
			(*c)[i].dy = (double*)calloc(iter,sizeof(double));
			(*c)[i].delta = (double*)calloc(iter,sizeof(double));
			(*c)[i].phi = (double*)calloc(iter,sizeof(double));
		}
	}
}

void freeCoefs(Coefs **c, int p){
	if((*c)[0].length > 0){
		int i;
		for(i=0;i<p;i++){
			free((*c)[i].x);
			free((*c)[i].dx);
			free((*c)[i].y);
			free((*c)[i].dy);
			free((*c)[i].delta);
			free((*c)[i].phi);
		}
		free(*c);
	}
}

void cudaMallocMap(Map **m, int p){
	*m = (Map*) malloc(sizeof(Map));
	//cudaMalloc((void**)m, sizeof(Map));
	(**m).length = p;
	//cudaMemset(&((**m).length), p, sizeof(int));
	if(p>0){
		//double **h_arr1 = (double**)malloc(sizeof(double*));
		//int **h_arr2 = (int**)malloc(sizeof(int*));
		//cudaMemcpy(h_arr1, &((**m).A), sizeof(double*), cudaMemcpyDeviceToHost);
		//cudaMalloc((void**)h_arr1, p*sizeof(double));
		//cudaMemcpy(h_arr2, &((**m).x), sizeof(int*), cudaMemcpyDeviceToHost);
		//cudaMalloc((void**)h_arr2, p*sizeof(int));
		//cudaMemcpy(h_arr2, &((**m).dx), sizeof(int*), cudaMemcpyDeviceToHost);
		//cudaMalloc((void**)h_arr2, p*sizeof(int));
		//cudaMemcpy(h_arr2, &((**m).y), sizeof(int*), cudaMemcpyDeviceToHost);
		//cudaMalloc((void**)h_arr2, p*sizeof(int));
		//cudaMemcpy(h_arr2, &((**m).dy), sizeof(int*), cudaMemcpyDeviceToHost);
		//cudaMalloc((void**)h_arr2, p*sizeof(int));
		//cudaMemcpy(h_arr2, &((**m).delta), sizeof(int*), cudaMemcpyDeviceToHost);
		//cudaMalloc((void**)h_arr2, p*sizeof(int));
		//cudaMemcpy(h_arr2, &((**m).phi), sizeof(int*), cudaMemcpyDeviceToHost);
		//cudaMalloc((void**)h_arr2, p*sizeof(int));
		//free(h_arr1);
		//free(h_arr2);
		cudaMalloc((void**)&((**m).A), p*sizeof(double));
		cudaMalloc((void**)&((**m).x), p*sizeof(int));
		cudaMalloc((void**)&((**m).dx), p*sizeof(int));
		cudaMalloc((void**)&((**m).y), p*sizeof(int));
		cudaMalloc((void**)&((**m).dy), p*sizeof(int));
		cudaMalloc((void**)&((**m).delta), p*sizeof(int));
		cudaMalloc((void**)&((**m).phi), p*sizeof(int));
	}
}

void cudaFreeMap(Map **m){
	//Map h_map;
	//cudaMemcpy(&h_map, *m, sizeof(Map), cudaMemcpyDeviceToHost);

	if((**m).length > 0){
		//cudaFree(h_map.A);
		//cudaFree(h_map.x);
		//cudaFree(h_map.dx);
		//cudaFree(h_map.y);
		//cudaFree(h_map.dy);
		//cudaFree(h_map.delta);
		//cudaFree(h_map.phi);
		cudaFree((**m).A);
		cudaFree((**m).x);
		cudaFree((**m).dx);
		cudaFree((**m).y);
		cudaFree((**m).dy);
		cudaFree((**m).delta);
		cudaFree((**m).phi);
	}
	free(*m);
}

void cudaMallocCoefs(Coefs **c, int iter, int p){
	if(iter>0){
		int i;
		*c = (Coefs*) malloc(p*sizeof(Coefs));
		//cudaMalloc((void**)c, p*sizeof(Coefs));
		(*c)[0].length = iter;
		for(i=0;i<p;i++){
			/*double **h_arr = (double**)malloc(sizeof(double*));
			cudaMemset(&((*c)[i].length), iter, sizeof(int));
			cudaMemcpy(h_arr, &((*c)[i].x), sizeof(double*), cudaMemcpyDeviceToHost);
			cudaMalloc((void**)h_arr, iter*sizeof(double));
			cudaMemcpy(h_arr, &((*c)[i].dx), sizeof(double*), cudaMemcpyDeviceToHost);
			cudaMalloc((void**)h_arr, iter*sizeof(double));
			cudaMemcpy(h_arr, &((*c)[i].y), sizeof(double*), cudaMemcpyDeviceToHost);
			cudaMalloc((void**)h_arr, iter*sizeof(double));
			cudaMemcpy(h_arr, &((*c)[i].dy), sizeof(double*), cudaMemcpyDeviceToHost);
			cudaMalloc((void**)h_arr, iter*sizeof(double));
			cudaMemcpy(h_arr, &((*c)[i].delta), sizeof(double*), cudaMemcpyDeviceToHost);
			cudaMalloc((void**)h_arr, iter*sizeof(double));
			cudaMemcpy(h_arr, &((*c)[i].phi), sizeof(double*), cudaMemcpyDeviceToHost);
			cudaMalloc((void**)h_arr, iter*sizeof(double));
			free(h_arr);*/
			cudaMalloc((void**)&((*c)[i].x), iter*sizeof(double));
			cudaMalloc((void**)&((*c)[i].dx), iter*sizeof(double));
			cudaMalloc((void**)&((*c)[i].y), iter*sizeof(double));
			cudaMalloc((void**)&((*c)[i].dy), iter*sizeof(double));
			cudaMalloc((void**)&((*c)[i].delta), iter*sizeof(double));
			cudaMalloc((void**)&((*c)[i].phi), iter*sizeof(double));
		}
	}
}

void cudaFreeCoefs(Coefs **c, int p){
	//Coefs h_coefs;
	int i;
	for(i=0;i<p;i++){
		//cudaMemcpy(&h_coefs, &((*c)[i]), sizeof(Coefs), cudaMemcpyDeviceToHost);
		cudaFree((*c)[i].x);
		cudaFree((*c)[i].dx);
		cudaFree((*c)[i].y);
		cudaFree((*c)[i].dy);
		cudaFree((*c)[i].delta);
		cudaFree((*c)[i].phi);
	}
	free(*c);
}

void readMap(FILE *fp, Map *m, int nr){
	char* line = (char*)malloc(200*sizeof(char));
	line = fgets(line, 200, fp);
	int dum1, dum2;
	if(strncmp(line, "     I  COEFFICIENT            ORDER EXPONENTS", 46)!=0){
		if(strncmp(line, "     ALL COMPONENTS ZERO ", 25)!=0){
			exit(EXIT_FAILURE);
		}
		else{
			(*m).A[0] = 1.0;
			(*m).x[0] = nr == 0?1:0;
			(*m).dx[0] = nr == 1?1:0;
			(*m).y[0] = nr == 2?1:0;
			(*m).dy[0] = nr == 3?1:0;
			(*m).delta[0] = nr == 4?1:0;
			(*m).phi[0] = nr == 5?1:0;
		}
	}
	for(int i=0;!strstr((line = fgets(line, 200, fp)), "------");i++){
		//TODO read chars ipv ints
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
	fprintf(stderr, "check and store muon mass\n");
	if (sscanf(line,"Muon Mass =   %lf MeV/c^2",&((*v).mass)) != 1) exit(EXIT_FAILURE);
	line = fgets(line, 200, fp);
	fprintf(stderr, "check and store muon momentum\n");
	if (sscanf(line,"Muon Momentum =   %lf MeV/c",&((*v).momentum)) != 1) exit(EXIT_FAILURE);
	line = fgets(line, 200, fp);
	fprintf(stderr, "check and store muon kin energy\n");
	if (sscanf(line,"Muon Kinetic Energy =   %lf MeV",&((*v).kinEn)) != 1) exit(EXIT_FAILURE);
	line = fgets(line, 200, fp);
	fprintf(stderr, "check and store Muon gamma\n");
	if (sscanf(line,"Muon gamma =   %lf",&((*v).gamma)) != 1) exit(EXIT_FAILURE);
	line = fgets(line, 200, fp);
	fprintf(stderr, "check and store Muon beta\n");
	if (sscanf(line,"Muon beta =  %lf",&((*v).beta)) != 1) exit(EXIT_FAILURE);
	line = fgets(line, 200, fp);
	fprintf(stderr, "check and store Muon Anomaly G\n");
	if (sscanf(line,"Muon Anomaly G =  %lf",&((*v).mAnomalyG)) != 1) exit(EXIT_FAILURE);
	line = fgets(line, 200, fp);
	fprintf(stderr, "check and store Muon Spin Tune G.gamma\n");
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

void scanCoefs(char *fileName, int *count){
	char* line = (char*)malloc(200*sizeof(char));
	FILE *fp = fopen(fileName, "r");
	if( fp == NULL ){
		fprintf(stderr, "Error while opening the coefficients file: %s\n", fileName);
		exit(EXIT_FAILURE);
	}
	for((*count)=0;fgets(line, 200, fp) != NULL;(*count)++){
		if( strncmp(line, "\n", 1)==0 || strncmp(line, "\0", 1)==0){
			(*count)--;
		}
	}
	if(!feof(fp)){
		fprintf(stderr, "Something was wrong with the coefficient file!");
		exit(EXIT_FAILURE);
	}
	fclose(fp);
	free(line);
}

void getCoefs(Coefs *c){
	fprintf(stderr, "Begin values of the 6 dimentions: ");

	scanf("%lf %lf %lf %lf %lf %lf",
		&((*c).x[0]),
		&((*c).dx[0]),
		&((*c).y[0]),
		&((*c).dy[0]),
		&((*c).delta[0]),
		&((*c).phi[0])
		);
}

void readCoefs(Coefs **c, char *fileName, int count){
	FILE *fp = fopen(fileName, "r");
	if( fp == NULL ){
		fprintf(stderr, "Error while opening the coefficients file: %s\n", fileName);
		exit(EXIT_FAILURE);
	}
	int i;
	for(i=0;i<count;i++){
		fscanf(fp, "%lf %lf %lf %lf %lf %lf",
			&((*c)[i].x[0]),
			&((*c)[i].dx[0]),
			&((*c)[i].y[0]),
			&((*c)[i].dy[0]),
			&((*c)[i].delta[0]),
			&((*c)[i].phi[0])
			);
	}
	fclose(fp);
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

void kernel(Coefs **c, Map *x, Map *dx, Map *y, Map *dy, Map *delta, Map *phi, int *particleCount, int *iter, char **outputFileName, int *separateFiles){
	for(int n = 0;n < *particleCount;n++){
		for(int i = 0;i < (*iter)-1;i++){
			calcCoefs(&((*c)[n]), i, x, &((*c)[n].x[i+1]));
			calcCoefs(&((*c)[n]), i, dx, &((*c)[n].dx[i+1]));
			calcCoefs(&((*c)[n]), i, y, &((*c)[n].y[i+1]));
			calcCoefs(&((*c)[n]), i, dy, &((*c)[n].dy[i+1]));
			calcCoefs(&((*c)[n]), i, delta, &((*c)[n].delta[i+1]));
			calcCoefs(&((*c)[n]), i, phi, &((*c)[n].phi[i+1]));
		}

		// show or save the coefficients 
		FILE* outputFile;
		char fullOutputFileName[200] = "";
		if(!strncmp(*outputFileName, "\0", 1)==0){
			if(*separateFiles==1){
				sprintf(fullOutputFileName, "part%09d.%s", n+1, *outputFileName);
			}else{
				sprintf(fullOutputFileName, "%s", *outputFileName);
			}
			if(n==0){
				outputFile = fopen(fullOutputFileName, "w");
			}else{
				if(*separateFiles==1){
					outputFile = fopen(fullOutputFileName, "w");
				}else{
					outputFile = fopen(fullOutputFileName, "a");
				}
			}
			if( outputFile == NULL ){
				fprintf(stderr, "Error while opening the output file: %s\n", fullOutputFileName);
				exit(EXIT_FAILURE);
			}
		}else{
			outputFile = stdout;
		}
		for(int i = 0;i < *iter;i++){
			fprintf(outputFile, "%10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n", (*c)[n].x[i], (*c)[n].dx[i], (*c)[n].y[i], (*c)[n].dy[i], (*c)[n].delta[i], (*c)[n].phi[i]);
		}
		fprintf(outputFile, "\n");
		if(!strncmp(*outputFileName, "\0", 1)==0){
			fclose(outputFile);
		}
	}
}

__global__ void testVals(Map *m, Coefs *c){
	m->length = 42;
	m->A[0] = 1.5;
	m->dx[20] = 5;
	c[0].delta[4] = 3.14159265;
}


int main(int argc, char **argv){
	char fileName[200] = "";
	char coefsFileName[200] = "";
	char *outputFileName = (char*)malloc(200*sizeof(char));
	outputFileName[0] = '\0';
	int separateFiles = 0;
	int xSize, dxSize, ySize, dySize, deltaSize, phiSize, argcCounter, particleCount = 1, iter = ITER;
	Map x, dx, y, dy, delta, phi;
	Map *dev_x, *dev_dx, *dev_y, *dev_dy, *dev_delta, *dev_phi;
	Coefs *c;
	Coefs *dev_c;
	Vars v;

	// Read the program arguments
	argcCounter = argc;
	while ((argcCounter > 1) && (argv[1][0] == '-'))
	{
		if(argv[1][2] == '=' && argv[1][3] != '\0'){
			switch (argv[1][1]){
			case 'm':
				sprintf(fileName, "%s",&argv[1][3]);
				break;

			case 'c':
				sprintf(coefsFileName ,"%s",&argv[1][3]);
				break;

			case 'o':
				sprintf(outputFileName ,"%s",&argv[1][3]);
				break;

			case 'i':
				sscanf(&argv[1][3] ,"%d", &iter);
				break;

			default:
				fprintf(stderr, "Wrong Argument: %s\n", argv[1]);
				exit(EXIT_FAILURE);
			}
		}else{
			switch (argv[1][1]){
			case 's':
				separateFiles = 1;
				break;
			case '-':
				if(!strstr(&argv[1][2], "help") || argv[1][6] != '\0'){
					fprintf(stderr, "Wrong Argument: %s\n", argv[1]);
					exit(EXIT_FAILURE);
				}
			case 'h':
				if(strstr(&argv[1][2], "help") ||  argv[1][2] == '\0'){
					printf("Calculates a certain amount of steps of a charged particle in a inhomogeneous\nmagnetic field.\n\n");
					printf("<executable> (-h|--help) | <executable> [-m=<mapFileName>] [-c=<coeffFileName>]\n[-o=<outputFileName> [-s]] [-i=<nr>]\n\n");
					printf("-h, --help\t\t Display help\n");
					printf("-m=<mapFileName>\t Set the map file to be <mapFileName>. If not set, it\n\t\t\t will be asked for in the program itself.\n");
					printf("-c=<coeffFileName>\t Set the coefficients file to be <coeffFileName>. If\n\t\t\t not set, it will be asked for in the program itself.\n\t\t\t Note that the coefficients file supports multiple\n\t\t\t particles, while if the program is run without this\n\t\t\t file, it supports only one particle.\n");
					printf("-o=<outputFileName>\t Set the output file to be <outputFileName>. If not\n\t\t\t set, it will default to stdout\n");
					printf("-s\t\t\t Choose if you want one output file or (if applicable)\n\t\t\t multiple output files. Note that this parameter can\n\t\t\t only be set if an output file is set.\n");
					printf("-i=<nr>\t\t\t Set the number of iterations to <nr>. If not set, it\n\t\t\t will default to 4000.\n");
					exit(EXIT_SUCCESS);
				}else{
					fprintf(stderr, "Wrong Argument: %s\n", argv[1]);
					exit(EXIT_FAILURE);
				}
				break;

			default:
				fprintf(stderr, "Wrong Argument: %s\n", argv[1]);
				exit(EXIT_FAILURE);
			}
		}
		++argv;
		--argcCounter;
	}
	if(separateFiles == 1 && strncmp(outputFileName, "\0", 1)==0){
		fprintf(stderr, "-s shouldn't be used without setting an output file");
		exit(EXIT_FAILURE);
	}

	// if not set in argument, ask for file name of the map file
	if(strncmp(fileName, "\0", 1)==0){
		fprintf(stderr, "Filename of the map: ");
		scanf("%s", fileName);
	}

	// use the map file to gather the sizes of the 6 coefficients
	fprintf(stderr, "open file\n");
	FILE *scanFileP = fopen(fileName, "r");
	fprintf(stderr, "check if file is NULL\n");
	if( scanFileP == NULL ){
		fprintf(stderr, "Error while opening the map file: %s\n", fileName);
		exit(EXIT_FAILURE);
	}
	fprintf(stderr, "Get map sizes\n");
	char* line = (char*)malloc(200*sizeof(char));
	do{
		line = fgets(line, 200, scanFileP);
	}while(!strstr(line, " A "));
	free(line);
	scanFile(scanFileP, &xSize);
	scanFile(scanFileP, &dxSize);
	scanFile(scanFileP, &ySize);
	scanFile(scanFileP, &dySize);
	scanFile(scanFileP, &deltaSize);
	scanFile(scanFileP, &phiSize);
	fclose(scanFileP);
	fprintf(stderr, "\nmap sizes: %d %d %d %d %d %d\n", xSize, dxSize, ySize, dySize, deltaSize, phiSize);

	// allocate memory for the map
	fprintf(stderr, "\nmap x\n");
	mallocMap(&x, xSize);
	cudaMallocMap(&dev_x, xSize);
	fprintf(stderr, "map dx\n");
	mallocMap(&dx, dxSize);
	cudaMallocMap(&dev_dx, dxSize);
	fprintf(stderr, "map y\n");
	mallocMap(&y, ySize);
	cudaMallocMap(&dev_y, ySize);
	fprintf(stderr, "map dy\n");
	mallocMap(&dy, dySize);
	cudaMallocMap(&dev_dy, dySize);
	fprintf(stderr, "map delta\n");
	mallocMap(&delta, deltaSize);
	cudaMallocMap(&dev_delta, deltaSize);
	fprintf(stderr, "map phi\n");
	mallocMap(&phi, phiSize);
	cudaMallocMap(&dev_phi, phiSize);

	// read some variables and the map lines from the map file
	fprintf(stderr, "open file\n");
	FILE *mapFileP = fopen(fileName, "r");
	fprintf(stderr, "check if file is NULL\n");
	if( mapFileP == NULL ){
		fprintf(stderr, "Error while opening the map file: %s\n", fileName);
		exit(EXIT_FAILURE);
	}
	fprintf(stderr, "read vars");
	readVars(mapFileP,&v);
	fprintf(stderr, "read x\n");
	readMap(mapFileP, &x, 0);
	fprintf(stderr, "read dx\n");
	readMap(mapFileP, &dx, 1);
	fprintf(stderr, "read y\n");
	readMap(mapFileP, &y, 2);
	fprintf(stderr, "read dy\n");
	readMap(mapFileP, &dy, 3);
	fprintf(stderr, "read delta\n");
	readMap(mapFileP, &delta, 4);
	fprintf(stderr, "read phi\n");
	readMap(mapFileP, &phi, 5);
	fclose(mapFileP);


	// allocate mempry for the coefficients

	// read the coefficients from user input
	if(strncmp(coefsFileName, "\0", 1)==0){
		fprintf(stderr, "map coefs\n");
		mallocCoefs(&c, iter, particleCount);
		cudaMallocCoefs(&dev_c, iter, particleCount);
		fprintf(stderr, "read coefs\n");
		getCoefs(c);
	}else{
		scanCoefs(coefsFileName, &particleCount);
		fprintf(stderr, "map coefs\n");
		mallocCoefs(&c, iter, particleCount);
		cudaMallocCoefs(&dev_c, iter, particleCount);
		fprintf(stderr, "read coefs\n");
		readCoefs(&c, coefsFileName, particleCount);
		fprintf(stderr, "Particle count: %d\n", particleCount);
	}


	// cuda test kernel
	testVals<<<1, 1>>>(dev_x, dev_c);
	int testval1;
	double testval2;
	double testval3;
	int length1;
	int length2;
	cudaMemcpy(&testval1, &(dev_x->dx[20]), sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&testval2, &(dev_x->A[0]), sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&testval3, &(dev_c[0].delta[4]), sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&length1, &(dev_x->length), sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&length2, &(dev_c->length), sizeof(int), cudaMemcpyDeviceToHost);

	// calculate the coefficients for 4000 iterations
	kernel(&c, &x, &dx, &y, &dy, &delta, &phi, &particleCount, &iter, &outputFileName, &separateFiles);


	// print test results
	fprintf(stderr, "Length map: %d\nLength coefs: %d\nValue map A[0]: %lf\nValue map dx[20]: %d\nValue coefs[0] delta[4]: %lf\n", length1, length2, testval2, testval1, testval3);


	// clean up the heap and tell that the computation is finished
	freeMap(&x);
	freeMap(&dx);
	freeMap(&y);
	freeMap(&dy);
	freeMap(&delta);
	freeMap(&phi);
	cudaFreeMap(&dev_x);
	cudaFreeMap(&dev_dx);
	cudaFreeMap(&dev_y);
	cudaFreeMap(&dev_dy);
	cudaFreeMap(&dev_delta);
	cudaFreeMap(&dev_phi);
	freeCoefs(&c, particleCount);
	cudaFreeCoefs(&dev_c, particleCount);
	free(outputFileName);
	fprintf(stderr, "Output is created. Press Enter to continue...\n");
	getchar();
	if(strncmp(coefsFileName, "\0", 1)==0){
		getchar();
	}
	return 0;
}
