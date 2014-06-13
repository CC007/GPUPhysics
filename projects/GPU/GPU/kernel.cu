#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
	char fileName[200] = "";
	char coefsFileName[200] = "";
	char outputFileName[200] = "";
	int separateFiles = 0;
	int xSize, dxSize, ySize, dySize, deltaSize, phiSize, argcCounter;
	Map x, dx, y, dy, delta, phi;
	Coefs c;
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

			case 's':
				sscanf(&argv[1][3] ,"%d", &separateFiles);
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
					printf("Calculates a certain amount of steps of a charged particle in a inhomogenus magnetic field.\n\n");
					printf("<executable> -h | <executable> [-m=<mapFileName>] [-c=<coeffFileName>]\n[-o=<outputFileName> [-s]]\n\n");
					printf("-h, --help\t\t Display help\n");
					printf("-m=<mapFileName>\t Set the map file to be <mapFileName>. If not set, it\n\t\t\t will be asked for in the program itself.\n");
					printf("-c=<coeffFileName>\t Set the coefficients file to be <coeffFileName>. If\n\t\t\t not set, it will be asked for in the program itself.\n\t\t\t Note that the coefficients file supports multiple\n\t\t\t particles, while if the program is run without this\n\t\t\t file, it supports only one particle.\n");
					printf("-o=<outputFileName>\t Set the output file to be <outputFileName>. If not\n\t\t\t set, it will default to stdout\n");
					printf("-s\t\t\t Choose if you want one output file or (if applicable)\n\t\t\t multiple output files. Note that this parameter can\n\t\t\t only be set if an output file is set.\n");
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
		fprintf(stderr, "Error while opening the file.\n");
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
	fprintf(stderr, "map dx\n");
	mallocMap(&dx, dxSize);
	fprintf(stderr, "map y\n");
	mallocMap(&y, ySize);
	fprintf(stderr, "map dy\n");
	mallocMap(&dy, dySize);
	fprintf(stderr, "map delta\n");
	mallocMap(&delta, deltaSize);
	fprintf(stderr, "map phi\n");
	mallocMap(&phi, phiSize);

	// allocate mempry for the coefficients
	fprintf(stderr, "map coefs\n");
	mallocCoefs(&c, ITER);

	// read some variables and the map lines from the map file
	fprintf(stderr, "open file\n");
	FILE *mapFileP = fopen(fileName, "r");
	fprintf(stderr, "check if file is NULL\n");
	if( mapFileP == NULL ){
		fprintf(stderr, "Error while opening the file.\n");
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

	// read the coefficients from user input
	fprintf(stderr, "read coefs\n");
	getCoefs(&c);

	// calculate the coefficients for 4000 iterations
	for(int i = 0;i < ITER-1;i++){
		calcCoefs(&c, i, &x, &(c.x[i+1]));
		calcCoefs(&c, i, &dx, &(c.dx[i+1]));
		calcCoefs(&c, i, &y, &(c.y[i+1]));
		calcCoefs(&c, i, &dy, &(c.dy[i+1]));
		calcCoefs(&c, i, &delta, &(c.delta[i+1]));
		calcCoefs(&c, i, &phi, &(c.phi[i+1]));
	}

	// show or save the coefficients 
	FILE* outputFile;
	if(!strncmp(outputFileName, "\0", 1)==0){
		outputFile = fopen(outputFileName, "w");
	}else{
		outputFile = stdout;
	}
	for(int i = 0;i < ITER;i++){
		fprintf(outputFile, "%10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n", c.x[i], c.dx[i], c.y[i], c.dy[i], c.delta[i], c.phi[i]);
	}
	fclose(outputFile);

	// clean up the heap and tell that the computation is finished
	freeMap(&x);
	freeMap(&dx);
	freeMap(&y);
	freeMap(&dy);
	freeMap(&delta);
	freeMap(&phi);
	freeCoefs(&c);
	fprintf(stderr, "Output is created. Press Enter to continue...\n");
	getchar();
	getchar();
	return 0;
}
