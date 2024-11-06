#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define OVR 10000
#define NN ((double)SECS * (double)HZ)
#define NNN HZ*WMIN*60

char SAVPATH[200] = "C:\\Users\\hakoha\\Documents\\";     /* raw data directory .csv, remember "\\"  */
char KALIBR[50] = "kalibr.txt";            /* file with calibration coefficients, (format: axivity id, scale x, scale y scale z, offset x, offset y, offset z) */
char RAWFILES[50] = "rawfiles.txt";        /* text file with the processed raw data files with paths */

double bfa[9] = { 1,-5.56637844,13.5626697,-19.07259356,17.08046126,-10.01575927,3.74451246,-0.81199192,0.07907978 }; /*  fourth-order ARMA coefficients of the butterworth band-pass filter a (0,2 - 15 Hz) */
double bfb[9] = { 0.0177713 ,0,-0.07108521 ,0,0.10662781 ,0,-0.07108521 ,0,0.0177713 };  /*   fourth-order ARMA coefficients of the band-pass filter b (0,2 - 15 Hz)  */



int SECS = 1;   /* used epoch sec (usually 1) */
int HZ=100;     /* used accelerometer HZ */
int UNS=5;      /* the length of the period used to decrease the wrist angle sec (recommended 5 sec) */
int STAX=2;    /* the number of axles /3 used to define the non-wear time (recommended 2) */
int WMIN = 15; /* length of period used for non-wear time min (recommended 15 min) */
double STDg = 13;  /* max std used to define the non-wear time (recommended 13 mG) */
double AxRan = 50;  /* acceleration range (mg) used to define the non-wear time (recommended 50 mG) */
int JAKO = 2;

int ENM=0;  /* ENMO calculation formula...  0=ENMO [MAX 0], 1=ENMO [abs], 2=ENMO [raw G] (usually 0) */
int mad=0;  /* including MADs (0 = no 1 = yes)  */
int bf=0;   /* the butterworth band-pass filter (0 = no 1 = yes)  (coefficients above) */ 
int sd=2;   /* 0 = only ENMO 1 ENMO+Non Wear 2 = ENMO+Non Wear+wrist angle */

int STEPHZ = 15;   /* STEPS: step frequency (recommended 15), 0 = no steps are counted (Windowed Peak Detection) */
double PEAKTR = 1.2;    /* STEPS: peak threshold g  (recommended 1.2 g) */ 
double VARTRE = 0.01; /* STEPS: threshold variance (recommended 0.01) */
int CONTRE = 6;    /* STEPS: continuity treshold per STEPHZ (recommended 6)  */
int CONWIN = 5;    /* STEPS: continuity window per STEPHZ (recommended 5) */
int SIMTRE = 1;    /* STEPS: similarity treshold g  (recommended 1)  */
int MINP = 3;      /* STEPS: minimum peak spacing per STEPHZ (recommended 3) */
int MST = 1;       /* STEPS: maximum time interval between steps s (recommended 1) */  

int sleep(int x);

typedef struct CsvList {
	char nm[500];
	struct CsvList* next;
} CsvList;

typedef struct Calibr_t {
	long id;
	double x_off, y_off, z_off, x_sc, y_sc, z_sc;
	struct Calibr_t* next;
} Calibr_t;

Calibr_t* calibrread(const char* calibr) {
	FILE* txtf;
	Calibr_t* cp, * cali, * start = NULL;
	char tl[100], * nt;

        txtf = fopen(calibr, "r");
	if (txtf) {
		while (fgets(tl, 499, txtf) != NULL) {
			cali = (Calibr_t*)malloc(sizeof(Calibr_t));
			cali->id = (long)atoi(strtok(tl, ","));
			cali->x_sc = (double)atof(strtok(NULL, ","));
			cali->y_sc = (double)atof(strtok(NULL, ","));
			cali->z_sc = (double)atof(strtok(NULL, ","));
			cali->x_off = (double)atof(strtok(NULL, ","));
			cali->y_off = (double)atof(strtok(NULL, ","));
			cali->z_off = (double)atof(strtok(NULL, ","));
			cali->next = NULL;
			if (start == NULL) {
				start = cali;
			}
			else {
				cp = start;
				while (cp->next != NULL) {
					cp = cp->next;
				}
				cp->next = cali;
			}
		}
		fclose(txtf);
	}
		
	return start;

}

CsvList * csvread(const char* csvs) {
	FILE* csv;
	char tl[500];
	CsvList* cp, *fname,*start=NULL;
	int i;

	csv = fopen(csvs, "r");
	if (csv) {
		while (fgets(tl, 499, csv) != NULL) {
			if (tl[0] == '\n');

			else {
				fname = (CsvList*)malloc(sizeof(CsvList));
				for (i = 0; tl[i] != 0; i++) if (tl[i] == '\n') tl[i] = 0;
				strcpy(fname->nm, tl);
				fname->next = NULL;
				if (start == NULL) {
					start = fname;
				}
				else {
					cp = start;
					while (cp->next != NULL) {
						cp = cp->next;
					}
					cp->next = fname;
				}
			}
		}
		fclose(csv);
	}
	
	return start;
}

long give_accid(char* pth) {
	int i,j=0;
	char n[500];
	strcpy(n, pth);
	for (i = 1; n[i] != 0; i++) if (n[i-1] == '_' && n[i]=='0') j = i;
	j=j-6;
	n[j + 5] = 0;
	return (long)atoi(n+j);
}




void give_output(char* pth, char *pth1) {

        char opth[200];
	int i,j=0;
        pth1[0]=0;
	strcpy(opth, pth);
	for (i = 0; opth[i] != 0; i++) if (opth[i] == 92) j = i;
	opth[j] = 0;
	strcat(pth1,SAVPATH);
        strcat(pth1,"E");
        strcat(pth1,opth+j+1); 
}

Calibr_t init_calibr(void) {
	Calibr_t c;

	c.id = 0;
	c.x_sc = 1;
	c.x_off = 0;
	c.y_sc = 1;
	c.y_off = 0;
	c.z_sc = 1;
	c.z_off = 0;
	c.next = NULL;
	return c;
}


Calibr_t give_coeff(long ide, Calibr_t *cals) {
	Calibr_t c;
	c=init_calibr();
	while (cals) {
		if (cals->id == ide) {
			c.x_sc = cals->x_sc;
			c.y_sc = cals->y_sc;
			c.z_sc = cals->z_sc;
			c.x_off = cals->x_off;
			c.y_off = cals->y_off;
			c.z_off = cals->z_off;
			break;
		}
		cals = cals->next;
	}
	
	return c;

}

void add_coeff(double* s, double x, int c) {
	int i;

	for (i = 0; i < c-1; i++) s[i] = s[i+1];
	s[c - 1] = x;
}

double arma(double x, double* ar, double* ma) {

	int i;
	double x1,y,z;


	y = 0; z = 0;

	for (i = 0; i < 8; i++) {
		y = y + bfa[i + 1] * -1 * ar[7 - i];
		z = z + bfb[i + 1] * ma[7 - i];
	}

	x1 = y + z + x * bfb[0];
	add_coeff(ar, x1,8);
	add_coeff(ma, x,8);

	return x1;
}

int peakvar(double *s, long p, double ka)
{
  int i,ee=0;
  double ms=0;

  if (p>=STEPHZ*MST) return 0;

  for (i=0;i<p;i++) {
    ms+=pow(s[i]-ka/p,2);   
  }
  ms = sqrt(ms/(p-1));
  if (ms >= VARTRE) ee=1;
  
  return ee;  
}

int stepscount(double vm, int *edv,double *pkg,long *p,long *epel, int *pea, long pos, long *peap, double *peah,double *epeh, double *mpk) 
{
        int st=0,ee,ce=0,i;

        if (vm>=PEAKTR) {
          *pea=1; 
          if (vm >= *peah) {
            *peap = pos;
            *peah = vm;  
          }
        }  
        else if(*pea==1) {

          if (*peap-*epel >=MINP && fabs(*peah - *epeh)<SIMTRE) {
            *epel = *peap;
            *epeh = *peah;
            ee = peakvar(pkg,*p,*mpk); 
            ce=0;
            for (i=1;i<CONTRE;i++) { 
              ce+=edv[i]; 
              edv[i-1] = edv[i];            
            }
            edv[CONTRE-1]=ee;
            ce+=ee;
            if (ce>=CONWIN) st=1;  
          } 
          *pea=0;
          *peap=0;
          *peah=0; 
          *p=-1; *mpk=0; 
        } 
        return st;
}


double mediaa(double* ka)
{
	double m = 0;
	int i, j,k,ix,f,x,*id;
	x = 3;
	id = (int*)calloc(x, sizeof(int));
	for (i = 0; i < x; i++) {
		m = 100000;
		for (j = 1; j <= 5; j++) {
			f = 0;
			if (ka[j] < m) {
				for (k = 0; id[k] > 0; k++) {
					if (id[k] == j) f = 1;
				}
				if (f == 0) { m = ka[j]; ix = j; }
			}
		}
		id[i] = ix;
	}
        free(id);
	return m;
}

double medi(double *m, int k)
{
  int i,j;
  for (i=1;i<=k;i++) {
    if (i%5==0) {
      j=(int)(i/5)-1;
      m[j]=mediaa(m+j*5);
    }
  }
  return j-1;  
}


double medofmed(double *s, int p)
{
  double x,*m;
  int k,i;
  k=p;
  
  m = (double*)calloc(p, sizeof(double));

 for(i=0;i<p;i++) m[i]=s[i];

  while (k>5) {
    k = medi(m,k);
  }
  x = mediaa(m);
  if (x == 0) x = 0.0001; 
  free(m);
  return x;
}


int main(int argc, char* argv[]) {

	char tl[200], * nt, fn[200],opr[300]="",num[50];
	char pth[200] = "";
	char pth1[200],pth2[200];
	char * dt, oname[500];
	double x[3], vm, enmo, men = 0, x1, x2, x3;
	FILE* csv, * wcsv;
        long i1, pc,sc;
	int i, c=0, sec = 0,sedt, sedpr,cad,b,sth=0,sps, *edv, pea;
	CsvList* csvs;
	Calibr_t* cals, clb;
	long accid,epel,peap;
	double ar1[8], ar2[8], ar3[8], ma1[8], ma2[8],ma3[8];
	double *madr,mka=0, madx=0, *sdx,*sdy,*sdz,mkax=0,mkay=0,mkaz=0;
	double ssdx = 0,ssdy=0,ssdz=0,*meax,*meay,*meaz,angl;
        double meamxx,meamxy,meamxz,meamnx,meamny,meamnz,stac;
        double bx1,bx2,bx3,*pkg,peah,epeh,mpk;

	if (argc > 1) strcpy(pth, argv[1]);

	strcpy(pth1, pth);
	strcpy(pth2, pth);

	strcat(pth1, RAWFILES);
	strcat(pth2, KALIBR);

	cals = calibrread(pth2);
	if (!(csvs = csvread(pth1))) return 0;
	
	madr = (double*)calloc(HZ*SECS, sizeof(double));
        sdx = (double*)calloc(NNN, sizeof(double));
        sdy = (double*)calloc(NNN, sizeof(double));
        sdz = (double*)calloc(NNN, sizeof(double));
        meax = (double*)calloc(UNS*HZ, sizeof(double));
        meay = (double*)calloc(UNS*HZ, sizeof(double));
        meaz = (double*)calloc(UNS*HZ, sizeof(double));
        pkg = (double*)calloc(MST*STEPHZ, sizeof(double));
        edv = (int*)calloc(CONTRE, sizeof(double));
     
  	dt = (char*)calloc(100, sizeof(char));

        if (STEPHZ>0) sth = (int)round(HZ/STEPHZ); 

	dt[0] = 0;


	while (csvs) {
		for (i = 0; i < 8; i++) { 
		 ar1[i] = 1; ar2[i] = 0; ar3[i] = 0; 
		 ma1[i] = 1; ma2[i] = 0; ma3[i] = 0;
		}

               for (i=0;i<UNS*HZ;i++) {
                   meax[i]=1; meay[i]=1; meaz[i]=1;   
                }
                strcpy(fn, csvs->nm);

		csv = fopen(fn, "r");
		if (csv) {
			accid = give_accid(csvs->nm);
			clb = give_coeff(accid, cals);
			give_output(csvs->nm, oname);
			wcsv = fopen(oname, "w");
 
			c = 0; b=0;
			i1 = 0;
			sec = 0;
			men = 0; mka = 0;
  		        mkax = 0; bx1=0;
           		mkay = 0; bx2=0;
         		mkaz = 0; bx3=0;
                        cad=0;
                        sedpr=0; mpk=0; peap=0;
                        sc=0; stac=0; sps=0; epel=1; pc=0; pea=0; peah=0; epeh=1;
                        for (i=0;i<CONTRE-2;i++) edv[i]=1;

			while (fgets(tl, 199, csv) != NULL) {
				if (tl[0]<48 || tl[0]>57) continue;
				for (i = 0; i < 4; i++) {
					if (i == 0) {
						dt = strtok(tl, ",");
						dt[19] = 0;

					}
					else {
						x[i - 1] = (double)atof(strtok(NULL, ","));

					}
				}

				if (mad) {
					madr[i1] = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
					mka += sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
				}

				x1 = clb.x_sc * x[0] + clb.x_off;
				x2 = clb.y_sc * x[1] + clb.y_off;
				x3 = clb.z_sc * x[2] + clb.z_off;

				if (bf) {  
					x1 = arma(x1, ar1, ma1);
					x2 = arma(x2, ar2, ma2);
					x3 = arma(x3, ar3, ma3);
				}

			        if (sd) {
                                    b++;
                                    bx1 += x1;bx2 += x2;bx3 += x3;
                                    if (b>=JAKO) {
                         	      sdx[i1] = bx1/JAKO;
				      sdy[i1] = bx2/JAKO;
				      sdz[i1] = bx3/JAKO;
				      mkax += bx1/JAKO;
				      mkay += bx2/JAKO;
				      mkaz += bx3/JAKO;      
                                       b=0; bx1=0; bx2=0; bx3=0;
                                    }
                                    i1++;
         			}  
                               
 
				vm = (double)sqrt(x1*x1 + x2*x2 + x3*x3);

                                if (ENM==0)
				  enmo = (vm - 1 > 0) ? vm - 1 : 0;
                                else if (ENM == 1)
                                  enmo = fabs(vm - 1); 
                                else
                                  enmo = vm;
                                
                                meax[cad]=x1; meay[cad]=x2; meaz[cad]=x3;
                                
                                cad++;
                                if (cad>=HZ*UNS) cad=0;


				men += enmo;
				c++;
                                

                                if (sth>0) {      /* askeleet */
                                  sc++; 
                                  stac += vm;
                                  if (sc%sth==0) {
                                    if (pc<MST*STEPHZ) { pkg[pc]=stac/sth; mpk+=stac/sth; }
                       
                                    sps += stepscount(stac/sth,edv,pkg,&pc,&epel,&pea,sc,&peap,&peah,&epeh, &mpk);
                                    stac=0; pc++;

                                  }        
                                }

				if (sd && i1 >=NNN) {   /* non-wear */
                                     meamxx = -10000; meamxy = -10000; meamxz = -10000;
                                     meamnx = 10000; meamny = 10000; meamnz = 10000;
                                     for (i=JAKO-1;i<NNN;i+=JAKO) {
						  ssdx += pow((1000 * mkax / NNN - (1000 * sdx[i])),2);
						  ssdy += pow((1000 * mkay / NNN - (1000 * sdy[i])),2);
						  ssdz += pow((1000 * mkaz / NNN - (1000 * sdz[i])),2);
                                                  if (meamxx < sdx[i]) meamxx = sdx[i];
                                                  if (meamxy < sdy[i]) meamxy = sdy[i];
                                                  if (meamxz < sdz[i]) meamxz = sdz[i];
                                                  if (meamnx > sdx[i]) meamnx = sdx[i];
                                                  if (meamny > sdy[i]) meamny = sdy[i];
                                                  if (meamnz > sdz[i]) meamnz = sdz[i];
	   			     }
                                     
			             sedt = 0; i1=0; 
			             if (sqrt(ssdx / NNN) < STDg || 1000*(meamxx - meamnx) <AxRan) sedt++;
				     if (sqrt(ssdy / NNN) < STDg || 1000*(meamxy - meamny) <AxRan) sedt++;
				     if (sqrt(ssdz / NNN) < STDg || 1000*(meamxz - meamnz) <AxRan) sedt++;
				     if (sedt >= STAX) sedpr +=WMIN;
				     else sedpr = 0;
                                       	ssdx=0; ssdy=0; ssdz=0; mkax=0; mkay=0; mkaz=0;  	
				} 


				if (c>=HZ) {
					c = 0;
					sec++;

	                                if (sd == 2) angl = atan(medofmed(meaz,HZ*UNS) / (sqrt(pow(medofmed(meax,HZ*UNS), 2) + pow(medofmed(meay,HZ*UNS), 2))))*180/ 3.141592653589793;

					if (sec >= SECS) {
					    if (mad) {
                                              for (i = 0; i < SECS * HZ; i++) madx += fabs(1000 * mka / NN - (1000 * madr[i]));
                                            }
					    men = men / NN * 1000;
					    if (men > OVR) men = 0;
                                            sprintf(opr,"%s,%4.3f",dt,men);
                                            if (mad) {
                                                sprintf(num,",%4.3f",madx / NN);
                                                strcat(opr,num);  
                                                madx=0;                                        
                                            }  
                                            if (sd) {
                                              sprintf(num,",%d",sedpr);
                                              
                                              strcat(opr,num);
                                              if (sd==2) {
                                                sprintf(num,",%2.3f",angl);
                                                strcat(opr,num);
                                              }
                                            }
                                            if (sth) {
                                              sprintf(num,",%d",sps);
                                              sps=0;
                                              strcat(opr,num);
   
                                            }

                                            fprintf(wcsv, "%s\n", opr);

					    men = 0;
					    sec = 0;
					    mka = 0;

                                           }

				}
				
			}
			fclose(csv);
			fclose(wcsv);
			printf("%s OK\n", csvs->nm);
                        sleep(10000);
		}
		
		csvs = csvs->next;

	}
	return 0;
}
