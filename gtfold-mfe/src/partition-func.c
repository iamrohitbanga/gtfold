#include <stdio.h>
#include "partition-func.h"
#include "energy.h"
#include "algorithms-partition.h"
#include "global.h"
#include "utils.h"

#ifdef _OPENMP 
#include "omp.h"
#endif

#define DEBUG_PF 1

#ifdef DEBUG_PF 
  #undef Ec
  #undef Eb
  #undef Ea
  #define Ec 0
  #define Eb 0
  #define Ea 0
  #define eS StackE
  #define eH HairE
  #define eL ILoopE
  #define Ed3 D3
  #define Ed5 D5
  #define auPenalty Etp
#endif

int StackE(int i, int j) { return 0; }
int HairE(int i, int j) { return 0; }
int ILoopE(int i, int j, int ip, int jp) { return 0; }
int D3(int i, int j, int k) { return 0; } 
int D5(int i, int j, int k) { return 0; }
int Etp(int i, int j) { return 0; }

double ** u;
double ** up;
double ** upm;
double ** ud;
double ** u1d;

double ** s1;
double ** s2;
double ** s3;
double ** u1;

int part_len;

static void create_partition_arrays();
static void init_partition_arrays();
static void fill_partition_arrays();
static void free_partition_arrays();

static void calc_u(int i, int j);
static void calc_ud(int i, int j);
static void calc_up(int i, int j);
static void calc_upm(int i, int j);
static void calc_u1(int i, int j);
static void calc_u1d(int i, int j);
static void calc_s1(int i, int j);
static void calc_s2(int i, int j);
static void calc_s3(int i, int j);
static double f(int j, int h, int l);

double f(int j, int h, int l){
	if(j - 1 == l)
		return 1;
	else 
		return exp(-Ed3(h,l,l+1)/RT);
}


void calculate_partition(int len) 
{
  int i, j;
  part_len = len;
  create_partition_arrays();
  init_partition_arrays();
  fill_partition_arrays();

  for (i = 0; i <= part_len+1; ++i) 
  {
    for (j = 0; j <= part_len+1; ++j)
      printf("%0.1f ",u[i][j]);
    printf("\n");
  }

  printf("%4.4f\n",u[1][part_len]);
}

void free_partition()
{
  free_partition_arrays();
}

void init_part_arrays_zeros(){
  int i,j,n;
  n = part_len+1;
  for(i=0; i<=n; ++i){
    for(j=0; j<=n; ++j){
      u[i][j]=0;
      up[i][j]=0;
      upm[i][j]=0;
      ud[i][j]=0;
      u1d[i][j]=0;

      s1[i][j]=0;
      s2[i][j]=0;
      s3[i][j]=0;
      u1[i][j]=0;
    }
  }
}

void init_partition_arrays()
{
  init_part_arrays_zeros();

  int i, j;							
  int n = part_len;
  for(i=1; i<=n-3; ++i){
    for(j=i; j<=i+3; ++j){
      u[i][j] = 1;
      up[i][j] = 0;
      ud[i][j] = 0;
      u1[i][j] = 0;
      u1d[i][j] =0 ;
      s1[i][j] = 0;
      s2[i][j] = 0;
      s3[i][j] = 0;
    }
  }	
  for(i=1; i<=n-4; ++i){
    s1[i][i+4] = 0;
    s2[i][i+4] = 0;
  }
  for(i=1; i<=n; ++i){
    u[i+1][i] = 1;
    u1[i+1][i] = 0;
    u1d[i+1][i] = 0;
  }
  for(i=1; i<=n-1; i++){
    u1[i+2][i] = 0;
  }
}

void fill_partition_arrays()
{
  int b,i,j;
  int n=part_len;
  for(b=TURN+1; b<n; ++b){
    for(i=1; i<=n-b; ++i){
      j=i+b;
         
      // Auxillary array
      calc_s1(i,j);
      calc_s2(i,j);
      calc_s3(i,j);
      
      if (j < n) { // handled boundary condition
        calc_u1d(i,j);
        calc_u1(i,j);
      }

      calc_upm(i,j);
      calc_up(i,j);
      
      calc_ud(i,j);
      calc_u(i,j);
    }
  }
}

void calc_s1(int h, int j)
{
		int l;						
		for (l = h; l < j; ++l)
		{
				s1[h][j] = up[h][l]*(exp(-(Ed5(h,l,h-1)+auPenalty(h,l))/RT))*
						(exp(-Ed3(h,l,l+1)/RT)*u[l+2][j]+(ud[l+1][j]+
																				up[l+1][j]*exp(-(auPenalty(l+1,j)/RT))));							
		}
}

void calc_s2(int h, int j)
{
	int l;							
	for (l = h; l < j; ++l)
	{
			s2[h][j] = up[h][l]*(exp(-(Ed5(h,l,h-1)+auPenalty(h,l))/RT))*
										(exp(-Ed3(h,l,l+1)/RT)*u1[l+2][j-1]+u1d[l+1][j-1]);							
	}
}

void calc_s3(int h, int j)
{
  int l;							
  for (l = h; l < j; ++l)
  {
    s3[h][j] = up[h][l]*(exp(-(Ed5(h,l,h-1)+auPenalty(h,l))/RT))*
      (f(j+1,h,l)*exp(-((j-l)*Eb)/RT) + 
       exp(-(Ed3(h,l,l+1)+Eb)/RT)*u1[l+2][j] + u1d[l+1][j]);
  }
}

void create_partition_arrays()
{
  int len = part_len + 2;						
  u = mallocTwoD(len,len);
  up = mallocTwoD(len,len);
  upm = mallocTwoD(len,len);
  ud = mallocTwoD(len,len);
  u1d = mallocTwoD(len,len);
  s1 = mallocTwoD(len,len);
  s2 = mallocTwoD(len,len);
  s3 = mallocTwoD(len,len);
  u1 = mallocTwoD(len,len);
}

void free_partition_arrays()
{
  int len = part_len + 2;
  freeTwoD(u,len,len);
  freeTwoD(up,len,len);
  freeTwoD(upm,len,len);
  freeTwoD(ud,len,len);
  freeTwoD(u1d,len,len);
  freeTwoD(s1,len,len);
  freeTwoD(s2,len,len);
  freeTwoD(s3,len,len);
  freeTwoD(u1,len,len);
}

void calc_upm(int i, int j){
	double a = Ea, b=Eb, c=Ec;
	double p_val = 0;
	int l,h;
	double quadraticSum = 0;
	
  if (canPair(RNA[i],RNA[j]) && j-i > TURN)
  {
   for(l=i+2; l<j; ++l){
      p_val += (up[i+1][l] * exp((-1)*(a+2*c+auPenalty(i+1,l))/RT) * 
          (exp((-1)*(Ed3(i+1,l,l+1)+b)/RT) * u1[l+2][j-1] + u1d[l+1][j-1]));
    }

    for(l=i+3; l<j; ++l){
      p_val += (up[i+2][l]*exp((-1)*(a+2*c+b+Ed3(i,j,i+1)+auPenalty(i+2,l))/RT) * 
          (exp((-1)*(Ed3(i+2,l,l+1)+b)/RT)*u1[l+2][j-1] + u1d[l+1][j-1]));
    }

    for(h=i+3; h<j-1; ++h){
      quadraticSum += (s2[h][j] * exp((-1)*(a+2*c+(h-i-1)*b)/RT));	
    }
    quadraticSum *= exp((-1)*Ed3(i,j,i+1)/RT);

    p_val += quadraticSum;
    upm[i][j] = p_val;
  }
  else {
    upm[i][j] = 0;   
  }
}

void calc_u1(int i, int j){
	double b=Eb, c=Ec;
	double p_val = u1d[i][j];
	int h;
	double quadraticSum = 0;
	
	for(h=i+1; h<j; ++h){
		quadraticSum += (s3[h][j] * exp((-1)*(c+(h-i)*b)/RT));
	}
	
	p_val += quadraticSum;
	
	u1[i][j] = p_val;
}

void calc_u1d(int i, int j){
	double b=Eb, c=Ec;
	double p_val = 0;
	int l;

	for(l=i+1; l<=j; ++l){
		p_val += (up[i][l]*exp((-1)*(c+auPenalty(i,l))/RT) * (f(j+1,i,l)*exp((-1)*(j-l)*b/RT) + exp((-1)*(Ed3(i,l,l+1)+b)/RT)*u1[l+2][j] + u1d[l+1][j]));
	}

	u1d[i][j] = p_val;
}

void calc_u(int i, int j)
{
	int uval = 1;
 int h;
	int ctr;

	for (h = i+1; h < j; ++h) {
		uval = up[h][j] * exp( -(Ed5(h,j,h-1) + auPenalty(h,j)) / RT );
	}

	uval = uval + ud[i][j];

	for (ctr = i+1; ctr < j-1; ++ctr) {
		uval = uval + s1[ctr][j];
	}
	u[i][j] = uval;
}

void calc_ud(int i, int j)
{
	int l;
  double udij = 0;

	for (l = i+1; l < j; ++l) 
  {
		double val1, val2, val3;
		val1 = up[i][l];
		val1 = val1 * exp(-auPenalty(i,l) / RT);

		val2 = u[l+2][j];
		val2 = val2 * exp(-Ed3(i,l,l+1)/RT);

		val3 = ud[l+1][j];
		val3 = val3 + up[l+1][j] * exp( -auPenalty(l+1,j) / RT );

		udij += (val1 * (val2 + val3));
	}
  ud[i][j]  = udij;
}

void calc_up(int i, int j)
{
  double up_val = 0.0;
  int p,q;

  if (canPair(RNA[i],RNA[j]) && j-i > TURN)
  {
    for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
      int minq = j-i+p-MAXLOOP-2;
      if (minq < p+1+TURN) minq = p+1+TURN;
      int maxq = (p==(i+1))?(j-2):(j-1);

      for (q = minq; q <= maxq; q++) {
        if (canPair(p,q)==0) continue;
        up_val += (up[p][q] * exp(-eL(i,j,p,q)/RT));
      }
    }

    up_val = up_val * exp(-Ed3(i,j,i+1)/RT);
    up_val = up_val + exp(-eH(i,j)/RT );
    up_val = up_val + exp(-eS(i,j)/RT ) * up[i+1][j-1];
    up_val = up_val + upm[i][j];

    up[i][j] = up_val;
  }
  else  {
    up[i][j] = 0;   
  }
}
