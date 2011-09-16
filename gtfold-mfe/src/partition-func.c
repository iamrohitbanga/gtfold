#include "partition-func.h"
#include "energy.h"
#include "algorithms-partition.h"


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
  part_len = len;
  create_partition_arrays();
  init_partition_arrays();
  fill_partition_arrays();
}

void free_partition()
{
  free_partition_arrays();
}

void init_partition_arrays()
{
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
  for(i=1; i<n; ++i){
    u[i+1][i] = 1;
    u1[i+1][i] = 0;
    u1d[i+1][i] = 0;
  }
  for(i=1; i<n-1; i++){
    u1[i+2][i] = 0;
  }
}

void fill_partition_arrays()
{
  int b,i;
  int n=part_len;
  for(b=1; b<n; ++b){
    for(i=1; i<=n-b; ++i){
      int j=i+b;
      calc_s1(i,j);
      calc_s2(i,j);
      calc_s3(i,j);
      calc_u1d(i,j);
      calc_u1(i,j);
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
  int len = part_len + 1;						
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
  int len = part_len + 1;
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
    if (l > part_len-2) continue;
		p_val += (up[i][l]*exp((-1)*(c+auPenalty(i,l))/RT) * (f(j+1,i,l)*exp((-1)*(j-l)*b/RT) + exp((-1)*(Ed3(i,l,l+1)+b)/RT)*u1[l+2][j] + u1d[l+1][j]));
	}

	u1d[l+1][j] = p_val;
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
}

void calc_ud(int i, int j)
{
	int l;

	for (l = i+1; l < j; ++l) {
		double val1, val2, val3;
		val1 = up[i][l];
		val1 = val1 * exp(-auPenalty(i,l) / RT);

		val2 = u[l+2][j];
		val2 = val2 * exp(-Ed3(i,l,l+1)/RT);

		val3 = ud[l+1][j];
		val2 = val2 + up[l+1][j] * exp( -auPenalty(l+1,j) / RT );

		ud[i][j] = ud[i][j] + (val1 * val2);
	}
}

void calc_up(int i, int j)
{
	double up_val = 0.0;
 int h,l;

	for (h = i+1; h < j-1; ++h) {
		for (l = h+1; l < j; ++l) {
			up_val = up[h][l] * exp(-eL(i,j,h,l)/RT);
		}
	}
	up_val = up_val * exp(-Ed3(i,j,i+1)/RT);
	up_val = up_val + exp(-eH(i,j)/RT );
	up_val = up_val + exp(-eS(i,j)/RT ) * up[i+1][j-1];
	up_val = up_val + upm[i][j];
}
