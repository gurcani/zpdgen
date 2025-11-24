#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gpdf.h"
int main(int argc, char *argv[]){
  unsigned long int t,tc;
  complex res;
  int lx,ly,l;
  int n,m;
  complex *F;
  complex *za,om;
  complex za_min, za_max, dza;
  int numx,numy;
  char buf[50];
  int nms[]={1,0,1,2,3,0};
  za_min=-6.0-6.0i;
  za_max=6.0+6.0i;
  dza=0.1+0.1i;
  numx=creal(za_max-za_min)/creal(dza);
  numy=cimag(za_max-za_min)/cimag(dza);
  F=malloc(sizeof(complex)*numx*numy);
  za=malloc(sizeof(complex)*numx*numy);
  for(l=0;l<3;l++){
    n=nms[l*2];
    m=nms[l*2+1];
    printf("computing I_%i%i...\n",n,m);
    tc=clock();
    for(lx=0;lx<numx;lx++){
      for(ly=0;ly<numy;ly++){
	za[lx*numy+ly]=za_min+creal(dza)*lx+I*cimag(dza)*ly;
	F[lx*numy+ly]=gpdf_inm(za[lx*numy+ly],0.0,0.09,n,m);
      }
    }
    printf("\n%lu msecs cpu time\n",clock()-tc);
  }
  exit(0);
}
