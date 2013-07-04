#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
//#include "dwt97.cpp"
#define WAVELET_TYPE int16_t
#define TIMESCALE 1
#define WW2D 1
#define WWW 1
#define WWH 1
#define WW_SCALE 1
//#define WWTIME 1
#define PACK 1
//32768
void d3fwt(WAVELET_TYPE* x, int sdiv, int d, int h, int w);
void d3iwt(WAVELET_TYPE* x, int sdiv, int d, int h, int w);
void fwt53(WAVELET_TYPE* x,int n, int stride);
void iwt53(WAVELET_TYPE* x,int n, int stride);
void main(){
	FILE *inf = fopen("input.raw", "rb");
	fseek(inf, 0L, SEEK_END);
	//size_t fsize = ftell(inf);
	size_t fsize = 32*128*256;
	fseek(inf, 0L, SEEK_SET);
	uint8_t *idata = (uint8_t*)malloc(fsize);
	WAVELET_TYPE *sdata = (WAVELET_TYPE*)calloc(fsize,sizeof(WAVELET_TYPE));
	WAVELET_TYPE *tdata = (WAVELET_TYPE*)calloc(fsize,sizeof(WAVELET_TYPE));
	fread(idata, 1, fsize, inf);
	fclose(inf);
	int i;
	for(i=0;i<fsize;i++){ // convert to short
		sdata[i] = idata[i]-127;
		tdata[i] = sdata[i];
	}
	//arrayimg* ydatao = (arrayimg*)idata;
	for(i=0;i<5;i++){
		d3fwt(sdata, 1<<i, 32, 128, 256);
	}
	//d3iwt(sdata, 32, 128, 256);
#ifdef DIFF
	for(i=0;i<fsize;i++){ 
		sdata[i] -= tdata[i];
	}
#endif
	FILE *outf = fopen("output.raw", "wb");
	fwrite(sdata, sizeof(WAVELET_TYPE), fsize, outf);
	fclose(outf);
	free(idata);
	free(sdata);
	free(tdata);
}

void d3fwt(WAVELET_TYPE* x, int sdiv, int d, int h, int w){
	int i,j;
#ifdef WW2D	
  #ifdef WWW
  for(i = 0;i<d;i++){
  	#pragma omp parallel for
    for(j = 0;j<w;j++){
      fwt53(&x[(i*w*h)+j], h, w);
    }
  }
  #endif
	#ifdef WWH
	
	for(i = 0;i<d;i++){
    #pragma omp parallel for
		for(j = 0;j<h;j++){
			fwt53(&x[((i*h)+j)*w], w, 1);
		}
	}
  #endif
#endif
#ifdef WWTIME
	for(i = 0;i<h;i++){
		#pragma omp parallel for    
		for(j = 0;j<w;j++){
			//fprintf(stderr,"I:%dJ:%d\n",i,j);
			fwt53(&x[(i*w)+j], d, w*h);
		}
	}
#endif
}

void d3iwt(WAVELET_TYPE* x, int sdiv, int d, int h, int w){
	int i,j;
  #ifdef WWTIME
	for(i = 0;i<h;i++){
    #pragma omp parallel for
		for(j = 0;j<w;j++){
			iwt53(&x[(i*w)+j], d, w*h);
		}
	}
#endif
#ifdef WW2D	
  #ifdef WWH
	for(i = 0;i<d;i++){
    #pragma omp parallel for
		for(j = 0;j<h;j++){
			iwt53(&x[((i*h)+j)*w], w, 1);
		}
	}
  #endif
  #ifdef WWW
  for(i = 0;i<d;i++){
    #pragma omp parallel for
    for(j = 0;j<w;j++){
      iwt53(&x[(i*w*h)+j], h, w);
    }
  }
  #endif
#endif
}
#define DF 1
void fwt53(WAVELET_TYPE* x,int n, int stride)
{
  WAVELET_TYPE a;
  int i;
  // Predict 1
  a=-2;
  //#pragma omp parallel for
  for (i=1;i<n-1;i+=2)
  {
  	//printf("P~%02d\n", i);
    x[i*stride]+=(x[(i-1)*stride]+x[(i+1)*stride])/(a);
  } 
  // Update 1
  a=+4;
  x[0*stride]+=(2*x[1*stride])/(a);
  x[(n-1)*stride]+=(2*x[(n-2)*stride])/(a);

  //#pragma omp parallel for
  for (i=2;i<n-DF;i+=2)
  {
   // printf("U:%02d\n", i);
    x[i*stride]+=(x[(i-1)*stride]+x[(i+1)*stride])/(a);
  }
  // Scale
  a=WW_SCALE; //Don't know the scaling factor
  //#pragma omp parallel for
  for (i=0;i<n;i+=2)
  {
    x[i*stride]*=a;
  }
  //#pragma omp parallel for
  for (i=1;i<n;i+=2)
  {
    x[i*stride]/=a;

  }
#ifdef PACK
  // Pack
  WAVELET_TYPE *tempbank = tempbank=(WAVELET_TYPE *)malloc(n*sizeof(WAVELET_TYPE));
  //#pragma omp parallel for
  for (i=0;i<n;i++)
  {
    if (i%2==0) tempbank[i/2]=x[i*stride];
    else tempbank[(n/2)+(i/2)]=x[i*stride];
  }
  //#pragma omp parallel for
  for (i=0;i<n;i++) x[i*stride]=tempbank[i];
  free(tempbank);
#endif
}

void iwt53(WAVELET_TYPE* x,int n, int stride)
{
  WAVELET_TYPE a;
  int i;
#ifdef PACK
  WAVELET_TYPE *tempbank = tempbank=(WAVELET_TYPE *)malloc(n*sizeof(WAVELET_TYPE));
  // Unpack
  //#pragma omp parallel for
  for (i=0;i<n/2;i++)
  {
    tempbank[i*2]=x[i*stride];
    tempbank[(i*2)+1]=x[(i+(n/2))*stride];
  }
  for (i=0;i<n;i++) x[i*stride]=tempbank[i];
  free(tempbank);
#endif
  // Undo scale
  a=WW_SCALE; //Don't know the scaling factor
  //#pragma omp parallel for
  for (i=0;i<n;i+=2)
  {
    x[i*stride]/=a;
  }
  //#pragma omp parallel for
  for (i=1;i<n;i+=2)
  {
    x[i*stride]*=a;
  }  
  // Undo Update 1
  a=-4;
  //#pragma omp parallel for
  for (i=2;i<n-DF;i+=2)
  {
    x[i*stride]+=(x[(i-1)*stride]+x[(i+1)*stride])/(a);
  }
  x[0*stride]+=((2*x[1*stride]))/(a);
  x[(n-1)*stride]+=(2*x[(n-2)*stride])/(a);

  // Undo Predict 1
  a=2;
  //#pragma omp parallel for
  for (i=1;i<n-1;i+=2)
  {
    x[i*stride]+=(x[(i-1)*stride]+x[(i+1)*stride])/(a);
  } 
  
}
