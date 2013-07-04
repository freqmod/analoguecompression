#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
//#include "dwt97.cpp"
#define WAVELET_TYPE int16_t
#define TIMESCALE 1
#define WREC 5
#define WW_SCALE 1
//enable stuff?
#define WW2D 1
#define WWW 1
#define WWH 1
//#define DIFF 1
#define REVERSE 1
#define WWTIME 1
#define PACK 1
//#define ABSDATA 1
//#define SIMPLENOISE 1
//32768
void d3fwt(WAVELET_TYPE* x, int sdiv, int d, int h, int w);
void d3iwt(WAVELET_TYPE* x, int sdiv, int d, int h, int w);
void fwt53(WAVELET_TYPE* x,int n, int stride);
void iwt53(WAVELET_TYPE* x,int n, int stride);
void main(){
	FILE *inf = fopen("input.raw", "rb");
	fseek(inf, 0L, SEEK_END);
	//size_t fsize = ftell(inf);
	size_t fsize = 32*128*256*3;
	fseek(inf, 0L, SEEK_SET);
	uint8_t *idata = (uint8_t*)malloc(fsize);
	WAVELET_TYPE *sdata = (WAVELET_TYPE*)calloc(fsize,sizeof(WAVELET_TYPE));
	WAVELET_TYPE *tdata = (WAVELET_TYPE*)calloc(fsize,sizeof(WAVELET_TYPE));
	int r = fread(idata, 1, fsize, inf);
	fclose(inf);
	int i,j,k;
	for(i=0;i<fsize;i++){ // convert to short
		sdata[i] = idata[i]-127;
		tdata[i] = sdata[i];
	}
	//for(k=0;k<100;k++){
	for(j=0;j<3;j++){
		for(i=0;i<WREC;i++){
			d3fwt(&sdata[j*32*128*256], i, 32, 128, 256);
		}
	}
	#if SIMPLENOISE
	for(j=0;j<3;j++){
		for(i=0;i<32*128*256;i++){ 
			if(abs(sdata[i])<10){
				sdata[(j*32*128*256)+i]=0;
			}
		}
	}
	#endif
	#ifdef REVERSE
	for(j=0;j<3;j++){
		for(i=WREC-1;i>=0;i--){
			d3iwt(&sdata[j*32*128*256], i, 32, 128, 256);
		}
	}
	#endif

	//}
	#ifdef DIFF
	for(i=0;i<fsize;i++){ 
		sdata[i] -= tdata[i];
	}
	#endif
	#ifdef ABSDATA
	for(i=0;i<3*32*128*256;i++){ 
		sdata[i] = abs(sdata[i]);
	}
	#endif
	FILE *outf = fopen("output.raw", "wb");
	fwrite(sdata, sizeof(WAVELET_TYPE), fsize, outf);
	fclose(outf);
	free(idata);
	free(sdata);
	free(tdata);
}

void d3fwt(WAVELET_TYPE* x, int sdivi, int d, int h, int w){
	int i,j;
	#define DDIV (d>>sdivi)
	#define HDIV (h>>sdivi)
	#define WDIV (w>>sdivi)
#ifdef WW2D	
  #ifdef WWW
  for(i = 0;i<DDIV;i++){
  	#pragma omp parallel for
    for(j = 0;j<WDIV;j++){
      fwt53(&x[(i*w*h)+j], HDIV, w);
    }
  }
  #endif
	#ifdef WWH
	for(i = 0;i<DDIV;i++){
    #pragma omp parallel for
		for(j = 0;j<HDIV;j++){
			fwt53(&x[((i*h)+j)*w], WDIV, 1);
		}
	}
  #endif
#endif
#ifdef WWTIME
	for(i = 0;i<HDIV;i++){
		#pragma omp parallel for    
		for(j = 0;j<WDIV;j++){
			//fprintf(stderr,"I:%dJ:%d\n",i,j);
			fwt53(&x[(i*w)+j], DDIV, w*h);
		}
	}
#endif
}

void d3iwt(WAVELET_TYPE* x, int sdivi, int d, int h, int w){
	int i,j,k;
	#define DDIV (d>>sdivi)
	#define HDIV (h>>sdivi)
	#define WDIV (w>>sdivi)

  #ifdef WWTIME
	for(i = 0;i<HDIV;i++){
    #pragma omp parallel for
		for(j = 0;j<WDIV;j++){
			iwt53(&x[(i*w)+j], DDIV, w*h);
		}
	}
#endif
#ifdef WW2D	
  #ifdef WWH
	for(i = 0;i<DDIV;i++){
    #pragma omp parallel for
		for(j = 0;j<HDIV;j++){
			iwt53(&x[((i*h)+j)*w], WDIV, 1);
		}
	}
  #endif
  #ifdef WWW
  for(i = 0;i<DDIV;i++){
    #pragma omp parallel for
    for(j = 0;j<WDIV;j++){
      iwt53(&x[(i*w*h)+j], HDIV, w);
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
  for (i=1;i<n-1;i+=2)
  {
    x[i*stride]+=(x[(i-1)*stride]+x[(i+1)*stride])/(a);
  } 
  // Update 1
  x[(n-1)*stride]+=(2*x[(n-2)*stride])/(a);
  a=+4;
  for (i=2;i<n-DF;i+=2)
  {
    x[i*stride]+=(x[(i-1)*stride]+x[(i+1)*stride])/(a);
  }
  x[0*stride]+=(2*x[1*stride])/(a);
#ifdef DOSCALE
  // Scale
  a=WW_SCALE; //Don't know the scaling factor
  for (i=0;i<n;i+=2)
  {
    x[i*stride]*=a;
  }
  for (i=1;i<n;i+=2)
  {
    x[i*stride]/=a;

  }
  #endif
#ifdef PACK
  // Pack
  WAVELET_TYPE *tempbank = tempbank=(WAVELET_TYPE *)malloc(n*sizeof(WAVELET_TYPE));
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
  for (i=0;i<n/2;i++)
  {
    tempbank[i*2]=x[i*stride];
    tempbank[(i*2)+1]=x[(i+(n/2))*stride];
  }
  for (i=0;i<n;i++) x[i*stride]=tempbank[i];
  free(tempbank);
#endif
  // Undo scale
  #ifdef DOSCALE
  a=WW_SCALE; //Don't know the scaling factor
  for (i=0;i<n;i+=2)
  {
    x[i*stride]/=a;
  }
  for (i=1;i<n;i+=2)
  {
    x[i*stride]*=a;
  }  
 #endif
  // Undo Update 1
  a=-4;
  for (i=2;i<n-DF;i+=2)
  {
    x[i*stride]+=(x[(i-1)*stride]+x[(i+1)*stride])/(a);
  }
  x[0*stride]+=((2*x[1*stride]))/(a);  
  // Undo Predict 1
  a=2;
  x[(n-1)*stride]+=(2*x[(n-2)*stride])/(a);
  for (i=1;i<n-1;i+=2)
  {
    x[i*stride]+=(x[(i-1)*stride]+x[(i+1)*stride])/(a);
  }   
}
