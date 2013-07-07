#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
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
#define NUMFRAMES 32
//#define ABSDATA 1
//#define SIMPLENOISE 1
//NUMFRAMES768
void d3fwt(WAVELET_TYPE* x, int sdiv, int d, int h, int w);
void d3iwt(WAVELET_TYPE* x, int sdiv, int d, int h, int w);
void fwt53(WAVELET_TYPE* x,int n, int stride);
void iwt53(WAVELET_TYPE* x,int n, int stride);
inline int getoffset(int level, int part, int recursion, int d, int h, int w){
	int res = 0;
	if((part & 1)==1){
		res += w >> (recursion - level);
	}
	if((part & 2)==2){
		res += (h >> (recursion - level))*w;
	}
	if((part & 4)==4){
		res += (d >> (recursion - level))*w*h;
	}
	return res;
}
inline int getlength(int level, int recursion, int original){
	return (original >> (recursion - level));
}

int main(){
	FILE *inf = fopen("input.raw", "rb");
	fseek(inf, 0L, SEEK_END);
	//size_t fsize = ftell(inf);
	size_t fsize = NUMFRAMES*128*256*3;
	fseek(inf, 0L, SEEK_SET);
	uint8_t *idata = (uint8_t*)malloc(fsize);
	WAVELET_TYPE *sdata = (WAVELET_TYPE*)calloc(fsize,sizeof(WAVELET_TYPE));
	WAVELET_TYPE *tdata = (WAVELET_TYPE*)calloc(fsize,sizeof(WAVELET_TYPE));
	uint8_t *mdata = (uint8_t*)calloc(fsize,sizeof(uint8_t));//bool
	int r = fread(idata, 1, fsize, inf);
	fclose(inf);
	int c,i,j,k,l,p, limit, count, reso, rsize, cutoff;
	for(i=0;i<fsize;i++){ // convert to short
		sdata[i] = idata[i]-127;
		tdata[i] = sdata[i];
	}
	//for(k=0;k<100;k++){
	for(j=0;j<3;j++){
		for(i=0;i<WREC;i++){
			d3fwt(&sdata[j*NUMFRAMES*128*256], i, NUMFRAMES, 128, 256);
		}
	}
	#if SIMPLENOISE
	for(j=0;j<3;j++){
		for(i=0;i<NUMFRAMES*128*256;i++){ 
			if(abs(sdata[i])<10){
				sdata[(j*NUMFRAMES*128*256)+i]=0;
			}
		}
	}
	#endif
	//remove unwanted colors
	for(c = 1;c<3;c++){
		for(l = 3;l<WREC;l++){
			for(p = 1;p<8;p++){
				int o = getoffset(l,p,WREC,NUMFRAMES,128,256);
				for(i = 0;i<getlength(l,WREC, NUMFRAMES);i++){
					for(j = 0;j<getlength(l,WREC, 128);j++){
						for(k = 0;k<getlength(l,WREC, 256);k++){
							sdata[(c*32*128*256)+o+((i*128+j)*256)+k] = 0;
						}
					}
				}
			}
		}
	}
	//SPIHT compression init
	uint8_t fracts[3][10] = {{1,1, 1,1, 2,3, 0,0, 0,0},
						     {1,8, 0,0, 2,3, 0,0, 0,0}, 
						     {1,8, 0,0, 2,3, 0,0, 0,0}};
#define COUNTS_SIZE (1<<((8*sizeof(WAVELET_TYPE)-1)))
	WAVELET_TYPE counts[COUNTS_SIZE]; 

	rsize = 0; //size of packed result
	reso = 0;
	for(c = 0;c<3;c++){
		for(l = 0;l<WREC;l++){
			for(p = 1;p<8;p++){
				if(fracts[c][l*2+1]>0){
					rsize += (getlength(l,WREC, NUMFRAMES)*getlength(l,WREC, 128)*getlength(l,WREC, 256)*fracts[c][l*2])/fracts[c][l*2+1];
				}
			}
		}
	}
	int preoffset = 0;
	//SPIHT compression, not correct yet
	WAVELET_TYPE *rdata = (WAVELET_TYPE*)calloc(rsize,sizeof(WAVELET_TYPE));
	for(c = 0;c<3;c++){
		for(l = 0;l<WREC;l++){
			for(p = 1;p<8;p++){
				int o = getoffset(l,p,WREC,NUMFRAMES,128,256);
				if(fracts[c][l*2+1]>0){
					//counting sort 
					memset(counts, 0, COUNTS_SIZE*sizeof(WAVELET_TYPE));
					for(i = 0;i<getlength(l,WREC, NUMFRAMES);i++){
						for(j = 0;j<getlength(l,WREC, 128);j++){
							for(k = 0;k<getlength(l,WREC, 256);k++){
								counts[abs(sdata[(c*32*128*256)+o+((i*128+j)*256)+k])]++;
								//if(sdata[(c*32*128*256)+o+((i*128+j)*256)+k]!=0)
									//printf("CO:%d=%d\n", sdata[(c*32*128*256)+o+((i*128+j)*256)+k], counts[abs(sdata[(c*32*128*256)+o+((i*128+j)*256)+k])]);
							}
						}
					}
					limit = (getlength(l,WREC, NUMFRAMES)*getlength(l,WREC, 128)*getlength(l,WREC, 256)*fracts[c][l*2])/fracts[c][l*2+1];
					count = 0;
					for(i = COUNTS_SIZE-1; i>=0; i--){
						count += counts[i];
						//printf("CA %i: %d\n", i, counts[i]);
						if(count > limit){
							printf("Breakatlimit %d: %d\n", i, counts[i]);
							count -= counts[i];
							break;
						}
					}
					cutoff = i;
					preoffset = reso;
					//Pack data into results
					for(i = 0;i<getlength(l,WREC, NUMFRAMES);i++){
						for(j = 0;j<getlength(l,WREC, 128);j++){
							for(k = 0;k<getlength(l,WREC, 256);k++){
								if(abs(sdata[(c*32*128*256)+o+((i*128+j)*256)+k]) > cutoff){
									mdata[(c*32*128*256)+o+((i*128+j)*256)+k] = 1;
									//printf("Testadd%d\n", reso);
									rdata[reso++] = sdata[(c*32*128*256)+o+((i*128+j)*256)+k];
									//printf("TestAdd%d\n", reso);
								} else if (abs(sdata[(c*32*128*256)+o+((i*128+j)*256)+k]) == limit &&
									count < limit){
									//this should idealy be distruted eavenly over the picture
									rdata[reso++] = sdata[(c*32*128*256)+o+((i*128+j)*256)+k];
									mdata[(c*32*128*256)+o+((i*128+j)*256)+k] = 1;
									//printf("TestBDD%d\n", reso);
									count++;
								} else{
									mdata[(c*32*128*256)+o+((i*128+j)*256)+k] = 0;
								}
							}
						}
					}
					//Add extra data to results if neccesary
					//printf("Offset after pack:Reso:%d Rem: %d (L:%d,P:%d;C:%d) [%d<%d]=%d*%d/%d\n", reso, rsize-reso, l, p, cutoff, count,limit, getlength(l,WREC, NUMFRAMES)*getlength(l,WREC, 128)*getlength(l,WREC, 256),fracts[c][l*2],fracts[c][l*2+1]);
					for(i=count;i<limit;i++){
						//printf("Padding%d\n",i);
						rdata[reso++] = 0;
					}
					if(count<limit){
						printf("Padded: R:%d, (C:%d,L:%d,P:%d) D:%d, retrieved %d (co:%d)\n", reso, c, l, p,  limit-count, count,cutoff);
					}
					if(reso-preoffset != limit){
						printf("Wrong progress %d!=%d (%d,%d) %d\n", reso-preoffset, limit, l,p, limit-count);	
					}
				}
			}
		}
	}
	if(reso != rsize){
		printf("Wrong progress at end %d\n", rsize-reso);
	}
	//SPIHT decompression

	printf("Offset:%d\n", getoffset(1,4,WREC,32,128,256));
	#ifdef REVERSE
	for(j=0;j<3;j++){
		for(i=WREC-1;i>=0;i--){
			d3iwt(&sdata[j*NUMFRAMES*128*256], i, NUMFRAMES, 128, 256);
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
	for(i=0;i<3*NUMFRAMES*128*256;i++){ 
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
