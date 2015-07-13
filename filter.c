/**
 *
 * ===========================================================================
 *
 * File: filter.c
 *
 * Abstract:
 *       generic filter routines.
 *
 * Henry A. Leinhos, June 2, 2009
 * Based on AMS fir.c code, 1997-2006
 *
 
 * ===========================================================================
 *
 *
 * ===========================================================================
 */


#include <math.h>
#include <stdlib.h>  /* for free() */
#include <stdio.h>
#include <string.h>
#include "filter.h"

/* usually, math.h provides INFINITY - just in case */
#ifndef INFINITY
# define INFINITY (HUGE_VAL)
#endif

#ifndef PI
# define PI 3.1415
#endif

#define DEBUG(X) X
#define DEBUG0 0

/* Filter Queue utility functions */

/* initialize a new queue 
 *
 * q - filter Queue pointer to Queue structure
 * buf - buffer pointer (if NULL, malloc new buf)
 * size - size of buffer in samples
 * 
 * returns pointer to first available buffer location
 * (TODO: constrain size of buf to power of 2 for FAST_QUEUE)
 */

Float *queue_init(Queue *q, Float *buf,int size) {

  unsigned int i,flag=0;

  /* check for required allocations, exit with NULL
   * if any failures
   */

#if FAST_QUEUE
  /* find next power of two larger than size */
  for(i=1;(i<=size)&&(i<(0x01<<31));i=i<<1); 
  i=i>>1; /* our buffer should be the next smallest power 
	   * to fit */
  if (i==(0x01<<31)) return(NULL);
#else
  i=size;
#endif 

  if(NULL == q) {
    flag=1;
    q = malloc(sizeof(Queue));
  }
  if(NULL == q) return(NULL); 
  if(NULL == buf) buf = malloc(sizeof(Float)*i);
  if(NULL == buf) {
    if (flag) free(q);
    return(NULL); 
  }
  
  q->buf = buf;  /* set buffer pointer */
  q->rp = 0;
  q->wp = 0;
  q->N = i;

#if FAST_QUEUE
  q->mask = i-1;
#else
  q->mask = ~0x0;
#endif

  DEBUG(printf("[queue_init] queue at 0x%x\n",q));
  DEBUG(printf("[queue_init] buf at 0x%x\n",q->buf));
  DEBUG(printf("[queue_init] rp,wp = %d, %d\n",q->rp,q->wp));
  DEBUG(printf("[queue_init] size = %d,actual size = %d\n",size,q->N));

  for (i=0;i<q->N;i++) q->buf[i] = 0.0;  /* clear buffer */

  
  return(q->buf);
}

/* push N data samples onto queue, return the number
 * of samples pushed
 */

int push_queue(Queue *q,Float *data,int N) {

  int n;

  /* add data to circ buffer up until buffer full */
  for (n=0;(n<N)&&(((q->wp+1)&q->mask)!=q->rp);n++) {
    q->buf[q->wp] = data[n];
    q->wp = (q->wp+1) & q->mask;  /* binary mask implements circular buffer */
  }

  return(n);
}

/* pop at most N data samples off queue, return the number
 * of samples poped
 */

int pop_queue(Queue *q,Float *data,int N) {

  int n;

  if(data != NULL) { /* if there is a destination, copy data */
    /* copy data off the circ buffer up until buff empty */
    for (n=0;n<N && (q->rp!=q->wp);n++) {
      data[n] = q->buf[q->rp];
      q->rp = (q->rp+1) & q->mask;
    }
  } else { /* otherwise, just update the rp */
    for (n=0;n<N && (q->rp!=q->wp);n++) {
      q->rp = (q->rp+1) & q->mask;
    }
   
  }
  return(n);
}

Filter *filter_init(  Filter *filt, /* actual filter structure */
		      Float *ff,   /* feedback filter coeffs */
		      Float *fb,   /* feedback filter coeffs */
		      Queue *q,    /* filter state history storage queue */
		      int N,       /* feedforward filter order */
		      int M,       /* feedback filter order */ 
		      int nch     /* number of channels in the data */) {


  int alloc_size = 0;
  int alloc_index = 0;

  /* disable automatic allocation  for now */
#if 0
  /* check for unalloc'd objects */
  if (NULL == filt) alloc_size += sizeof(Filter);
  if (NULL == ff)   alloc_size += sizeof(Float)*N;
  if (NULL == fb)   alloc_size += sizeof(Float)*M;
  if (NULL == q)    alloc_size += sizeof(Queue);

  if(NULL == filt) {
    if(NULL == (filt = malloc(alloc_size))) return (NULL);
    filt->ptr = (void *) filt;  /* save alloc pointer */
  } else if (alloc_size > 0)
    if(NULL == (filt->ptr=malloc(alloc_size))) return(NULL);
#else
  
  if ((filt == NULL) ||
      (ff == NULL)   ||
      (fb == NULL)   ||
      (q == NULL) ) return (NULL);
#endif  

  /* assuming all the objects are properly allocated,
   * set up the filter structure and initialize objects
   */

  filt->ff = ff;  /* point to feed forward filter coeffs */
  filt->fb = fb;  /* feed back filer coeffs */
  filt->q = q;    /* filter state history */
  filt->N = N;    /* feedforward filter order */
  filt->M = M;    /* feedback filter order */
  filt->nch = nch; /* number of channels to filter */

  DEBUG(printf("[filter_init] filt struct at 0x%x\n",filt));
  DEBUG(printf("[filter_init] queue struct at 0x%x\n",filt->q));
  DEBUG(printf("[filter_init] ff buf at 0x%x\n",filt->ff));
  DEBUG(printf("[filter_init] fb buf at 0x%x\n",filt->fb));
  DEBUG(printf("[filter_init] ff/fb size is %d/%d\n",filt->N,filt->M));
  DEBUG(printf("[filter_init] nch is %d\n",filt->nch));

  return(filt);

}


/* Filter Functions */

/* Direct Form II IIR filter routine */

/* 
 *             v(n)
 * x(n) --->+---------b0-->+---> y(n)
 *          ^      |       ^
 *          |     z-1      |
 *          | -    |       |
 *          +<--a1----b1-->+
 *          ^      |       ^
 *          |     z-1      |
 *          | -    |       |
 *          +<--a2----b2-->+
 *          ^      |       ^
 *          |     z-1      |
 *          |      |       |
 *
 *             etc...
 *
 */

/* multi-channel organization:
 *
 * the circular buffer contains nch data samples at-a-time so that
 * each channel output is calculated as a block
 *
 *  
 */
/* note:  this version directly manipulates the circular buffer queue */

Float *filter(Filter *f, Float *in, Float *out, int nsamps) {
  int i,n,m,k,nch;
  int v_ind1,v_ind2;
  int size;
  Float *v,v_n=0.0,a,b;

  nch = f->nch;  /* get number of channels */
  size = f->q->N;   /* get buffer size */

  for (i=0;i<nsamps;i++) {
    for(k=0;k<nch;k++) { /* start with input value */
      v_ind2 = f->q->wp-nch+k;                  /* point to new history sample */
      v_ind2 &= f->q->mask;
      f->q->buf[v_ind2] = in[i*nch+k];     /* add input value */
    }
    for (m=1;m<f->M;m++) { /* calc next intermediate value */
      a = f->fb[m];  /* get the mth filter coeff */
      for(k=0;k<nch;k++) { /* for each kth channel */
	v_ind1 = f->q->wp-(m+1)*nch+k;  /* index into data history */
	v_ind1 &= f->q->mask;
	v_ind2 = f->q->wp-nch+k;        /* new history sample */
	v_ind2 &= f->q->mask;
	v_n = f->q->buf[v_ind1];        /* get value from history */
                                              /* (note the binary mask trick) */
	f->q->buf[v_ind2] -= a*v_n;      /* accumulate fb values */
      }
    }
 
    /* calc contribution by feedforward filter */
    for(k=0;k<nch;k++) {
      v_ind2 = f->q->wp-nch+k;                  /* point to new history sample */
      v_ind2 &= f->q->mask;
      out[i*nch+k] = f->ff[0] * f->q->buf[v_ind2] ;
    }
    for(n=1;n<f->M;n++) { /* apply ff filter */
      b = f->ff[n];
      for(k=0;k<nch;k++) {
	v_ind1 = f->q->wp-(n+1)*nch+k;  /* index back into history */
	v_n = f->q->buf[(v_ind1) & f->q->mask]; /* get value from history */
	out[i*nch+k] += b * v_n;
      }
    }
 
    for(k=0;k<nch;k++) { /* update history write pointer */
      f->q->wp =  (f->q->wp+1) & f->q->mask;
      /* bump read pointer if buffer is full */
      if(f->q->wp == f->q->rp) f->q->rp =  (f->q->rp+1) & f->q->mask;
    }
  }

  return(out);
}
/* Filter Design Functions */

/* fir filter design using windowing method (and Hamming window) */

Float *fir1(Float *coeffs,Float bw, Float fc, int size) {
  
  int k;
  Float sum=0.0;

  for (k=0;k<size;k++) {
    coeffs[k] = sin(2*PI*bw/2*(k-size/2+0.00001))/(2*PI*bw/2*(k-size/2+0.00001)); /* sinc function */
    coeffs[k] *= 0.54-0.46*cos(2*PI*k/size);  /* Hamming window */
#if 1
    sum+= coeffs[k];     /* normalize to dc gain */
#else
    sum+= coeffs[k]*coeffs[k];  /* normalize to power */
#endif    
    coeffs[k] *= cos(2*PI*k*fc);
  }

#if DEBUG0
  printf("coeffs = [\n");
#endif  
  for (k=0;k<size;k++){
    coeffs[k]/=sum;  /* normalize filter */
#if DEBUG0
    printf("%f;\n",coeffs[k]);
#endif
  }
#if DEBUG0
  printf("];\n");
#endif
  return(coeffs);
}
/* butterworth filter design routine */
/* returns 2x(size+1) polynomial coeffs
 * for a discrete-time Butterworth filter of order size
 *
 * first size+1 elements are the numerator polynomial coeffs,
 * next size+1 elements are the denominator polynomial coeffs.
 */

/* outstanding issues: 
 * 
 * 
 * filter center frequency is forced to be zero for now (need to rotate poles/zeros about origin)
 * filter bw is 2x that used in Octave/Matlab
 *
 * HAL 6/1/06
 */

Float *butter(Float *coeffs,float bw, float fc, int size) {

  Float T;
  Float denom;
  Complex roots[MAX_BUTTERWORTH_SIZE*2],rootz[MAX_BUTTERWORTH_SIZE*2];
  Complex poly[(MAX_BUTTERWORTH_SIZE+1)*2],temp[MAX_BUTTERWORTH_SIZE+1],cur,prev,gain;
  int k,j,N;


  /* prewarp frequency specification to take into account bilinear trasformation */
  T = 2.0;
  bw = 2/T*tan(PI*bw/T);

  DEBUG(printf("[FIR_butter] warped cuttoff freq is %f\n",bw));
  /* analog butterworth filters have roots at (wc^2)*(-1)^(1/n) */
  /* or (wc)*exp(j*PI/2)*exp(j*PI/N/2)*exp(j*PI/N*k) with k=0:N-1 */

  /* allocate space for roots 2*size (zeros and poles)*/
 
  if(size>MAX_BUTTERWORTH_SIZE) return(NULL);

  for(k=0;k<4*size;k++) { /* zero out roots */
    roots[k].real = 0.0;
    roots[k].imag = 0.0;
  }

  for(k=0;k<size;k++) { /* calc S-plane zeros/poles */
    roots[k].real = -INFINITY;  /* zeros at -Inf */
    roots[k].imag = -INFINITY;
    DEBUG(printf("[FIR_butter]   S-zeros(%d)= %f + %fj\n",k,roots[k].real,roots[k].imag));

    roots[size+k].real =  bw*cos((size+(2*k+1))*PI/size/2); /* pole locations */
    roots[size+k].imag =  bw*sin((size+(2*k+1))*PI/size/2);
    DEBUG(printf("[FIR_butter]   S-poles(%d)= %f + %fj\n",k,roots[k+size].real,roots[k+size].imag));

  }

  /* transform S-Plane roots to Z-plane roots via bilinear transform */

  /* from Octave-forge bilinear.m:  
   * ## ----------------  -------------------------  ------------------------
   * ## Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
   * ##      2 z-1        pole: -1                   zero: -1
   * ## S -> - ---        gain: T/(2-xT) (NO?)       gain: (2-xT)/T (NO?)
   * ##      T z+1             -> (2/T+x)               ->  1/(2/T+x)
   * ## ----------------  -------------------------  ------------------------
   */

  /* the real part of the poles will be at (4-T^2)/((4+T^2)+cos(Xk))
   * and the imaginary part of the poles will be at 4*T*sin(Xk)/(.) 
   *
   * where Xk = (size+(2*k+1))*PI/size/2
   */
  /* note: all the zeros at -Inf map to z=-1 */

  /* this should really be implemented via complex arithmetic functions */

  gain.real = 1.0;  /* initialize transform gain */
  gain.imag = 0.0;
  T *= bw;  /* scale sample time for desired bandwidth */
  for (k=0;k<size;k++) {
    /* new zero locations */
    rootz[k].real = -1.0;
    rootz[k].imag = 0.0;

    /* new pole locations */
    denom = (4+T*T)-4*T*cos((size+(2*k+1))*PI/size/2);
    rootz[k+size].real = (4-T*T)/denom;
    rootz[k+size].imag = 4*T*sin((size+(2*k+1))*PI/size/2)/denom;
    /* multilply by zero gain */
    DEBUG(printf("[FIR_butter]   zero gain(%d)= %f + j*%f\n",k,(2*bw/T-rootz[k].real),-rootz[k].imag));

    DEBUG(printf("[FIR_butter]   pole gain(%d)=1/(%f + j*%f)\n",k,(2*bw/T-rootz[k+size].real),-rootz[k+size].imag));

    cur.real = 2*bw/T-rootz[k].real;
    cur.imag = -rootz[k].imag;

    denom = ((2*bw/T)-rootz[k+size].real)*((2*bw/T)-rootz[k+size].real)+(rootz[k+size].imag*rootz[k+size].imag);
    
    /* gain due to zero = (2/T)+z */
    /*    gain.real = (gain.real*(2*bw/T-rootz[k].real)-gain.imag*(-rootz[k].imag)); */
    /*    gain.imag = (gain.real*(-rootz[k].imag)+gain.imag*(2*bw/T-rootz[k].real));  */
    complex_mult(&gain,&cur,&gain);

    cur.real = 2*bw/T-rootz[k+size].real;
    cur.imag = -rootz[k+size].imag;
    DEBUG(printf("[FIR_butter]   gain(%d)= %f + %fj\n",k,gain.real,gain.imag));
    
    /* gain due to pole 1/[(2/T)+p]*/
    /*   gain.real = (gain.real*(2*bw/T-rootz[k+size].real)+gain.imag*(-rootz[k+size].imag))/denom; */
    /*    gain.imag = (-gain.real*(-rootz[k+size].imag)+gain.imag*(2*bw/T-rootz[k+size].real))/denom; */

    complex_div(&gain,&cur,&gain);
    
    DEBUG(printf("[FIR_butter]   gain(%d)= %f + %fj\n",k,gain.real,gain.imag));

    DEBUG(printf("[FIR_butter]   Z-zeros(%d)= %f + %fj\n",k,rootz[k].real,rootz[k].imag));
    DEBUG(printf("[FIR_butter]   Z-poles(%d)= %f + %fj\n",k,rootz[k+size].real,rootz[k+size].imag));
  }

  N = size+1;
  for(k=0;k<(size+1);k++) { /* clear output polynomial */
    poly[k].real = 0.0;  /* numerator */
    poly[k].imag = 0.0;
    poly[k+N].real = 0.0;  /* denominator */
    poly[k+N].imag = 0.0;
  }
  poly[0].real = 1.0/gain.real; /* numerator */
  poly[0].imag = 0.0; 
  poly[N].real = 1.0; /* denominator */
  poly[N].imag = 0.0;

  /*    y(2:(j+1)) = y(2:(j+1)) - v(j) .* y(1:j); */

  /* do numerator first */
  for(j=0;j<(size);j++) {
    /* copy current poly to temp space */
    memcpy(temp,poly,(size+1)*sizeof(Complex));
    for (k=1;k<=(j+1);k++) {
      /* y(k)=y(k) + z(j)*y(k-1) */
      cur.real = temp[k].real;
      cur.imag = temp[k].imag;
      prev.real = temp[k-1].real;
      prev.imag = temp[k-1].imag;

      DEBUG(printf("[FIR_butter] (%d) zero(%d) = %f + %fj\n",k,j,rootz[j].real,rootz[j].imag));
      poly[k].real = cur.real - (rootz[j].real*prev.real - rootz[j].imag*prev.imag);
      poly[k].imag = cur.imag - (rootz[j].real*prev.imag + rootz[j].imag*prev.real);
      DEBUG(printf("[FIR_butter]      num(%d) = %f + %fj\n",k,poly[k].real,poly[k].imag));
    }

  }
  
  /* now do denominator */
  for(j=0;j<(size);j++) {
    /* copy current poly to temp space */
    memcpy(temp,&(poly[size+1]),(size+1)*sizeof(Complex));
     for (k=1;k<=(j+1);k++) {
      /* y(k)=y(k) + z(j)*y(k-1) */
      cur.real = temp[k].real;
      cur.imag = temp[k].imag;
      prev.real = temp[k-1].real;
      prev.imag = temp[k-1].imag;

      DEBUG(printf("[FIR_butter] (%d) pole(%d) = %f + %fj\n",k,j,rootz[j+size].real,rootz[j+size].imag));
      poly[k+N].real = cur.real - (rootz[j+size].real*prev.real - rootz[j+size].imag*prev.imag);
      poly[k+N].imag = cur.imag - (rootz[j+size].real*prev.imag + rootz[j+size].imag*prev.real);
      DEBUG(printf("[FIR_butter]      denom(%d) = %f + %fj\n",k,poly[k+N].real,poly[k+N].imag));
    }

  }

  /* check for any imagnary parts greater than 0.01 */
  for (k=0;k<(size+1)*2;k++) if((poly[k].imag*poly[k].imag)>0.0001) return(NULL);

  /* copy real part of num/denom polynomials to output Float vector */

  for (k=0;k<(size+1)*2;k++) coeffs[k] = poly[k].real;
  return(coeffs);
} /* end butter() */


/* besself - Bessel filter design 
 *
 */

/* 
 * this function  needs a root finder, followed by
 * a polynomial multiplication, or a direct bilinear
 * transform of polynomial coeffs
 *
 * poly multiplicatoin is simply convolution of coeffs
 * so that we will need a conv function implementation
 *
 */

/* Implementation */

/* bilinear mapping is S-> 2 z-1
 *                         - ---
 *                         T z+1
 * such that
 *
 *                k
 *   k     /2 z-1\
 *  S   -> |- ---|
 *         \T z+1/
 *
 * mult Num and Denom by (z+1)^N gives numerator coeffs of
 *
 * a_k *(2/T)^k * (z+1)^(N-k) * (z+1)^k
 *
 * so we need to find powers of polynomials (z+/-1)^k
 */




/*% Octave/Matlab code
 *
 * [b,a] = besself(n,Wo)
 *
 * TF = B_n(0)/B_n(s/Wo)
 *
 * generate bessel polynomial
 *
 * Wo = w/Ts;  % (fo/fs) * fs
 *
 * l = (N:-1:0); % reverse bessel polynomial;
 *
 * %yn =(1/(2*Wo)).^k .* factorial(N+k)./(factorial(N-k).*factorial(k));
 * yn = (1./(Wo.^l)) .* (1./(2.^(N-l))) .* factorial(2*N-l)./ (factorial(N-l).*factorial(l))
 * y0 = (1./(2.^(N))) .* factorial(2*N)./ (factorial(N).*factorial(0));
 *
 *  get roots of the of the polynomials
 *
 * sp = roots(yn);
 * [zz zp zg] = bilinear(zeros(1,N),sp,y0,Ts)
 *
 * [b a] = bilinear(y0,yn,Ts);
 *  transform roots to z domain
 *
 */


/* calc complex out=in1*in2 */
Complex *complex_mult(Complex *in1,Complex *in2, Complex *out) {

  Complex temp1, temp2;
  temp1.real = in1->real;
  temp1.imag = in1->imag;
  temp2.real = in2->real;
  temp2.imag = in2->imag;

  out->real = temp1.real*temp2.real-temp1.imag*temp2.imag;
  out->imag = temp1.real*temp2.imag + temp1.imag*temp2.real;

  return(out);
}

/* calc complex out=in1/in2
 * 
 * 1/(c+jd) = (c-jd)/(c^2 + d^2)\
 * 
 * so that (a+jb)/(c+jd) = (ac+bd) + j(-ad+bc)
 */
Complex *complex_div(Complex *in1,Complex *in2, Complex *out) {  

  Complex temp1, temp2;
  Float denom;

  temp1.real = in1->real;
  temp1.imag = in1->imag;
  temp2.real = in2->real;
  temp2.imag = in2->imag;

  denom = temp2.real*temp2.real + temp2.imag*temp2.imag;

  out->real = (temp1.real*temp2.real + temp1.imag*temp2.imag)/denom;
  out->imag = (-(temp1.real*temp2.imag)+(temp1.imag*temp2.real))/denom;

  return(out);
}

/* linearly interpolate value from (Nx2) table X,Y for input 
 * value Xi
 */
Float interp1(double *X, double *Y, double Xi, int N) {

  int k=0;
  Float val=0.0;

  if(!X || !Y) return 0.0;

  /* find nearest neighbors for XI */
  
  while((X[k]<Xi) && (k<N)) k++;

  /* check for input value too low */
  if(k==0 && N>1) {
    /* interpolate (backward) based on the slope of the
     * first two points
     */
    val = Y[0] - (X[0]-Xi) * (Y[1]-Y[0])/(X[1]-X[0]);
    return val;
  }
  /* check for input value too high */
  if(k==N && N>1) {
    /* interpolate (forward) base on the solop of the 
     * last two points 
     */
    val = Y[N-1] + (Xi-X[N-1]) * (Y[N-1]-Y[N-2])/(X[N-1]-X[N-2]);
    return val;
  }

  /* otherwise, interpolate between data points */

  val = Y[k-1] + (Xi-X[k-1]) * (Y[k]-Y[k-1])/(X[k]-X[k-1]);
  return val;
}

#if defined(TEST_FILTER)


#include <stdio.h>
#include "matrix.h"
#include "model.h"
#include "filter.h"
#include "string.h"

#include "autopilot.h"
/* run a simple test for this file  */
/* compile this with 
 * gcc -g -DTEST_FILTER -lm -lautopilot filter.c -o filter 
 */

/* TODO:  implement some sort of unit test */

#define FILTER_ORDER 4
#define FB_SIZE 4
#define NCH 6
#define TEST_SIZE (32*NCH)+1

main()

{
  Filter filt,*f;
  Float buf[TEST_SIZE];
  Float coeffs[(FILTER_ORDER+1)*2];
  Float *result;
  static Float in[256*NCH],out[256*NCH];

  int N=FILTER_ORDER+1,M=FILTER_ORDER+1;
  int nch = NCH;
  int n,k,i;

  Queue q;

  Float X[]={-1.0,3.0,4.5,10.1,25.0};
  Float Y[]={1.0,2.0,-1.0,-2.5, 10.0};
  Float Xi = 3.1;
  Float Yo = 0.0;
  int len = sizeof(X)/sizeof(Float);

  /* initialize multichannel input test data */
  for (k=0;k<nch;k++)
    in[k]=1.0;out[k] = 0.0;
  for(n=1;n<256;n++) {
    for(k=0;k<nch;k++) {
      out[n*nch+k] = in[n*nch+k] = 1.0*sin(k*0.5/400);
    }
  }

  for (i=0;i<1;i++) {

    printf("Initializing Queue\n");

    if((result=queue_init(&q,buf,TEST_SIZE))!=NULL)
      printf("success!\n");

    printf("Calculating filter coeffs\n");

    if((result=butter(coeffs,0.05,0.0,FILTER_ORDER))!=NULL){
      printf("FF/FB coeffs = [\n");
      for(n=0;n<N;n++) printf("     %2.6f  %2.6f;\n",coeffs[n],coeffs[n+N]);
      printf("  ]\n");
    } else {
      printf("Failed\n");
      return;
    }
    printf("Initializing Filter\n");
  
    if((f=filter_init(&filt,coeffs,&coeffs[FILTER_ORDER+1],&q,N,M,nch))!=NULL)
      printf("success!\n");

    printf("Running filter\n");
    
    if((result = filter(&filt,&in[256/2*i],out,256/2))!=NULL) {
      printf("success\n");
      printf("in/out = [\n");
      for(n=256/2*i;n<256/2*(i+1);n++) {
	printf("   ");
	for(k=0;k<nch;k++) 
	  printf("%2.4f, ",in[n*nch*(i+1)+k]);
	printf("  ");
	for(k=0;k<nch;k++) 
	  printf("%2.4f, ",out[n*nch*(i+1)+k]);
	printf("\n");
      }
      printf("   ];\n");
    }
  }  /* end i loop */


  /* test interp1 function too */

  printf("testing interp1()\n");
  
  printf("Y = [ ");
  for (k=0;k<len;k++) printf("%2.5f ", Y[k]);
  printf("]\n");
  printf("X = [ ");
  for (k=0;k<len;k++) printf("%2.5f ", X[k]);
  printf("]\n");


  Yo = interp1(X,Y,Xi,len);
  
  printf("input was %2.5f, output is %2.5f\n",Xi,Yo);

}

#endif
