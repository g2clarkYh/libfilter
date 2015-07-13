/* filter.h - header file for filter utilities */

/* HAL 6-9-09 */
#ifndef _FILT_H
#define _FILT_H

/* misc AMS compatibility defs */

#define ERR_LEV0 0

/* misc limits and constants */
#define MAX_BUTTERWORTH_SIZE 8
#define FAST_QUEUE 1 
#define DOUBLE_FLOATS 1

#if DOUBLE_FLOATS
typedef double Float; /* allows to play precission games */
#else
typedef float Float;  
#endif

/* basic complex number type */

typedef struct { 
  Float real;
  Float imag;
} Complex;

/* filter state history queue  */

typedef struct {
  Float *buf;    /* history storage */
  int rp;        /* read pointer */
  int wp;        /* write pointer */
  int N;         /* size of buffer */
  int nch;       /* number of channels */
  unsigned int mask;  /* buffer counter mask */
} Queue;

/* FIR/IIR filter definition TODO: make multi-channel
 * (for now, this will be a homogeneous multi-channel 
 *  filter.  TODO: add multiple filter coeff vectors 
 */
typedef struct {
  Float *ff;   /* feedforward filter coeffs */
  Float *fb;   /* feedback filter coeffs */
  Queue *q;    /* filter state history storage queue */
  int N;       /* feedforward filter order */
  int M;       /* feedback filter order */
  int nch;     /* number of channels to filter */
  void *ptr;  /* malloc pointer for alloc/dealloc mem */
} Filter;

/* function prototypes */
Complex *complex_mult(Complex *in1,Complex *in2, Complex *out);
Complex *complex_div(Complex *in1,Complex *in2, Complex *out);

Float *queue_init(Queue *q, Float *buf,int size);
int push_queue(Queue *q,Float *data,int N);
int pop_queue(Queue *q,Float *data,int N);

Filter *filter_init(  Filter *filt, /* actual filter structure */
		      Float *ff,   /* feedback filter coeffs */
		      Float *fb,   /* feedback filter coeffs */
		      Queue *q,    /* filter state history storage queue */
		      int N,       /* feedforward filter order */
		      int M,       /* feedback filter order */ 
		      int nch     /* number of channels in the data */);

Float *filter(Filter *f, Float *in, Float *out, int nsamps);
Float *butter(Float *coeffs,float bw, float fc, int size);
Float interp1(Float *X, Float *Y, Float Xi, int N);
#endif
