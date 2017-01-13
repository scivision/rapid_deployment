#include <stdio.h>
#include <stdint.h>

/* gcc -fPIC -shared -o liblp.so lp.c -lfftw3f */

#include <fftw3.h>

// static variables
static fftwf_complex *in, *out;
static fftwf_plan p;
static int zdec=0;
  
int lp(float *z, int len, float *pfb_out, int dec, float *win)
{
  int i,j,nwin;
  float *inf, *outf;

  if ( zdec != dec )
  {
    if(zdec != 0)
    {
      printf("free\n");
      fftwf_free(in);
      fftwf_free(out);
      fftwf_destroy_plan(p);
    }
    printf("fftw alloc\n");
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * dec);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * dec);
    p = fftwf_plan_dft_1d(dec, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  }
  inf=(float *)in;
  outf=(float *)out;
  nwin = len/dec;
  for(i=0 ; i < nwin ; i++)
  {
    for(j=0 ; j < dec ; j++)
    {
      in[j][0] = z[(i*dec + j)*2]*win[j];
      in[j][1] = z[(i*dec + j)*2+1]*win[j];
    }
    fftwf_execute(p); /* repeat as needed */
    pfb_out[2*i]=outf[0];
    pfb_out[2*i+1]=outf[1];
  }
}

/* 
   polyphase filterbank 
   z = complex floating point array to be processed (float32 re, float32 im)
   len = length of z. must be n*dec
   pfb_out = output of polyphase filterbank (complex float32). 
             it is two dimensional with shape (len/dec, n_channels).
   dec = decimation length, or number of points in fft of pfb
   *win = complex float array defining the filter coefficients.
   pl = how many times longer is window than fft length dec. real float32 window coefficients
   *channels = array of channels to output (fft indices)
   n_channels = number of channels to output
 */
int pfb(float *z, int len, float *pfb_out, int dec, float *win, int pl, int64_t *channels, int n_channels)
{
  int i,j,k,nwin;
  float *inf;
  float *outf;
  for(i=0 ; i < n_channels ; i++)
  {
    printf("channel %d %d\n",i,(int)channels[i]);
  }
  if ( zdec != dec )
  {
    if(zdec != 0)
    {
      fftwf_free(in);
      fftwf_free(out);
      fftwf_destroy_plan(p);
    }
    in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * dec);
    out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * dec);
    p = fftwf_plan_dft_1d(dec, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  }
  inf=(float *)in;
  outf=(float *)out;
  nwin = len/dec;
  for(i=0 ; i < nwin ; i++)
  { 
    for(j=0 ; j < dec ; j++)
    {
      in[j][0]=0.0;
      in[j][1]=0.0;

      for(k=0 ; k < pl ; k++)
      {
        float zre = z[(i*dec + k*dec + j)*2];
        float zim = z[(i*dec + k*dec + j)*2 + 1];        
        float wre = win[(k*dec + j)*2];
        float wim = win[(k*dec + j)*2+1];
        in[j][0] += zre*wre - zim*wim;
        in[j][1] += zre*wim + zim*wre;
      }
    }
    fftwf_execute(p); /* repeat as needed */
    for(j=0 ; j<n_channels ; j++)
    {
      pfb_out[(j*nwin + i)*2]=out[channels[j]][0];
      pfb_out[(j*nwin + i)*2+1]=out[channels[j]][1];
    }
  }  
  fftwf_free(in);
  fftwf_free(out);
  fftwf_destroy_plan(p);
  zdec=0;
}

