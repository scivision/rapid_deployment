#ifndef FLIGHTS_H
#define FLIGHTS_H

#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "sgo_beacon_opts.h"
#include <math.h>

#define SGO_PI 3.14159265358979323846264338327950288

typedef struct complex_float_str {
  float re;
  float im;
} complex_float;

typedef struct utc_day_t_str {
  int y;
  int m;
  int d;
} utc_day_t;

void complex_mul(complex_float *a, complex_float *res);
bool is_utc_day_equal(utc_day_t *d1, utc_day_t *d2);
void get_utc_day(utc_day_t *d);
void get_utc_day_f(utc_day_t *d, int secs);
double get_unix_time_now();


/*!
  A class that contains information about one satellite
 */
class satellite
{
 public:
  /*!
    The starting and ending time of pass.
   */
  double start_time, end_time;

  /*!
    Center frequency of channels 1 and 2
   */
  double freq[2];
  
  /*!
    Which usrp channels are to be used
   */
  int usrp_channel[2];

  /*!
    The requested bandwidth that this satellite is stored with
   */
  double bw[2];
  
  /*!
    How much additional decimation is needed
   */
  int dec[2];
  
  /*!
    phase increment
   */
  complex_float dp[2];
  
  /*!
    cordic oscillator 
   */
  complex_float osc[2];

  /*!
     Peak elevation of the satellite on this pass
   */
  double peak_elevation;

  /*!
     Satellite name 
  */
  char *name;
  
  

  /*! 
    File name and handle for channels 1 & 2
  */
  char *fname[2];
  FILE *file[2];
  
  /*!
    real and imaginary accumulator
   */
  float re_acc[2];
  float im_acc[2];
  
  /*!
    accumulator index
   */
  int acc[2];
  
  /*!
    time in sample units.
   */
  double t[2];
  
  int noverruns;
  
  double tstart;
};


/*!
  A class that encapsulates the collection of all satellite overpasses.
 */
class flights
{
 public:
  satellite *satellites;
  utc_day_t day;

  int *active;            /* which satellites are active */
  int *active_idx;        /* which satellites are active */
  int nsatellites;
  int nactive;
  sgo_beacon_options *o;
  /*!
    data directory prefix.
    e.g., 
    /data/beacon/data/2011.06.28/
   */
  char dirname[8192];
  FILE *logfile;
  bool ok;
  
  
  /*
    Read in data file.
   */
  flights(sgo_beacon_options *o, utc_day_t *day);
  
  /* 
     Destructor. Delete all allocated memory.
   */
  ~flights();

  /*!
    Determine which satellites are active.

    \return number of active satellites
   */
  int find_active_satellites(); 

  
  /*!
    reset the cordics for satellite s
   */
  void reset_cordic(satellite *s, double t, int chan);

};

#endif
