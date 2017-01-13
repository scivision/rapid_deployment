#include "flights.h"

#define DEBUG_FLIGHTS
// #define DEBUG_FLIGHTS_2

/*
  res = a*res
 */
void complex_mul(complex_float *a, complex_float *res){
  float tmp;
  tmp = res->re;
  res->re = a->re*tmp - a->im*res->im;
  res->im = a->im*tmp + a->re*res->im;
}

bool is_utc_day_equal(utc_day_t *d1, utc_day_t *d2)
{
  if(d1->y == d2->y && d1->m == d2->m && d1->d == d2->d)
    return(true);
  else
    return(false);
}

void get_utc_day(utc_day_t *d){
  time_t t;
  struct tm tmptime;
  t = time(NULL);
  localtime_r(&t,&tmptime);
  d->y = tmptime.tm_year+1900;
  d->m = tmptime.tm_mon+1;
  d->d = tmptime.tm_mday;
}

void get_utc_day_f(utc_day_t *d, int secs){
  time_t t;
  struct tm tmptime;
  t = time(NULL);
  t = t + (time_t)secs;
  localtime_r(&t,&tmptime);
  d->y = tmptime.tm_year+1900;
  d->m = tmptime.tm_mon+1;
  d->d = tmptime.tm_mday;
}

double get_unix_time_now(){
  struct timeval tv;
  double tnow;  

  gettimeofday(&tv, NULL);
  tnow = ((double)tv.tv_sec)+(double)tv.tv_usec*1.0/1e6;
  return(tnow);
}

/*
  calculate cordic coefficients for different satellites.
 */
void flights::reset_cordic(satellite *s, double t, int chan){
  double dphase, dt;

  dt = (((double)o->decimation)/((double)o->sr));
  dphase = -2.0*SGO_PI*((double)(1e6*s->freq[chan]) - (double)o->usrp_freq[s->usrp_channel[chan]]);
#ifdef DEBUG_FLIGHTS_2
  printf("freqdiff %s %1.5f %1.5f\n",s->name,1e6*s->freq[0] - o->usrp_freq[s->usrp_channel[0]],1e6*s->freq[1] - o->usrp_freq[s->usrp_channel[1]]);
#endif

  s->dp[chan].re = (float)cos(dphase*dt);
  s->dp[chan].im = (float)sin(dphase*dt);

  s->osc[chan].re = (float)cos(dphase*t);
  s->osc[chan].im = (float)sin(dphase*t);
}



/* 
   Read in a data file. A flights object is always associated with one UTC day. 

   Use opt->datadir + "/" + YYYY.MM.DD format
   e.g., /data/beacon/data/2011.06.28/
*/
flights::flights(sgo_beacon_options *opt, utc_day_t *dday)
{
  FILE *df;
  char line[4096];
  double startt, endt, freq1, freq2, bw1, bw2, peak_el;
  char satname[1024];
  char fname[4096];
  char tmp[8192];
  char *c;
  int n;
  int fi=0;
  satellite *s;
  int len;

  ok = true;
  day.y = dday->y; day.m = dday->m; day.d = dday->d; 
  sprintf(dirname,"%s/%04d.%02d.%02d",opt->datadir.c_str(),day.y, day.m, day.d);

  sprintf(fname,"%s/passes.txt",dirname);
#ifdef DEBUG_FLIGHTS
  printf("dirname %s\n",dirname);
  printf("opening %s\n",fname);
#endif
  
  /* open log files */
  sprintf(tmp,"%s/flights.log",dirname);
  logfile = (FILE *)fopen(tmp, "a");
  if(logfile == NULL){
    ok = false;
    return;
  }
  
  nactive=0;
  o=opt;

  df = fopen(fname,"r");
  if(df == NULL){
    fclose(logfile);
    ok = false;
    return;
  }
  
  /* 
     first line contains header 
     one the first slurp, we only count the number of lines
  */
  c = fgets(line,4096,df);
  n=8;
  while(n==8){
    n = fscanf(df,"%lg %lg %lg %lg %lg %lg %lg %s",
	       &startt, &endt, &freq1, &freq2, &bw1, &bw2, &peak_el, satname);
    if(n==8)
      fi++;
  }
  fclose(df);

  nsatellites = fi;

  /* 
     and then we store them to the objects
  */
  satellites = new satellite[nsatellites];

  df = fopen(fname,"r");
  c = fgets(line,4096,df);

  for(fi=0 ; fi<nsatellites ; fi++){
    s = &satellites[fi];
    
    n = fscanf(df,"%lg %lg %lg %lg %lg %lg %lg %s",
	       &s->start_time, 
	       &s->end_time, 
	       &s->freq[0], 
	       &s->freq[1], 
	       &s->bw[0], 
	       &s->bw[1], 
	       &s->peak_elevation, 
	       satname);

#ifdef DEBUG_FLIGHTS_2
    printf("%lg %lg %lg %lg %lg %lg %lg %s\n",
	   s->start_time, 
	   s->end_time, 
	   s->freq[0], 
	   s->freq[1], 
	   s->bw[0], 
	   s->bw[1], 
	   s->peak_elevation, 
	   satname);
#endif    
    s->dec[0] = (int)floor( o->sr/o->decimation/(1000.0*s->bw[0]) );
    s->dec[1] = (int)floor( o->sr/o->decimation/(1000.0*s->bw[1]) );

    s->bw[0] = o->sr/((double)o->decimation)/((double)s->dec[0]);
    s->bw[1] = o->sr/((double)o->decimation)/((double)s->dec[1]);

    /* reset accumulators */
    s->acc[0] = 0; s->acc[1] = 0;
    s->re_acc[0] = 0.0; s->re_acc[1] = 0.0;
    s->im_acc[0] = 0.0; s->im_acc[1] = 0.0;

    s->t[0] = 0.0;
    s->t[1] = 0.0;

    s->noverruns = 0;

    /*
      Determine which USRP channel corresponds to which logical channel
     */
    if( fabs(s->freq[0]*1e6 - o->usrp_freq[0]) < fabs(s->freq[0]*1e6 - o->usrp_freq[1]) ){
      s->usrp_channel[0] = 0;
      s->usrp_channel[1] = 1;
      if( fabs(s->freq[0]*1e6 - o->usrp_freq[0]) > o->sr/2.0/o->decimation )
	printf("Warning, beacon %s with freq %1.10lf does not fit band %1.10lf! Decrease decimation.\n",satname, s->freq[0]*1e6, o->usrp_freq[0]);
      if( fabs(s->freq[1]*1e6 - o->usrp_freq[1]) > o->sr/2.0/o->decimation )
	printf("Warning, beacon %s with freq %1.10lf does not fit band %1.10lf! Decrease decimation.\n",satname, s->freq[1]*1e6, o->usrp_freq[1]);
    }else{
      s->usrp_channel[0] = 1;
      s->usrp_channel[1] = 0;

      if( fabs(s->freq[0]*1e6 - o->usrp_freq[1]) > o->sr/2.0/o->decimation )
	printf("Warning, beacon %s with freq %1.10lf does not fit band %1.10lf! Decrease decimation.\n",satname, s->freq[0]*1e6, o->usrp_freq[1]);
      if( fabs(s->freq[1]*1e6 - o->usrp_freq[0]) > o->sr/2.0/o->decimation )
	printf("Warning, beacon %s with freq %1.10lf does not fit band %1.10lf! Decrease decimation.\n",satname, s->freq[1]*1e6, o->usrp_freq[0]);

    }

    /* 
       if valid line in conf file
       - initialize file names
       - setup cordic coefficients
     */
    if(n==8){
      len = strlen(satname);
      s->name = (char *)malloc(sizeof(char)*(len+10));
      strcpy(s->name, satname);
      
      sprintf(tmp,"%s/flight-%06d.001.bin",dirname,fi);
      len = strlen(tmp);

      s->fname[0] = (char *)malloc(sizeof(char)*(len+10));
      strcpy(s->fname[0], tmp);

      sprintf(tmp,"%s/flight-%06d.002.bin",dirname,fi);
      len = strlen(tmp);
      s->fname[1] = (char *)malloc(sizeof(char)*(len+10));
      strcpy(s->fname[1], tmp);

      reset_cordic(s,0.0,0);
      reset_cordic(s,0.0,1);
    }else{
      break;
    }
  }
  
  /* 
     Active flag array for all passes. 
     set to 0 when inactive and 1 when active 
   */
  active = (int *)malloc(sizeof(int)*nsatellites);
  for(fi=0 ; fi < nsatellites ; fi++)
  {
    active[fi] = 0;
  }

  /* 
     Indices of active satellites. Only the first nactive are used.
     But we allocate everything so we don't have to realloc on the fly.
   */
  active_idx = (int *)malloc(sizeof(int)*nsatellites);

#ifdef DEBUG_FLIGHTS
  printf("number of flights %d\n", nsatellites);
#endif
  fclose(df);
}



/* 
   Go through all passes and figure out which are currently active
 */
int flights::find_active_satellites(){

  int i, bestidx;
  satellite *s;
  double tnow;
  double best, dt;
  
  if(!ok)
  {
    nactive = 0;
    return(0);
  }

  tnow = get_unix_time_now();
  
  best = 24.0*60.0*60.0;
  bestidx = 0;
  
  nactive = 0;
  for(i=0 ; i<nsatellites ; i++)
  {
    s = &satellites[i];
    
#ifdef DEBUG_FLIGHTS
    /* which satellite is next */
    dt = s->start_time - tnow;
    if(dt > 0.0 && dt < best){
      best = dt;
      bestidx = i;
    }
#endif
    
    if(tnow > s->start_time && tnow < s->end_time)
    {
        active_idx[nactive++] = i;
        
        if(active[i] == 0)
        {
          /* 
             if this is a new satellite:
             - open file handles
             - reset cordic
          */ 
          active[i] = 1;
#ifdef DEBUG_FLIGHTS
          printf("Starting to acquire satellite %s\n",s->name);
#endif
          s->tstart = tnow;
          s->file[0] = fopen(s->fname[0],"wb");
          s->file[1] = fopen(s->fname[1],"wb");
          
          reset_cordic(s, 0.0, 0);
          reset_cordic(s, 0.0, 1);
        }
    }
    else
    {
      if(active[i]==1){
        /* 
           This ends the dump:
           - we close the file 
           - write entry in logfile 
        */
        fprintf(logfile, "%06d %1.3lf %1.3f %s %s %1.3lf %1.3lf %d %d %1.3lf %1.3lf %d %s\n",i, s->tstart, tnow, s->fname[0], s->fname[1], s->freq[0], s->freq[1], s->dec[0], s->dec[1], s->bw[0], s->bw[1], s->noverruns, s->name);
        fflush(logfile);

        fclose(s->file[0]);
        fclose(s->file[1]);
      }
      active[i] = 0;
    }
  }
#ifdef DEBUG_FLIGHTS
  s = &satellites[bestidx];
  if(best < (60*60*24-1))
  {
    printf("Acquiring %d, next satellite %s in %f s\n", nactive, s->name, best);
  }
#endif
  return(nactive);
}

/*
  Destructor, remove all allocated memory.
 */
flights::~flights()
{
  int i;
  satellite *s;
  if(!ok){
    return;
  }
  for(i=0 ; i<nsatellites ; i++){
    s = &satellites[i];
    free(s->fname[0]);
    free(s->fname[1]);
    free(s->name);
  }
  delete(satellites);
  free(active);
  free(active_idx);
  fclose(logfile);
}

/*
  test creation and deletion.
 */
#ifdef TEST_FLIGHTS
int main(int argc, char **argv)
{
  flights *f;
  FILE *log;
  int n;
  sgo_beacon_options *o;

  o = new sgo_beacon_options(150.0, 400.0, 64, 3600, 1.0, 1.0, 64e6, ".");

  log = fopen("flights.log","w");
  f = new flights((char *)"passes.txt",o);
  n = f->find_active_satellites(log);
  printf("found %d active satellites\n",n);
  delete(f);
  fclose(log);
}

#endif
