#include "sgo_beacon_sink.h"

#define DEBUG_BEACON_SINK

#define datatype std::complex<float>

sgo_beacon_sink_sptr make_sgo_beacon_sink(sgo_beacon_options *o, uhd_usrp_source_sptr u)
{
  return sgo_beacon_sink_sptr(new sgo_beacon_sink(o,u));
}

/*!
  what day is it
  do we have a directory for this day?
  - if not
    - fetch beacon.tle
    - predict future passes
  do we have a directory for now + 30 minutes?
  - if not
    - fetch beacon.tle
    - predict future passes
 */
void sgo_beacon_sink::update_flights()
{
  utc_day_t d;
  utc_day_t d2;
  char tmp[512];

  get_utc_day(&d);
  
  if(f[flight_idx] == NULL)
  {
    fprintf(rlogfile,"%1.3f Updating flights.\n", get_unix_time_now());
    f[flight_idx] = new flights(o,&d);
    fprintf(rlogfile,"%1.3f Number of satellite passes detected %d %04d.%02d.%02d ok? %d\n", get_unix_time_now(), f[flight_idx]->nsatellites,d.y,d.m,d.d,f[flight_idx]->ok);
    fflush(rlogfile);
  }
  
  /* if day changes 30 minutes from now, create new flights */
  get_utc_day_f(&d2,1800);
  if(!is_utc_day_equal(&d, &d2))
  {
    //fprintf(rlogfile,"%1.3f Updating flights 2.\n", get_unix_time_now());
    if(f[(flight_idx+1) % 2] != NULL)
    {
      /* is the other flight object on the new or old day? */
      if(!is_utc_day_equal(&d2, &f[flight_idx]->day))
      {
	fprintf(rlogfile,"%1.3f Deleting old flights.\n", get_unix_time_now());
	flight_idx = (flight_idx + 1) % 2;
	delete(f[flight_idx]);
	f[flight_idx] = new flights(o,&d2);
	fprintf(rlogfile,"%1.3f Number of satellite passes detected %d %04d.%02d.%02d ok? %d\n", get_unix_time_now(), f[flight_idx]->nsatellites,d2.y,d2.m,d2.d,f[flight_idx]->ok);
	fflush(rlogfile);
      }
    }
    else
    {
      fprintf(rlogfile,"%1.3f Creating new flights.\n", get_unix_time_now());
      flight_idx = (flight_idx + 1) % 2;
      f[flight_idx] = new flights(o,&d2);
      fprintf(rlogfile,"%1.3f Number of satellite passes detected %d %04d.%02d.%02d ok? %d\n", get_unix_time_now(), f[flight_idx]->nsatellites,d2.y,d2.m,d2.d,f[flight_idx]->ok);
      fflush(rlogfile);
    }
  }
}

/*!
  This is a sink that eats a two-channel data stream, which has the following format:
  I1Q1I2Q2... 

  The data is processed in the following way:
  - These two streams are first separated
  - We shift the center frequency to baseband
  - Integrate and decimate to the closest bandwidth that is more or equal to the requested bandwidth 
 */
sgo_beacon_sink::sgo_beacon_sink(sgo_beacon_options *opt, uhd_usrp_source_sptr u) :
  gr_sync_block ("sgo_beacon_sink",
                 gr_make_io_signature(1, 1, sizeof(std::complex<float>)),
                 gr_make_io_signature(0, 0, 0))
{
  int ret;
  char dirstr[512];
  char tmp[8192];
  
  usrp = u;
  o = opt;
  oldn = 0;
  
  flight_idx = 0;
  f[0] = NULL; f[1] = NULL;

  sprintf(tmp,"%s/receiver.log",opt->datadir.c_str());
  rlogfile = (FILE *)fopen(tmp, "a");
  write_rlogpreamble();

  /* 
     initialize flights()
   */
  update_flights();
  prev_update = 0;
  
  /* allocate temporary output vector */
  tmpvec = (float *)malloc(sizeof(float)*16384);
  tmplen=16384;
}

void sgo_beacon_sink::write_rlogpreamble()
{
  double tnow = get_unix_time_now();
  fprintf(rlogfile,"%1.3f Beacon receiver started\n",tnow);
  fprintf(rlogfile,"%1.3f Subdev %d %s\n",tnow, 0,usrp->get_device()->get_rx_subdev_name(0).c_str());
  fprintf(rlogfile,"%1.3f Subdev %d %s\n",tnow, 1,usrp->get_device()->get_rx_subdev_name(1).c_str());
  fprintf(rlogfile,"%1.3f Gain %1.2f %1.2f\n", tnow, o->gain[0], o->gain[1]);
  fprintf(rlogfile,"%1.3f Freq %1.5f %1.5f\n", tnow, o->usrp_freq[0], o->usrp_freq[1]);
  fprintf(rlogfile,"%1.3f Spectrum sign %f %f\n", tnow, o->sign[0], o->sign[1]);
  fprintf(rlogfile,"%1.3f Locked LO %d %d\n",tnow, usrp->get_dboard_sensor("lo_locked", 0).to_bool(),usrp->get_dboard_sensor("lo_locked", 1).to_bool());
  fflush(rlogfile);
}

void sgo_beacon_sink::handle_usrp_overrun(int n)
{
  int i;
  double tnow = get_unix_time_now();
  satellite *s;

  fprintf(rlogfile,"%1.3f %d overruns/underruns detected. Discontiniuty in data.\n", tnow, n);
  fflush(rlogfile);
  oldn = n;
}

/*
  The work function eats the usrp data vector by vector.
 */
int sgo_beacon_sink::work(int noutput_items,
			  gr_vector_const_void_star &  	input_items,
			  gr_vector_void_star &  	output_items)
{
  float *host_items;
  double tnow;
  int j, i, no;

  tnow = get_unix_time_now();

  /* every five seconds, see if new satellites are visible */
  if(tnow - prev_update > 5)
  {
    /* check if we need to get ready for a change in day */
    update_flights();
    if(!f[flight_idx]->ok)
      printf("ERROR: Could not find satellites! Check configuration.\n");
    for(i=0 ; i<2 ; i++)
    {
      if(f[i] != NULL)
	f[i]->find_active_satellites();
    }
    prev_update = tnow;
    fflush(stdout);
  }

  /*
    For all active satellites, shift in frequency, decimate and dump to file
    We have two different days when there can be active flights.
    When the day changes, there can be flights on both days that are active.
   */
  for(j=0 ; j<2 ; j++){
    if(f[j] != NULL)
    {
      for(i=0 ; i<f[j]->nactive ; i++){
	shift_and_decimate(0, f[j], f[j]->active_idx[i], (float*)input_items[0], noutput_items);
	shift_and_decimate(1, f[j], f[j]->active_idx[i], (float*)input_items[0], noutput_items);
      }
    }
  }
  /* eat all items */
  return(noutput_items);
}

/*
  Shift frequency to baseband and decimate items
  
  @param n       Which channel (0 = freq1, 1 = freq2)
  @param fi      Which flight
  @param items   Interleaved items
  @param nitems  Number of items
 */
void sgo_beacon_sink::shift_and_decimate(int chan, flights *fo, int fi, float *items, int nitems)
{
  int k, oi, u_chan;
  float re_in, im_in, re_osc, im_osc;
  satellite *s;
  
  s = &fo->satellites[fi];
  u_chan = s->usrp_channel[chan];
    
  /* if our temp vector is too short, increase length */
  if( nitems/2/s->dec[chan] > (2*tmplen) ){
    free(tmpvec);
    tmplen = (nitems/2/s->dec[chan])*2 + 1024;
    tmpvec = (float *)malloc(sizeof(float)*tmplen);
  }
  
  /*
    Eat all elements. 
    there are two channels interleaved, so there are actually 
    nitems/2 items in the array.
   */
  oi = 0;
  k = 0;
  s->t[chan] += ((double)o->decimation)*((double)nitems)/2.0/((double)o->sr);
  while(k < nitems/2)
  {
    /* 
       de-interleave
     */
    re_in = items[4*k + u_chan*2];
    im_in = items[4*k + u_chan*2 + 1];
    
    /* 
       Advance cordic by one tick 
    */
    complex_mul(&s->dp[chan],&s->osc[chan]);
    
    re_osc = s->osc[chan].re;
    im_osc = s->osc[chan].im;

    /*
      output down-converted signal:
    */
    s->re_acc[chan] += re_in*re_osc - im_in*im_osc;
    s->im_acc[chan] += im_in*re_osc + re_in*im_osc;

    s->acc[chan]++;
    
    if(s->acc[chan] == s->dec[chan])
    {
      tmpvec[2*oi] = s->re_acc[chan];
      tmpvec[2*oi+1] = s->im_acc[chan];
      s->acc[chan] = 0;
      s->re_acc[chan] = 0.0;
      s->im_acc[chan] = 0.0;
      oi++;
    }
    k++;
  }
  
  /* 
     Increase stability by using precise cos and sin updates
     Does this really work?
  */
  fo->reset_cordic(s, s->t[chan], chan);
  int outc = fwrite(tmpvec,sizeof(float),oi*2,s->file[chan]); 
}
