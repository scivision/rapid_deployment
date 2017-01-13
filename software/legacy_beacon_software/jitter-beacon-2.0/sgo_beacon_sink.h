#ifndef SGO_BEACON_SINK_H
#define SGO_BEACON_SINK_H

#define datatype std::complex<float>

#include <gr_top_block.h>
#include <gr_uhd_usrp_source.h>
#include <gr_uhd_api.h>


#include <sys/time.h>
#include <sys/stat.h>
#include "flights.h"
#include <gr_null_sink.h> 
#include <gr_io_signature.h>
#include <gr_head.h>
#include <stdexcept>
#include <unistd.h>
#include <iostream>
#include <complex>
#include <stdio.h>
#include <boost/program_options.hpp>
#include <math.h>

#include <gr_top_block.h>
#include <gr_file_sink.h>

#include "sgo_beacon_opts.h"
#include "flights.h"

class sgo_beacon_sink;
typedef boost::shared_ptr<sgo_beacon_sink> sgo_beacon_sink_sptr;
typedef boost::shared_ptr<uhd_usrp_source> uhd_usrp_source_sptr;

sgo_beacon_sink_sptr make_sgo_beacon_sink(sgo_beacon_options *o, uhd_usrp_source_sptr usrp);

/*!
  Our own custom sink that eats dual channel single precision floating point IF data.
  We avoid interleaving to ensure that both channels have equal number of samples that are 
  coherent.
*/
class sgo_beacon_sink :public gr_sync_block
{
private:
  sgo_beacon_sink(sgo_beacon_options *o, uhd_usrp_source_sptr usrp);
  friend sgo_beacon_sink_sptr make_sgo_beacon_sink(sgo_beacon_options *o, uhd_usrp_source_sptr usrp);
  
  /*!
    All receiver related configuration
    \see sgo_beacon_opts.h
   */
  sgo_beacon_options *o;	
  
  int oldn;
    
  /*! 
     This contains information about the satellite passes
     we need two of these. one to hold the current day, and 
     one to hold the next.
   */
  flights *f[2];
  int flight_idx;

  /*
    receiver log file
   */
  FILE *rlogfile;
  
  /*!
     tmp vector for frequency shifting and decimation
   */
  float *tmpvec;
  int tmplen;

  /*!
     Time of previous satellite update 
     Used for updating active satellites every five seconds. 
  */
  double prev_update;
  
  uhd_usrp_source_sptr usrp;
  
  /*!
     This eats up samples from the usrp
     
     \param noutput_items Number of items going in to work()
     \param input_items   Pointer to data vector
     \param output_items  (unused)
     \return              The number of items consumed
  */
  int work(int noutput_items,
	   gr_vector_const_void_star &input_items,
	   gr_vector_void_star &output_items);
  
  /*!
     Shift in frequency, decimate and dump to file.
     \param chan       Requested channel {0,1}
     \param flight_idx Which satellite pass is in question
     \param vec        IF data, IQ and channels interleaved (I1Q1I2Q2...)
  */
  void shift_and_decimate(int chan, flights *fo, int flight_idx, float *vec, int nitems);
  
  /*!
    Write all startup information to receiver log.
   */
  void write_rlogpreamble();
  
  /*!
    Handle USRP overrun. Log to receiver.log, quit if necessary (it usually isn't).
   */
  void handle_usrp_overrun(int n);
  
  /*!
    Read in flights.log files, update flights if day is about to change.
   */
  void update_flights();
  
};

#endif
