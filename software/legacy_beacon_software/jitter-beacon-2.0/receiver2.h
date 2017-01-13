#ifndef RECEIVER_H
#define RECEIVER_H
#endif

#include "sgo_beacon_sink.h"
#include <gr_top_block.h>
#include <gr_uhd_usrp_source.h>
#include <gr_null_sink.h>
#include <gr_uhd_api.h>
#include <gr_interleave.h>

#include <gr_file_sink.h>
#include <gr_realtime.h>


class usrp_beacon_receiver;
typedef boost::shared_ptr<usrp_beacon_receiver> usrp_beacon_receiver_sptr;


usrp_beacon_receiver_sptr make_usrp_beacon_receiver_cfile(int which, 
							  int decim, 
							  double freq1, 
							  double freq2, 
							  float gain1,
							  float gain2,
							  const std::string &outputdir,
							  double offset,
							  double offset2,
							  int useDriverFreq);

/*!
  Top block, which configures the USRP to acquire two 
  channels at specified frequencies, and then forwards
  them to the sgo_beacon_sink block
  \see sgo_beacon_sink
 */
class usrp_beacon_receiver : public gr_top_block
{
private:
  usrp_beacon_receiver(int which, 
		       long refclock, 
		       std::string addr,
		       std::string subdev,
		       int decim, 
		       double freq1, 
		       double freq2, 
		       float gain1,
		       float gain2, 
		       const std::string &outputdir,
		       double offset,
		       double offset2,
		       int useDriverFreq);
    
  friend usrp_beacon_receiver_sptr make_usrp_beacon_receiver(int which, 
							     long refclock,
							     std::string addr,
							     std::string subdev,
							     int decim, 
							     double freq1, 
							     double freq2, 
							     float gain1,
							     float gain2, 
							     const std::string &outputdir,
							     double offset,
							     double offset2,
							     int useDriverFreq);
  
  int d_which;
  int d_decim;
  double d_freq1;
  double d_freq2;
  float d_gain1;
  float d_gain2;
  sgo_beacon_options* o;
  double d_offset;
  double d_offset2;
  int d_useDriverFreq;
  
 public:
  gr_block_sptr d_beacon_sink;
};
