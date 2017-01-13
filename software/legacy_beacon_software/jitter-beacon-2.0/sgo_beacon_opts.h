#ifndef _SGO_BEACON_OPTS_H_
#define _SGO_BEACON_OPTS_H_

#include <stdio.h>
#include <gr_top_block.h>
//#include <usrp_source_base.h>
//#include <usrp_source_c.h>
//#include <usrp_source_s.h>
#include <gr_file_sink.h>

/*!
  Store all beacon satellite receiver related documentation in this class.
 */
class sgo_beacon_options
{
public:
  /*! Center frequency of channel 1 & 2 */
  double usrp_freq[2];
  
  /*! Original USRP sample rate before decimation */
  double sr;
  
  /*! Decimation */
  double decimation;
  
  /*! How often do we update out schedule file (seconds) */
  int updateInterval;
  
  /*! Is spectrum flipped? Channel 1 & 2 */
  double sign[2];
  
  /*! DDC frequency */
  double dxc_freq[2];
  
  /*! Directory for data */
  std::string datadir;
  
  double gain[2];
  
  /* handles to the subdevices (e.g., for AGC) */
  //db_base_sptr subdev[2];
  
  /* do we stop on an overrun? */
  int stop_on_overrun;
  
  sgo_beacon_options(double f1, double f2, int d, int ui, double s1, double s2, double r, const std::string &outputdir);
  
};

#endif
