#include "sgo_beacon_opts.h"

/*
  Nothing fancy here yet. Just store options.
 */
sgo_beacon_options::sgo_beacon_options(double f1, 
				       double f2, 
				       int d, 
                                       int ui, 
				       double s1, 
				       double s2, 
				       double sample_rate,
				       const std::string &data_dir)
{
  usrp_freq[0] = f1;
  usrp_freq[1] = f2;
  sr = sample_rate;
  sign[0] = s1;
  sign[1] = s2;
  updateInterval = ui;
  decimation = d;
  stop_on_overrun=0;
  datadir = data_dir;
}
