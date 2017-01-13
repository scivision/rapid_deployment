/*!
 * (c) 2010,2011 Juha Vierinen
 *  
 * USRP2
 *
 * Dual channel beacon satellite receiver.
 * Every now and then read in a schedule file, 
 * and record the timeslots and bandwidths specified by
 * the schedule file passes.txt.
 */

#include "receiver2.h"

namespace po = boost::program_options;

// Shared pointer constructor
usrp_beacon_receiver_sptr make_usrp_beacon_receiver(int which, 
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
						    int useDriverFreq)
{
  return gnuradio::get_initial_sptr(new usrp_beacon_receiver(which, 
							     refclock,
							     addr,
							     subdev,
							     decim, 
							     freq1, 
							     freq2, 
							     gain1, 
							     gain2,
							     outputdir,
							     offset,
							     offset2,
							     useDriverFreq));
}

// Hierarchical block constructor, with no inputs or outputs
usrp_beacon_receiver::usrp_beacon_receiver(int which, long refclock,
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
					   int useDriverFreq) : 
  gr_top_block("usrp_beacon_receiver"),
  d_which(which), d_decim(decim), d_freq1(freq1), d_freq2(freq2), 
  d_gain1(gain1), d_gain2(gain2), d_offset(offset), d_offset2(offset2), d_useDriverFreq(useDriverFreq)
{
  uhd_usrp_source_sptr usrp;

  //  usrp = uhd_make_usrp_source(uhd::device_addr_t("addr0=192.168.10.2,recv_buff_size=100000000"),uhd::io_type_t::COMPLEX_FLOAT32, 2);
  usrp = uhd_make_usrp_source(addr,uhd::io_type_t::COMPLEX_FLOAT32, 2);
  usrp->set_clock_config(uhd::clock_config_t::external(), 0);
  usrp->set_subdev_spec(subdev,0);
  // use external reference. gpsdo+tvrx2

  //usrp->set_time_unknown_pps(0.0);

  //usrp->set_subdev_spec("B:0 A:0",0);
  usrp->set_samp_rate((double)refclock/decim);
  usrp->set_bandwidth((double)refclock/decim,0);
  usrp->set_bandwidth((double)refclock/decim,1);
  
  printf("Decim %d refclock %ld\n",decim, refclock);

  uhd::tune_result_t tf1 = usrp->set_center_freq(freq1,0);
  uhd::tune_result_t tf2 = usrp->set_center_freq(freq2,1);
  
  printf("Tune results 1 rf %1.10f ddc %1.10f\n",tf1.actual_rf_freq,tf1.actual_dsp_freq);
  printf("Tune results 2 rf %1.10f ddc %1.10f\n",tf2.actual_rf_freq,tf2.actual_dsp_freq);
  
  printf("Subdev %d %s\n",0,usrp->get_device()->get_rx_subdev_name(0).c_str());
  printf("Subdev %d %s\n",1,usrp->get_device()->get_rx_subdev_name(1).c_str());
  
  /* what is the actual center frequency */
  double actual_freq1;
  double actual_freq2;
  if(useDriverFreq){
    actual_freq1 = usrp->get_center_freq(0);
    actual_freq2 = usrp->get_center_freq(1);
  }
  else
  {
    printf("Using custom freq\n");
    actual_freq1 = freq1 + offset;
    actual_freq2 = freq2 + offset2;
  }
  /*
    with TVRX2 you cannot trust the tuning result, as it is plainly wrong. 
    we have to use a calibrated frequency offset.
    400 MHz = 400 MHz + 33.47291594699999706108 Hz
    150 MHz = 150 MHz - 1.70262008222810568014 Hz
   */


  usrp->set_gain(gain1,0);
  usrp->set_gain(gain2,1);
  
  printf("Actual sample rate %1.10f\n",usrp->get_samp_rate());
  printf("%d %s\n",0,usrp->get_dboard_sensor("lo_locked", 0).to_pp_string().c_str());
  printf("%d %s\n",1,usrp->get_dboard_sensor("lo_locked", 1).to_pp_string().c_str());
  printf("%s\n",usrp->get_mboard_sensor("ref_locked", 0).to_pp_string().c_str());
  //  printf("%s\n",usrp->get_mboard_sensor("gps_locked", 0).to_pp_string().c_str());

  printf("%s\n",usrp->get_dboard_sensor("rssi", 0).to_pp_string().c_str());
  printf("%s\n",usrp->get_dboard_sensor("rssi", 1).to_pp_string().c_str());

  gr_rt_status_t rts = gr_enable_realtime_scheduling();

  printf("Realtime scheduling %d.\n",rts);
  
  o=new sgo_beacon_options(actual_freq1, 
                           actual_freq2, 
                           d_decim, 
			   3600,
			   1.0,
			   1.0,
			   refclock,
			   outputdir);
  o->gain[0] = (double)gain1;
  o->gain[1] = (double)gain2;
  
  d_beacon_sink = make_sgo_beacon_sink(o, usrp);
  
  gr_interleave_sptr interleave;
  interleave = gr_make_interleave(2*sizeof(float));
  connect(usrp, 0, interleave, 0);
  connect(usrp, 1, interleave, 1);
  connect(interleave, 0, d_beacon_sink, 0);
}

int main(int argc, char *argv[])
{
  int which = 0;                       // specify which USRP board
  //usrp_subdev_spec spec1(0,0);         // specify the d'board side
  //usrp_subdev_spec spec2(0,1);         // specify the d'board side
  int decim = 64;                      // set the decimation rate
  double freq1 = 150e6;                // set the frequency
  double freq2 = 400e6;                // set the frequency
  double offset = 0.0;                 // sample clock offset
  double offset2 = 0.0;                 // sample clock offset
  int useDriverFreq = 0;
  float gain1 = 1;                    // set the gain; -1 will set the mid-point gain
  float gain2 = 1;                    // set the gain; -1 will set the mid-point gain
  char *dirname;                       // output directory name
  long refclock = 100000000;
  std::string outputdir = ".";
  std::string addr = "recv_frame_size=4096,num_recv_frames=4096";
  std::string subdev = "A:RX1 A:RX2";

  po::options_description cmdconfig("Program options: usrp_text_rx [options] filename");
  cmdconfig.add_options()
    ("help,h", "produce help message")
    ("which,W", po::value<int>(&which), "select which USRP board")
    ("subdev-spec,r", po::value<std::string>(), "Subdev spec (e.g., A:RX1 A:RX2, or B:0 A:0)")
    ("addr,a", po::value<std::string>(), "Device address string (e.g., recv_frame_size=4096,num_recv_frames=4096, or addr0=192.168.10.2,recv_buff_size=100000000")
    ("decim,d", po::value<int>(&decim), "set fgpa decimation rate to DECIM")
    ("driverFreq,u", po::value<int>(&useDriverFreq), "use driver frequency.")
    ("freq1,f", po::value<double>(), "set frequency to FREQ1")
    ("freq2,F", po::value<double>(), "set frequency to FREQ2")
    ("gain1,g", po::value<float>(), "set gain 1 in dB (default is midpoint)")
    ("gain2,G", po::value<float>(), "set gain 2 in dB (default is midpoint)")
    ("refclock,e", po::value<long>(), "reference clock frequency (default 64 MHz)")
    ("dirname,o", po::value<std::string>(), "Ouput directory name (default=.)")
    ("offset,s", po::value<double>(), "Frequency offset at 150 MHz (default=0.0)")
    ("offset2,t", po::value<double>(), "Frequency offset at 400 MHz (default=0.0)")
    ;

  po::options_description fileconfig("Input file options");
  fileconfig.add_options()
    ("filename", po::value<std::string>(), "input file")
    ;

  po::positional_options_description inputfile;
  inputfile.add("filename", -1);

  po::options_description config;
  config.add(cmdconfig).add(fileconfig);
  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
	    options(config).positional(inputfile).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << cmdconfig << "\n";
    return 1;
  }

  if(vm.count("dirname")) {
    outputdir = vm["dirname"].as<std::string>();
  }
  
  if(vm.count("freq1")) {
    freq1 = vm["freq1"].as<double>();
  }
  if(vm.count("freq2")) {
    freq2 = vm["freq2"].as<double>();
  }

  if(vm.count("gain1")) {
    gain1 = vm["gain1"].as<float>();
  }
  if(vm.count("gain2")) {
    gain2 = vm["gain2"].as<float>();
  }

  if(vm.count("refclock")) {
    refclock = vm["refclock"].as<long>();
  }

  if(vm.count("offset")) {
    offset = vm["offset"].as<double>();
  }
  if(vm.count("offset2")) {
    offset2 = vm["offset2"].as<double>();
  }

  if(vm.count("subdev-spec")) {
    subdev = vm["subdev-spec"].as<std::string>();
  }
  if(vm.count("addr")) {
    addr = vm["addr"].as<std::string>();
  }

  std::cout << "which:   " << which << std::endl;
  std::cout << "addr:   " << addr << std::endl;
  std::cout << "subdev:   " << subdev << std::endl;
  std::cout << "decim:   " << decim << std::endl;
  std::cout << "freq1:    " << freq1 << std::endl;
  std::cout << "freq2:    " << freq2 << std::endl;
  std::cout << "gain1:    " << gain1 << std::endl;
  std::cout << "gain2:    " << gain2 << std::endl;
  std::cout << "outputdir:    " << outputdir << std::endl;

  usrp_beacon_receiver_sptr top_block = make_usrp_beacon_receiver(which, refclock, addr, subdev, 
								  decim, freq1, freq2, 
								  gain1, gain2, outputdir, offset,offset2, useDriverFreq);
  top_block->run();
  
  return 0;
}
