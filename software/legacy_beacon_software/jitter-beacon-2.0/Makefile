
INCS = -I/usr/local/include -I/usr/local/include/gnuradio  #-I/usr/include/gnuradio/swig -I/usr/include/python2.6  -I/usr/include/gnuradio -DOMNITHREAD_POSIX=1 -Wno-strict-aliasing -Wno-parentheses -fPIC -DPIC 
LIBS =  -L/usr/lib -L/usr/local/lib -lboost_program_options -lm -luhd -lgnuradio-uhd -lgnuradio-core
CXXFLAGS=-O3 $(INCS) 
CFLAGS=-O3 -lm 

all:
	g++ $(CXXFLAGS) -o beacon receiver2.cc sgo_beacon_sink.cc flights.cc sgo_beacon_opts.cc $(LIBS)

#LIBS =  -lstdc++ -lboost_thread-gcc43-mt-1_35 -lrt -lboost_date_time-gcc43-1_35 -lgnuradio-core -lgruel -lfftw3f -lgsl -lusrp -lgslcblas -lm -lgromnithread -lboost_program_options -lgnuradio-usrp
