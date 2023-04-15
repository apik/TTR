# add directory containing "firefly.pc" to PKG_CONFIG_PATH 

FF_CFLAGS=$(shell pkg-config --cflags firefly)
FF_LIBS=$(shell pkg-config --libs firefly)

TTR_CFLAGS= -Ofast 

TTR:TTR.cpp
	g++ $(TTR_CFLAGS) $(FF_CFLAGS) TTR.cpp $(FF_LIBS) -o TTR

