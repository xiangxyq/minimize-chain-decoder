
#
# This file is used for build
#
# author: xiangxyq
# E-mail: xiangxyq@foxmail.com
#

CC = gcc
CXX = g++
STRIP = strip


OPENBLASINC = ../../tools/openblas/install/include
OPENBLASLIBS = ../../tools/openblas/install/lib/libopenblas.a

ALSAINC = ../../tools/alsa/install/include
ALSALIBS = ../../tools/alsa/install/lib

CC_OPTIONS = -std=c99 -I.. -I$(OPENBLASINC) -I$(ALSAINC) -c -pthread -O3 -Wall -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated-declarations -Winit-self
CXX_OPTIONS = -std=c++11 -I.. -I$(OPENBLASINC) -I$(ALSAINC) -c -pthread -O3 -Wall -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated-declarations -Winit-self


LIBS = -ldl  -lm $(OPENBLASLIBS) -L$(ALSALIBS) -lasound -lpthread

OBJ = audio-capture.o asr-math.o mel-banks.o srfft.o compute-feat.o main.o

decoder.bin : $(OBJ)
	$(CXX) -o decoder.bin $(OBJ) $(LIBS)
	$(STRIP) decoder.bin
	rm *.o

#audio
audio-capture.o : ../audio/audio-capture.cc
	$(CXX) $(CXX_OPTIONS) ../audio/audio-capture.cc

#utils
asr-math.o : ../utils/asr-math.cc
	$(CXX) $(CXX_OPTIONS) ../utils/asr-math.cc

# feat
mel-banks.o : ../feat/mel-banks.cc
	$(CXX) $(CXX_OPTIONS) ../feat/mel-banks.cc

srfft.o : ../feat/srfft.cc
	$(CXX) $(CXX_OPTIONS) ../feat/srfft.cc

compute-feat.o : ../feat/compute-feat.cc
	$(CXX) $(CXX_OPTIONS) ../feat/compute-feat.cc

#main
main.o : ../main/chain-decoder.cc
	$(CXX) $(CXX_OPTIONS) ../main/chain-decoder.cc -o main.o

clean :
	rm -rf *.o
