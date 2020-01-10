
#
# This file is used for build
#
# author: xiangxyq
# E-mail: xiangxyq@foxmail.com
#

CC = gcc
CXX = g++
STRIP = strip

USE_OPEN_BLAS=true

ifeq ($(USE_OPEN_BLAS) , true)
OPENBLASINC = ../../tools/openblas/install/include
OPENBLASLIBS = ../../tools/openblas/install/lib/libopenblas.a
endif

ALSAINC = ../../tools/alsa/install/include
ALSALIBS = ../../tools/alsa/install/lib

ifeq ($(USE_OPEN_BLAS), true)
CC_OPTIONS = -std=c99 -I.. -I$(OPENBLASINC) -I$(ALSAINC) -c -pthread -O3 -Wall -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated-declarations -Winit-self -D USE_OPENBLAS
CXX_OPTIONS = -std=c++11 -I.. -I$(OPENBLASINC) -I$(ALSAINC) -c -pthread -O3 -Wall -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated-declarations -Winit-self -D USE_OPENBLAS
LIBS = -ldl  -lm $(OPENBLASLIBS) -L$(ALSALIBS) -lasound -lpthread
else
CC_OPTIONS = -std=c99 -I.. -I$(ALSAINC) -c -pthread -O3 -Wall -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated-declarations -Winit-self
CXX_OPTIONS = -std=c++11 -I.. -I$(ALSAINC) -c -pthread -O3 -Wall -Wno-sign-compare -Wno-unused-local-typedefs -Wno-deprecated-declarations -Winit-self
LIBS = -ldl  -lm -L$(ALSALIBS) -lasound -lpthread
endif

OBJ = asr-math.o hash-list.o matrix.o tokens.o \
	audio-capture.o  mel-banks.o srfft.o compute-feat.o \
	fst.o word-symbols.o compents.o chain-compute.o decode.o \
	asr-apis.o main.o

decoder.bin : $(OBJ)
	$(CXX) -o decoder.bin $(OBJ) $(LIBS)
	$(STRIP) decoder.bin
	rm *.o

#utils
asr-math.o : ../utils/asr-math.cc
	$(CXX) $(CXX_OPTIONS) ../utils/asr-math.cc

hash-list.o : ../utils/hash-list.cc
	$(CXX) $(CXX_OPTIONS) ../utils/hash-list.cc

matrix.o : ../utils/matrix.cc
	$(CXX) $(CXX_OPTIONS) ../utils/matrix.cc

tokens.o : ../utils/tokens.cc
	$(CXX) $(CXX_OPTIONS) ../utils/tokens.cc

#audio
audio-capture.o : ../audio/audio-capture.cc
	$(CXX) $(CXX_OPTIONS) ../audio/audio-capture.cc

# feat
mel-banks.o : ../feat/mel-banks.cc
	$(CXX) $(CXX_OPTIONS) ../feat/mel-banks.cc

srfft.o : ../feat/srfft.cc
	$(CXX) $(CXX_OPTIONS) ../feat/srfft.cc

compute-feat.o : ../feat/compute-feat.cc
	$(CXX) $(CXX_OPTIONS) ../feat/compute-feat.cc

#models
fst.o : ../models/fst.cc
	$(CXX) $(CXX_OPTIONS) ../models/fst.cc

word-symbols.o : ../models/word-symbols.cc
	$(CXX) $(CXX_OPTIONS) ../models/word-symbols.cc

#nnet3
compents.o : ../nnet3/compents.cc
	$(CXX) $(CXX_OPTIONS) ../nnet3/compents.cc

chain-compute.o : ../nnet3/chain-compute.cc
	$(CXX) $(CXX_OPTIONS) ../nnet3/chain-compute.cc

#decode
decode.o : ../decode/decode.cc
	$(CXX) $(CXX_OPTIONS) ../decode/decode.cc

#asr-itf
asr-apis.o : ../asr-itf/asr-apis.cc
	$(CXX) $(CXX_OPTIONS) ../asr-itf/asr-apis.cc

#main
main.o : ../main/chain-decoder.cc
	$(CXX) $(CXX_OPTIONS) ../main/chain-decoder.cc -o main.o

clean :
	rm -rf *.o
