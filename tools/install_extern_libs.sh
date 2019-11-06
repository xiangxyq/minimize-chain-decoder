#!/bin/bash

GCC=gcc

# Install alsa, portaudio and openblas libs
echo "*********************** Delete last time compile files ***********************"
rm -rf alsa/ portaudio/ openblas/

################################################## alsa libs download && compile ##################################################
echo "*********************** Download alsa-lib-1.1.7.tar.bz2 ***********************"
if [ ! -f "./alsa-lib-1.1.7.tar.bz2" ];then
	wget -T 10 -t 3 http://www.mirrorservice.org/sites/ftp.alsa-project.org/pub/lib/alsa-lib-1.1.7.tar.bz2
fi

tar -jxvf alsa-lib-1.1.7.tar.bz2
mv alsa-lib-1.1.7 alsa

echo "*********************** Compile alsa ***********************"
cd alsa/
CC=${GCC} CFLAGS="-O3" ./configure "--enable-shared"
make
mkdir -p install/lib
cp src/.libs/* install/lib -R
cp include install/include -R
cd ../

################################################## portaudio libs download && compile ##################################################
./install_portaudio.sh


################################################## OpenBLAS libs download && compile ##################################################
echo "*********************** Download OpenBLAS v0.3.6.tar.gz ***********************"
if [ ! -f "./v0.3.6.tar.gz" ];then
	wget -T 10 -t 3 https://github.com/xianyi/OpenBLAS/archive/v0.3.6.tar.gz
	#git clone https://github.com/xianyi/OpenBLAS.git
fi

tar zxvf v0.3.6.tar.gz
mv OpenBLAS-0.3.6 openblas

echo "*********************** Compile openblas ***********************"
cd openblas/
make CC=${GCC} ONLY_CBLAS=1 USE_THREAD=0 COMMON_OPT=" -O3"
make install PREFIX=`pwd`/install

echo "*********************** All libs Compiled Done!!! ***********************"
