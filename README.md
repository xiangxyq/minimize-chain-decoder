# minimize-chain-decoder
## Introduction

Minimize kaldi nnet3 chain decoder

Chinese introduce: [README_CN]()

## 1. Install external libs: alsa, openblas

Moved to tools folder, then run below command:

```
./install_extern_libs.sh
```

if message "All libs Compiled Done!!!" shows that compile alsa, openblas successfully!

Notice: openblas just run a general compile, if you need optimize it, should specify the CPU target.

You can know more about in this [github](https://github.com/xianyi/OpenBLAS)


## 2. Audio Capture Introduction

audio capture code is in src/audio/ path, there is only one api void Read(std::vector<BaseFloat> &data)

Default capture audio data is 2 seconds, if need modify, you can change AUDIO_LEN macro;

## 3. MFCC Introduction

Please refer my blog: https://blog.csdn.net/cj1989111/article/details/102954071

