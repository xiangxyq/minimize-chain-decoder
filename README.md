# minimize-chain-decoder
## Introduction

Minimize kaldi nnet3 chain decoder, it can recognize domain words in embedded platform.

I used it in CPU: 1GHz, ram:16m, flash can ignore, because this decoder executable file is very small, contains models(based on your model), my demo is 484K.

About 2s sound wave, nosiy envirenment time cost is 300~400ms.

I also try porting it to STM32, STM32 memory resource is crisis, chain model decode need a lot of memory(Token,ForwardLinks), so you can refer online-wav-gmm-decode-faster.cc

where ProcessEmitting() caculate LogLikelihood() instead of chain neural networks compute result. I finished this and can work on STM32.

Others: This decoder not include vad , ns , agc ..., if these done, The recognize result will be better.

## 1. Install external libs: alsa, openblas

Moved to tools folder, then run below command:

```
./install_extern_libs.sh
```

if message "All libs Compiled Done!!!" shows that compile alsa, openblas successfully!

Notice: openblas just run a general compile, if you need optimize it, should specify the CPU target.

You can know more about in this [github](https://github.com/xianyi/OpenBLAS)

**sooner I will add a generic matrix operation for Optional**


## 2. Audio Capture Introduction

audio capture code is in src/audio/ path, there is only one api void Read(std::vector<BaseFloat> &data)

Default capture audio data is 2 seconds, if need modify, you can change AUDIO_LEN macro;

## 3. MFCC Introduction

Please refer my blog: https://blog.csdn.net/cj1989111/article/details/102954071


## 4. NNET3 Introduction

// TODO

## 5. HCLG search Introduction

// TODO

## 6. How to compile and run

Before step 1 need done.

cd src/build/

run make

then you can see decoder.bin executable file, in this demo, you can say wakeup word "智能管家", it can recognized and very few misidentifications.

## 7. How to use your models

**Notice**: models should follows as these ---> features has no ivector , no pitch, and mfcc 40 dims or 13 dims

The models I trained under [aishell2](https://github.com/kaldi-asr/kaldi/tree/master/egs/aishell2)

The trained scripts currently not open source.

Please refer to [convert_models](https://github.com/xiangxyq/minimize-chain-decoder/blob/master/tools/convert_models/README.md)


