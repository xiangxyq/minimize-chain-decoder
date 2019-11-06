# minimize-chain-decoder
## Introduction

Minimize kaldi nnet3 chain decoder, and speech front add simple vad

Chinese introduce: [README_CN]()

## 1. Install external libs: alsa, portaudio, openblas

Moved to tools folder, then run below command:

```
./install_extern_libs.sh
```

if message "All libs Compiled Done!!!" shows that compile alsa, portaudio, openblas successfully!

Notice: openblas just run a general compile, if you need optimize it, should specify the CPU target.

You can know more about in this [github](https://github.com/xianyi/OpenBLAS)

