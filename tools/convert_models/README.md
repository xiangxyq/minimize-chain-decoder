# Introduction
1. This script need run in **Ubuntu 16.04**, other version not verified.

2. Need install libatlas-base-dev

```
if you meet 'libatlas.so.3: cannot open shared object file: No such file or directory'

Please install libatlas-base-dev

Unbuntu install command: sudo apt-get install libatlas-base-dev
```

3. openfst install

You can enter kaldi tools dir, and make, then openfst will download automatic and compiled.

You just need configure it to your PATH environment, do as this:

```
vim ~/.bashrc
export OPEN_FST=/xxx/kaldi/tools/openfst/bin/i686-m64 #xxx is your path
export PATH=$OPEN_FST:$PATH

Modifed finished, then run this:
source ~/.bashrc
```

4. Generate Headers

Sript and test models packed in convert_models.tar.gz, In case the model file is damaged when upload, so need unzip it first.

run tar zxvf convert_models.tar.gz

Place your models(features has no ivector , no pitch, and mfcc 40 dims or 13 dims) in the existing format, then run convert_models.sh scripts

Will get headers in ./models/models_headers/ dir.

5. copy headers in src/models/ then re-compile source.