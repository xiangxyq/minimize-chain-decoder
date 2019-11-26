#include "../feat/compute-feat.h"
#include <stdio.h>
#include <iostream>
#include "../audio/audio-capture.h"
#include "../asr-itf/asr-apis.h"

int main(int argc, char *argv[])
{
  InitAsrModels();

  while (true)
  {
    std::vector<BaseFloat> waveform;
    std::vector<std::string> result;
    BaseFloat likelihood = 0.0f;

    Read(waveform);
    AsrDecodingBegin(waveform);

    AsrDecodingResult(result, &likelihood);
    for(int32 i = 0; i < result.size(); i++)
      std::cout << result[i] << " ";
    std::cout << std::endl;
  }

  AsrEnd();

  return 0;
}
