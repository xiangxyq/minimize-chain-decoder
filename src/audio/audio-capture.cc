/*
* This file is used for audio capture
*
* author: xiangxyq
* E-mail: xiangxyq@foxmail.com
*/

/* Use the newer ALSA API */
#define ALSA_PCM_NEW_HW_PARAMS_API
#include <alsa/asoundlib.h>
#include <iostream>
#include "audio-capture.h"

#define CHNNAL_NUM    (1)
#define SAMPLE_RATE   (16000)
#define DATA_BASE_LEN (1000000)
#define AUDIO_LEN     (2)    //if set 2, the capture audio data length is 2 second

void Read(std::vector<BaseFloat> &data)
{
  long loops;
  int32 rc;
  int32 size;
  snd_pcm_t *handle;
  snd_pcm_hw_params_t *params;
  unsigned int val;
  int32 dir;
  snd_pcm_uframes_t frames;
  char *buffer;

  /* Open PCM device for recording (capture). */
  rc = snd_pcm_open(&handle, "default", SND_PCM_STREAM_CAPTURE, 0);
  if (rc < 0) 
  {
    std::cout << "unable to open pcm device: " << snd_strerror(rc) << std::endl;
    exit(1);
  }

  /* Allocate a hardware parameters object. */
  snd_pcm_hw_params_alloca(&params);
  /* Fill it in with default values. */
  snd_pcm_hw_params_any(handle, params);
  /* Set the desired hardware parameters. */
  /* Interleaved mode */
  snd_pcm_hw_params_set_access(handle, params, SND_PCM_ACCESS_RW_INTERLEAVED);
  /* Signed 16-bit little-endian format */
  snd_pcm_hw_params_set_format(handle, params, SND_PCM_FORMAT_S16_LE);
  /* Mono channels (Mono) */
  snd_pcm_hw_params_set_channels(handle, params, CHNNAL_NUM);
  /* 16000 bits/second sampling rate */
  val = SAMPLE_RATE;
  snd_pcm_hw_params_set_rate_near(handle, params, &val, &dir);
  /* Set period size to 32 frames. */
  frames = 32;
  snd_pcm_hw_params_set_period_size_near(handle, params, &frames, &dir);

  /* Write the parameters to the driver */
  rc = snd_pcm_hw_params(handle, params);
  if (rc < 0) 
  {
      std::cout << "unable to set hw parameters: " << snd_strerror(rc) << std::endl;
      exit(1);
  }

  /* Use a buffer large enough to hold one period */
  snd_pcm_hw_params_get_period_size(params, &frames, &dir);
  size = frames * 2 * CHNNAL_NUM; /* 2 bytes/sample, 1 channels */
  buffer = (char *)malloc(size);

  /* We want to loop for 2 seconds */
  snd_pcm_hw_params_get_period_time(params, &val, &dir);
  loops = AUDIO_LEN * DATA_BASE_LEN / val;

  while (loops > 0)
  {
    loops--;
    rc = snd_pcm_readi(handle, buffer, frames);
    if (rc == -EPIPE) 
    {
      /* EPIPE means overrun */
      std::cout << "overrun occurred" << std::endl;
      snd_pcm_prepare(handle);
    } 
    else if (rc < 0) 
      std::cout << "error from read: " << snd_strerror(rc) << std::endl;
    else if (rc != (int)frames) 
      std::cout << "short read, read "<< rc << " frames" << std::endl;
    
    int16 *temp = (int16 *) buffer;
    for (int i = 0; i < size / 2; ++i)
      data.push_back(static_cast<BaseFloat> (temp[i]));
  }

  snd_pcm_drain(handle);
  snd_pcm_close(handle);
  free(buffer);
}
