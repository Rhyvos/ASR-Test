#pragma once
#include "Structs.h"
#include <vector>
#include "Config.h"
class Audio
{
public:
	Audio(Config * cf = nullptr);
	~Audio(void);
	Config * cf;
	int trace;
	int frame_size;
	int frame_overlap;
	int low_freq;
	int high_freq;
	int cep_number;
	int regression_window;
	int fft_frame_size;
	ParamAudio * ParamAudioFile(const char * audio_src, const char * label_src);
	
	std::vector<Label*> ExtractLabels(const char * label_src);
//	void LoadHmm(HMM * hmm, std::string hmm_src);
	
};

