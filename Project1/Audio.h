#pragma once
#include "Structs.h"
#include <vector>
class Audio
{
public:
	Audio(void);
	~Audio(void);
	int frame_size;
	int frame_overlap;
	int low_freq;
	int high_freq;
	int cep_number;
	int regression_window;
	int fft_frame_size;
	ParamAudio * ParamAudioFile(const char * audio_src, const char * label_src);
	
	std::vector<Label*> ExtractLabels(const char * label_src);
};

