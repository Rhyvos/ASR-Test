#pragma once
#include <string>
struct Wav												
	{
		char ChunkID[4];				// "RIFF"				
		int ChunkSize;					// 4 + (8 + SubChunk1Size) + (8 + SubChunk2Size)
		char Format[4];					// "WAVE"
		char Subchunk1ID[4];			// "fmt "
		int Subchunk1Size;				// bytes remaining in subchunk, 16 if uncompressed
		short AudioFormat;				// 1 = uncompressed
		short NumChannels;				// mono or stereo
		int SampleRate;					
		int ByteRate;					// == SampleRate * NumChannels * BitsPerSample/8
		short BlockAlign;				// == NumChannels * BitsPerSample/8
		short BitsPerSample;
		char Subchunk2ID[4];			// "data"
		int Subchunk2Size;				// == NumSamples * NumChannels * BitsPerSample/8
	};


struct Label
{
public:
	int start;
	int end;
	std::string name;
	Label(std::string name, long int start, long int end): start(start), end(end), name(name) {}
	Label(Label && other) : start(other.start), end(other.end), name(other.name) {}
	Label& operator= (Label && other)  {start=other.start; end=other.end; name=other.name; }
private:
	Label();
};

struct ObservationSegment
{
	float ** coef;
	float ** delta;
	float ** acc;
	int frame_lenght;
	int frames;
	Label * l;
};

struct ParamAudio
{
	ObservationSegment * os;
	int segments;
	int frame_size;
	int frame_overlap;
	std::string audio_src;
	std::string laber_scr;
	Wav audio_header;
	float ** coef_first;			//ptr on first element just for correct deleting
	float ** delta_first;
	float ** acc_first;
	int param_frames;
	~ParamAudio()
	{
	for(int i=0;i<segments;i++){
		if(os[i].l)
			delete os[i].l;
	}
	for (int i = 0; i < param_frames; i++)
	{
		if(coef_first[i])
			delete[] coef_first[i];
		if(delta_first[i])
			delete[] delta_first[i];
		if(acc_first[i])
			delete[] acc_first[i];
	}
	if(coef_first)
		delete[] coef_first;
	if(delta_first)
		delete[] delta_first;	
	if(acc_first)
		delete[] acc_first;
	if(os)
		delete[] os;
	}
};

struct State
{
	float * mean;
	float * var;

	float * old_mean;
	float * old_var;
	float obs_count;

	int state_nr;
	float g_const;
};



