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

struct FrameParam						
{
	float * mfcc;
	float * delta;
	float * acc;
	int n;
};

struct Label
{
public:
	long long int start;
	long long int end;
	std::string name;
	Label(std::string name, long int start, long int end): start(start), end(end), name(name) {}
	Label(Label && other) : start(other.start), end(other.end), name(other.name) {}
private:
	Label();
};

struct ObservationSegment
{
	FrameParam * fr_array;
	int frames;
	Label l;
};

struct ParamAudio
{
	ObservationSegment * os;
	int segments;
	int frame_size;
	int frame_overlap;
	char * audio_src;
	char * laber_scr;
	Wav audio_header;
};


