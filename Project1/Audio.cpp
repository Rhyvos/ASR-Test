#include "Audio.h"
#include <fstream>
#include "MFCC.h"
#include <vector>
#include <stdlib.h> 
#pragma warning (disable : 4996)
#define BitsPerByte 8
Audio::Audio(void)
{
	frame_size = 400;
	frame_overlap = 160;
	low_freq = 0;
	high_freq = 0;
	cep_number = 12;
	regression_window = 2;
	fft_frame_size = 2;
	while(fft_frame_size<frame_size) fft_frame_size<<=1;
}


Audio::~Audio(void)
{
}


ParamAudio * Audio::ParamAudioFile(const char * audio_src, const char * label_src)
{
	FILE * audio;
	audio = fopen(audio_src,"rb");
    FILE * label;
	label = fopen(label_src,"rb");
	if(!audio)
	{
		if(label)
			fclose(label);
		fprintf(stderr,"Can't open file: %s",audio_src);
		return NULL;
	}
	
		
	
	ParamAudio * pm;
	pm = new ParamAudio;
	
	int bytes;

	bytes = fread(&pm->audio_header,1,sizeof(pm->audio_header),audio);
	if(bytes != sizeof(pm->audio_header))
	{
		fprintf(stderr,"Wrong audio file: %s",audio_src);
		return NULL;
	}
	high_freq = high_freq ? high_freq : (pm->audio_header.SampleRate/2);

	MFCC *mfcc = new MFCC(fft_frame_size,low_freq,high_freq,pm->audio_header.SampleRate);

	short buffer;
	float *buffer_array = new float[frame_size];
	int index = 0;
	int frame_index = 0;
	int frames = floor((pm->audio_header.Subchunk2Size/(pm->audio_header.BitsPerSample/BitsPerByte)) / (frame_size-frame_overlap));

	while(frames*(frame_size-frame_overlap)+frame_overlap <( pm->audio_header.Subchunk2Size/(pm->audio_header.BitsPerSample/BitsPerByte)))
		frames++;

	float ** coef, **delta, ** acc;

	coef = new float*[frames];
	delta = new float*[frames];
	acc = new float*[frames];


	for(int i=0 ; i<frames ;i++ )
	{
		coef[i] = new float[cep_number];
		delta[i] = new float[cep_number];
		acc[i] = new float[cep_number];
	}


	while(fread(&buffer,1,sizeof(buffer),audio))									// audio parametrizing loop
	{
		buffer_array[index++] = (float)buffer;

		if(index == frame_size)
		{
			mfcc->Compute(index,cep_number,buffer_array,coef[frame_index++]);			// get frame parametrs
			for(int i=0; i<frame_overlap; i++)
			{
				buffer_array[i] = buffer_array[frame_size - frame_overlap + i];		//copy overlaping data for each frame
			}
			index = frame_overlap;
		}
	}

	if(index > frame_overlap)
	{
		mfcc->Compute(index,cep_number,buffer_array,coef[frame_index++]);	
	}
	
	mfcc->AddRegression(coef,delta,frame_index,cep_number,regression_window);
	mfcc->AddRegression(delta,acc,frame_index,cep_number,regression_window);

	std::vector<Label> Labels = ExtractLabels(label_src);
	


	if(audio)
		fclose(audio);
	if(label)
		fclose(label);
	delete[] buffer_array;
	delete mfcc;
	return NULL;
}

std::vector<std::string> split(char * str,const char * delimiter) {
  std::vector<std::string> internal;
  char * buff;
  buff = strtok(str,delimiter);
  while(buff != NULL) {
    internal.push_back(buff);
	buff = strtok(NULL,delimiter);
  }
  
  return internal;
}

std::vector<Label> Audio::ExtractLabels(const char * label_src)
{
	std::ifstream input(label_src);
	std::vector<std::string> split_buffer;
	std::vector<Label> output;
	int index=0;
	char buffer[256];
	while(input.good())
	{
		input.getline(buffer,256);
		split_buffer = split(buffer," \n\t");
		if(split_buffer.size() != 3)								// 3 because startpoint endpoint and label name per line
		{
			fprintf(stderr,"Label file:%s Corrupted at line %d: %s",label_src,index+1,buffer);
			return output;
		}
		else
		{
			output.emplace_back(split_buffer[2],std::atol(split_buffer[0].c_str()),atol(split_buffer[1].c_str()));
			index++;
		}
	}
	return output;
}

