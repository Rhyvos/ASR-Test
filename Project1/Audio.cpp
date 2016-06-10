#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif


#include "Audio.h"
#include <fstream>
#include "MFCC.h"
#include <vector>
#include <stdlib.h> 
#pragma warning (disable : 4996)
#define BitsPerByte 8
#define Second 10000000						// 1s = 1000000000ns but we have time variable in 100ns units
Audio::Audio(void)
{
	frame_size = 400;
	frame_overlap = 240;
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
	
	if(!audio)
	{
		fprintf(stderr,"Can't open file: %s\n",audio_src);
		return NULL;
	}
	
		
	
	ParamAudio * pm;
	pm = new ParamAudio;
	
	int bytes;

	bytes = fread(&pm->audio_header,1,sizeof(pm->audio_header),audio);
	if(bytes != sizeof(pm->audio_header))
	{
		fprintf(stderr,"Wrong audio file: %s\n",audio_src);
		return NULL;
	}
	high_freq = high_freq ? high_freq : (pm->audio_header.SampleRate/2);

	MFCC *mfcc = new MFCC(fft_frame_size,low_freq,high_freq,pm->audio_header.SampleRate);

	short buffer;
	float *buffer_array = new float[frame_size];
	int index = 0;
	int frame_index = 0;
	int frames = 0;// floor((pm->audio_header.Subchunk2Size/(pm->audio_header.BitsPerSample/BitsPerByte)) / (frame_size-frame_overlap));

	while(frames*(frame_size-frame_overlap)+frame_overlap <( pm->audio_header.Subchunk2Size/(pm->audio_header.BitsPerSample/BitsPerByte)))
		frames++;
	frames--;

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

	int c=0;
	while(fread(&buffer,1,sizeof(buffer),audio))									// audio parametrizing loop
	{
		buffer_array[index++] = (float)buffer;
		c++;
		if(index == frame_size)
		{
			mfcc->Compute(index,cep_number,buffer_array,coef[frame_index++]);			// get frame parametrs
			for(int i=0; i<frame_overlap; i++)
				buffer_array[i] = buffer_array[frame_size - frame_overlap + i];		//copy overlaping data for each frame
			index = frame_overlap;
		}
	}

	/*if(index > frame_overlap)
	{
		mfcc->Compute(index,cep_number,buffer_array,coef[frame_index++]);	
	}*/
	
	mfcc->AddRegression(coef,delta,frame_index,cep_number,regression_window);
	mfcc->AddRegression(delta,acc,frame_index,cep_number,regression_window);
	
	std::vector<Label*> Labels = ExtractLabels(label_src);
	
	if(Labels.size() <= 0)
	{
		fprintf(stderr,"Bad label file: %s, proceed segmentation without label file\n",label_src);
		if(label_src != NULL)
			pm->laber_scr = label_src;
		else 
			pm->laber_scr = "";
		pm->audio_src = audio_src;

		pm->segments = 1;											//because only one segment with whole audio params
		pm->frame_size = frame_size;
		pm->frame_overlap = frame_overlap;
		pm->coef_first = coef;
		pm->delta_first = delta;
		pm->acc_first = acc;
		pm->param_frames = frames;
		pm->coef_num = cep_number;
		int sample_delta = Second/pm->audio_header.SampleRate;
		int i=0;
		pm->os = new ObservationSegment[pm->segments];

		pm->os[0].coef = coef;
		pm->os[0].delta = delta;
		pm->os[0].acc = acc;
		pm->os[0].frame_lenght = cep_number;
		pm->os[0].frames = frames;
		pm->os[0].l = new Label("",0,0);
	}
	else
	{
		pm->laber_scr = label_src;
		pm->audio_src = audio_src;
		pm->segments = Labels.size();
		pm->frame_size = frame_size;
		pm->frame_overlap = frame_overlap;
		pm->coef_first = coef;
		pm->delta_first = delta;
		pm->acc_first = acc;
		pm->param_frames = frames;
		pm->coef_num = cep_number;
		int sample_delta = Second/pm->audio_header.SampleRate;
		int i=0;
		pm->os = new ObservationSegment[pm->segments];

		for(std::vector<Label*>::iterator l = Labels.begin(); l != Labels.end(); l++)
		{
			int start_frame = (*l)->start/sample_delta;
			int end_frame = (*l)->end/sample_delta;
			start_frame = floor(start_frame/(frame_size-frame_overlap));
			end_frame = floor(end_frame/(frame_size-frame_overlap));
			if(end_frame - start_frame <= 0 || start_frame>=frames)
			{
				fprintf(stderr,"Label[%d]: wrong segment lenght: %d or start frame: %d out of bound: %d, skipping that label\n",i,end_frame - start_frame,start_frame,frames);
			}
			else
			{
				if(start_frame>=0 && start_frame<frames)
				{
					pm->os[i].coef = coef+start_frame;
					pm->os[i].delta = delta+start_frame;
					pm->os[i].acc = acc+start_frame;
				}
				else
				{
					fprintf(stderr,"Label[%d]: wrong start frame\n",i,start_frame);
					continue;
				}

				
				pm->os[i].frame_lenght = cep_number;
				if(end_frame > frames)
					pm->os[i].frames = frames - start_frame;
				else
					pm->os[i].frames = end_frame - start_frame;
				pm->os[i].l = *l;
				i++;
			}
			
		}
		if(i != Labels.size())
		{
			for(int j=i;j<Labels.size();j++)
				delete Labels.at(j);
			pm->segments = i;
		}
	}



	if(audio)
		fclose(audio);



	delete[] buffer_array;
	delete mfcc;
	return pm;
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

std::vector<Label*> Audio::ExtractLabels(const char * label_src)
{
	std::vector<Label*> output;
	if(label_src == NULL)
		return output;

	std::ifstream input(label_src);
	std::vector<std::string> split_buffer;
	
	int start,end;
	int index=0;
	char buffer[256];
	if(input)
		while(input.good())
		{
			input.getline(buffer,256);
			split_buffer = split(buffer," \n\t");
			if(split_buffer.size() != 3)								// 3 because startpoint endpoint and label name per line
			{
				fprintf(stderr,"Label file:%s Corrupted at line %d: %s\n",label_src,index+1,buffer);
			}
			else
			{
				start = atoi(split_buffer[0].c_str());
				end = atoi(split_buffer[1].c_str());
				if(start>=end)
				{
					fprintf(stderr,"Label file:%s Corrupted at line %d: %s (start time higher then end time)\n",label_src,index+1,buffer);
				}
				else
					output.push_back(new Label(split_buffer[2],start,end));
				index++;
			}
		}
	return output;
}



