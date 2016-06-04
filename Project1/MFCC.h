#pragma once
#include "FFT.h"
class MFCC
{
public:
	MFCC(int window_lenght = 512, float low_freq = 0, float high_freq = 8000, int sample_rate = 16000);
	~MFCC(void);
	void Compute(int num_samples, int num_cep, float * input, float * output);
	void Preemphasis(int num_samples, float a_param, float * input);
	void AddRegression(float ** input, float ** output, int n, int num_cep, int win_length);
private:
	void DCT(int num_samples, bool inverse, float * input, float * output);
	float ** filter_bank;
	int * freq_res;
	inline float MFCC::hz_to_mel(float hz);
	inline float MFCC::mel_to_hz(float mel);
	
	FFT * fft;

	const int window_lenght;

	int Filters_Number;
	float preemphasis_param;
	int CepLifter;
	float * CepLift;
	int CepLiftSize;

	
public:
	void WeightCepstrum(float * input , int n);
	void GenCepLift(int CepLift, int n);
//	void SetMinDuration(void);
};
