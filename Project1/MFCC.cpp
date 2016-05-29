#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif

#define _USE_MATH_DEFINES
#include "MFCC.h"
#include "FFT.h"
#include <cmath>
#include <stdio.h>

#define TEST
#ifdef TEST
float Mel(int k,float fres);
float WarpFreq (float fcl, float fcu, float freq, float minFreq, float maxFreq , float alpha);
void InitFBank(int frameSize, long sampPeriod, int numChans,
			   float lopass, float hipass, bool usePower, bool takeLogs,
			   bool doubleFFT,
			   float alpha, float warpLowCut, float warpUpCut,short *loChan,float*loWt );  
#endif // TEST


MFCC::MFCC(int window_lenght, float low_freq, float high_freq, int sample_rate): window_lenght(window_lenght)
{
	preemphasis_param = 0.97;
	CepLifter = 22;
	Filters_Number = 24;

	float tmp_low = hz_to_mel(low_freq);
	float tmp_high = hz_to_mel(high_freq);
	float step = (tmp_high-tmp_low)/(float)(Filters_Number+1);
	float * freq_res = new float[Filters_Number+2];
	filter_bank = new float *[Filters_Number];

	

	for(int i=0; i < Filters_Number+2 ; i++)
	{
		freq_res[i] = mel_to_hz(tmp_low + (i * step));
		freq_res[i] = ((this->window_lenght+1)*freq_res[i])/sample_rate;
	}

	for (int i = 0; i < Filters_Number; i++)
	{
		filter_bank[i] = new float[this->window_lenght/2];
		for (int j = 0; j < this->window_lenght/2; j++)
		{
			if(j<freq_res[i])
				filter_bank[i][j]=0;

			if(j>=freq_res[i] && j<=freq_res[i+1])
				filter_bank[i][j]=(j-freq_res[i])/(freq_res[i+1]-freq_res[i]);

			if(j>=freq_res[i+1] && j<=freq_res[i+2])
				filter_bank[i][j]=(freq_res[i+2]-j)/(freq_res[i+2]-freq_res[i+1]);;

			if(j>freq_res[i+2])
				filter_bank[i][j]=0;
		}
	}

	fft = new FFT();

	

	delete[] freq_res;
}


MFCC::~MFCC(void)
{
	if (filter_bank)
	{
		for(int i=0;i<Filters_Number;i++)
				delete[] filter_bank[i];
			delete[] filter_bank;
	}
	delete fft;
	delete[] CepLift;
}





void MFCC::Compute(int num_samples, int num_cep, float * input, float * output)
{
	
	if(window_lenght < num_samples || num_samples <= 0)
	{
		fprintf(stderr, "[MFCC::Compute] Wrong audio window size: %d\n",window_lenght);
		return;
	}

	if(Filters_Number<num_cep || num_cep <= 0)
	{
		fprintf(stderr, "[MFCC::Compute] Wrong coefficients number: %d\n",num_cep);
		return;
	}





	float *tmp_data = new float[window_lenght];
	float *real_data = new float[window_lenght];
	float *imag_data = new float[window_lenght];
	float ek;
	float *filter_energy = new float[Filters_Number];

	for(int i=0 ; i<num_samples; i++)
	{
		tmp_data[i]=input[i];
	}

	for(int i=num_samples ; i<window_lenght; i++)
	{
		tmp_data[i]=0.0;
	}

	for (int i = 0; i < Filters_Number; i++)
	{
		filter_energy[i]=0;
	}


	Preemphasis(num_samples,preemphasis_param,tmp_data);					//only audio data
	fft->Window_Func(FFT::Hamming,num_samples,tmp_data);					//only audio data

	fft->Transform(window_lenght,false,tmp_data,NULL,real_data,imag_data);

#ifndef TEST
	for (int i = 0; i < Filters_Number; i++)
	{
		for (int j = 0; j < window_lenght/2; j++)
		{
			ek=sqrt(real_data[j]*real_data[j]+imag_data[j]*imag_data[j]);
			filter_energy[i]+=ek*filter_bank[i][j]; 
		}

		filter_energy[i]=log(filter_energy[i]);
	}  
#else
	 float t1,t2;
	 const float melfloor = 1.0;
	 int bin;
#pragma error(disable : 4700)
	 short *loChan = new short[window_lenght/2+1];
	 float *loWt= new float[window_lenght/2+1];
	 InitFBank(num_samples,625,Filters_Number,-1,-1,false,true,false,1,0,0,loChan,loWt);
	 for (int k = 2; k <= window_lenght/2; k++) {             /* fill bins */
      t1 = real_data[k-1]; t2 = imag_data[k-1];
      
         ek = sqrt(t1*t1 + t2*t2);
	  
      bin = loChan[k];
      t1 = loWt[k]*ek;
	  
      if (bin>0) filter_energy[bin-1] += t1;
	  if (bin<Filters_Number) filter_energy[bin] += ek - t1;
   }
   //for (bin=1; bin<=info.numChans; bin++) { printf("fbank[%d] = %f\n",bin,fbank[bin]);}getchar();
   /* Take logs */
   
	 for (bin=0; bin<Filters_Number; bin++) { 
         t1 = filter_energy[bin];
		 if (t1<melfloor){ t1 = melfloor; }
		 filter_energy[bin] = log(t1);
      }
	 delete[] loChan;
	 delete[] loWt;
#endif // !test


	float *tmp = new float[Filters_Number];
	DCT(Filters_Number,false,filter_energy,tmp);

	for (int i = 0; i < num_cep; i++)
	{
		output[i] = tmp[i];
	}
	WeightCepstrum(output,num_cep);
	
	delete[] tmp_data;
	delete[] real_data;
	delete[] imag_data;
	delete[] tmp;
	delete[] filter_energy;
}


void MFCC::DCT(int num_samples, bool inverse, float * input, float * output)
{
   int j,k,numChan;
   float mfnorm,pi_factor,x;
   mfnorm = sqrt(2.0/(float)num_samples);
   pi_factor = M_PI/(float)num_samples;

 

   for (j=0; j<num_samples; j++)  {
	   output[j] = 0.0; x = (float)(j + 1)* pi_factor;
      for (k=0; k<num_samples; k++)
	    output[j] += input[k] * cos(x*(k+0.5)); 
      output[j] *= mfnorm;
   }  
}


inline float MFCC::hz_to_mel(float hz)
{
	return 1127 * log(1+hz/700);
}

inline float MFCC::mel_to_hz(float mel)
{
	return 700 * (exp(mel/1127)-1);
}

void MFCC::Preemphasis(int num_samples, float a_param, float * input)
{
   for (int i=num_samples-1; i>=1 ;i--){
      input[i] -= input[i-1] * a_param;
   }

   input[0] *= (1.0-a_param);

}


void MFCC::AddRegression(float ** input, float ** output, int n, int num_cep, int win_length)
{
	int sigma = 0;
	int back, forw;
	float sum = 0.0;

	for(int i=1; i<=win_length; i++)
		sigma += i*i;
	sigma *= 2;

	for(int i=0; i<n; i++)
	{
		for (int j = 0; j < num_cep; j++)
		{
			back = forw = i;
			sum = 0.0;
			for (int k = 1; k <= win_length; k++)
			{
				if(back-1 >= 0) back --;
				if(forw+1 <  n) forw ++;
				sum += k * (input[forw][j] - input[back][j]);  
			}
			output[i][j] = sum / sigma;
		}
	}
}


void MFCC::WeightCepstrum(float * input , int n)
{
   int i,j;
   
   if ( n > CepLiftSize)
      GenCepLift(CepLifter,n);
   
   for (i=0;i<n;i++)
		input[i] *= CepLift[i];
}


void MFCC::GenCepLift(int CepLiftCoef, int n)
{
   int i;
   float a, Lby2;
 
   if (CepLift==NULL ||CepLiftSize < n)
      CepLift = new float[n];
   a = M_PI/CepLiftCoef;
   Lby2 = CepLiftCoef/2.0;
   for (i=0;i<n;i++)
		CepLift[i] = 1.0 + Lby2*sin((i+1) * a);
   CepLiftSize = n;
}


#ifdef TEST
float Mel(int k,float fres)
{
	return 1127 * log(1 + (k-1)*fres);
}

float WarpFreq (float fcl, float fcu, float freq, float minFreq, float maxFreq , float alpha)
{
	if (alpha == 1.0)
		return freq;
	else {
		float scale = 1.0 / alpha;
		float cu = fcu * 2 / (1 + scale);
		float cl = fcl * 2 / (1 + scale);

		float au = (maxFreq - cu * scale) / (maxFreq - cu);
		float al = (cl * scale - minFreq) / (cl - minFreq);

		if (freq > cu)
			return  au * (freq - cu) + scale * cu ;
		else if (freq < cl)
			return al * (freq - minFreq) + minFreq ;
		else
			return scale * freq ;
	}
}

void InitFBank(int frameSize, long sampPeriod, int numChans,
			   float lopass, float hipass, bool usePower, bool takeLogs,
			   bool doubleFFT,
			   float alpha, float warpLowCut, float warpUpCut,short *loChan,float*loWt )
{
	float mlo,mhi,ms,melk;
	int k,chan,maxChan,Nby2;
	int fftN;
	float fres;
	int klo,khi;
	/* Save sizes to cross-check subsequent usage */
	frameSize = frameSize; numChans = numChans;
	sampPeriod = sampPeriod; 
	usePower = usePower; takeLogs = takeLogs;

	/* Calculate required FFT size */
	fftN = 2;   
	while (frameSize>fftN) fftN *= 2;
	if (doubleFFT) 
		fftN *= 2;
	Nby2 = fftN / 2;
	fres = 1.0E7/(sampPeriod * fftN * 700.0);
	maxChan = numChans+1;

	/* set lo and hi pass cut offs if any */
	klo = 2; khi = Nby2;       /* apply lo/hi pass filtering */
	mlo = 0; mhi = Mel(Nby2+1,fres);

	if (lopass>=0.0) {
		mlo = 1127*log(1+lopass/700.0);
		klo = (int) ((lopass * sampPeriod * 1.0e-7 * fftN) + 2.5);
		if (klo<2) klo = 2;
	}
	if (hipass>=0.0) {
		mhi = 1127*log(1+hipass/700.0);
		khi = (int) ((hipass * sampPeriod * 1.0e-7 * fftN) + 0.5);
		if (khi>Nby2) khi = Nby2;
	}


	/* Create vector of fbank centre frequencies */
	float * cf = new float[maxChan+1];
	ms = mhi - mlo;

	for (chan=1; chan <= maxChan; chan++) {
		if (alpha == 1.0) {
			cf[chan] = ((float)chan/(float)maxChan)*ms + mlo;

		}
		else {
			/* scale assuming scaling starts at lopass */
			float minFreq = 700.0 * (exp (mlo / 1127.0) - 1.0 );
			float maxFreq = 700.0 * (exp (mhi / 1127.0) - 1.0 );
			float cf1 = ((float)chan / (float) maxChan) * ms + mlo;

			cf1 = 700 * (exp (cf1 / 1127.0) - 1.0);

			cf[chan] = 1127.0 * log (1.0 + WarpFreq (warpLowCut, warpUpCut, cf1, minFreq, maxFreq, alpha) / 700.0);
		}

	}

	/* Create loChan map, loChan[fftindex] -> lower channel index */

	for (k=1,chan=1; k<=Nby2; k++){
		melk = Mel(k,fres);
		if (k<klo || k>khi) loChan[k]=-1;
		else {
			while (cf[chan] < melk  && chan<=maxChan) ++chan;
			loChan[k] = chan-1;
		}

	}

	/* Create vector of lower channel weights */   

	for (k=1; k<=Nby2; k++) {
		chan = loChan[k];
		if (k<klo || k>khi) loWt[k]=0.0;
		else {
			if (chan>0) 
				loWt[k] = ((cf[chan+1] - Mel(k,fres)) / 
				(cf[chan+1] - cf[chan]));
			else
				loWt[k] = (cf[1]-Mel(k,fres))/(cf[1] - mlo);
		}
	}
	/* Create workspace for fft */

	delete[] cf;
}  
#endif // TEST
