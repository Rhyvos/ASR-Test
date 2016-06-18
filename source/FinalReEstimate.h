#pragma once
#include <vector>
#include "HMM.h"
class FinalReEstimate
{
public:
	FinalReEstimate(Config * cf);
	~FinalReEstimate(void);
	std::vector<HMM *> hmms;
	HMM ** hmm_seq;							//list of hmms ordered by apperence at label file
	Config *cf;
	int trace;
	float *** outprob;
	double ** alphat;
	double ** alphat1;
	double *** beta;
	float * occt;
	int *qHi;
	int *qLo;
	int frames;
	int seq_num;
	double pruneThresh;
	double pr;
	float minFrwdP; 
	
	float *** GetProbability(ParamAudio * pa);
	float OutP(ParamAudio * pa, int fr_number, int segment, int state);
	double Beta(void);
	double LAdd(double x, double y);
	double MaxModelProb(int q, int t, int minq);
	void AddHmm(HMM *  hmm);
	void LoadHmm(std::string hmm_src);
	void ForwardBackward(ParamAudio * pa);
	void ListHmms(ParamAudio * pa);
	void SetBeamTaper(int q, int t);
	void Alpha(ParamAudio * pa);
	void StepAlpha(int t, int * start, int * end, int Q, int T, double pr);
	void ZeroAlpha(int qlo, int qhi);
	void SetOcct(HMM * hmm , int q, double * aqt, double * bqt, double * bq1t, double pr);
	void UpdateTransitionParams(HMM *  hmm, int t, int q, double * aqt , double * bqt , double * bqt1 , double * bq1t , double pr);
	void UpdateParms(HMM *  hmm, ParamAudio *  pa, int t, double * aqt, double * bqt, double pr);
	void UpdateModels(void);
	void FreeMemory(void);
};

