#include "FinalReEstimate.h"


FinalReEstimate::FinalReEstimate(void)
{
}


FinalReEstimate::~FinalReEstimate(void)
{
}


void FinalReEstimate::AddHmm(HMM *  hmm)
{
	hmms.push_back(hmm);
	if(hmm->minimum_duration < 0)
		hmm->SetMinDuration();
}


void FinalReEstimate::LoadHmm(std::string hmm_src)
{
	HMM * h = new HMM(hmm_src);
	h->SetMinDuration();
	AddHmm(h);
}
