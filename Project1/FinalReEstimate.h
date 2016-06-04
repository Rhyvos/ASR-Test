#pragma once
#include <vector>
#include "HMM.h"
class FinalReEstimate
{
public:
	FinalReEstimate(void);
	~FinalReEstimate(void);
	std::vector<HMM *> hmms;
	void AddHmm(HMM *  hmm);
	void LoadHmm(std::string hmm_src);
};

