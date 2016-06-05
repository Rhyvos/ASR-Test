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


void FinalReEstimate::ForwardBackward(ParamAudio * pa)
{
	ListHmms(pa);
	SetBeamTaper(pa->segments,pa->param_frames);
}


void FinalReEstimate::ListHmms(ParamAudio * pa)
{
	hmm_seq = new HMM*[pa->segments];

	for (int i = 0; i < pa->segments; i++)
	{
		
		for (int j = 0; j < hmms.size(); j++)
		{
			
			if(hmms[j]->name.compare(pa->os[i].l->name)==0)
			{
				hmm_seq[i]=hmms[j];
				break;
			}
			fprintf(stderr,"Can't find hmm:%s",pa->os[i].l->name);
		}
	}
}


void FinalReEstimate::SetBeamTaper(int Q, int T)
{
	int q,dq,i,t;
	qHi = new int[T];
	qLo = new int[T];
   /* Set leading taper */
	q=0;dq=hmm_seq[q]->minimum_duration;i=0;
   for (t=0;t<T;t++) {
      while (i==dq) {
         i=0;
         if (q<Q-1)
		 {
			q++;
			dq=hmm_seq[q]->minimum_duration;
		 }
         else dq=-1;
      }
      qHi[t]=q;
      i++;
   }
   q=Q-1;dq=hmm_seq[q]->minimum_duration;i=0;
   for (t=T-1;t>=0;t--) {
      while (i==dq) {
         i=0;
         if (q>0) q--,dq=hmm_seq[q]->minimum_duration;
         else dq=-1;
      }
      qLo[t]=q;
      i++;
   }
      
         for (t=0;t<T;t++)
         printf("%d: %d to %d\n",t,qLo[t],qHi[t]);
         
}
