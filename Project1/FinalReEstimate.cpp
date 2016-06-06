#include "FinalReEstimate.h"


FinalReEstimate::FinalReEstimate(void)
{
	pruneThresh = 2000.0;
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
	outprob = GetProbability(pa);
	printf("Beta=%f\n",Beta());
}


void FinalReEstimate::ListHmms(ParamAudio * pa)
{
	hmm_seq = new HMM*[pa->segments];
	seq_num = pa->segments;
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
      
}


float *** FinalReEstimate::GetProbability(ParamAudio * pa)
{
	float *** tmp;
	tmp= new float**[pa->param_frames];
	frames = pa->param_frames;
	int startq,endq;
	for (int t = 0; t < pa->param_frames; t++)
	{
		tmp[t]= new float*[seq_num];
		if(t<pa->param_frames-1)
		{
			startq = qHi[t+1];
			endq = (qLo[t+1]==0?0:((qLo[t]>=qLo[t+1])?qLo[t]:qLo[t+1]-1)); 
		}
		else
		{
			startq = seq_num-1;
			endq = qLo[pa->param_frames-1];
		}

		if(endq>0) endq--;
		for (int i = 0; i < seq_num; i++)
		{
			if(i<=startq && i>=endq)
			{
				tmp[t][i] = new float[hmm_seq[i]->states-2];
				for (int k = 0; k < hmm_seq[i]->states-2; k++)
				{
					tmp[t][i][k] = OutP(pa,t,i,k);
				}
			}
		} 
	}
	return tmp;
}


float FinalReEstimate::OutP(ParamAudio * pa, int fr_number, int segment, int state_nr)
{
	if(pa->coef_num*3!=hmm_seq[segment]->vector_size)
	{
		fprintf(stderr,"FinalReEstimate::OutP:Audio param vector size diffrent then state vector size: %d!=%s \n",pa->coef_num,hmm_seq[segment]->vector_size);
		return 0;
	}

	if(hmm_seq[segment]->states<state_nr-2)
	{
		fprintf(stderr,"FinalReEstimate::OutP: Can't find state: %d \n",state_nr);
		return 0;
	}

	float sum;
	float x;
	State * state = hmm_seq[segment]->state;
	sum = state[state_nr].g_const;


	for (int i = 0; i < pa->coef_num; i++)
	{
		if(pa->coef_first != NULL)
		{
			x= pa->coef_first[fr_number][i] - state[state_nr].mean[i];
			sum += (x*x)/state[state_nr].var[i];
		}

		if(pa->delta_first != NULL)
		{
			x= pa->delta_first[fr_number][i] - state[state_nr].mean[i+pa->coef_num];
			sum += (x*x)/state[state_nr].var[i+pa->coef_num];
		}

		if(pa->acc_first != NULL)
		{
			x= pa->acc_first[fr_number][i] - state[state_nr].mean[i+(pa->coef_num*2)];
			sum += (x*x)/state[state_nr].var[i+(pa->coef_num*2)];
		}
	}
	return -0.5*sum;
}


double FinalReEstimate::Beta(void)
{
	int Q,endq,startq,st,lst=0,q_at_gMax,t,q;
	double *maxP,*bq1t1,*bqt = NULL;
	maxP = new double[seq_num];
	HMM * hmm;
	double a1N=0.0,a,y,x,gMax,lMax;
	beta = new double**[frames];
	for (int i = 0; i < frames; i++)
	{
		beta[i]= new double*[seq_num];
		for (int j = 0; j < seq_num; j++)
		{
			beta[i][j] = new double[hmm_seq[j]->states];
		}
	}
	Q = seq_num-1;
	qHi[frames-1] = Q; endq = qLo[frames-1];
	gMax = LZERO;   q_at_gMax = 0;
	for (q=Q; q>=endq; q--)
	{
		hmm = hmm_seq[q]; st = hmm->states-1;

		beta[frames-1][q][st] = (q==Q)?0.0:beta[frames-1][q+1][lst]+a1N;
		
		for (int i=1;i<st;i++) 
			beta[frames-1][q][i] = hmm->transition[i][st]+beta[frames-1][q][st]; 

		x = LZERO;
		for (int j=1; j<st; j++){
			a = hmm->transition[0][j]; y = beta[frames-1][q][j];
			if (a>LSMALL && y > LSMALL)
				x = LAdd(x,a+outprob[frames-1][q][j-1]+y);
		}
		beta[frames-1][q][1] = x;
		lst = st; a1N = hmm->transition[0][st];
		if (x>gMax) {
			gMax = x; q_at_gMax = q;
		}
	}

	for (t=frames-2;t>=0;t--) {      

      gMax = LZERO;   q_at_gMax = 0;    /* max value of beta at time t */
      startq = qHi[t+1];
      endq = (qLo[t+1]==0)?0:((qLo[t]>=qLo[t+1])?qLo[t]:qLo[t+1]-1);
	  while (endq>0 && hmm_seq[endq-1]->minimum_duration==0) endq--;
      /* start end-point at top of beta beam at t+1  */
      /*  unless this is outside the beam taper.     */
      /*  + 1 to allow for state q+1[1] -> q[N]      */
      /*  + 1 for each tee model preceding endq.     */
    
      
      for (q=startq;q>=endq;q--) {
         lMax = LZERO;                 /* max value of beta in model q */
		 hmm = hmm_seq[q]; 
         st = hmm->states-1;
         bqt = beta[t][q];
         bq1t1 = (q==Q)?NULL:beta[t+1][q+1];
         
         bqt[st] = (bq1t1==NULL)?LZERO:bq1t1[0];
         if (q<startq && a1N>LSMALL)
            bqt[st]=LAdd(bqt[st],beta[t][q+1][lst]+a1N);
         for (int i=st-1;i>0;i--){
            x = hmm->transition[i][st] + bqt[st];
            if (q>=qLo[t+1]&&q<=qHi[t+1])
               for (int j=1;j<st;j++) {
                  a = hmm->transition[i][j]; y = beta[t+1][q][j];
                  if (a>LSMALL && y>LSMALL)
                     x = LAdd(x,a+outprob[t+1][q][j-1]+y);
               }
            bqt[i] = x;
            if (x>lMax) lMax = x;
            if (x>gMax) {
               gMax = x; q_at_gMax = q;
            }
         }
        
         x = LZERO;
         for (int j=1; j<st; j++){
            a = hmm->transition[0][j];
            y = bqt[j];
            if (a>LSMALL && y>LSMALL)
               x = LAdd(x,a+outprob[t][q][j-1]+y);
         }
         bqt[0] = x;
         maxP[q] = lMax;
         lst = st; a1N = hmm->transition[1][st];
      }
      while (gMax-maxP[startq] > pruneThresh) {
         beta[t][startq] = NULL;
         --startq;                   /* lower startq till thresh reached */
         if (startq<0) fprintf(stderr,"SetBeta: Beta prune failed sq < 1\n");
      }
      while(qHi[t]<startq) {        /* On taper */
         beta[t][startq] = NULL;
         --startq;                   /* lower startq till thresh reached */
         if (startq<0) fprintf(stderr,"SetBeta: Beta prune failed on taper sq < 1\n");
      }
     qHi[t] = startq;
      while (gMax-maxP[endq]>pruneThresh){
         beta[t][endq] = NULL;
         ++endq;                   /* raise endq till thresh reached */
         if (endq>startq) {
            return(LZERO);
         }
      }
      qLo[t] = endq;
   
   }

	return bqt[0];
}


double FinalReEstimate::Alpha(void)
{
	return 0;
}


double FinalReEstimate::LAdd(double x, double y)
{
   double temp,diff,z;
   double minLogExp = -log(-LZERO);
   if (x<y) {
      temp = x; x = y; y = temp;
   }
   diff = y-x;
   if (diff<minLogExp) 
      return  (x<LSMALL)?LZERO:x;
   else {
      z = exp(diff);
      return x+log(1.0+z);
   }
}
