#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
#include "FinalReEstimate.h"


FinalReEstimate::FinalReEstimate(Config * cf)
{
	this->cf=cf;
	pruneThresh = 2000.0;
    minFrwdP = 10.0; 
	if (cf != nullptr)
	{
		if(cf->Exist("PRUNETHRESH"))
			pruneThresh = cf->GetConfig("PRUNETHRESH");

		if(cf->Exist("MINFRWDP"))
			minFrwdP = cf->GetConfig("MINFRWDP");

		if(cf->Exist("TRACE"))
			trace = cf->GetConfig("TRACE"); 
	}
}


FinalReEstimate::~FinalReEstimate(void)
{
}


void FinalReEstimate::AddHmm(HMM *  hmm)
{
	hmms.push_back(hmm);
	if(hmm->minimum_duration < 0)
		hmm->SetMinDuration();
	
	hmm->ResetOldParams();
}


void FinalReEstimate::LoadHmm(std::string hmm_src)
{
	HMM * h = new HMM(hmm_src, cf);
	h->SetMinDuration();
	AddHmm(h);
}


void FinalReEstimate::ForwardBackward(ParamAudio * pa)
{
	if (trace&TRACE::TOP) 
		printf("Final ReEstimate %d HMM's with data: %s\n", hmms.size(), std::string(pa->audio_src).c_str());
	ListHmms(pa);
	SetBeamTaper(pa->segments,pa->param_frames);
	outprob = GetProbability(pa);
	pr = Beta();
	Alpha(pa);
	FreeMemory();
}


void FinalReEstimate::ListHmms(ParamAudio * pa)
{
	hmm_seq = new HMM*[pa->segments];
	seq_num = pa->segments;
	int j;
	for (int i = 0; i < pa->segments; i++)
	{
		
		for (j = 0; j < hmms.size(); j++)
		{
			
			if(hmms[j]->name.compare(pa->os[i].l->name)==0)
			{
				hmm_seq[i]=hmms[j];
				break;
			}
			
			
		}
		if(j == hmms.size())
			fprintf(stderr,"Can't find hmm:%s\n",std::string(pa->os[i].l->name).c_str());
	}
}


void FinalReEstimate::SetBeamTaper(int Q, int T)
{
	int q,dq,i,t;
	qHi = new int[T];
	qLo = new int[T];
   /* Set leading taper */
	q=0;
	dq=hmm_seq[q]->minimum_duration;
	i=0;
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
			else
			{
				tmp[t][i] = NULL;
			}
		} 
	}
	return tmp;
}


float FinalReEstimate::OutP(ParamAudio * pa, int fr_number, int segment, int state_nr)
{
	if(pa->coef_num*3!=hmm_seq[segment]->vector_size)
	{
		fprintf(stderr,"FinalReEstimate::OutP:Audio param vector size diffrent then state vector size: %d!=%d \n",pa->coef_num,hmm_seq[segment]->vector_size);
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
	for (int i = 0; i < seq_num; i++)
	{
		maxP[i]=0.0;
	}
	HMM * hmm;
	double a1N=0.0,a,y,x,gMax,lMax;
	beta = new double**[frames];
	for (int i = 0; i < frames; i++)
	{
		beta[i]= new double*[seq_num];
		for (int j = 0; j < seq_num; j++)
		{
			beta[i][j] = new double[hmm_seq[j]->states];
			for (int k = 0; k < hmm_seq[j]->states; k++)
			{
				beta[i][j][k] = LZERO;
			}
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
		 delete[] beta[t][startq];
         beta[t][startq] = NULL;
         --startq;                   /* lower startq till thresh reached */
         if (startq<0) fprintf(stderr,"SetBeta: Beta prune failed sq < 1\n");
      }
      while(qHi[t]<startq) {        /* On taper */
         delete[] beta[t][startq];
		 beta[t][startq] = NULL;
         --startq;                   /* lower startq till thresh reached */
         if (startq<0) fprintf(stderr,"SetBeta: Beta prune failed on taper sq < 1\n");
      }
     qHi[t] = startq;
      while (gMax-maxP[endq]>pruneThresh){

		 delete[] beta[t][endq];
         beta[t][endq] = NULL;
         ++endq;                   /* raise endq till thresh reached */
         if (endq>startq) {
            return(LZERO);
         }
      }
      qLo[t] = endq;
   
   }
	delete[] maxP;
	if (trace&TRACE::DEEP) {
        printf(" Utterance prob per frame = %e\n",bqt[0]/frames);
   }
	return bqt[0];
}


void FinalReEstimate::Alpha(ParamAudio * pa)
{
	int endq,st,startq;
	HMM * hmm;
	double x,a,a1N=0.0;
	double * aq, ** tmp, *laq;
	int sq,eq,i,j,q,Nq,lNq,maxn=0;
	double y;
	double *aqt, * bqt, *bq1t, *bqt1, *aqt1;

	
	
	alphat= new double*[seq_num];
	for (int j = 0; j < seq_num; j++)
	{
		alphat[j] = new double[hmm_seq[j]->states];
	}
	alphat1= new double*[seq_num];
	for (int j = 0; j < seq_num; j++)
	{
		alphat1[j] = new double[hmm_seq[j]->states];
	}
	for (int i = 0; i < hmms.size(); i++)
	{
		if (hmms[i]->states > maxn)
		{
			maxn = hmms[i]->states;
		}
	}
	occt = new float[maxn];


	startq = 0;
	endq = qHi[0];
   for (int q=0; q<=endq; q++)
   {
	  hmm = hmm_seq[q]; 
	  st = hmm->states-1;
      aq = alphat[q];
      aq[0] = (q==0)?0.0:alphat[q-1][0]+a1N;
      
      for (int j=1;j<st;j++) {
		  a = hmm->transition[0][j];
         aq[j] = (a>LSMALL)?aq[0]+a+outprob[0][q][j-1]:LZERO;
      }
      x = LZERO;
      for (int i=1;i<st;i++) {
		  a = hmm->transition[i][st];
         if (a>LSMALL)
            x = LAdd(x,aq[i]+a);
      }
      aq[st] = x;
      a1N = hmm->transition[0][st];
   }


	for (int t=0;t<frames;t++) {


      

      if (t>0)
		  StepAlpha(t,&startq,&endq,seq_num-1,frames-1,pr);
    
      for (int q=startq;q<=endq;q++) { 
         /* increment accs for each active model */
         hmm = hmm_seq[q];
         //up_hmm = up_qList[q];
         aqt = alphat[q];
         bqt = beta[t][q];

         bqt1 = (t==frames-1) ? NULL:beta[t+1][q];
         aqt1 = (t==1)      ? NULL:alphat1[q];
		 bq1t = (q==seq_num-1) ? NULL:beta[t][q+1];
         SetOcct(hmm,q,aqt,bqt,bq1t,pr);
         
         UpdateParms(hmm,pa,t,aqt,bqt,pr);
   
         UpdateTransitionParams(hmm,t,q,aqt,bqt,bqt1,bq1t,pr);
      } 
    
   }
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


double FinalReEstimate::MaxModelProb(int q, int t, int minq)
{
	double *aq,*bq,*bq1;
   double maxP,x;
   int Nq1,Nq,i,qx,qx1;
   HMM * hmm;
   
   if (q==0)
      maxP = LZERO;
   else {
	   bq1 = beta[t][q-1]; Nq1 = hmm_seq[q-1]->states-1;
      maxP = (bq1==NULL)?LZERO:alphat[q-1][Nq1] + bq1[Nq1];
	  for (qx=q-1;qx>minq && hmm_seq[qx]->transition[1][Nq1] > LSMALL;qx--){
         qx1 = qx-1;
		 bq1 = beta[t][qx1]; Nq1 = hmm_seq[qx1]->states-1;
         x=(bq1==NULL)?LZERO:alphat[qx1][Nq1]+bq1[Nq1];
         if (x > maxP) maxP = x;
      }
   }
   hmm = hmm_seq[q]; Nq = hmm->states-1;   
   bq=beta[t][q];
   if (bq != NULL) {
      aq = alphat[q]; 
      for (i=0;i<Nq;i++)
         if ((x=aq[i]+bq[i]) > maxP) maxP = x;
   }
   return maxP;
}


void FinalReEstimate::StepAlpha(int t, int * start, int * end, int Q, int T, double pr)
{
   double *aq,*laq,**tmp, **alphat,**alphat1;
   int sq,eq,i,j,q,Nq,lNq;
   double x=0.0,y,a,a1N=0.0;
   HMM * hmm;
   
   alphat  = this->alphat;
   alphat1 = this->alphat1;

   /* First prune beta beam further to get alpha beam */
   
   sq = qLo[t-1];    /* start start-point at bottom of beta beam at t-1 */

   while (pr-MaxModelProb(sq,t-1,sq)>minFrwdP){
      ++sq;                /* raise start point */
      if (sq>qHi[t]) 
         fprintf(stderr,"StepAlpha: Alpha prune failed sq(%d) > qHi(%d)\n",sq,qHi[t]);
   }
   if (sq<qLo[t])       /* start-point below beta beam so pull it back */
      sq = qLo[t];
   
   eq = qHi[t-1]<Q?qHi[t-1]+1:qHi[t-1];
   /* start end-point at top of beta beam at t-1  */
   /* JJO : + 1 to allow for state q-1[N] -> q[1] */
   /*       + 1 for each tee model following eq.  */
   while (pr-MaxModelProb(eq,t-1,sq)>minFrwdP){
      --eq;             /* lower end-point */
      if (eq<sq) 
         fprintf(stderr,"StepAlpha: Alpha prune failed eq(%d) < sq(%d)\n",eq,sq);
   }
   while (eq<Q && hmm_seq[eq]->minimum_duration==0) eq++;
   if (eq>qHi[t])  /* end point above beta beam so pull it back */
      eq = qHi[t]; 
      
   
   
   /* Now compute current alpha column */
   tmp = this->alphat1; this->alphat1 = this->alphat; this->alphat = tmp;
   alphat  = this->alphat;
   alphat1 = this->alphat1;

   if (sq>0) ZeroAlpha(0,sq-1);

   Nq = (sq == 0) ? 0:hmm_seq[sq-1]->states-1;

   for (q = sq; q <= eq; q++) {
      lNq = Nq; hmm = hmm_seq[q]; Nq = hmm->states-1; 
      aq = alphat[q]; 
      laq = alphat1[q];
      if (q==0)
         aq[0] = LZERO;
      else{
         aq[0] = alphat1[q-1][lNq];
         if (q>sq && a1N>LSMALL) /* tee Model */
            aq[0] = LAdd(aq[0],alphat[q-1][0]+a1N);
      }
      for (j=1;j<Nq;j++) {
         a = hmm->transition[0][j];
         x = (a>LSMALL)?a+aq[0]:LZERO;
         for (i=1;i<Nq;i++){
            a = hmm->transition[i][j]; y = laq[i];
            if (a>LSMALL && y>LSMALL)
               x = LAdd(x,y+a);
         }
         aq[j] = x + outprob[t][q][j-1];
      }
      x = LZERO;
      for (i=2;i<Nq;i++){
         a = hmm->transition[i][Nq]; y = aq[i];
         if (a>LSMALL && y>LSMALL)
            x = LAdd(x,y+a);
      }
      aq[Nq] = x; a1N = hmm->transition[0][Nq];
   }
   if (eq<Q) ZeroAlpha(eq+1,Q);

   
   if (t==T){
      if (fabs((x-pr)/T) > 0.001)
         fprintf(stderr,"StepAlpha: Forward/Backward Disagree %f/%f\n",x,pr);
      
   }

   *start=sq; *end=eq;
}


void FinalReEstimate::ZeroAlpha(int qlo, int qhi)
{
	HMM * hmm;
   int Nq,j,q;
   double *aq;
   
   for (q=qlo;q<=qhi;q++) {   
	   hmm = hmm_seq[q]; 
      Nq = hmm->states-1; 
      aq = alphat[q];
      for (j=0;j<=Nq;j++)
         aq[j] = LZERO;
   }
}


void FinalReEstimate::SetOcct(HMM * hmm , int q, double * aqt, double * bqt, double * bq1t, double pr)
{
	int i,N;
   double x;
 
   
   N=hmm->states-1;
   for (i=0;i<=N;i++) {
		x =aqt[i] + bqt[i];
      if (i==1 && bq1t != NULL && hmm->transition[0][N] > LSMALL)
		  x = LAdd(x,aqt[0]+bq1t[0]+hmm->transition[0][N]);
      x -= pr;
      occt[i] = (x>MINEARG) ? exp(x) : 0.0;
   }
   
}


void FinalReEstimate::UpdateTransitionParams(HMM *  hmm, int t, int q, double * aqt , double * bqt , double * bqt1 , double * bq1t , double pr)
{
   int i,j,N;
   float *ti, *ai;
   float *outprob1;
   double sum,x;
   

   N = hmm->states-1;
   if (bqt1!=NULL) outprob1 = outprob[t+1][q];  /* Bug fix */
   else outprob1 = NULL;
   for (i=0;i<N;i++)
      hmm->obs_count[i] += occt[i];
   for (i=0;i<N;i++) {
	   ti = hmm->trans_count[i]; ai = hmm->transition[i];
      for (j=1;j<=N;j++) {
         if (i==0 && j<N) {                  /* entry transition */
            x = aqt[0]+ai[j]+outprob[t][q][j-1]+bqt[j]-pr;
            if (x>MINEARG) ti[j] += exp(x);
         } else
            if (i>0 && j<N && bqt1!=NULL) {     /* internal transition */
               x = aqt[i]+ai[j]+outprob1[j-1]+bqt1[j]-pr;
               if (x>MINEARG) ti[j] += exp(x);
            } else
               if (i>0 && j==N) {                  /* exit transition */
                  x = aqt[i]+ai[N]+bqt[N]-pr;
                  if (x>MINEARG) ti[N] += exp(x);
               }
         if (i==0 && j==N && ai[N]>LSMALL && bq1t != NULL){ /* tee transition */
            x = aqt[0]+ai[N]+bq1t[0]-pr;
            if (x>MINEARG) ti[N] += exp(x);
         }
      }
   }
}


void FinalReEstimate::UpdateParms(HMM *  hmm, ParamAudio *  pa, int t, double * aqt, double * bqt, double pr)
{
   int N;
   double x, Lr,z ;
   float *var, *mn;
   N = hmm->states-1;
   for (int j=1;j<N;j++) 
   {
 	  if(bqt == NULL)
		x = aqt[j]-pr;
	  else if(aqt == NULL)
		x = bqt[j]-pr;
	  else
		x =aqt[j] + bqt[j] - pr;
         if (-x<minFrwdP) 
		 {
                 Lr = exp(x);
                 hmm->state[j-1].obs_count += Lr;
				 mn = hmm->state[j-1].old_mean;
				 var = hmm->state[j-1].old_var;
				 for (int k=0;k<pa->coef_num;k++) 
				 {
					x=pa->coef_first[t][k]-hmm->state[j-1].mean[k];
					mn[k] += x;
                    var[k] += x*x*Lr;

					x=pa->delta_first[t][k]-hmm->state[j-1].mean[k+pa->coef_num];
					mn[k+pa->coef_num] += x;
                    var[k+pa->coef_num] += x*x*Lr;

					x=pa->acc_first[t][k]-hmm->state[j-1].mean[k+(pa->coef_num*2)];
					mn[k+(pa->coef_num*2)] += x;
                    var[k+(pa->coef_num*2)] += x*x*Lr;
				 }
          }
      }
}



void FinalReEstimate::UpdateModels(void)
{
	for (int i = 0; i < hmms.size(); i++)
	{
		hmms[i]->UpdateMean();
		hmms[i]->UpdateVar();
		hmms[i]->UpdateTransition();
	}
}


void FinalReEstimate::FreeMemory(void)
{

	for (int j = 0; j < seq_num; j++)
	{
		delete[] alphat[j];
	}

	for (int j = 0; j < seq_num; j++)
	{
		delete[] alphat1[j];
	}

	for (int i = 0; i < frames; i++)
	{
		for (int j = 0; j < seq_num; j++)
		{
			if (outprob[i][j]!=NULL)
			{
				delete[] outprob[i][j];
			}
			
		}
		delete[] outprob[i];
	}

	
	for (int i = 0; i < frames; i++)
	{
		
		for (int j = 0; j < seq_num; j++)
		{
			delete[] beta[i][j];
		}
		delete beta[i];
	}

	delete[] beta;
	delete[] outprob;
	delete[] alphat;
	delete[] alphat1;
	delete[] occt;
	delete[] qHi;
	delete[] qLo;
	delete[] hmm_seq;
}
