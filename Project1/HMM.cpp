#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif

#define _USE_MATH_DEFINES
#include "HMM.h"
#include <cmath>
#include <algorithm> 
#include <fstream>
#include <iomanip>
#include <sstream>
HMM::HMM(void)
{

}


HMM::~HMM(void)
{
	for(int i=0;i<states-2;i++)
	{
		delete[] state[i].mean;
		delete[] state[i].var;
		delete[] state[i].old_mean;
		delete[] state[i].old_var;
	}
	delete[] state;
	for(int i=0;i<states;i++)
	{
		delete[] transition[i];
		delete[] trans_count[i];
	}
	delete[] transition;
	delete[] trans_count;
	delete[] obs_count;
}


HMM::HMM(std::string hmm_src) //load from file
{
}


HMM::HMM(int vector_size, int states, std::string l_name): vector_size(vector_size), states(states), name(l_name)
{

	if(states<3 || vector_size<1)
	{
		fprintf(stderr, "HMM: Wrong initialisaton parametrs: vector_size: %d or states: %d\n",vector_size,states);
		return;
	}
	
	state = new State[states-2];				// we don't need first and last state mean and variance
	minVar  = 1.0E-2;
	epsilon = 1.0E-4;
	minLogExp = -log(-LZERO);
	minimum_duration =-1;
	for(int i=0 ;i<states-2;i++)
	{
		state[i].mean = new float[vector_size];
		state[i].var = new float[vector_size];
		state[i].old_mean = new float[vector_size];
		state[i].old_var = new float[vector_size];
		for(int j=0;j<vector_size;j++)
		{
			state[i].mean[j]=0.0;
			state[i].var[j]=1.0;
		}
		state[i].state_nr=i+1;
	}

	transition = new float*[states];
	trans_count = new float*[states];
	obs_count = new float[states];
	for(int i=0 ;i<states;i++)
	{
		transition[i] = new float[states];
		trans_count[i] = new float[states];
		for(int j=0;j<states;j++)
			transition[i][j] = LZERO;
		if(i!=0 && i!=(states-1))
		{
			transition[i][i] = log(0.6);
			transition[i][i+1] = log(0.4);
		}
		if(i==0)
			transition[i][i+1] = log(1.0); 

	}




}


void HMM::Initialise(ParamAudio * pa, int iteration)
{
	float totalP,newP,delta;
	bool converged = false;   
	int * states_vec;
	int * mixes;
	int j;
	GetMean(pa);
	GetVariance(pa);
	FindGConst();
	totalP = LZERO;
	for (int i = 0; !converged && i < iteration; i++)
	{
		ResetOldParams();
		for (j = 0, newP = 0; j < pa->segments; j++)
		{
			if(!pa->os[j].l->name.compare(name) && pa->os[j].frames >= states-2)
			{
				states_vec = new int[pa->os[j].frames];
				mixes = new int[pa->os[j].frames];

				newP += ViterbiAlg(pa->os+j,states_vec,mixes);
				UpdateCounts(pa->os+j,states_vec);
				delete[] mixes;
				delete[] states_vec;
			}

			if(pa->os[j].frames < states-2)
				fprintf(stderr,"HMM::Initialise():Segment(%d) dont have enough frames(%d), minimum frames=%d\n",j,pa->os[j].frames,states-2);
		}
		newP /= pa->segments;
		delta = newP - totalP;
		converged = (i>0) && (fabs(delta) < epsilon);
		if (!converged)
		{	
			UpdateMean();
			UpdateVar();
			UpdateTransition();
		}
		totalP = newP;

	}

}


void HMM::GetMean(ParamAudio * pa)
{
	int statenumber=0;
	float frames_per_state;
	int *state_fr_counter=new int[states-2];
	for (int i = 0; i < states-2; i++)
	{
		state_fr_counter[i]=0;
	}


	for(int i=0;i<pa->segments;i++)
	{
		if(pa->os[i].l == NULL || !pa->os[i].l->name.compare(name))
		{
			frames_per_state = (float)pa->os[i].frames/(float)(states-2);
			for (int j= 0; j < pa->os[i].frames; j++)
			{
				statenumber = (int)((float)j/(float)(frames_per_state));
				for (int k = 0; k < pa->os[i].frame_lenght; k++)
				{
					if(pa->os[i].coef != NULL)
						state[statenumber].mean[k] += (pa->os[i].coef[j][k]);
					if(pa->os[i].delta != NULL)
						state[statenumber].mean[k+pa->os[i].frame_lenght] += (pa->os[i].delta[j][k]);
					if(pa->os[i].acc != NULL)
						state[statenumber].mean[k+(2*pa->os[i].frame_lenght)] += (pa->os[i].acc[j][k]);
				}
				state_fr_counter[statenumber]++;
			}

		}

	}

	for (int i = 0; i < states - 2; i++)
	{
		for (int j = 0; j < vector_size; j++)
		{
			state[i].mean[j] /= state_fr_counter[i];
		}
	}

	delete[] state_fr_counter;
	
}




void HMM::GetVariance(ParamAudio * pa)
{
	int statenumber=0;
	float frames_per_state;
	int *state_fr_counter=new int[states-2];
	double ** temp_var = new double*[states-2];
	double x;
	for (int i = 0; i < states-2; i++)
	{
		state_fr_counter[i]=0;
		temp_var[i] = new double[vector_size];
		for (int j = 0; j < vector_size; j++)
		{
			temp_var[i][j]=0.0;
		}
	}

	

	for(int i=0;i<pa->segments;i++)
	{
		if(pa->os[i].l == NULL || !pa->os[i].l->name.compare(name))
		{
			frames_per_state = (float)pa->os[i].frames/(float)(states-2);
			for (int j= 0; j < pa->os[i].frames; j++)
			{
				statenumber = (int)((float)j/(float)(frames_per_state));
				for (int k = 0; k < pa->os[i].frame_lenght; k++)
				{
					if(pa->os[i].coef != NULL)
					{
						x=(pa->os[i].coef[j][k])-state[statenumber].mean[k];
						temp_var[statenumber][k] += x*x;
					}

						
					if(pa->os[i].delta != NULL)
					{
						x=(pa->os[i].delta[j][k])-state[statenumber].mean[k+pa->os[i].frame_lenght] ;
						temp_var[statenumber][k+pa->os[i].frame_lenght] += x*x;
					}

						
					if(pa->os[i].acc != NULL)
					{
						x=(pa->os[i].acc[j][k])-state[statenumber].mean[k+(2*pa->os[i].frame_lenght)];
						temp_var[statenumber][k+(2*pa->os[i].frame_lenght)] += x*x;
					}
				}
				state_fr_counter[statenumber]++;
			}

		}

	}

	for (int i = 0; i < states - 2; i++)
	{
		for (int j = 0; j < vector_size; j++)
		{
			temp_var[i][j]/=state_fr_counter[i];
			state[i].var[j] = temp_var[i][j]<minVar?minVar:temp_var[i][j];
		}
	}

	for (int i = 0; i < states-2; i++)
	{
			delete[] temp_var[i];
	}
	delete[] state_fr_counter;
	delete[] temp_var;

}


void HMM::FindGConst(void)
{
	float z;
	float sum;
	for(int i=0 ;i<states-2;i++)
	{
		sum = vector_size*log(2*M_PI);
		for (int j = 0; j < vector_size; j++)
		{
			z = (state[i].var[j]<=MINLARG)?LZERO:log(state[i].var[j]);
			sum += z;
		}
		state[i].g_const=sum;
	}

}


float HMM::OutP(ObservationSegment * os, int fr_number, int state_nr)
{
	if(os->frame_lenght*3!=vector_size)
	{
		fprintf(stderr,"HMM.OpuP:Audio param vector size diffrent then state vector size: %d!=%s \n",os->frame_lenght*3,vector_size);
		return 0;
	}

	float sum;
	float x;
	sum = state[state_nr-1].g_const;


	for (int i = 0; i < os->frame_lenght; i++)
	{
		if(os->coef != NULL)
		{
			x= os->coef[fr_number][i] - state[state_nr-1].mean[i];
			sum += (x*x)/state[state_nr-1].var[i];
		}

		if(os->delta != NULL)
		{
			x= os->delta[fr_number][i] - state[state_nr-1].mean[i+os->frame_lenght];
			sum += (x*x)/state[state_nr-1].var[i+os->frame_lenght];
		}

		if(os->acc != NULL)
		{
			x= os->acc[fr_number][i] - state[state_nr-1].mean[i+(os->frame_lenght*2)];
			sum += (x*x)/state[state_nr-1].var[i+(os->frame_lenght*2)];
		}
	}
	return -0.5*sum;
}


float HMM::ViterbiAlg(ObservationSegment * os, int * states_vec, int * mixes)
{

   int currState,prevState,bestPrevState;
   int frame_id;
   float  bestP,currP,tranP,prevP;
   short ** trace_back;
   float *lastP;
   float *thisP;

   trace_back = new short*[os->frames];

   for (int i=0; i < os->frames; i++)
   {
		trace_back[i] = new short[states];
   }   

   lastP = new float[states];
   thisP = new float[states];

 
   
   for (currState=1; currState < states-1; currState++) {
      tranP = transition[0][currState];
	 
      if (tranP<LSMALL) 
         lastP[currState] = LZERO;
      else
         lastP[currState] = tranP + OutP(os,0,currState);
      trace_back[0][currState] = 0;
   }
   
  
  
   for (frame_id=1; frame_id<os->frames; frame_id++) {
      
         
      for (currState=1;currState<states-1;currState++) {

         bestPrevState=1;
         tranP = transition[1][currState]; 
		 prevP = lastP[1];
         bestP = (tranP<LSMALL) ? LZERO : tranP+prevP;

         for (prevState=2;prevState<states-1;prevState++) {
            tranP = transition[prevState][currState];
		
            prevP = lastP[prevState];
            currP = (tranP<LSMALL) ? LZERO : tranP+prevP;
            if (currP > bestP) {
               bestPrevState=prevState; bestP=currP;
            }
         }
         if (bestP<LSMALL)
            currP = thisP[currState] = LZERO;
         else {
            currP = OutP(os,frame_id,currState);
            thisP[currState] = bestP+currP;
         }
        
         trace_back[frame_id][currState]=bestPrevState;
      }
	  for (int i = 0; i < states; i++)
	  {
		  lastP[i]=thisP[i];
	  }
   }
   

   bestPrevState=1;
   tranP = transition[1][states-1]; prevP = lastP[1];
   
   bestP=(tranP<LSMALL) ? LZERO : tranP+prevP;
   for (prevState=2;prevState<states-1;prevState++) {
      tranP = transition[prevState][states-1]; prevP = lastP[prevState];
	  
      currP = (tranP<LSMALL) ? LZERO : tranP+prevP;
      if (currP > bestP) {
         bestPrevState=prevState; bestP=currP;
      }
   }  

   
   ComputeTraceBack(os->frames,bestPrevState,states_vec,trace_back);
	for (int i = 0; i < os->frames; i++)
	{
		delete[] trace_back[i];
	}
   delete[] trace_back;
   delete[] lastP;
   delete[] thisP;
  
   return bestP;
}


void HMM::ComputeTraceBack(int frames_in_seg, int state, int * states_vec, short ** trace_back)
{
   for (int i=frames_in_seg-1; i>=0; i--) {
      states_vec[i] = state;
      state=trace_back[i][state];
   }
}




void HMM::UpdateCounts(ObservationSegment * os, int * state_vec)
{
	int state_nr, last;
	float x;
	last = 0;
	for (int i = 0; i < os->frames; i++)
	{
		state_nr = state_vec[i];
		state_nr--;
		state[state_nr].obs_count++;
		for (int j = 0; j < os->frame_lenght; j++)
		{
			x = os->coef[i][j] - state[state_nr].mean[j];
			state[state_nr].old_mean[j] += x;
			state[state_nr].old_var[j] += x*x;

			x = os->delta[i][j] - state[state_nr].mean[j+os->frame_lenght];
			state[state_nr].old_mean[j+os->frame_lenght] += x;
			state[state_nr].old_var[j+os->frame_lenght] += x*x;

			x = os->acc[i][j] - state[state_nr].mean[j+(os->frame_lenght*2)];
			state[state_nr].old_mean[j+(os->frame_lenght*2)] += x;
			state[state_nr].old_var[j+(os->frame_lenght*2)] += x*x;
		}

		obs_count[last] += 1.0;
		trans_count[last][state_nr+1] += 1.0;
		last = state_nr+1;
        if (i==os->frames-1){ 
           obs_count[state_nr+1] += 1.0;
           trans_count[state_nr+1][states-1] += 1.0;
        }
	}
}


void HMM::ResetOldParams(void)
{

	for(int i=0 ;i<states-2;i++)
	{
		for(int j=0;j<vector_size;j++)
		{
			state[i].old_mean[j]=0.0;
			state[i].old_var[j]=0.0;
		}
		state[i].obs_count=0;
	}

	for(int i=0 ;i<states;i++)
	{
		obs_count[i] = 0.0;
		for(int j=0;j<states;j++)
			trans_count[i][j] = 0;
	}

}



void HMM::UpdateMean(void)
{
	float x;
	for (int i = 0; i < states-2; i++)
	{
		for (int j = 0; j < vector_size; j++)
		{
			x= state[i].mean[j] + state[i].old_mean[j] / state[i].obs_count;
			state[i].old_mean[j] = state[i].mean[j];
			state[i].mean[j] = x;
		}
	}
}


void HMM::UpdateVar(void)
{
	float x,z;
	for (int i = 0; i < states-2; i++)
	{
		for (int j = 0; j < vector_size; j++)
		{
			x = state[i].mean[j]-state[i].old_mean[j];
			z = state[i].old_var[j]/state[i].obs_count - x*x;
			state[i].var[j] = (z<minVar)?minVar:z;
		}
	}
	FindGConst();
}


void HMM::UpdateTransition(void)
{
	float sum,x;
	for (int i = 0; i < states-1; i++)
	{
		transition[i][0] = LZERO;
		sum = 0;
		for (int j = 1; j < states; j++)
		{
			x=trans_count[i][j]/obs_count[i];
			transition[i][j] = x;
			sum+=x;
		}
		for (int j = 1; j < states; j++)
		{
			x = transition[i][j]/sum;
			transition[i][j] = (x<MINLARG) ? LZERO : log(x);
		}

	}
}


void HMM::ReEstimate(ParamAudio * pa, int iterations)
{
	float ** prob = NULL;
	int max_obs=0,ntu;
	double ap,bp;
	float spr,newP,delta,totalP;
	bool converged = false;
	for (int i = 0; i < pa->segments; i++)
	{
		if(max_obs<pa->os[i].frames)
			max_obs = pa->os[i].frames;
	}

	alpha = new double*[states];
	beta = new double*[states];

	for (int i = 0; i < states; i++)
	{
		alpha[i] = new double[max_obs];
		beta[i] = new double[max_obs];
	}

	totalP = LZERO;

	for (int i = 0; !converged && i < iterations; i++)
	{
		ResetOldParams();
		ntu=0;
		newP=0.0;
		for (int j = 0; j < pa->segments; j++)
		{
			if(!pa->os[j].l->name.compare(name))
			{
				prob = GetProbability(pa->os+j);
				
				if ((ap = GetAlpha(prob,pa->os[j].frames))>LSMALL)
				{
					bp = GetBeta(prob,pa->os[j].frames);
					spr= (ap + bp) / 2.0;
					newP += spr; 
					ntu++;
					UpdateRestCount(prob,pa->os+j,spr);
				}
				if(prob != NULL)
				{
					for (int i = 0; i < states-2; i++)
					{
						if(prob[i])
							delete[] prob[i];
					}

					delete[] prob; 
				} 
			}
		}
	
		UpdateMean();
		UpdateVar();
		UpdateTransition();

		newP /= ntu;
		delta = newP - totalP;
		converged = (i>0) && (fabs(delta) < epsilon);
		
		totalP = newP;

		printf("Ave LogProb at iter %d = %10.5f using %d examples",
			i,totalP,ntu);
         if (i >= 1)
            printf("  change = %10.5f",delta);
         printf("\n");

	}

	for (int i = 0; i < states; i++)
	{
		delete[] alpha[i];
		delete[] beta[i];
	}
	delete[] alpha;
	delete[] beta;
	

}


float ** HMM::GetProbability(ObservationSegment * os)
{
	float ** prob;
	prob = new float*[states-2];
	for (int i = 1; i < states-1; i++)
	{
		prob[i-1] = new float[os->frames];
		for (int j = 0; j < os->frames; j++)
		{
			prob[i-1][j] = OutP(os,j,i);
		}
	}

	return prob;
}


double HMM::GetAlpha(float ** prob, int frames)
{

   
   double x,a;

   alpha[0][0] = 0.0;
   for (int j=1;j<states-1;j++) {             
      a=transition[0][j];
      if (a<LSMALL)
         alpha[j][0] = LZERO;
      else
         alpha[j][0] = a+prob[j-1][0];
   }
   alpha[states-1][0] = LZERO;
   
   for (int t=1;t<frames;t++) {           
      for (int j=1;j<states-1;j++) {
         x=LZERO ;
         for (int i=1;i<states-1;i++) {
            a=transition[i][j];
            if (a>LSMALL)
               x = LAdd(x,alpha[i][t-1]+a);
         }
         alpha[j][t]=x+prob[j-1][t];
      }
      alpha[0][t] = alpha[states-1][t] = LZERO;
   }
   x = LZERO ;                      
   for (int i=1;i<states-1;i++) {
      a=transition[i][states-1];
      if (a>LSMALL)
         x=LAdd(x,alpha[i][frames-1]+a); 
   }  
   alpha[states-1][frames-1] = x;
   
   
   return x;
}


double HMM::LAdd(double x, double y)
{
   double temp,diff,z;
   
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


double HMM::GetBeta(float ** prob, int frames)
{
   double x,a;

   beta[states-1][frames-1] = 0.0;
   for (int i=1;i<states-1;i++)                
      beta[i][frames-1]=transition[i][states-1];
   beta[0][frames-1] = LZERO;
   for (int t=frames-2;t>=0;t--) {           
      for (int i=0;i<states;i++)
         beta[i][t]=LZERO ;
      for (int j=1;j<states-1;j++) {
         x=prob[j-1][t+1]+beta[j][t+1];
         for (int i=1;i<states-1;i++) {
            a=transition[i][j];
            if (a>LSMALL)
               beta[i][t]=LAdd(beta[i][t],x+a);
         }
      }
   } 
   x=LZERO ;
   for (int j=1;j<states-1;j++) {
      a=transition[0][j];
      if (a>LSMALL)
         x=LAdd(x,beta[j][0]+a+prob[j-1][0]); 
   }
   beta[0][0] = x;
  
   return x;
}


void HMM::UpdateRestCount(float ** prob,ObservationSegment * os,double pr)
{
	
   double x,y,Lr;
   float * occr = new float[states-1];
   occr[0] = 1.0;

   for (int i=1;i<states-1;i++) {
      x=LZERO ;
      for (int t=0;t<os->frames;t++)
         x=LAdd(x,alpha[i][t]+beta[i][t]);
      x -= pr;
      if (x>MINEARG) 
         occr[i] = exp(x);
      else
         occr[i] = 0.0;
   }

   for (int i=1; i<states-1; i++)
	   obs_count[i] += occr[i];

   for (int j=1;j<states-1;j++) {

	  if (transition[0][j]>LSMALL) {
         x = transition[0][j] + prob[j-1][0] + beta[j][0] - pr;
         if (x>MINEARG) {
            y =  exp(x);
            trans_count[0][j] += y; obs_count[0] += y;
         }
      }
   }
   for (int i=1;i<states-1;i++) {       
      for (int j=1;j<states-1;j++) {
         if (transition[i][j]>LSMALL) {
            x=LZERO; 
            for (int t=0;t<os->frames-1;t++)
               x=LAdd(x,alpha[i][t]+transition[i][j]+prob[j-1][t+1]+beta[j][t+1]);
            x -= pr;
            if (x>MINEARG)
               trans_count[i][j] += exp(x);
         }
      }
   }
   for (int i=1; i<states-1; i++) {  
	  if (transition[i][states-1]>LSMALL) {
         x = transition[i][states-1] + alpha[i][os->frames-1] - pr;
         if (x>MINEARG)
            trans_count[i][states-1] += exp(x);
      }     
   }

   for (int i = 1; i < states -1; i++)
   {
	   for (int j = 0; j < os->frames; j++)
	   {
		   Lr = alpha[i][j]+beta[i][j] - pr;
		   y = exp(Lr);
		   state[i-1].obs_count+=y;
		   for (int k = 0; k < os->frame_lenght; k++)
		   {
			    x = os->coef[j][k] - state[i-1].mean[k];
				state[i-1].old_mean[k] += x * y;
				state[i-1].old_var[k] += x*x * y;

				x = os->delta[j][k] - state[i-1].mean[k+os->frame_lenght];
				state[i-1].old_mean[k+os->frame_lenght] += x * y;
				state[i-1].old_var[k+os->frame_lenght] += x*x * y;

				x = os->acc[j][k] - state[i-1].mean[k+(os->frame_lenght*2)];
				state[i-1].old_mean[k+(os->frame_lenght*2)] += x * y;
				state[i-1].old_var[k+(os->frame_lenght*2)] += x*x * y;
		   }

	   }
   }




   delete[] occr;
}



void HMM::SetMinDuration(void)
{
	int *md;
	int *so;
	int nds,k,i,d;

	md = new int[states];
	so = new int[states];

	for (i=0,nds=0;i<states;i++) so[i]=md[i]=-1;
	FindSO(md,&nds,states-1);
	for (i=0;i<nds;i++) so[md[i]]=i;
    for (i=0;i<states;i++) md[i]=states-1;
    for (k=0,md[0]=0;k<nds;k++) {
         i=so[k];
         if (i<0 || i>=states)  continue;
         /* Find minimum duration to state i */
         for (int j=0;j<states-1;j++)
			 if (transition[j][i]>LSMALL) {
               d=md[j]+((i==states-1)?0:1);
               if (d<md[i]) md[i]=d;
            }
    }
	if (md[states-1]<0 || md[states-1]>=states) {
		fprintf(stderr,"SetMinDurs: Transition matrix with discontinuity");
		minimum_duration = (transition[0][states-1]>LSMALL ? 0 : 1 ); 
    }
    else
		minimum_duration = md[states-1];

	delete[] md;
	delete[] so;
}


void HMM::FindSO(int * so, int * d, int s)
{
   so[s]=0; 
   for (int p=0;p<states-1;p++) { 
      if (transition[p][s]>LSMALL && p!=s)
         if (so[p]<0) 
            FindSO(so,d,p);
   }
   so[s]=++(*d)-1; 
}


std::string FormatFloat(float f)
{
	int e=0;
	bool m = (f<1&&f>-1)?false:true;
	if(f!= 0.0)
		while((f<1 || f>= 10) && (f>-1 || f<=-10))
		{
			if(f<1 && f>-1)
				f*=10;
			if(f>=10 || f<=-10)
				f/=10;

			e++;
		}
	if(e == 0)
		m=true;
	
	std::stringstream stream;
	stream << std::fixed << std::setprecision(6) << f <<"e"<<(m?"+":"-")<< std::setfill('0') << std::setw(3) << e;


	return stream.str();

}

void HMM::SaveHmm(std::string out_dir)
{
	std::ofstream output(out_dir,std::ofstream::out);
	
	for (int i = 0; i < states-2; i++)
	{
		output<<"\n<STATE> "<<state[i].state_nr+1<<std::endl;
		output<<"<MEAN> "<<vector_size<<std::endl;
		for (int j = 0; j < vector_size; j++)
		{
			output<<FormatFloat(state[i].mean[j])<<" ";

		}
		output<<"\n<VARIANCE> "<<vector_size<<std::endl;
		for (int j = 0; j < vector_size; j++)
		{
			output<<FormatFloat(state[i].var[j])<<" ";

		}
		output<<"\n<GCONST> "<<FormatFloat(state[i].g_const);

	}
	output<<"\n<TRANSP> "<<states<<std::endl;
	for (int i = 0; i <states; i++)
	{
		for (int j = 0; j < states; j++)
		{
			output<<FormatFloat(exp(transition[i][j]))<<" ";
		}
		output<<"\n";
	}

	output.close();

}
