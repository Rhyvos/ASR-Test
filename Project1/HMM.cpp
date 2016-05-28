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



HMM::HMM(void)
{

}


HMM::~HMM(void)
{
	for(int i=0;i<states-2;i++)
	{
		delete[] state[i].mean;
		delete[] state[i].var;
	}
	delete[] state;
	for(int i=0;i<states;i++)
	{
		delete[] transition[i];
	}
	delete[] transition;
}


HMM::HMM(char * hmm_src) //load from file
{
}


HMM::HMM(int vector_size, int states): vector_size(vector_size), states(states)
{

	if(states<3 || vector_size<1)
	{
		fprintf(stderr, "HMM: Wrong initialisaton parametrs: vector_size: %d or states: %d\n",vector_size,states);
		return;
	}
	
	state = new State[states-2];				// we don't need first and last state mean and variance
	minVar  = 1.0E-2;
	for(int i=0 ;i<states-2;i++)
	{
		state[i].mean = new float[vector_size];
		state[i].var = new float[vector_size];
		for(int j=0;j<vector_size;j++)
		{
			state[i].mean[j]=0.0;
			state[i].var[j]=1.0;
		}

		state[i].state_nr=i+1;
	}

	transition = new float*[states];
	for(int i=0 ;i<states;i++)
	{
		transition[i] = new float[states];
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


void HMM::Initialise(ParamAudio * pa, int iteration, std::string label_name)
{
	float totalP,newP,delta;
	bool converged = false;   
	int * states_vec;
	int * mixes;

	GetMean(pa,label_name);
	GetVariance(pa,label_name);
	FindGConst();

	for (int i = 0; !converged && i < iteration; i++)
	{
		for (int j = 0, newP = 0; j < pa->segments; j++)
		{
			if(!pa->os[j].l->name.compare(label_name))
			{
				states_vec = new int[pa->os[j].frames];
				mixes = new int[pa->os[j].frames];

				newP += ViterbiAlg(pa->os+j,states_vec,mixes);
				for(int i=0;i<pa->os[j].frames;i++)
					printf("%d\t",states_vec[i]);
				delete[] mixes;
				delete[] states_vec;
			}
		}

	}

}


void HMM::GetMean(ParamAudio * pa, std::string label_name)
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
		if(pa->os[i].l == NULL || !pa->os[i].l->name.compare(label_name))
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




void HMM::GetVariance(ParamAudio * pa, std::string label_name)
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
		if(pa->os[i].l == NULL || !pa->os[i].l->name.compare(label_name))
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

	for(int i=0;i<os->frame_lenght;i++)
	   printf("%f\t%f\t%f\n",os->coef[fr_number][i],state[state_nr-1].mean[i],state[state_nr-1].var[i]);
	for(int i=0;i<os->frame_lenght;i++)
	   printf("%f\t%f\t%f\n",os->delta[fr_number][i],state[state_nr-1].mean[i+os->frame_lenght],state[state_nr-1].var[i+os->frame_lenght]);
	for(int i=0;i<os->frame_lenght;i++)
	   printf("%f\t%f\t%f\n",os->acc[fr_number][i],state[state_nr-1].mean[i+(os->frame_lenght*2)],state[state_nr-1].var[i+(os->frame_lenght*2)]);

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

   
   for (int i = 0; i < os->frames; i++)
   {
	   for (int j = 0; j < states; j++)
	   {
		   printf("%d\t",trace_back[i][j]);

	   }
	   printf("\n");

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


