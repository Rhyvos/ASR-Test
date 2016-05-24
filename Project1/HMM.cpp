#define _USE_MATH_DEFINES
#include "HMM.h"
#include <cmath>
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
			transition[i][j] = 0.0;
		if(i!=0 && i!=(states-1))
		{
			transition[i][i] = 0.6;
			transition[i][i+1] = 0.4;
		}
		if(i==0)
			transition[i][i+1] = 1.0; 

	}




}


void HMM::Initialise(ParamAudio * pa, int iteration, std::string label_name)
{
	GetMean(pa,label_name);
	GetVariance(pa,label_name);
	FindGConst();
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
						temp_var[statenumber][k] = x*x;
					}

						
					if(pa->os[i].delta != NULL)
					{
						x=(pa->os[i].delta[j][k])-state[statenumber].mean[k+pa->os[i].frame_lenght] ;
						temp_var[statenumber][k+pa->os[i].frame_lenght] = x*x;
					}

						
					if(pa->os[i].acc != NULL)
					{
						x=(pa->os[i].acc[j][k])-state[statenumber].mean[k+(2*pa->os[i].frame_lenght)];
						temp_var[statenumber][k+(2*pa->os[i].frame_lenght)] = x*x;
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
			state[i].var[j] = temp_var[i][j]/state_fr_counter[i];
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
