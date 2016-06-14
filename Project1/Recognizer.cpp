#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif

#include "Recognizer.h"
#include <fstream>

Recognizer::Recognizer(void)
{
	null_token.like=LZERO;
	null_token.path=NULL;
	genThresh = LSMALL;
	wordpen = 0.0;
	genBeam = 300;
}


Recognizer::~Recognizer(void)
{
	for (int i = 0; i < hmm_nodes.size(); i++)
	{
		for (int j = 0; j < hmm_nodes[i]->hmm->states-1; j++)
		{
			delete[] hmm_nodes[i]->seIndexes[j];
		}
		delete[] hmm_nodes[i]->seIndexes;
		delete[] hmm_nodes[i]->states;
		delete hmm_nodes[i];
	}
}


void Recognizer::AddHmm(HMM * hmm)
{
	Node * tmp;
	int i,min,max,N;
	float ** trP;
	tmp= new Node();
	tmp->hmm = hmm;
	tmp->states = new TokenSet[hmm->states-1];
	for (int i = 0; i < hmm->states-1; i++)
	{
		tmp->states[i].tok.like = LZERO;
		tmp->states[i].tok.path = NULL;

	}
	
	tmp->seIndexes = new short*[hmm->states-1];
	N=hmm->states-1;
	trP=hmm->transition;
	for (int i = 1; i <= hmm->states-1; i++)
	{
		tmp->seIndexes[i-1] = new short[2];
		for (min=(i==N-1)?1:0;min<N;min++) /* Don't want tee transitions */
            if (trP[min][i]>LSMALL) break;
         for (max=N-1;max>1;max--)
            if (trP[max][i]>LSMALL) break;
		 tmp->seIndexes[i-1][0]=min;
         tmp->seIndexes[i-1][1]=max;
	}


	hmm_nodes.push_back(tmp);
}


void Recognizer::LoadHmm(std::string hmm_src)
{
	HMM * h = new HMM(hmm_src);
	AddHmm(h);
}



void Recognizer::DoRecognition(ParamAudio * pa)
{

	//resettokens()
	int max=0;
	std::string tmp="";
	for (int i = 0; i < hmm_nodes.size(); i++)
	{
		hmm_nodes[i]->states[0].tok.like = 0;
		hmm_nodes[i]->exit = new TokenSet();
		if(hmm_nodes[i]->hmm->states > max)
			max = hmm_nodes[i]->hmm->states;
	}
	max--;
	
	sBuf=new TokenSet[max];
	for (int i = 0; i < max; i++)
	{
		sBuf[i].tok=null_token;
	}

	for (int i = 0; i < pa->os[0].frames; i++)
	{

		ProccesObservation(pa->os,i);
		if(tmp != std::string(genMaxNode->hmm->name))
		{
			tmp = std::string(genMaxNode->hmm->name);
			printf("Observation @%d most like node %s\n",i,tmp.c_str());
			
		}
		
		//tokenpropagation(bestwordtoken)
	}

	if(wordMaxTok.path != NULL)
	{
		ReadPath(wordMaxTok.path);
		float start, end;
		start = (wordMaxTok.path->prev==NULL)?0:wordMaxTok.path->prev->frame;
		end = pa->os[0].frames;
		start++;
		end++;
		start = (start * 160.0 * 226.0) / 10000000.0;
		end = (end * 160.0 * 226.0) / 10000000.0;
		printf("%f %f %s \n",start,end,std::string(wordMaxNode->hmm->name).c_str());
	}
	
	FreeMemory();

	for (int i = 0; i < hmm_nodes.size(); i++)
	{
		delete hmm_nodes[i]->exit;

	}

	delete[] sBuf;

}


void Recognizer::ProccesObservation(ObservationSegment * os, int t)
{

	genMaxTok=wordMaxTok=null_token;
    genMaxNode=wordMaxNode=NULL;
	for (int i = 0; i < hmm_nodes.size(); i++)
	{
		StepNode(os,t,hmm_nodes[i]);
	}
	genThresh=genMaxTok.like-genBeam;
	if (genThresh<LSMALL) genThresh=LSMALL;
	

	if(wordMaxTok.like>genThresh)
	{
		for (int i = 0; i < hmm_nodes.size(); i++)
		{
			TokenPropagation(hmm_nodes[i],&wordMaxTok,wordMaxNode,t);
		}
	}


}


void Recognizer::StepNode(ObservationSegment * os, int t, Node * node)
{

   HMM *hmm;
   Token tok,max;
   TokenSet *res,cmp,*cur;
   int i,j,k,N,endi;
   float outp;
   float ** trP;

   max=null_token;
   hmm=node->hmm; 
   N=hmm->states-1;
   trP=hmm->transition;
   for (j=1,res=sBuf+1;j<N;j++,res++) 
   { 
	  i=node->seIndexes[j-1][0]; 
      endi=node->seIndexes[j-1][1];
      cur=node->states+i;
      res->tok=cur->tok; 
      res->tok.like+=trP[i][j];
      for (i++,cur++;i<=endi;i++,cur++)
	  {
         cmp.tok=cur->tok;
         cmp.tok.like+=trP[i][j];
         if (cmp.tok.like > res->tok.like)
               res->tok=cmp.tok;
      }
      if (res->tok.like>genThresh) { 
		  outp=node->hmm->OutP(os,t,j);
         res->tok.like+=outp;
         if (res->tok.like>max.like)
            max=res->tok;
      } 
      else {
         res->tok=null_token;
      }
   }

   


   for (i=0,res=sBuf,cur=node->states; i<N;i++,res++,cur++) {
      cur->tok=res->tok; 
   }
   if (max.like>genMaxTok.like) {
      genMaxTok=max;
      genMaxNode=node;
   }
   node->max=max.like;

   i=node->seIndexes[N-1][0];
   endi=node->seIndexes[N-1][1];
   
   res=node->exit;
   cur=node->states+i;

 
   res->tok=cur->tok; 
   

   res->tok.like+=trP[i][N];

   for (i++,cur++;i<=endi;i++,cur++) {
      cmp.tok=cur->tok; 
      cmp.tok.like+=trP[i][N];

         if (cmp.tok.like > res->tok.like) 
            res->tok=cmp.tok;
     
   }
   if (res->tok.like>LSMALL){
      tok=res->tok;
      if (tok.like > wordMaxTok.like) {
         wordMaxTok=tok;
         wordMaxNode=node;
      }
      
    
      node->exit->tok=null_token;
   }


}


void Recognizer::TokenPropagation(Node * node, Token * tok,Node * prevn, int frame)
{
	Path * tmp;
	tmp = new Path();
	paths.push_back(tmp);
	tmp->frame = frame;
	tmp->like = tok->like + wordpen;
	tmp->prev = tok->path;
	tmp->node = prevn;
	node->states[0].tok = *tok;
	node->states[0].tok.like += wordpen;
	node->states[0].tok.path = tmp;
}


void Recognizer::ReadPath(Path *  path)
{
	if(path->prev != NULL)
		ReadPath(path->prev);
	
	float start, end;
	start = (path->prev==NULL)?0:path->prev->frame;
	end = path->frame;
	start++;
	end++;
	start = (start * 160.0 * 226.0) / 10000000.0;
	end = (end * 160.0 * 226.0) / 10000000.0;
	printf("%f %f %s \n",start,end,std::string(path->node->hmm->name).c_str());
}


void Recognizer::FreeMemory(void)
{
	for (int i = 0; i < paths.size(); i++)
	{
		delete paths[i];
	}
}
