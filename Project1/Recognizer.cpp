#include "Recognizer.h"
#include <fstream>

Recognizer::Recognizer(void)
{
	null_token.like=LZERO;
	null_token.path=NULL;
	genThresh = LSMALL;
	wordpen = 0;
	genBeam = 300;
}


Recognizer::~Recognizer(void)
{
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
			TokenPropagation(hmm_nodes[i],&wordMaxTok,t);
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
   for (j=1,res=sBuf+1;j<N;j++,res++) {  /* Emitting states first */
	  i=node->seIndexes[j-1][0]; 
      endi=node->seIndexes[j-1][1];
      cur=node->states+i;

      res->tok=cur->tok; 
	  //printf("res=cur =  %f\n",res->tok.like);

      res->tok.like+=trP[i][j];

      for (i++,cur++;i<=endi;i++,cur++) {
		 //printf("cur++\n");
         cmp.tok=cur->tok;
		 //printf("trP[i][j] =  %f\n",exp(trP[i][j]));
         cmp.tok.like+=trP[i][j];
         
         if (cmp.tok.like > res->tok.like)
               res->tok=cmp.tok;

      }
      if (res->tok.like>genThresh) { /* State pruning */
		  outp=node->hmm->OutP(os,t,j);
		 //printf("res.tok.like =  %f\n",res->tok.like);
         res->tok.like+=outp;
		 //printf("state:%d\n",j);
		 //printf("res.tok.like + outp =  %f\n",res->tok.like);
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

   /* Set up pruning limits */
   if (max.like>genMaxTok.like) {
      genMaxTok=max;
      genMaxNode=node;
   }
   node->max=max.like;
   //printf("max.like=%f\n",max.like);
   i=node->seIndexes[N-1][0]; /* Exit state (ignoring tee trP) */
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
      tok.like=res->tok.like;
      if (tok.like > wordMaxTok.like) {
         wordMaxTok=tok;
         wordMaxNode=node;
      }
      
    
      node->exit->tok=null_token;
   }


}


void Recognizer::TokenPropagation(Node * node, Token * tok, int frame)
{
	Path * tmp;
	tmp = new Path();
	tmp->frame = frame;
	tmp->like = tok->like + wordpen;
	tmp->prev = tok->path;
	node->states[0].tok = *tok;
	node->states[0].tok.like += wordpen;
	node->states[0].tok.path = tmp;
}
