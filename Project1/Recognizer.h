#pragma once
#include "HMM.h"
class Recognizer
{
public:
	Recognizer(void);
	~Recognizer(void);
	std::vector<Node *> hmm_nodes;
	Node *genMaxNode;     /* Most likely node in network */
	Node *wordMaxNode;    /* Most likely word end node in network */

	Token genMaxTok;      /* Most likely token */
	Token wordMaxTok;     /* Most likely word end token */
	Token null_token;
	float genThresh;
	float wordpen;
	float genBeam;
	Token *tBuf;             /* Buffer Array[2..N-1] of tok for StepHMM1 */
	TokenSet *sBuf;
	void AddHmm(HMM * hmm);
	void LoadHmm(std::string hmm_src);
	void DoRecognition(ParamAudio * pa);
	void ProccesObservation(ObservationSegment * os, int t);
	void StepNode(ObservationSegment * os, int t, Node * node);
	void TokenPropagation(Node * node, Token * tok,Node * prevn, int frame);
	void ReadPath(Path *  path);
};

