ECHO ASR Initialise
ASR -c conf.ini -h hmm_list.txt -s train.txt -i
PAUSE
ECHO ASR Training
ASR -c conf.ini -h hmm_list.txt -s train.txt -t
PAUSE
ECHO ASR Final ReEstimation Iteration = 1
ASR -c conf.ini -h hmm_list.txt -s train.txt -f
PAUSE
ECHO ASR Final ReEstimation Iteration = 1
ASR -c conf.ini -h hmm_list.txt -s train.txt -f
PAUSE
ECHO ASR Final ReEstimation Iteration = 1
ASR -c conf.ini -h hmm_list.txt -s train.txt -f
PAUSE