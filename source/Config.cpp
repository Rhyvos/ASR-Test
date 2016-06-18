#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif
#include "Config.h"
#include <fstream>
#include <stdlib.h> 
#include <algorithm>

std::vector<std::string> split1(char * str,const char * delimiter) {
  std::vector<std::string> internal;
  char * buff;
  buff = strtok(str,delimiter);
  while(buff != NULL) {
    internal.push_back(buff);
	buff = strtok(NULL,delimiter);
  }
  
  return internal;
}

Config::Config(std::string src)
{
	std::ifstream input(src,std::ifstream::in);
	std::vector<std::string> tmp;
	char buffer[256];
	while (input.good()) {
		input.getline(buffer,256);
		tmp = split1(buffer," \t");
		if(tmp.size() == 2)
		{
			params pa(tmp[0],atof(tmp[1].c_str()));
			p.push_back(pa);
		}
  }
	input.close();
}


Config::~Config(void)
{
}


float Config::GetConfig(std::string param)
{
	for (int i = 0; i < p.size(); i++)
	{
		if(p[i].param_name.compare(param) == 0)
			return p[i].value;
	}
	return 0;
}


bool Config::Exist(std::string param)
{
	for (int i = 0; i < p.size(); i++)
	{
		if(p[i].param_name.compare(param) == 0)
			return true;
	}
	return false;
}
