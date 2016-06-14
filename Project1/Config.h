#pragma once
#include <string>
#include <vector>
#include "Config.h"
struct params
{
	std::string param_name;
	float value;
	params(std::string s, float f) : param_name(s), value(f) {}
	params(const params & copy) :  param_name(copy.param_name), value(copy.value) {}
	params(params && move) :  param_name(move.param_name), value(move.value) {}
};

class Config
{
public:
	Config(std::string);
	~Config(void);

	std::vector<params> p;

	float GetConfig(std::string param);
	bool Exist(std::string param);
};

