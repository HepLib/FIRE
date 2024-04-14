#pragma once

#include "common.h"
#include "point.h"
#include "equation.h"
#include <sys/stat.h>
#include <dirent.h>

int parse_config(const string& filename, set<point, indirect_more> &points, string &output);
vector<COEFF> split_coeff(const string &s);
int s2i(const char *digit, int &result);
int s2u(const char *digit, unsigned int &result);
int s2v(const char *digit, vector<t_index> &result);
int s2sf(const char *digit, SECTOR &result);  // to sector fast directly
int s2vv(const char *digit, vector<vector<t_index> > &result);
int s2vvv(const char *digit, vector<vector<vector<t_index> > > &result);
void parseArgcArgv(int argc, char *argv[]);
