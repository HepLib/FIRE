#pragma once

#include "common.h"
#include <cstdio>
#include <iostream>
#include <string>

int iniCalc(const char *path, const string & vars);
string calc(const string &buf);
void closeCalc();
