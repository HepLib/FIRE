/**
 * @file gateToFermat.h
 * @author Mikhail Tentyukov <tentukov@physik.uni-bielefeld.de>
 *
 * This file is a part of the program gateToFermat.
 * Copyright (C) Mikhail Tentyukov <tentukov@physik.uni-bielefeld.de>
 *
 * This is the wrapper to the program FERMAT
 * http://www.bway.net/~lewis

 * There are the following environment variables recognized by the program:
 *  * GTF_THEPATH - full path to the Fermat executable.
 */

#ifndef GATETOFERMAT_H
#define GATETOFERMAT_H 1

#include "common.h"
#include <cstdio>
#include <iostream>
#include <string>

/**
 * Stop Fermat and worker threads.
 */
void closeCalc();

/**
 * Send job to fermat worker thread.
 * @param buf buffer to be evaluated.
 * @param thread_number number of thread that'll receive the job.
 * @return pointer to the result of calculations.
 */
char *calc(const string &buf, int thread_number);

/**
 * Make some initializations and start Fermat.
 * @param path path to Fermat executable. Is overridden by GTF_THEPATH environment variable.
 * @param pvars list of polynomial variables separated by '\\n'.
 * @param threads number of threads.
 * @return 0 if everything is alright, error status otherwise.
 */
int iniCalc(const char *path, const char *pvars, int threads);

#endif
