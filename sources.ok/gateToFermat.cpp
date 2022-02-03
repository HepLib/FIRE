/**
 * @file gateToFermat.cpp
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

/*
 *  This file is a part of the program gateToFermat.
 *  Copyright (C) Mikhail Tentyukov <tentukov@physik.uni-bielefeld.de>
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License version 2 as
 *  published by the Free Software Foundation.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
*/

#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <cstring>
#include <sys/wait.h>
#include <sys/stat.h>
#include <cerrno>
#include <csignal>

#include "gateToFermat.h"

/**
 * Minimum memory allocation step for buffers.
 */
constexpr size_t DELTA_OUT = {128};
/**
 * Default read buffer size.
 */
constexpr int FROMFERMATBUFSIZE = {1024};

static const char *g_thepath = nullptr;
static int g_fermatExecNum = 0;

/*The output buffer (baseout) and related control p_refers.
stopout p_refs to the end of allocated space, fullout p_refs to
the end of space filled by an actual content:*/
static char *g_baseout[MAX_THREADS];
static char *g_fullout[MAX_THREADS];
static char *g_stopout[MAX_THREADS];
static FILE *g_to[MAX_THREADS];
static FILE *g_from[MAX_THREADS];
static pid_t g_childpid[MAX_THREADS];
static char *g_fbuf[MAX_THREADS];


/*The inline function places one char to the output buffer with possible
expansion of the buffer:*/
static void addtoout(char ch, int thread_number) {
    if (g_fullout[thread_number] >= g_stopout[thread_number]) {
        size_t l = g_stopout[thread_number] - g_baseout[thread_number] + DELTA_OUT;
        char *ptr = static_cast<char *>(realloc(g_baseout[thread_number], l));
        if (ptr == nullptr) {
            exit(2);
        }
        g_fullout[thread_number] = ptr + (g_fullout[thread_number] - g_baseout[thread_number]);
        g_stopout[thread_number] = ptr + l;
        g_baseout[thread_number] = ptr;
    }
    *g_fullout[thread_number]++ = ch;
}

/*Starts Fermat and swallows its stdin/stdout:*/
static pid_t openprogram(FILE **to, FILE **from) {
    int fdin[2], fdout[2];
    pid_t childpid;
    if (pipe(fdin) < 0) {
        perror("pipe");
        return (0);
    }
    if (pipe(fdout) < 0) {
        perror("pipe");
        return (0);
    }

    if ((childpid = fork()) == -1) {
        perror("fork");
        return (0);
    }
    if (childpid == 0) {
        /* Child process closes up input side of pipe */
        close(fdin[1]);
        close(fdout[0]);
        /*redirect stdin*/
        if (dup2(fdin[0], STDIN_FILENO) < 0) {
            perror("dup2 stdin");
            return 0;
        }
        /*redirect stdout*/
        if (dup2(fdout[1], STDOUT_FILENO) < 0) {
            perror("dup2 stdout");
            return 0;
        }
        /*stderr > /dev/null*/
        dup2(open("/dev/null", O_WRONLY), STDERR_FILENO);
        /*argc[0] = must be the same as g_thepath:
        Fermat uses it to find the full path:*/
        execlp(g_thepath, g_thepath, NULL);
        exit(0);
    } else { /* Parent process closes up output side of pipe */
        close(fdin[0]);
        close(fdout[1]);
        if ((*to = fdopen(fdin[1], "w")) == nullptr) {
            kill(childpid, SIGKILL);
            return 0;
        }
        if ((*from = fdopen(fdout[0], "r")) == nullptr) {
            fclose(*to);
            kill(childpid, SIGKILL);
            return 0;
        }
    }
    return childpid;
}

/*reads the stream 'from' up to the line 'terminator' (only 'thesize' first characters are compared):*/
static void readup(const char *terminator, size_t thesize, int thread_number) {
    char *c;
    for (;;) {
        do {
            for (c = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, g_from[thread_number]); *c <= ' '; c++)
                if (*c == '\0') break;
        } while (*c <= ' ');
        if (strncmp(terminator, c, thesize) == 0) {
            return;
        }
    }
}

/*Starts Fermat and makes some initializations:*/
static pid_t initFermat(FILE **to, FILE **from, const char *pvars, int thread_number) {
    pid_t childpid;
    char *ch;
    childpid = openprogram(to, from);
    if (childpid == 0) {
        return -1;
    }
    /*Fermat is running*/

    /*Switch off floating point representation:*/
    fputs("&d\n0\n", *to);

    /*Set prompt as '\n':*/
    fputs("&(M=' ')\n", *to);
    fflush(*to);

    readup("> Prompt", 8, thread_number);
    readup("Elapsed", 7, thread_number);
    char *temp = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, *from);
    if (!temp) {
        printf("Fermat init error\n");
        abort();
    }

    /*Switch off timing:*/
    fputs("&(t=0)\n", *to);
    fflush(*to);
    readup("0", 1, thread_number);
    temp = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, *from);
    if (!temp) {
        printf("Fermat init error\n");
        abort();
    }

    /*Switch on "ugly" printing: no spaces in int integers;
    do not omit '*' for multiplication:*/
    fputs("&U\n", *to);
    fflush(*to);
    readup("0", 1, thread_number);
    temp = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, *from);
    if (!temp) {
        printf("Fermat init error\n");
        abort();
    }

    /*NEW*/
    /*Switch off suppression of output to terminal of int polys.:*/
    fputs("&(_s=0)\n", *to);
    fflush(*to);
    readup("0", 1, thread_number);
    temp = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, *from);
    if (!temp) {
        printf("Fermat init error\n");
        abort();
    }

    /*Set polymomial variables:*/
    /*pvars looks like "a\nb\nc\n\n":*/
    while (*pvars > '\n') {
        /*Copy the variable up to '\n' into g_fbuf[thread_number]:*/
        /* Modification by A.Smirnov: for compatibility with Mathematica*/
        for (*(ch = g_fbuf[thread_number]) = '\0'; ((*ch = *pvars) > '\n' && (*ch = *pvars) != ','); ch++, pvars++) {
            if (ch - g_fbuf[thread_number] >= FROMFERMATBUFSIZE) {
                return -2;
            }
        }
        *ch = '\0';
        if (*pvars != '\0')pvars++;
        /*Set g_fbuf[thread_number] as a polymomial variable:*/
        fprintf(*to, "&(J=%s)\n", g_fbuf[thread_number]);
        fflush(*to);
        readup("0", 1, thread_number);
        temp = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, *from);
        if (!temp) {
            printf("Fermat init error\n");
            abort();
        }
    }

    /*Switch off fastest shortcut polygcd method (it is buggy in fermat-3.6.4!):*/
    fputs("&(_t=0)\n", *to);
    fflush(*to);
    readup("0", 1, thread_number);

    /*Set the number of bits for each symbol power. Fermat provides 128 bits in total.
          We choose 16 bits per symbol here, for a maximum of 128/16 = 8 symbols.
          If the problem exceeds this, Fermat reverts to the slow method.
          The default is 9 bits, for a maximum power of 511. This is *not* sufficient!
          16 b`its allows a maximum symbol power of 65535, which "should be" sufficient.*/
    fprintf(*to, "&(_o=16)\n");
    fflush(*to);
    readup("0", 1, thread_number);
    temp = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, *from);
    if (!temp) {
        printf("Fermat init error\n");
        abort();
    }
    return childpid;
}

/*Public function:*/
/* Make some initializations and start Fermat. pvars is a
list of polynomial variables separated by '\n', like: "a\nb\nc"*/
int iniCalc(const char *path, const char *pvars, int threads) {
    char *ch;
    if (((ch = getenv("GTF_THEPATH")) != nullptr) && (*ch != '\0')) {
        g_thepath = ch;
    }
    if (*path != '\0') {
        g_thepath = path;
    }
    if (g_thepath == nullptr) {
        perror("g_thepath is not set in iniCalc!");
        abort();
    }
    g_fermatExecNum = threads;
    for (int thread_number = 0; thread_number != g_fermatExecNum; ++thread_number) {
        g_baseout[thread_number] = nullptr;
        g_fullout[thread_number] = nullptr;
        g_stopout[thread_number] = nullptr;
        g_to[thread_number] = nullptr;
        g_from[thread_number] = nullptr;
        g_childpid[thread_number] = 0;
        g_fbuf[thread_number] = static_cast<char *>(malloc(FROMFERMATBUFSIZE));

        if ((g_childpid[thread_number] = initFermat(&g_to[thread_number], &g_from[thread_number], pvars, thread_number)) < 0) {/*Error?*/
            return int(g_childpid[thread_number]);
        }
        /*Allocate output buffer:*/
        if ((g_baseout[thread_number] = static_cast<char *>(malloc(DELTA_OUT))) != nullptr) {
            g_fullout[thread_number] = g_baseout[thread_number];
            g_stopout[thread_number] = g_baseout[thread_number] + DELTA_OUT;
        } else {/*Memory allocation error, kill started program and exit:*/
            kill(g_childpid[thread_number], SIGKILL);
            waitpid(g_childpid[thread_number], nullptr, 0);
            return -1;
        }
    }
    return 0;
}

char *calc(const string &buf, int thread_number) {
    const char *c;
    /*Feed the buffer to Fermat:*/
    for (c = buf.c_str(); *c != '\0'; c++) {
        if (putc(*c, g_to[thread_number]) != *c) {
            cout << endl << endl << "Fermat crash on fputc!" << endl << buf <<endl;
            abort();
        }
    }
    if (putc('\n', g_to[thread_number]) != '\n') {
        cout << endl << endl << "Fermat crash on fputc2!" << endl << buf <<endl;
        abort();
    };/*stroke the line*/
    fflush(g_to[thread_number]);/*Now Fermat starts to work*/

    /*Read the Fermat answer (up to "\n\n") and collect everything into the output buffer:*/
    /*Ignore leading empty lines:*/
    do {
        c = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, g_from[thread_number]);
        if (c == nullptr) {
            cout << endl << endl << "Fermat crash on fgets!" << endl << buf <<endl;
            abort();
        };
        while (*c <= ' ') {
            if (*c == '\0') break;
            c++;
        }
    } while (*c <= ' ');

    /*initialize the output buffer:*/
    g_fullout[thread_number] = g_baseout[thread_number];
    do {
        if (*c == '\0') break;

        bool previous_line_had_no_new_line_symbol = false;

        while (*c != '\n') {
            if (*c == '\0') {
                // line ended with \0 at FROMFERMATBUFSIZE
                previous_line_had_no_new_line_symbol = true;
                break;
            }
            /*ignore '`' and spaces:*/
            if ((*c != ' ') && (*c != '`')) { addtoout(*c, thread_number); }
            c++;
        }

        c = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, g_from[thread_number]);
        if (c == nullptr) {
            cout << endl << endl << "Fermat crash on fgets2!" << endl << buf <<endl;
            abort();
        };

        if ((*c == '\n') && previous_line_had_no_new_line_symbol) {
            // that is the first \n, not second
            c = fgets(g_fbuf[thread_number], FROMFERMATBUFSIZE, g_from[thread_number]);
            if (c == nullptr) {
                cout << endl << endl << "Fermat crash on fgets3!" << endl << buf <<endl;
                abort();
            };
        }
    } while (*c != '\n');
    /*the empty line is the Fermat prompt*/

    addtoout('\0', thread_number);/*Complete the line*/
    return (g_baseout[thread_number]);
}

/*Public function:*/
void closeCalc() {
    for (int thread_number = 0; thread_number != g_fermatExecNum; ++thread_number) {
        fputs("&q\n", g_to[thread_number]);
        fflush(g_to[thread_number]);
        fclose(g_to[thread_number]);
        fclose(g_from[thread_number]);
        kill(g_childpid[thread_number], SIGKILL);
        waitpid(g_childpid[thread_number], nullptr, 0);
        free(g_baseout[thread_number]);
        g_baseout[thread_number] = nullptr;
        g_fullout[thread_number] = nullptr;
        g_stopout[thread_number] = nullptr;
        free(g_fbuf[thread_number]);
    }
}
