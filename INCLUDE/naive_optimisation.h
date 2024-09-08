#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include "util.h"
#include "ttp_heuristic.h"
#include "ttp_util.h"
#include "math_util.h"
#include "logistic_regression.h"

int NaiveTTPSolver(CCdatagroup *dat, CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol),
    IteratedTTPSolver(CCdatagroup *dat, CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol, int nIters);