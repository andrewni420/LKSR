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

ttpSolution bestTourPacking(ttpInstance *instance, int** tours, int nTours, bool isCityPacking, float* capacities, int nCapacities);

int initializedHillClimber(ttpInstance *instance, ttpSolution *sol, double strength, int durationWithoutImprovement, int maxGens),
    classifierHillClimber(ttpInstance *instance, ttpSolution *sol, logisticRegressor *classifier, double strength, double factor, double weight, int durationWithoutImprovement, int maxGens),
    linearHillClimber(ttpInstance *instance, ttpSolution *sol, int nonImprovementDuration, int generations, double *initLinear, double *initPercent),
    tourOrReverse(ttpInstance *instance, ttpSolution *sol, float *capacities, int nCapacities),
    initializedLinearClimber(ttpInstance *instance, ttpSolution *sol, int nonImprovementDuration, int generations, double *linear, double *percent),
    gradientClimber(ttpInstance *instance, ttpSolution *sol, int nSteps, double lr), 
    quadraticHillClimber(ttpInstance *instance, ttpSolution *sol, logisticRegressor *quadratic_heuristic, int nonImprovementDuration, int generations, double *capacity),
    initializedQuadraticClimber(ttpInstance *instance, ttpSolution *sol, logisticRegressor *quadratic_heuristic, int nonImprovementDuration, int generations, double *capacity);