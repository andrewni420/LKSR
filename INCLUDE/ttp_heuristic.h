
#include "ttp_util.h"
#include "logistic_regression.h"

void combinedHeuristic(ttpInstance *instance, ttpSolution *sol, double rDistCapacity, double profitRatioCapacity),
    greedyProfitRatio(ttpInstance *instance, ttpSolution *sol), greedyRDist(ttpInstance *instance, ttpSolution *sol),
    greedyRDistCapped(ttpInstance *instance, ttpSolution *sol, double capacity), greedyProfitRatioCapped(ttpInstance *instance, ttpSolution *sol, int capacity),
    linearHeuristicCapped(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, double capacity, bool sorted), 
    linearHeuristic(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, bool sorted),
    quadraticHeuristicCapped(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, double capacity, bool sorted),
    quadraticHeuristic(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, bool sorted);

ttpObjective quadraticHeuristicIterated(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, double *capacity);
ttpObjective dh(ttpInstance *instance, ttpSolution *sol);

double estimateLinear(ttpInstance *instance, int estimator), estimatePercent(ttpInstance *instance, int estimator);
int packIterative(ttpInstance *instance, ttpSolution *sol, double C, double D, int Q);

extern int heuristic_evaluations;