
#ifndef LOGISTIC_REGRESSION_H

#define LOGISTIC_REGRESSION_H

typedef struct logisticRegressor{
    double *weight;
    double bias;
    int nWeights;
} logisticRegressor;

typedef struct trainingOutput{
    double loss; 
    double accuracy;
    int numWrong;
    int total;
} trainingOutput;

void initRegressor(logisticRegressor *regressor), freeRegressor (logisticRegressor *regressor), 
    backwardCE(logisticRegressor *regressor, double *input, double output, int label, double lr),
    initWeights(logisticRegressor *regressor, int nWeights, double std), print(logisticRegressor *regressor), normalize(logisticRegressor *regressor);

double activation(logisticRegressor *regressor, double *input), forward(logisticRegressor *regressor, double *input),
    score(logisticRegressor *regressor, double **inputs, int *label, int nPoints), CELoss(double output, int label);

trainingOutput fit(logisticRegressor *regressor, double **inputs, const bool *label, const double lr, const int epochs, const int batchSize, const int nPoints);

#endif 