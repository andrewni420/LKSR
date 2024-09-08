#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "util.h"
#include "math_util.h"

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


void initRegressor(logisticRegressor *regressor){
    regressor->weight = (double *)NULL;
}

void initWeights(logisticRegressor *regressor, int nWeights, double std){
    regressor->weight = CC_SAFE_MALLOC(nWeights,double);
    regressor->nWeights=nWeights;
    for (int i=0;i<nWeights;i++){
        regressor->weight[i] = randNormal(0,std);
    }
    regressor->bias = 0;
}

void freeRegressor (logisticRegressor *regressor)
{
    CC_IFFREE (regressor->weight, double);
}


double activation(const logisticRegressor *regressor, const double *input){
    double output = 0;
    for (int i=0;i<regressor->nWeights;i++){
        output+=regressor->weight[i]*input[i];
    }
    return regressor->bias+output;
}

double forward(logisticRegressor *regressor, double *input){

    // regressor->weight[0] = 3.338779;
    // regressor->weight[1] = -0.418764;
    // regressor->bias = -1.759471;
    // input[0] = 1620.401198;
    // input[1] = 1.140187;
    // double act = activation(regressor, input);
    // printf("activation %lf\n",act);
    // double out = sigmoid(act);
    // printf("output %lf\n",out);
    // exit(0);

    // output 0.000000 loss -1.000000 lr 0.111111 input [1620.401198, 1.140187]pre bias -1.759471 weight [3.338779, -0.418764]
    // bias -1.648360 weight [183.383357, -0.292077]
    return sigmoid(activation(regressor, input));
}

double score(const logisticRegressor *regressor, const double **inputs, const int *label, const int nPoints){
    double accuracy = 0;
    for (int i=0;i<nPoints;i++){
        double output = activation(regressor, inputs[i]);
        accuracy+=output>0 ? label[i] : 1-label[i];
    }
    return accuracy/nPoints;
}

void backwardCE(logisticRegressor *regressor, const double *input, const double output, const int label, const double lr){
    double loss = output-label;
    // printf("output %lf loss %lf lr %lf input ", output, loss, lr);
    // printDoubleArray(input, regressor->nWeights);
    // printf("pre bias %lf weight ",regressor->bias);
    // printDoubleArray(regressor->weight, regressor->nWeights);
    // printf("\n");
    
    for (int i=0;i<regressor->nWeights;i++) regressor->weight[i] -= lr*input[i]*loss;
    regressor->bias-=lr*loss;
    // printf("bias %lf weight ", regressor->bias);
    // printDoubleArray(regressor->weight, regressor->nWeights);
    // printf("\n\n");
}

double CELoss(double output, int label){
    return -(label==1? log(output+1E-9) : log(1 - output + 1E-9))/2;
}


trainingOutput fit(logisticRegressor *regressor, double **inputs, const bool *label, double lr, const int epochs, const int batchSize, const int nPoints){
    double loss = 0;
    double accuracy = 0;
    int *indices = CC_SAFE_MALLOC(nPoints, int);
    double *output = CC_SAFE_MALLOC(batchSize, double);
    range(indices,nPoints);
    int pointsPerEpoch = batchSize*(nPoints/batchSize);

    for (int epoch=0;epoch<epochs;epoch++){
        loss = 0;
        accuracy = 0;
        
        shuffle(indices,nPoints);
        for (int i=0;i<nPoints/batchSize;i++){
            for (int j=0;j<batchSize;j++) {
                // printf("weights ");
                // printDoubleArray(regressor->weight, regressor->nWeights);
                // printf(" bias %lf\n",regressor->bias);
                output[j] = forward(regressor, inputs[indices[i*batchSize + j]]);
                // printf("output %lf, correct %d\n",output[j], label[indices[i*batchSize + j]]);
                loss +=CELoss(output[j],label[indices[i*batchSize + j]]);
                int dAcc = output[j]>0.5 ? label[indices[i*batchSize + j]] : 1-label[indices[i*batchSize + j]];
                accuracy+=dAcc;
                // printf("j %d input [%lf, %lf] output %lf activation %lf label %d\n", j, inputs[indices[i*batchSize + j]][0], inputs[indices[i*batchSize + j]][1], output[j], activation(regressor, inputs[indices[i*batchSize + j]]), label[indices[i*batchSize + j]]);
            }
            for (int j=0;j<batchSize;j++) backwardCE(regressor, inputs[indices[i*batchSize + j]],output[j], label[indices[i*batchSize + j]] ? 1 : 0, lr/(double) batchSize);
        }
        loss /= pointsPerEpoch;
        accuracy /= pointsPerEpoch;
        // printf("loss %lf, accuracy %lf\n",loss,accuracy);
        // double arr[] = {regressor->weight[0],regressor->weight[1],regressor->bias};
        // printDoubleArray(&arr, 3);
        // printf("\n");
    }
    CC_IFFREE(indices, int);
    CC_IFFREE(output, double);

    loss = 0;
    accuracy = 0;

    for (int i=0;i<nPoints;i++){
        double output = forward(regressor, inputs[i]);
        loss +=CELoss(output,label[i]);
        int dAcc = output>0.5 ? label[i] : 1-label[i];
        accuracy+=dAcc;
    }

    trainingOutput result;
    result.loss = loss/nPoints;
    result.accuracy = accuracy/nPoints;
    result.total = nPoints;
    result.numWrong = (int)(nPoints-accuracy);

    return result;
}

void printRegressor(logisticRegressor *regressor){
    printf("Regressor(");

}

void normalize(logisticRegressor *regressor){
    double sumsq = 0;
    for (int i=0;i<regressor->nWeights;i++) sumsq+=regressor->weight[i]*regressor->weight[i];
    sumsq = sqrt(sumsq)+1E-9;
    for (int i=0;i<regressor->nWeights;i++) regressor->weight[i]/=sumsq;
    regressor->bias/=sumsq;
}

