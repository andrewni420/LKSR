#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include "util.h"
#include "math_util.h"
#include "ttp_util.h"
#include "logistic_regression.h"
#include "ttp_heuristic.h"
#include "kdtree.h"



ttpSolution bestTourPacking(ttpInstance *instance, int** tours, int nTours, bool isCityPacking, float* capacities, int nCapacities){
    if (instance==NULL) {printf("best Tour Packing instance null\n"); fflush(stdout);}
    if (capacities==NULL) {printf("best Tour Packing null capacities\n"); fflush(stdout);}
    ttpSolution bestSol;
    double bestObj = -1*DBL_MAX;
    bestSol.packing = (bool *)NULL;
    bestSol.isCityPacking = isCityPacking;
    for (int i=0;i<nTours;i++){
        if (tours[i]==NULL) {printf("best Tour Packing null Tour\n"); fflush(stdout);}
        ttpSolution sol;
        sol.tour = tours[i];
        sol.isCityPacking = isCityPacking;
        sol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
        combinedHeuristic(instance, &sol, capacities[0], capacities[1]);
        ttpObjective obj = evaluate(instance, &sol);
        printf("bestTourPacking i %d obj %lf\n",i,obj.objective);
        if (obj.objective>bestObj){
            CC_IFFREE(bestSol.packing, bool);
            bestObj = obj.objective;
            bestSol = sol;
        } else {
            CC_IFFREE(sol.packing, bool);
        }
    }

    int *bestTour = bestSol.tour;
    if (bestTour==NULL) {printf("best Tour Packing null bestTour\n"); fflush(stdout);}
    bestSol.tour = CC_SAFE_MALLOC(instance->numberOfNodes+1, int);
    for (int i=0;i<=instance->numberOfNodes;i++) bestSol.tour[i] = bestTour[i];

    return bestSol;
}

int tourOrReverse(ttpInstance *instance, ttpSolution *sol, float *capacities, int nCapacities){
    int **tourPtrs = CC_SAFE_MALLOC(2, int *);
    tourPtrs[0] = CC_SAFE_MALLOC(instance->numberOfNodes+1, int);
    tourPtrs[1] = CC_SAFE_MALLOC(instance->numberOfNodes+1, int);
    copyInt(sol->tour, tourPtrs[0], instance->numberOfNodes);
    copyInt(sol->tour, tourPtrs[1], instance->numberOfNodes);
    tourPtrs[0][instance->numberOfNodes]=0;
    tourPtrs[1][instance->numberOfNodes]=0;

    reverseInt(tourPtrs[1], instance->numberOfNodes+1);
    ttpSolution bestTourSol = bestTourPacking(instance, tourPtrs, 2, sol->isCityPacking, capacities, nCapacities);

    CC_IFFREE(tourPtrs[0], int);
    CC_IFFREE(tourPtrs[1], int);
    CC_IFFREE(tourPtrs, int*);

    return 0;
}


int initializedHillClimber(ttpInstance *instance, ttpSolution *sol, double strength, int durationWithoutImprovement, int maxGens){
    ttpObjective obj = evaluate(instance, sol);
    double bestObjective = obj.objective;
    ttpSolution newSol;
    newSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    newSol.tour = sol->tour;
    newSol.isCityPacking = sol->isCityPacking;

    int counter = 0;
    for (int i=0;i<maxGens;i++){
        if (counter>durationWithoutImprovement) break;
        if (i%10==0 && CHECK_TIME()) break;
        for (int j=0;j<instance->numberOfItems;j++) {
            if (randUnif(0,1)<strength/(double)instance->numberOfItems){
                newSol.packing[j]=rand()%2;
            } else {
                newSol.packing[j]=sol->packing[j];
            }
        }
                
        ttpObjective newObj = evaluate(instance, &newSol);

        if (newObj.objective >= bestObjective && newObj.finalweight<=instance->capacityOfKnapsack) {

            // for the stopping criterion: check if there was an actual improvement
            if (newObj.objective > bestObjective) counter = 0;
            bestObjective = newObj.objective;
            for (int j=0;j<instance->numberOfItems;j++) sol->packing[j]=newSol.packing[j];
        } else counter ++;
        // if (i%1000==0) printf("generation %d, best %lf\n", i, bestObjective);

    }

    CC_IFFREE(newSol.packing, bool);
    return 0;
}


int classifierHillClimber(ttpInstance *instance, ttpSolution *sol, logisticRegressor *classifier, double strength, double factor, double weight, int durationWithoutImprovement, int maxGens) {
    double **inputs = CC_SAFE_MALLOC(instance->numberOfItems,double*);
    for (int i = 1; i < instance->numberOfNodes; i++) {
        for (int j = 0; j < instance->itemsPerCity; j++) {
            int itemIndex = getItemIndex(instance, sol->tour[i], j);
            int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, i, j);
            inputs[packingIndex] = CC_SAFE_MALLOC(3, double);
            inputs[packingIndex][0] = instance->normalizedProfitRatios[itemIndex];
            inputs[packingIndex][1] = instance->normalizedRDist[sol->tour[i]];
            inputs[packingIndex][2] = inputs[packingIndex][0]*inputs[packingIndex][1];
        }
    }

    int *indicesToChange = CC_SAFE_MALLOC(instance->numberOfItems, int);
    double *activations = CC_SAFE_MALLOC(instance->numberOfItems, double);

    // logisticRegressor classifier;
    // initRegressor(&classifier);
    // initWeights(&classifier, 2, 1);

    trainingOutput results = fit(classifier, inputs, sol->packing, 0.1, 10, 9, instance->numberOfItems);
    normalize(classifier);

    int nIndices = 0;
    for (int i=0;i<instance->numberOfItems;i++){
        activations[i] = activation(classifier, inputs[i]);
    }
    scaleQuartile(activations,instance->numberOfItems,factor);
    for (int i=0;i<instance->numberOfItems;i++){
        if (fabs(activations[i])<1 || (activations[i]>0 != sol->packing[i]))  
        indicesToChange[nIndices++]=i;
    }
    
    ttpObjective obj = evaluate(instance, sol);
    double bestObjective = obj.objective;
    ttpSolution newSol;
    newSol.tour = sol->tour;
    newSol.isCityPacking = sol->isCityPacking;
    newSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);

    double totalWeightInit = obj.finalweight;
    // printf("init weight %lf objective %lf capacity %lf\n", totalWeightInit, obj.objective, instance->capacityOfKnapsack);
    double nIterReduction = 1;
    double nIterOffset = 0;

    int counter = 0;
    for (int i = 0; i < maxGens; i++) {
        if (counter>durationWithoutImprovement) break;
        if (i%10==0 && CHECK_TIME()) break;
        
        copyBool(sol->packing, newSol.packing, instance->numberOfItems);
        double totalWeight = totalWeightInit;

        for (int j_ = 0; j_ < nIndices; j_++) {
            int j = indicesToChange[j_];
            if (randUnif(0,1) < strength/(double) nIndices){
                // double output = activation(classifier, inputs[j]);
                double output = activations[j];
                int changed = (randUnif(0,1) > sigmoid(output*weight)) ? 0 : 1;
                // if (fabs(output)>5 && changed!=newSol.packing[j]) 
                // printf("iteration %d large output and changed output %lf prob %lf prev %d actual %d\n",i, output, sigmoid(output*factor), newSol.packing[j], changed);
                newSol.packing[j] = changed;
            }
            if (sol->isCityPacking){
                if (sol->packing[j]) totalWeight-=instance->weight[j];
                if (newSol.packing[j]) totalWeight+=instance->weight[j];
            } else {
                int city = sol->tour[j/instance->itemsPerCity+1];
                int itemNumber = j%instance->itemsPerCity;
                int itemIndex = getItemIndex(instance, city,itemNumber);
                if (sol->packing[j]) totalWeight-=instance->weight[itemIndex];
                if (newSol.packing[j]) totalWeight+=instance->weight[itemIndex];
            }
        }

        // ttpObjective newObj = evaluate(instance, &newSol);

        // if (totalWeight!=newObj.finalweight){
        //     printf("totalWeight construction has unequal weights total %lf final %lf\n", totalWeight, newObj.finalweight);
        //     exit(0);
        // }

        int nIters = 0;
        
        // printf("total weight reduction from %lf\n", totalWeight);
        while (totalWeight>instance->capacityOfKnapsack){
            int j_ = (int)randUnif(0,nIndices);
            int j = indicesToChange[j_];
            if (newSol.packing[j] && randUnif(0,1) <  strength/(double) nIndices){
                // printf("total weight %lf capacity %lf\n", totalWeight, instance->capacityOfKnapsack);
                // double output = activation(classifier, inputs[j]);
                double output = activations[j];
                // printf("index %d output %lf prob %lf\n", j, output, sigmoid(output*factor*nIterReduction));
                output = sigmoid(output*weight*nIterReduction);
                output =nIterOffset + (1-2*nIterOffset)*output;
                int changed = (randUnif(0,1) > output) ? 0 : 1;
                if (sol->isCityPacking){
                    totalWeight+=(changed - (newSol.packing[j] ? 1 : 0)) * instance->weight[j];
                } else {
                    int city = sol->tour[j/instance->itemsPerCity+1];
                    int itemIndex = j%instance->itemsPerCity;
                    totalWeight+=(changed - (newSol.packing[j] ? 1 : 0)) * instance->weight[getItemIndex(instance, city,itemIndex)];
                }
                newSol.packing[j] = changed;
                ++nIters;
                if (nIters%nIndices==nIndices-1) {
                    // printf("reducing factor to %lf offset to %lf\n", factor*nIterReduction/10, nIterOffset);
                    // printDoubleArray(classifier->weight, 2);
                    // printf("\n");
                    nIterReduction/=10;
                    nIterOffset += (0.5-nIterOffset)/5;
                    }
            }
        }
        // if (nIters<100) nIterReduction*=2;
        // printf("total weight reduction to %lf\n", totalWeight);
        
        
        ttpObjective newObj = evaluate(instance, &newSol);

        if (totalWeight!=newObj.finalweight){
            printf("totalWeight constraint enforcement has unequal weights total %lf final %lf\n", totalWeight, newObj.finalweight);
            exit(0);
        }

        if (newObj.objective >= bestObjective && newObj.finalweight <= instance->capacityOfKnapsack) {
            for (int k=0;k<instance->numberOfItems;k++) sol->packing[k] = newSol.packing[k];

            // for the stopping criterion: check if there was an actual improvement
            if (newObj.objective > bestObjective) {
                // printf("iteration %d, best %lf prev %lf\n", i, newObj.objective, bestObjective);

                trainingOutput results = fit(classifier, inputs, sol->packing, 0.1, 10, 9, instance->numberOfItems);
                normalize(classifier);
                nIndices = 0;
                for (int k=0;k<instance->numberOfItems;k++){
                    activations[k] = activation(classifier, inputs[k]);
                }
                scaleQuartile(activations,instance->numberOfItems,factor);
                for (int i=0;i<instance->numberOfItems;i++){
                    if (fabs(activations[i])<1 || (activations[i]>0 != sol->packing[i]))  
                    indicesToChange[nIndices++]=i;
                }
                // printf("nIndices %d\n", nIndices);
                // printf("update loss %lf accuracy %lf\n",results.loss, results.accuracy);
                counter = 0;
            }

            bestObjective = newObj.objective;
            totalWeightInit = newObj.finalweight;
        } else counter++;

    }

    CC_IFFREE(newSol.packing, bool);
    // CC_IFFREE(profitRatios, double);
    // CC_IFFREE(rDist, double);
    for (int i=0;i<instance->numberOfItems;i++) CC_IFFREE(inputs[i], double);
    CC_IFFREE(inputs, double*);
    CC_IFFREE(indicesToChange, int);
    CC_IFFREE(activations, double);

    // printf("CLASSIFIER OVER\n");

    return 0;
}

// int classifierHillClimber(ttpInstance *instance, ttpSolution *sol, logisticRegressor *classifier, int numChanged, int durationWithoutImprovement, int maxGens) {
//     // printf("CLASSIFIER HILL CLIMBER\n");
//     // numChanged = 5;
//     // double *rDist = reverseDistances(instance, sol->tour);
//     // normalizeMAD(&(rDist[1]), instance->numberOfNodes-1);
//     // double *profitRatios =  CC_SAFE_MALLOC(instance->numberOfItems, double);
//     // for (int i = 0; i < instance->numberOfItems; i++)
//     //     profitRatios[i] = instance->profit[i] / instance->weight[i];
//     // normalizeMAD(profitRatios, instance->numberOfItems);


//     double **inputs = CC_SAFE_MALLOC(instance->numberOfItems,double*);
//     for (int i = 1; i < instance->numberOfNodes; i++) {
//         for (int j = 0; j < instance->itemsPerCity; j++) {
//             int itemIndex = getItemIndex(instance, sol->tour[i], j);
//             int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, i, j);
//             inputs[packingIndex] = CC_SAFE_MALLOC(3, double);
//             inputs[packingIndex][0] = instance->normalizedProfitRatios[itemIndex];
//             inputs[packingIndex][1] = instance->normalizedRDist[sol->tour[i]];
//             inputs[packingIndex][2] = inputs[packingIndex][0]*inputs[packingIndex][1];
//         }
//     }

//     int *indicesToChange = CC_SAFE_MALLOC(instance->numberOfItems, int);
//     double *activations = CC_SAFE_MALLOC(instance->numberOfItems, double);

//     // logisticRegressor classifier;
//     // initRegressor(&classifier);
//     // initWeights(&classifier, 2, 1);

//     trainingOutput results = fit(classifier, inputs, sol->packing, 0.1, 10, 9, instance->numberOfItems);
//     normalize(classifier);

//     int nIndices = 0;
//     for (int i=0;i<instance->numberOfItems;i++){
//         activations[i] = activation(classifier, inputs[i]);
//         if (fabs(activations[i])<1 || (activations[i]>0 != sol->packing[i]))  
//         indicesToChange[nIndices++]=i;
//     }
//     // printf("nIndices %d\n", nIndices);

//     // printf("initial loss %lf accuracy %lf\n",results.loss, results.accuracy);

//     ttpObjective obj = evaluate(instance, sol);
//     double bestObjective = obj.objective;
//     ttpSolution newSol;
//     newSol.tour = sol->tour;
//     newSol.isCityPacking = sol->isCityPacking;
//     newSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);

//     double totalWeightInit = obj.finalweight;
//     // printf("init weight %lf objective %lf capacity %lf\n", totalWeightInit, obj.objective, instance->capacityOfKnapsack);
//     double nIterReduction = 1;
//     double nIterOffset = 0;

//     int counter = 0;
//     for (int i = 0; i < maxGens; i++) {
//         if (counter>durationWithoutImprovement) break;
//         if (i%10==0 && CHECK_TIME()) break;
        
//         copyBool(sol->packing, newSol.packing, instance->numberOfItems);
//         double totalWeight = totalWeightInit;
//         double factor = 200 - 190*i/(double)maxGens;

//         for (int j_ = 0; j_ < nIndices; j_++) {
//             int j = indicesToChange[j_];
//             if (randUnif(0,1) < (double) numChanged/(double) nIndices){
//                 // double output = activation(classifier, inputs[j]);
//                 double output = activations[j];
//                 int changed = (randUnif(0,1) > sigmoid(output*factor)) ? 0 : 1;
//                 // if (fabs(output)>5 && changed!=newSol.packing[j]) 
//                 // printf("iteration %d large output and changed output %lf prob %lf prev %d actual %d\n",i, output, sigmoid(output*factor), newSol.packing[j], changed);
//                 newSol.packing[j] = changed;
//             }
//             if (sol->isCityPacking){
//                 if (sol->packing[j]) totalWeight-=instance->weight[j];
//                 if (newSol.packing[j]) totalWeight+=instance->weight[j];
//             } else {
//                 int city = sol->tour[j/instance->itemsPerCity+1];
//                 int itemNumber = j%instance->itemsPerCity;
//                 int itemIndex = getItemIndex(instance, city,itemNumber);
//                 if (sol->packing[j]) totalWeight-=instance->weight[itemIndex];
//                 if (newSol.packing[j]) totalWeight+=instance->weight[itemIndex];
//             }
//         }

//         ttpObjective newObj = evaluate(instance, &newSol);

//         if (totalWeight!=newObj.finalweight){
//             printf("totalWeight construction has unequal weights total %lf final %lf\n", totalWeight, newObj.finalweight);
//             exit(0);
//         }

//         int nIters = 0;
        
//         // printf("total weight reduction from %lf\n", totalWeight);
//         while (totalWeight>instance->capacityOfKnapsack){
//             int j_ = (int)randUnif(0,nIndices);
//             int j = indicesToChange[j_];
//             if (newSol.packing[j] && randUnif(0,1) < (double) numChanged/(double) nIndices){
//                 // printf("total weight %lf capacity %lf\n", totalWeight, instance->capacityOfKnapsack);
//                 // double output = activation(classifier, inputs[j]);
//                 double output = activations[j];
//                 // printf("index %d output %lf prob %lf\n", j, output, sigmoid(output*factor*nIterReduction));
//                 output = sigmoid(output*factor*nIterReduction);
//                 output =nIterOffset + (1-2*nIterOffset)*output;
//                 int changed = (randUnif(0,1) > output) ? 0 : 1;
//                 if (sol->isCityPacking){
//                     totalWeight+=(changed - (newSol.packing[j] ? 1 : 0)) * instance->weight[j];
//                 } else {
//                     int city = sol->tour[j/instance->itemsPerCity+1];
//                     int itemIndex = j%instance->itemsPerCity;
//                     totalWeight+=(changed - (newSol.packing[j] ? 1 : 0)) * instance->weight[getItemIndex(instance, city,itemIndex)];
//                 }
//                 newSol.packing[j] = changed;
//                 ++nIters;
//                 if (nIters%nIndices==nIndices-1) {
//                     // printf("reducing factor to %lf offset to %lf\n", factor*nIterReduction/10, nIterOffset);
//                     // printDoubleArray(classifier->weight, 2);
//                     // printf("\n");
//                     nIterReduction/=10;
//                     nIterOffset += (0.5-nIterOffset)/5;
//                     }
//             }
//         }
//         // if (nIters<100) nIterReduction*=2;
//         // printf("total weight reduction to %lf\n", totalWeight);
        
        
//         newObj = evaluate(instance, &newSol);

//         if (totalWeight!=newObj.finalweight){
//             printf("totalWeight constraint enforcement has unequal weights total %lf final %lf\n", totalWeight, newObj.finalweight);
//             exit(0);
//         }

//         if (newObj.objective >= bestObjective && newObj.finalweight <= instance->capacityOfKnapsack) {
//             for (int k=0;k<instance->numberOfItems;k++) sol->packing[k] = newSol.packing[k];

//             // for the stopping criterion: check if there was an actual improvement
//             if (newObj.objective > bestObjective) {
//                 // printf("iteration %d, best %lf prev %lf\n", i, newObj.objective, bestObjective);

//                 trainingOutput results = fit(classifier, inputs, sol->packing, 0.1, 10, 9, instance->numberOfItems);
//                 normalize(classifier);
//                 nIndices = 0;
//                 for (int k=0;k<instance->numberOfItems;k++){
//                     activations[k] = activation(classifier, inputs[k]);
//                     if (fabs(activations[k])<1 || (activations[k]>0 != sol->packing[k]))  
//                     indicesToChange[nIndices++]=k;
//                 }
//                 // printf("nIndices %d\n", nIndices);
//                 // printf("update loss %lf accuracy %lf\n",results.loss, results.accuracy);
//                 counter = 0;
//             }

//             bestObjective = newObj.objective;
//             totalWeightInit = newObj.finalweight;
//         } else counter++;

//     }

//     CC_IFFREE(newSol.packing, bool);
//     // CC_IFFREE(profitRatios, double);
//     // CC_IFFREE(rDist, double);
//     for (int i=0;i<instance->numberOfItems;i++) CC_IFFREE(inputs[i], double);
//     CC_IFFREE(inputs, double*);
//     CC_IFFREE(indicesToChange, int);
//     CC_IFFREE(activations, double);

//     // printf("CLASSIFIER OVER\n");

//     return 0;
// }

int quadraticHillClimber(ttpInstance *instance, ttpSolution *sol, logisticRegressor *quadratic_heuristic, int nonImprovementDuration, int generations, double *capacity){
    logisticRegressor copy;
    copy.weight = CC_SAFE_MALLOC(3, double);
    copyDouble(quadratic_heuristic->weight, copy.weight, 3);
    copy.bias = quadratic_heuristic->bias;
    copy.nWeights = quadratic_heuristic->nWeights;


    ttpSolution newSol;
    newSol.tour = sol->tour;
    newSol.isCityPacking = sol->isCityPacking;
    newSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    double initPercent = estimatePercent(instance, 0);

    quadraticHeuristicIterated(instance, sol, &copy, capacity);
    ttpObjective newObj = evaluate(instance, sol);
    double bestObj = newObj.objective;

    // printf("best initial obj %lf weight %lf percent %lf \n", bestObj, newObj.finalweight, initPercent);
    int counter = 0;

    for (int i=0;i<generations;i++){
        // printf("QHILLCLIMBER GENERATION %d\n",i);
        if (CHECK_TIME()) break;
        // for (int j=0;j<copy.nWeights;j++) copy.weight[j]=quadratic_heuristic->weight[j]*exp(randNormal(0,0.1*log(stdvsmean(instance, 3))));
        // copy.bias = quadratic_heuristic->bias*exp(randNormal(0,0.1*log(stdvsmean(instance, 3))));
        // for (int j=0;j<copy.nWeights;j++) copy.weight[j]=quadratic_heuristic->weight[j]*exp(randNormal(0,0.1));
        for (int j=0;j<copy.nWeights;j++) copy.weight[j]=quadratic_heuristic->weight[j]+randNormal(0,0.05);
        // copy.bias = quadratic_heuristic->bias*exp(randNormal(0,0.1));

        // for (int j=0;j<copy.nWeights;j++) copy.weight[j]=quadratic_heuristic->weight[j]+randNormal(0,0.2);
        // copy.bias = quadratic_heuristic->bias+randNormal(0,0.2);
        
        // printBoolArray(newSol.packing, instance->numberOfItems);
        // printf("\n");
        // quadraticHeuristicCapped(instance, &newSol, &copy, percent*instance->capacityOfKnapsack, true);
        newObj = quadraticHeuristicIterated(instance, &newSol, &copy, capacity);
        // printf("capacity %lf\n",capacity);
        // printBoolArray(newSol.packing, instance->numberOfItems);
        // printf("\n");

        // newObj = evaluate(instance, &newSol);
        // printf("new objective %lf\n",newObj.objective);

        if (newObj.objective>bestObj){
            bestObj = newObj.objective;
            normalize(&copy);
            // copyBool(newSol.packing, sol->packing, instance->numberOfItems);
            copyDouble(copy.weight,quadratic_heuristic->weight, 3);
            quadratic_heuristic->bias = copy.bias;
            counter = 0;
            // printf("(ga) generation %d best objective %lf weight %lf capacity %lf parameters: w ", i, bestObj, newObj.finalweight, *capacity);
            // printDoubleArray(copy.weight, 3);
            // printf(" b %lf\n", copy.bias);
        } else counter++;

        if (counter>nonImprovementDuration) break;
    }
    CC_IFFREE(newSol.packing, bool);
    freeRegressor(&copy);

    return 0;

}


int initializedQuadraticClimber(ttpInstance *instance, ttpSolution *sol, logisticRegressor *quadratic_heuristic, int nonImprovementDuration, int generations, double *capacity){
    ttpObjective initialObj = evaluate(instance, sol);
    bool *packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    for (int i=0;i<instance->numberOfItems;i++) packing[i] = sol->packing[i];
    quadraticHillClimber(instance, sol, quadratic_heuristic, nonImprovementDuration, generations, capacity);
    ttpObjective finalObj = evaluate(instance, sol);

    if (finalObj.objective<=initialObj.objective){
        for (int i=0;i<instance->numberOfItems;i++) sol->packing[i] = packing[i];
    }
    CC_IFFREE(packing, bool);
    return 0;
}


int linearHillClimber(ttpInstance *instance, ttpSolution *sol, int nonImprovementDuration, int generations, double *linear, double *percent){
    // Runs through 1.5*generations to determine the optimal value. 
    logisticRegressor regressor;
    regressor.weight = CC_SAFE_MALLOC(2, double);
    regressor.weight[0] = (double)(*linear);
    regressor.weight[1] = -1;
    regressor.nWeights=2;
    regressor.bias =  -1;

    linearHeuristicCapped(instance, sol, &regressor, *percent*instance->capacityOfKnapsack, true);
    ttpObjective obj = evaluate(instance, sol);
    double bestObj = obj.objective;

    ttpSolution newSol;
    newSol.tour = sol->tour;
    newSol.isCityPacking = sol->isCityPacking;
    newSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    // printf("iteration 0 best objective %.3e weight %lf linear %lf percent %lf\n", bestObj, obj.finalweight, *linear, *percent);

    int counter = 0;

    for (int i=0;i<generations;i++){
        if (CHECK_TIME()) break;
        double newLinear = *linear*exp(randNormal(0,0.5));
        double newPercent = *percent+randNormal(0,0.1);
        
        if (newPercent<0) newPercent = 0;
        if (newPercent>1) newPercent = 1;

        regressor.weight[0] = newLinear;
        
        linearHeuristicCapped(instance, &newSol, &regressor, newPercent*instance->capacityOfKnapsack, true);

        ttpObjective newObj = evaluate(instance, &newSol);
        // printf("newLinear %lf newPercent %lf newObj %.3e weight %lf\n", newLinear, newPercent, newObj.objective, newObj.finalweight);

        if (newObj.objective>bestObj){
            bestObj = newObj.objective;
            for (int i=0;i<instance->numberOfItems;i++) sol->packing[i] = newSol.packing[i];
            counter = 0;
            *linear = newLinear;
            *percent = newPercent;
            // printf("(ga) generation %d best objective %lf weight %lf linear %lf percent %lf\n", i, bestObj, newObj.finalweight, *linear, *percent);
        } else counter++;

        if (counter>nonImprovementDuration) break;
    }

    double lowLinear = log(0.5);
    double highLinear = log(2);
    double lowPct = *percent-0.15<0 ? 0 : *percent-0.15;
    double highPct = *percent+0.15>1 ? 1 : *percent+0.15;
    for (int climbIter = 0; climbIter<2; climbIter++){
        for (int i=0;i<generations/4;i++){
            if (CHECK_TIME()) break;
            double newLinear, newPercent;
            if (climbIter==0){
                newLinear = *linear*exp(lowLinear+(highLinear-lowLinear)*i*2/(double) generations);
                newPercent = *percent;
            } else {
                newLinear = *linear;
                newPercent = lowPct+(highPct-lowPct)*i*2/(double) generations;
            }
            regressor.weight[0] = newLinear;
            linearHeuristicCapped(instance, &newSol, &regressor, newPercent*instance->capacityOfKnapsack, true);
            ttpObjective newObj = evaluate(instance, &newSol);
            if (newObj.objective>bestObj){
                bestObj = newObj.objective;
                for (int i=0;i<instance->numberOfItems;i++) sol->packing[i] = newSol.packing[i];
                *linear = newLinear;
                *percent = newPercent;
                // printf("(local search) generation %d best objective %lf weight %lf profit %lf linear %lf percent %lf\n", i, bestObj, newObj.finalweight, newObj.profit, *linear, *percent);
            } 
        }
    }
    
    CC_IFFREE(newSol.packing, bool);
    freeRegressor(&regressor);

    // printf("linear heuristic best Obj %lf\n", bestObj);

    return 0;
}

int initializedLinearClimber(ttpInstance *instance, ttpSolution *sol, int nonImprovementDuration, int generations, double *linear, double *percent){
    ttpObjective initialObj = evaluate(instance, sol);
    bool *packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    for (int i=0;i<instance->numberOfItems;i++) packing[i] = sol->packing[i];
    linearHillClimber(instance, sol, nonImprovementDuration, generations, linear, percent);
    ttpObjective finalObj = evaluate(instance, sol);

    if (finalObj.objective<=initialObj.objective){
        for (int i=0;i<instance->numberOfItems;i++) sol->packing[i] = packing[i];
    }
    CC_IFFREE(packing, bool);
    return 0;
}

double singleGradient(ttpInstance *instance, ttpSolution *sol, int packingIdx){
    bool temp = sol->packing[packingIdx];
    double tempW = instance->weight[packingIdx];
    // instance->weight[packingIdx]=1;

    sol->packing[packingIdx] = true;
    ttpObjective withObj = evaluate(instance, sol);
    sol->packing[packingIdx] = false;
    ttpObjective withoutObj = evaluate(instance, sol);

    sol->packing[packingIdx] = temp;
    instance->weight[packingIdx]= tempW;
    // printf("dt %lf\n", withObj.time*instance->rentingRatio-withoutObj.time*instance->rentingRatio);
    return withObj.objective-withoutObj.objective;
}

int accumulateGradients(ttpInstance *instance, ttpSolution *sol, double *theta, double lr){
    double *weights = CC_SAFE_MALLOC(instance->numberOfNodes, double);
    double weight = 0;
    for (int i=0;i<instance->numberOfNodes;i++){
        if (i>0) {
            for (int j=0;j<instance->itemsPerCity;j++){
                int packingIdx = sol->isCityPacking ? getItemIndex(instance, sol->tour[i],j) : getPackingIndex(instance, i,j);

                if (sol->packing[packingIdx]) {
                    int itemIdx = getItemIndex(instance, sol->tour[i],j);
                    weight+=instance->weight[itemIdx];
                }
            }
        }
        weights[i]=weight;
    }

    double *wgrad = CC_SAFE_MALLOC(instance->numberOfNodes, double);
    double testTime = 0;
    for (int i=instance->numberOfNodes-1;i>0;i--){
        double dist = getDistance(instance, sol->tour[i], sol->tour[(i+1)%instance->numberOfNodes]);
        double speed = getSpeed(instance, weights[i]);
        speed = weights[i]>instance->capacityOfKnapsack ? instance->minSpeed : speed;

        int itemIdx = getItemIndex(instance, sol->tour[i],0);
        double curWeight = 1;
        double curProfit = instance->profit[itemIdx];

        double spd2 = getSpeed(instance, weights[i]-curWeight);
        spd2 = (weights[i]-curWeight)>instance->capacityOfKnapsack ? instance->minSpeed : spd2;
        testTime += instance->rentingRatio*dist*(1/speed-1/spd2);

        int packingIdx = sol->isCityPacking ? getItemIndex(instance, sol->tour[i],0) : getPackingIndex(instance, i,0);

        // double sg = singleGradient(instance, sol, packingIdx);
        
        wgrad[i] = ((i==instance->numberOfNodes-1) ? 0 : wgrad[i+1])+instance->rentingRatio*dist*(instance->maxSpeed-instance->minSpeed)/(instance->capacityOfKnapsack*speed*speed);
        // printf("i %d w %lf tt %lf profit1 %lf profit2 %lf sg %lf\n", i, wgrad[i]*curWeight, testTime,curProfit-wgrad[i]*curWeight,curProfit-testTime, sg);
    }


    
    int nChanged = 0;
    for (int i=1;i<instance->numberOfNodes;i++){
        for (int j=0;j<instance->itemsPerCity;j++){
            int packingIdx = sol->isCityPacking ? getItemIndex(instance, sol->tour[i],j) : getPackingIndex(instance, i,j);
            int itemIdx = getItemIndex(instance, sol->tour[i],j);
            double grad = instance->profit[itemIdx]-wgrad[i]*instance->weight[itemIdx];
            // printf("gradient %lf\n", grad);
            // double actualGrad = singleGradient(instance, sol, packingIdx);
            // printf("gradient %lf actual %lf profit %lf wpct %lf idx %d packed? %d\n", grad, actualGrad, instance->profit[itemIdx], instance->weight[itemIdx]/(double)instance->capacityOfKnapsack, i, sol->packing[packingIdx]);
            
            theta[packingIdx]+=lr*grad;
            if (theta[packingIdx]>1) theta[packingIdx]=1;
            if (theta[packingIdx]<0) theta[packingIdx]=0;
        }
    }

    printDoubleArray(theta, instance->numberOfItems);
    printf("\n");
    return 0;

}

int gradientClimber(ttpInstance *instance, ttpSolution *sol, int nSteps, double lr){
    double *theta = CC_SAFE_MALLOC(instance->numberOfItems, double);
    for (int i=0;i<instance->numberOfItems;i++) theta[i] = sol->packing[i] ? 1 : 0;
    for (int i=0;i<nSteps;i++){
        accumulateGradients(instance, sol, theta, lr);
        int nChanged=0;
        for (int j=0;j<instance->numberOfItems;j++){
            bool tempPacking = sol->packing[j];
            sol->packing[j] = theta[j]>0.5 ? true : false;
            if (sol->packing[j]!=tempPacking) nChanged++;
        }
        printf("number changed: %d\n", nChanged);
        ttpObjective obj = evaluate(instance, sol);
        printf("gradient iteration %d objective %lf\n", i, obj.objective);
        // exit(0);
    }
    return 0;
}
