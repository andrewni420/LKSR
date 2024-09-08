#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "util.h"
#include "math_util.h"
#include "ttp_util.h"
#include "logistic_regression.h"

int heuristic_evaluations = 0;


typedef struct val_i_j{
    double val;
    int i;
    int j;
} val_i_j;

void greedyRDistCapped(ttpInstance *instance, ttpSolution *sol, double capacity){
    double weight = 0;
    for (int i=0;i<instance->numberOfItems;i++) sol->packing[i]=false;
    // sol->packing = CC_SAFE_MALLOC(instance->numberOfNodes, bool);
    for (int nodeIndex=instance->numberOfNodes-1;nodeIndex>0;nodeIndex--){
        for (int itemNumber=instance->itemsPerCity-1;itemNumber>=0;itemNumber--){
            int itemIndex = getItemIndex(instance, sol->tour[nodeIndex], itemNumber);
            if (weight+instance->weight[itemIndex]<capacity) {
                weight+=instance->weight[itemIndex];
                int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, nodeIndex,itemNumber);
                sol->packing[packingIndex]=true;
            } else {
                return;
            }
        }
    }
} 

void greedyRDist(ttpInstance *instance, ttpSolution *sol){greedyRDistCapped(instance, sol, instance->capacityOfKnapsack);}


void greedyProfitRatioCapped(ttpInstance *instance, ttpSolution *sol, int capacity){
    int weight = 0;

    bool *packing =  CC_SAFE_MALLOC(instance->numberOfItems, bool);
    for (int i=0;i<instance->numberOfItems;i++) packing[i]=false;
    doubleArrayEntry *indices = CC_SAFE_MALLOC(instance->numberOfItems, doubleArrayEntry);
    for (int i=0;i<instance->numberOfItems;i++){
        indices[i].value = instance->profit[i]/instance->weight[i];
        indices[i].index = i;
    }

    qsort(indices, (size_t) instance->numberOfItems, sizeof(indices[0]), comparingDouble);

    for (int i=0;i<instance->numberOfItems;i++){
        if (weight+instance->weight[indices[i].index]<capacity){
            weight+=instance->weight[indices[i].index];
            packing[indices[i].index]=true;
        } else break;
    }

    // sol->packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    updateWithIndexPacking(instance, sol, packing);

    CC_IFFREE(packing, bool);
    CC_IFFREE(indices, doubleArrayEntry);
}
void greedyProfitRatio(ttpInstance *instance, ttpSolution *sol){ greedyProfitRatioCapped(instance,sol,instance->capacityOfKnapsack);}

void combinedHeuristic(ttpInstance *instance, ttpSolution *sol, double rDistCapacity, double profitRatioCapacity){
    ttpSolution rDistSol;
    rDistSol.tour = sol->tour;
    rDistSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    rDistSol.isCityPacking = sol->isCityPacking;
    greedyRDistCapped(instance, &rDistSol, instance->capacityOfKnapsack*rDistCapacity);

    ttpSolution profitRatioSol;
    profitRatioSol.tour = sol->tour;
    profitRatioSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    profitRatioSol.isCityPacking = sol->isCityPacking;
    greedyProfitRatioCapped(instance, &profitRatioSol, instance->capacityOfKnapsack*profitRatioCapacity);

    // sol->packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    for (int i=0;i<instance->numberOfItems;i++) sol->packing[i]=rDistSol.packing[i] || profitRatioSol.packing[i];
    CC_IFFREE(rDistSol.packing, bool);
    CC_IFFREE(profitRatioSol.packing, bool);
}

// void quadraticHeuristic(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor){
//     // printf("requested capacity %lf\n", capacity);
//     double weight = 0;
//     double distance = 0;
//     // bool *packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
//     // for (int i=0;i<instance->numberOfItems;i++) packing[i]=false;
//     double *input = CC_SAFE_MALLOC(regressor->nWeights, double);

//     for (int i = instance->numberOfNodes-1; i >0; i--) {
//         for (int j = 0; j < instance->itemsPerCity; j++) {
//             int itemIndex = getItemIndex(instance, sol->tour[i], j);
//             input[0] = instance->normalizedProfitRatios[itemIndex];
//             input[1] = instance->normalizedRDist[sol->tour[i]];
//             input[2] = input[0]*input[1];
//             int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, i,j);
//             if (weight+instance->weight[itemIndex]>instance->capacityOfKnapsack) sol->packing[packingIndex]=false;
//             else sol->packing[packingIndex]=activation(regressor, input)>0;
//             // sol->packing[packingIndex]=activation(regressor, input)>0;
//             weight+=sol->packing[packingIndex] ? instance->weight[itemIndex] : 0;
//             // printf("index %d packing idx %d weight %lf x %lf y %lf val %lf\n", itemIndex, packingIndex, weight, input[0], input[1], activation(regressor, input));
//         }
//     }

//     // printIntArray(sol->tour, instance->numberOfNodes+1);
//     // printf("\n");

//     // CC_IFFREE(packing, bool);
//     CC_IFFREE(input, double);
// }

int comparingVal_i_j(const void *a, const void *b){
    val_i_j *a_ = (val_i_j *) a;
    val_i_j *b_ = (val_i_j *) b;
    if (a_->val > b_->val) return -1;
    else if (a_->val < b_->val) return 1;
    else return 0;
}

ttpObjective quadraticHeuristicIterated(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, double *capacity){
    // printf("initial capacity %lf\n", *capacity);
    // static int consecutive_shrink_counter=0;
    static int step_count=1;  
    // static double first_shrink_fraction=0;
    // static int max_consecutive_shrinks=0;
    // printf("consecutive shrink counter %d max %d\n", consecutive_shrink_counter, max_consecutive_shrinks);

    double step_multiplier = 2;
    double step_divider = 8;
    double start_multiplier=2;

    // int shrinks_to_divide=20;
    // double fraction_to_multiply = 0.5;

    double weight = 0;
    double distance = 0;
    double *input = CC_SAFE_MALLOC(regressor->nWeights, double);
    ttpObjective bestObj;
    double bestCapacity = *capacity;

    val_i_j *indices =  CC_SAFE_MALLOC(instance->numberOfItems, val_i_j);
    for (int i = 1; i < instance->numberOfNodes; i++) {
        for (int j = 0; j < instance->itemsPerCity; j++) {
            int itemIndex = getItemIndex(instance, sol->tour[i], j);
            input[0] = instance->normalizedProfitRatios[itemIndex];
            input[1] = instance->normalizedRDist[sol->tour[i]];
            input[2] = input[0]*input[1];
            indices[itemIndex].val = activation(regressor,input);
            // printf("item index %d input %lf, %lf, value %lf\n",itemIndex, input[0], input[1], indices[itemIndex].value);
            indices[itemIndex].i = i;
            indices[itemIndex].j = j;
        }
    }

    qsort(indices, (size_t) instance->numberOfItems, sizeof(indices[0]), comparingVal_i_j);
    for (int i=0;i<instance->numberOfItems;i++) sol->packing[i]=false;
    
    int i=0;
    // int *itemIndices = CC_SAFE_MALLOC(instance->numberOfItems,int);
    // int ctr=0;
    
    for (;i<instance->numberOfItems;i++){
        int itemIndex = getItemIndex(instance, sol->tour[indices[i].i], indices[i].j);
        int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, indices[i].i,indices[i].j);
        if (weight+instance->weight[itemIndex]<=*capacity){
            // printf("pack item %d weight %lf\n",itemIndex, instance->weight[itemIndex]);
            weight+=instance->weight[itemIndex];
            sol->packing[packingIndex]=true;
            // itemIndices[ctr++]=itemIndex;
        } else break;
    }

    // qsort(itemIndices,ctr,sizeof(int),compareInt);
    // printf("quadratic heuristic items:\n");
    // printIntArray(itemIndices,ctr);
    // printf("\n");

    bestObj = evaluate(instance, sol);
    double bestObjective = bestObj.objective;
    
    //minimum unpacked index
    int minIndex = i;
    int very_first_index=i;
    heuristic_evaluations++;
    bestCapacity = bestObj.finalweight;

    int stepSize = step_count;
    int num_consecutive_shrinks = 0;
    bool start=true;
    int prevMinIndex = minIndex;
    while (stepSize>0){
        // printf("UPSTEP SIZE %d\n",stepSize);
        
        minIndex+=stepSize;
        if (minIndex>instance->numberOfItems) minIndex = instance->numberOfItems;
        for (;i<minIndex;i++){
            int itemIndex = getItemIndex(instance, sol->tour[indices[i].i], indices[i].j);
            int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, indices[i].i,indices[i].j);
            sol->packing[packingIndex]=true;
        }
        ttpObjective obj = evaluate(instance,sol);
        // printf("up step index %d size %d objective %lf bestObjective %lf\n", minIndex,stepSize,obj.objective,bestObjective);
        heuristic_evaluations++;
        if (obj.objective>bestObjective){
            bestObjective = obj.objective;
            bestObj = obj;
            bestCapacity = obj.finalweight;
            stepSize*=start?start_multiplier:step_multiplier;
            prevMinIndex = minIndex;
            // printf("improvement at %d shrinks, counter = %d, shrink_frac=%lf, step_count=%d\n",num_consecutive_shrinks, consecutive_shrink_counter, first_shrink_fraction, step_count);
            // if (num_consecutive_shrinks<2 && num_consecutive_shrinks>-1) {
            //     consecutive_shrink_counter=0;
            //     max_consecutive_shrinks=0;
            // } 
            // if (num_consecutive_shrinks>=2) max_consecutive_shrinks = (num_consecutive_shrinks<max_consecutive_shrinks || max_consecutive_shrinks==0) ? num_consecutive_shrinks : max_consecutive_shrinks;
            // if (num_consecutive_shrinks==0) first_shrink_fraction=0.9*first_shrink_fraction+0.1;
            // if (first_shrink_fraction>fraction_to_multiply) {
            //     step_count*=1.2;
            //     first_shrink_fraction=0;
            //     printf("multiply step by 1.5\n");
            // }
            // num_consecutive_shrinks=-1;
        } else {
            // minIndex-=stepSize;
            minIndex = prevMinIndex;
            if (stepSize<step_divider && stepSize>1) stepSize=1;
            else stepSize/=step_divider;
            for (;i>=minIndex;i--){
                int itemIndex = getItemIndex(instance, sol->tour[indices[i].i], indices[i].j);
                int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, indices[i].i,indices[i].j);
                sol->packing[packingIndex]=false;
            }
            i = minIndex;
            // if (num_consecutive_shrinks==0)first_shrink_fraction*=0.9;
            // if (num_consecutive_shrinks>-1) num_consecutive_shrinks++;
            
            // if (num_consecutive_shrinks==2) consecutive_shrink_counter++;
            // if (consecutive_shrink_counter*max_consecutive_shrinks>shrinks_to_divide) {
            //     step_count*=0.8;
            //     consecutive_shrink_counter=0;
            //     max_consecutive_shrinks=0;
            // }
        }
        start = false;
    }
    // if (num_consecutive_shrinks>=2) max_consecutive_shrinks = (num_consecutive_shrinks<max_consecutive_shrinks || max_consecutive_shrinks==0) ? num_consecutive_shrinks : max_consecutive_shrinks;
            

    stepSize = step_count;
    // num_consecutive_shrinks=0;
    start=true;
    prevMinIndex = minIndex;
    while (stepSize>0){
        // printf("DOWNSTEP SIZE %d\n",stepSize);
        minIndex-=stepSize;
        if (minIndex<0) minIndex=0;
        
        for (;i>=minIndex;i--){
            int itemIndex = getItemIndex(instance, sol->tour[indices[i].i], indices[i].j);
            int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, indices[i].i,indices[i].j);
            sol->packing[packingIndex]=false;
        }
        i=minIndex;
        ttpObjective obj = evaluate(instance,sol);
        // printf("down step index %d step %d obj %lf bestObj %lf\n", minIndex,stepSize, obj.objective,bestObjective);
        heuristic_evaluations++;
        // printf("i %d minIndex %d\n",i,minIndex);
        if (obj.objective>bestObjective){
            bestObjective = obj.objective;
            bestObj = obj;
            bestCapacity = obj.finalweight;
            stepSize*=start?start_multiplier:step_multiplier;
            prevMinIndex = minIndex;
            // printf("improvement at %d shrinks, counter = %d, shrink_frac=%lf, step_count=%d\n",num_consecutive_shrinks, consecutive_shrink_counter, first_shrink_fraction, step_count);
            // if (num_consecutive_shrinks<2 && num_consecutive_shrinks>-1) {
            //     consecutive_shrink_counter=0;
            //     max_consecutive_shrinks = 0;
            // }
            // if (num_consecutive_shrinks>=2) max_consecutive_shrinks = (num_consecutive_shrinks<max_consecutive_shrinks || max_consecutive_shrinks==0) ? num_consecutive_shrinks : max_consecutive_shrinks;
            // if (num_consecutive_shrinks==0) first_shrink_fraction=0.9*first_shrink_fraction+0.1;
            // if (first_shrink_fraction>fraction_to_multiply) {
            //     step_count*=1.2;
            //     first_shrink_fraction=0;
            //     printf("multiply step by 1.5\n");
            // }
            // num_consecutive_shrinks=-1;
        } else {
            if (minIndex==0 && prevMinIndex==0) stepSize=0;
            minIndex = prevMinIndex;
            if (stepSize<step_divider && stepSize>1) stepSize=1;
            else stepSize/=step_divider;
            for (;i<minIndex;i++){
                int itemIndex = getItemIndex(instance, sol->tour[indices[i].i], indices[i].j);
                int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, indices[i].i,indices[i].j);
                sol->packing[packingIndex]=true;
            }
            // if (num_consecutive_shrinks==0)first_shrink_fraction*=0.9;
            // if (num_consecutive_shrinks>-1) num_consecutive_shrinks++;
            // if (num_consecutive_shrinks==2) consecutive_shrink_counter++;
            // if (consecutive_shrink_counter*max_consecutive_shrinks>shrinks_to_divide) {
            //     step_count*=0.8;
            //     consecutive_shrink_counter=0;
            //     max_consecutive_shrinks=0;
            // }
        }
        start=false;
    }
    // if (num_consecutive_shrinks>=2) max_consecutive_shrinks = (num_consecutive_shrinks<max_consecutive_shrinks || max_consecutive_shrinks==0) ? num_consecutive_shrinks : max_consecutive_shrinks;
            
    // obj = evaluate(instance, sol);
    // printf("evaluations %d weight %lf bestObjective %lf\n",evaluations, *capacity, obj.objective);
    if (fabs((*capacity-bestCapacity)/(*capacity))<0.1){
        *capacity = bestCapacity;
    } else{
        *capacity = 0.5*(*capacity) + 0.5*bestCapacity;
    }
    
    // *capacity = bestCapacity;
    // printf("total difference %d\n",i-very_first_index);

    CC_IFFREE(indices, val_i_j);
    CC_IFFREE(input, double);
    return bestObj;
}


void quadraticHeuristicCapped(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, double capacity, bool sorted){
    // printf("requested capacity %lf\n", capacity);
    double weight = 0;
    double distance = 0;
    double *input = CC_SAFE_MALLOC(regressor->nWeights, double);

    val_i_j *indices = NULL;
    if (sorted) indices = CC_SAFE_MALLOC(instance->numberOfItems, val_i_j);
    for (int i = 1; i < instance->numberOfNodes; i++) {
        for (int j = 0; j < instance->itemsPerCity; j++) {
            int itemIndex = getItemIndex(instance, sol->tour[i], j);
            input[0] = instance->normalizedProfitRatios[itemIndex];
            input[1] = instance->normalizedRDist[sol->tour[i]];
            input[2] = input[0]*input[1];
            if (sorted){
                indices[itemIndex].val = activation(regressor,input);
                // printf("item index %d input %lf, %lf, value %lf\n",itemIndex, input[0], input[1], indices[itemIndex].value);
                indices[itemIndex].i = i;
                indices[itemIndex].j = j;
            } else {
                int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, i,j);
                sol->packing[packingIndex]=activation(regressor, input)>0;
            }
        }
    }

    if (sorted) {
        qsort(indices, (size_t) instance->numberOfItems, sizeof(indices[0]), comparingDouble);
        for (int i=0;i<instance->numberOfItems;i++){
            int itemIndex = getItemIndex(instance, sol->tour[indices[i].i], indices[i].j);
            int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, indices[i].i,indices[i].j);
            if (weight+instance->weight[itemIndex]<capacity){
                // printf("pack item %d weight %lf\n",indices[i].index, instance->weight[indices[i].index]);
                weight+=instance->weight[itemIndex];
                sol->packing[packingIndex]=true;
            } else break;
        }
    }

    CC_IFFREE(indices, val_i_j);
    CC_IFFREE(input, double);
    // CC_IFFREE(rDist, double);
    // CC_IFFREE(profitRatios, double);
}
void quadraticHeuristic(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, bool sorted){quadraticHeuristicCapped(instance,sol,regressor, instance->capacityOfKnapsack, sorted);}


void linearHeuristicCapped(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, double capacity, bool sorted){
    // printf("requested capacity %lf\n", capacity);
    double weight = 0;
    double distance = 0;
    bool *packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    for (int i=0;i<instance->numberOfItems;i++) packing[i]=false;
    double *input = CC_SAFE_MALLOC(regressor->nWeights, double);

    // double *rDist = reverseDistances(instance, sol->tour);
    // double *profitRatios =  CC_SAFE_MALLOC(instance->numberOfItems, double);
    // for (int i = 0; i < instance->numberOfItems; i++)
    //     profitRatios[i] = instance->profit[i] / instance->weight[i];

    doubleArrayEntry *indices = NULL;
    if (sorted) indices = CC_SAFE_MALLOC(instance->numberOfItems, doubleArrayEntry);
    for (int i = 1; i < instance->numberOfNodes; i++) {
        for (int j = 0; j < instance->itemsPerCity; j++) {
            int itemIndex = getItemIndex(instance, sol->tour[i], j);
            input[0] = instance->profitRatios[itemIndex];
            input[1] = instance->rDist[sol->tour[i]];
            if (sorted){
                indices[itemIndex].value = activation(regressor,input);
                // printf("item index %d input %lf, %lf, value %lf\n",itemIndex, input[0], input[1], indices[itemIndex].value);
                indices[itemIndex].index = itemIndex;
            } else {
                int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, i,j);
                sol->packing[packingIndex]=activation(regressor, input)>0;
            }
        }
    }

    if (sorted) {
        qsort(indices, (size_t) instance->numberOfItems, sizeof(indices[0]), comparingDouble);
        for (int i=0;i<instance->numberOfItems;i++){
            if (weight+instance->weight[indices[i].index]<capacity){
                // printf("pack item %d weight %lf\n",indices[i].index, instance->weight[indices[i].index]);
                weight+=instance->weight[indices[i].index];
                packing[indices[i].index]=true;
            } else break;
        }
        updateWithIndexPacking(instance, sol, packing);
    }

    CC_IFFREE(indices, doubleArrayEntry);
    CC_IFFREE(packing, bool);
    CC_IFFREE(input, double);
    // CC_IFFREE(rDist, double);
    // CC_IFFREE(profitRatios, double);
}
void linearHeuristic(ttpInstance *instance, ttpSolution *sol, logisticRegressor *regressor, bool sorted){linearHeuristicCapped(instance,sol,regressor, instance->capacityOfKnapsack, sorted);}


double estimateWeightRatio(ttpInstance *instance, int mode){
    if (mode ==0){
        // add(add(sub(mul(div(X3, X2), 0.314), X0), sub(div(div(X3, X1), 0.248), X0)), mul(X7, 0.314))
        double x3 = instance->capacityOfKnapsack;
        double x2 = instance->rentingRatio;
        double x1 = instance->numberOfItems;
        double x6 = 0.1;
        double x0 = instance->numberOfNodes;
        double temp = x3*0.537*0.537/(x6*x1);
        double temp2 = x3*0.537*0.675;
        temp2 = x3/sqrt(x0)+temp2;
        return temp+temp2;
    }
    return 0;
}

double estimateLinear(ttpInstance *instance, int estimator){
    double x0 = instance->numberOfNodes;
    double x2 = instance->rentingRatio;
    double x3 = instance->capacityOfKnapsack;
    double x4 = instance->itemsPerCity;
    double x5 = instance->maxSpeed;
    double x6 = instance->minSpeed;

    if (estimator==0){
        //capacity / renting ratio / 1.75
        return x3/x2/1.75;
    } else if (estimator==1){
        // sub(sub(mul(mul(div(X3, X2), -0.706), -0.691), mul(div(X3, div(X4, mul(-0.691, X6))), 0.200)), X0)
        double temp = x3/x2*(-0.706)*(-0.691);
        double temp2 = 0.2*x3/(x4/(-0.691*x6));
        return temp - temp2 - x0;
    } else if (estimator==2){
        // div(div(X3, X2), add(0.993, add(div(X5, add(0.993, X5)), div(X5, X7))))
        double totalWeight = 0;
        for (int i=0;i<instance->numberOfItems;i++) totalWeight+=instance->weight[i];
        int x7 = (int)((11*instance->capacityOfKnapsack/totalWeight)+1E-3);
        double temp = x5/(x5+0.993)+x5/x7+0.993;
        return x3/x2/temp;
    } else {
        return 0;
    }
    
}

double estimatePercent(ttpInstance *instance, int estimator){
    double totalWeight = 0;
    for (int i=0;i<instance->numberOfItems;i++) totalWeight+=instance->weight[i];
    int capacityFactor = (int)((11*instance->capacityOfKnapsack/totalWeight)+1E-3);

    if (estimator==0){
        return 0.997576762823811-0.02626174*capacityFactor;
    } else if (estimator==1){
        double percents[] = {0.974709344085561, 0.95237689413325, 0.9203905922557376, 0.8862318850248149, 0.8592531312054973, 0.8344101398181115, 0.8101255521970334, 0.7895009140494049, 0.7629576476369783, 0.7414157699418626};
        return percents[capacityFactor-1];
    } else {
        return 0;
    }
    
}

ttpObjective dh(ttpInstance *instance, ttpSolution *sol){
    calcRDist(instance, sol->tour);
    val_i_j *indices = CC_SAFE_MALLOC(instance->numberOfItems, val_i_j);
    for (int i = 1; i < instance->numberOfNodes; i++) {
        for (int j = 0; j < instance->itemsPerCity; j++) {
            int itemIndex = getItemIndex(instance, sol->tour[i], j);
            double speed = instance->maxSpeed-(instance->maxSpeed-instance->minSpeed)*instance->weight[itemIndex]/instance->capacityOfKnapsack;
            indices[itemIndex].val = instance->profit[itemIndex]-instance->rentingRatio*instance->rDist[sol->tour[i]]/speed;
            indices[itemIndex].i = i;
            indices[itemIndex].j = j;
        }
    }

    for (int i=0;i<instance->numberOfItems;i++) sol->packing[i]=false;

    double weight =0;
    ttpObjective bestObj = evaluate(instance,sol);
    double bestObjective = bestObj.objective;
    int num_evals=0;
    qsort(indices, (size_t) instance->numberOfItems, sizeof(indices[0]), comparingVal_i_j);
    for (int i=0;i<instance->numberOfItems;i++){
        printf("i %d\n",i);
        int itemIndex = getItemIndex(instance, sol->tour[indices[i].i], indices[i].j);
        int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, indices[i].i,indices[i].j);
        if (weight+instance->weight[itemIndex]<instance->capacityOfKnapsack){
            
            sol->packing[packingIndex]=true;
            ttpObjective newObj = evaluate(instance,sol);
            num_evals++;
            if (newObj.objective>bestObj.objective){
                bestObj = newObj;
                weight+=instance->weight[itemIndex];
            } else {
                sol->packing[packingIndex] = false;
            }
        }
    }

    printf("{\"objective: %lf, \"num_evals\":%d\"}\n",bestObj.objective,num_evals);

    return bestObj;
}

void packIterative_inner(ttpInstance *instance, ttpSolution *sol, double Alpha) {
    sol->isCityPacking = true;
    int         numItems = instance->numberOfItems;
    int         muItems  = numItems/100;
    bool        jFlag;
    if (muItems > 1)
        jFlag   = true;
    else {    
        jFlag   = false;
        muItems = 2;
        
    }
    long    remSpace = instance->capacityOfKnapsack;
    doubleArrayEntry *scores;
    scores = CC_SAFE_MALLOC(numItems,doubleArrayEntry);

    for (int itemId  = 0; itemId < numItems; itemId++) {
        sol->packing[itemId] = false;
        int  cityId  = getCityFromItemIndex(instance, itemId);
        //if (Len[cityId]>0)    
        scores[itemId].value=pow(instance->profitRatios[itemId],Alpha) / instance->rDist[cityId];
        scores[itemId].index = itemId;
        //else
        //    scores[itemId] = std::numeric_limits<double>::infinity();
    }

    
    ttpObjective bestObj = evaluate(instance, sol);
    qsort(scores,numItems,sizeof(doubleArrayEntry),comparingDouble);

    // for (int i=0;i<instance->numberOfItems;i++){
    //     printf("(%d, %lf)\n",scores[i].index, scores[i].value);
    // }
    // exit(0);
    
    int     index = 0, numPicked = 0;
    while (index<numItems && remSpace>0 && muItems>=2) {
        int itemId = scores[index].index;
        if (remSpace >= (long) instance->weight[itemId]) {
            sol->packing[itemId] = true;
            remSpace    -= (long) instance->weight[itemId];
        }
        if (index%muItems==0 || !jFlag) {
            ttpObjective newObj = evaluate(instance,sol);
            if (newObj.objective >= bestObj.objective) {
                numPicked  = index + 1;
                bestObj = newObj;
            }
            else {
                while (index >= numPicked)
                {
                    itemId = scores[index].index;
                    if (sol->packing[itemId]) {
                        sol->packing[itemId] = false;
                        remSpace    += (long) instance->weight[itemId];
                    }
                    index--;
                }
                muItems = ceil(muItems/2.0);
            }
        }    
        index++;
    }
    CC_IFFREE(scores,doubleArrayEntry);
}


int packIterative(ttpInstance *instance, ttpSolution *sol, double C, double D, int Q) { /// c=5, d=2.5, q=20
    
    double   epsilon   = 0.1;
    int      numCities = instance->numberOfNodes;
    int      numItems  = instance->numberOfItems;

    ttpSolution     leftPlan;
    leftPlan.packing = CC_SAFE_MALLOC(numItems, bool);
    for (int i=0;i<numItems;i++) leftPlan.packing[i]=false;
    leftPlan.tour=  sol->tour;
    leftPlan.isCityPacking = sol->isCityPacking;
    
    ttpSolution     middlePlan;
    middlePlan.packing = CC_SAFE_MALLOC(numItems, bool);
    for (int i=0;i<numItems;i++) middlePlan.packing[i]=false;
    middlePlan.tour=  sol->tour;
    middlePlan.isCityPacking = sol->isCityPacking;

    ttpSolution     rightPlan;
    rightPlan.packing = CC_SAFE_MALLOC(numItems, bool);
    for (int i=0;i<numItems;i++) rightPlan.packing[i]=false;
    rightPlan.tour=  sol->tour;
    rightPlan.isCityPacking = sol->isCityPacking;

    // printf("C %lf  D %lf Q %lf\n", C, D, Q);


    packIterative_inner(instance, &leftPlan, C-D); 
    ttpObjective leftTotalGain  = evaluate(instance,&leftPlan);  
    // printf("eval alpha %lf obj %lf\n",C-D, leftTotalGain.objective);  
    
    packIterative_inner(instance, &middlePlan,C); 
    ttpObjective middleTotalGain = evaluate(instance,&middlePlan);  
    // printf("eval alpha %lf obj %lf\n",C, middleTotalGain.objective);  
    
    packIterative_inner(instance, &rightPlan,C+D); 
    ttpObjective rightTotalGain = evaluate(instance,&rightPlan);  
    // printf("eval alpha %lf obj %lf\n",C+D, rightTotalGain.objective);  
    
    ttpSolution   bestPlan = middlePlan;
    
    for (int i = 0; i < Q; i++) {
        if (leftTotalGain.objective>middleTotalGain.objective && rightTotalGain.objective>middleTotalGain.objective) {
            if (leftTotalGain.objective > rightTotalGain.objective) {
                middleTotalGain = leftTotalGain;  
                bestPlan = leftPlan;
                C = C - D;
            }
            else {
                middleTotalGain = rightTotalGain;  
                bestPlan = rightPlan;
                C = C + D;
            }
        }
        else if (leftTotalGain.objective > middleTotalGain.objective) {
            middleTotalGain = leftTotalGain;  
            bestPlan = leftPlan;
            C = C - D;
        } 
        else if (rightTotalGain.objective > middleTotalGain.objective) {
            middleTotalGain = rightTotalGain;  
            bestPlan = rightPlan;
            C = C + D;
        }
        
        D /= 2;

        packIterative_inner(instance, &leftPlan, C-D);
        leftTotalGain  = evaluate(instance, &leftPlan);  
        // printf("eval alpha %lf obj %lf\n",C-D, leftTotalGain.objective);  
        
        packIterative_inner(instance, &rightPlan, C+D);
        rightTotalGain = evaluate(instance, &rightPlan); 
        // printf("eval alpha %lf obj %lf\n",C+D, rightTotalGain.objective);   
        
        if (abs(leftTotalGain.objective  - middleTotalGain.objective) < epsilon && 
            abs(rightTotalGain.objective - middleTotalGain.objective) < epsilon)
            break;
    }    
    
    copyBool(bestPlan.packing,sol->packing,instance->numberOfItems);
    
    CC_IFFREE(leftPlan.packing, bool);
    CC_IFFREE(rightPlan.packing, bool);
    CC_IFFREE(middlePlan.packing, bool);
    
    return 0;
} 