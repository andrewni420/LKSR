#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "util.h"
#include "math_util.h"
#include "ttp_util.h"
#include "kdtree.h"
#include "linkern.h"
#include "delaunay.h"

double zeit_initial = -1;
double zeit_bound = 10;

double getSpeed(ttpInstance *instance, double weight){
    return (instance->maxSpeed-weight*(instance->maxSpeed-instance->minSpeed)/instance->capacityOfKnapsack);
}

int getItemIndex(ttpInstance *instance, int city, int itemNumber){return (city-1)+itemNumber*(instance->numberOfNodes-1);}
int getCityFromItemIndex(ttpInstance *instance, int itemIndex){return itemIndex%(instance->numberOfNodes-1)+1;}
int itemNumberFromItemIndex(ttpInstance *instance, int itemIndex){return itemIndex/(instance->numberOfNodes-1);}
int getPackingIndex(ttpInstance *instance, int tourIdx, int itemNumber){return (tourIdx-1)*instance->itemsPerCity+itemNumber;}
int tourIndexFromPackingIndex(ttpInstance *instance, int packingIndex){return packingIndex/instance->itemsPerCity+1;}
int itemNumberFromPackingIndex(ttpInstance *instance, int packingIndex){return packingIndex%instance->itemsPerCity;}
int getCityPackingIndex(ttpInstance *instance, int city, int itemNumber){return (city-1)*instance->itemsPerCity+itemNumber;}
bool isPacked(ttpInstance *instance, ttpSolution *sol, int tourIdx, int itemNumber){
  return sol->isCityPacking ? getItemIndex(instance, sol->tour[tourIdx], itemNumber) : getPackingIndex(instance, tourIdx, itemNumber);
}

void objective_print(ttpObjective *obj){
  printf("(d %lf w %lf t %lf p %lf obj %lf)",obj->distance,obj->finalweight,obj->time,obj->profit,obj->objective);
}

void sol_print(ttpInstance *instance, ttpSolution *sol){
  printf("Sol(\tisCityPacking %d\n\ttour: ",sol->isCityPacking);
  printIntArray(sol->tour,instance->numberOfNodes);
  printf("\n\tpacking: ");
  printBoolArray(sol->packing, instance->numberOfItems);
  printf(")\n");
}

int copySolution(ttpInstance *instance, ttpSolution *dest, ttpSolution *src){
  dest->isCityPacking = src->isCityPacking;
  for (int i=0;i<instance->numberOfItems;i++) dest->packing[i]=src->packing[i];
  for (int i=0;i<instance->numberOfNodes;i++) dest->tour[i]=src->tour[i];
  return 0;
}

void init_params(adaptive_params *params){
  params->probImprove = 0;
  params->Gmax_threshold = 0;
  params->check_objective_frequency = 0;
  params->upper_quadratic_max_nonimprovement = 0;
  params->upper_quadratic_max = 0;
  params->upper_mbfs_max = 0;
  params->quadratic_max_nonimprovement = 0;
  params->quadratic_max = 0;
  params->mbfs_max = 0;
  params->numOutputs = 0;
  params->outIters = NULL;
  params->pgch_round_counter =0;
  params->pgch_nonimprovement_counter=0;
  params->normal_nonimprovement_counter=0;
}


void set_params(ttpInstance *instance, adaptive_params *params, double time){
    params->probImprove = 1;
    params->normal_rounds=0;
    params->ttp_rounds=1;
    params->pgch_rounds=0;
    params->mbfs_sorted = false;
    params->Gmax_threshold = -CCutil_MAXDOUBLE;
    params->check_objective_frequency = 1;
    params->upper_quadratic_max_nonimprovement = 10;
    params->upper_quadratic_max = 100;
    params->upper_mbfs_max = -1; 
    params->upper_mbfs_sorted = false;
    params->quadratic_max_nonimprovement = 10;
    params->quadratic_max = 100;
    params->mbfs_max = -1; 
    params->kicktype = CC_LK_WALK_KICK;
    params->fast_pgch = false;
    params->initial_pgch_max = -1;
    params->upper_improvement_prob = 1;
    // int maxIter = instance->numberOfNodes>MAXLKITERS ? MAXLKITERS : instance->numberOfNodes;
    int maxIter = NUMLKITERS(instance->numberOfNodes);
    int numIncreasing = (int)log2(maxIter/2);
    int numDecreasing = numIncreasing/2;
    int num_out = numIncreasing+numDecreasing+3;
    
    CC_IFFREE(params->outIters, int);
    params->numOutputs = num_out; 
    params->outIters = CC_SAFE_MALLOC(num_out, int); 
    params->outIters[0]=0;
    for (int i=0;i<=numIncreasing;i++){
      int iter = exp2(i);
      params->outIters[i+1]=iter;
    } 

    for (int i=0;i<=numDecreasing;i++){
      int iter = maxIter-exp2(2*(numDecreasing-i));
      params->outIters[i+numIncreasing+2]=iter;
    }

    int pgch_time_limit = 0;

    if (instance->numberOfNodes<1000){
      pgch_time_limit = 150;
      if (time>150){
        params->pgch_rounds=1;
        params->ttp_rounds=3;
        params->kicktype = CC_LK_RANDOM_KICK;
      }
    }
    if (instance->numberOfNodes>1000){
      pgch_time_limit=100;
      if (instance->itemsPerCity>5) pgch_time_limit = 150;
      params->normal_rounds = 1;
      params->ttp_rounds = 0;
      params->upper_improvement_prob = 0.1;
      // if (time>1000){
      int normal_time_limit = 100;
      if (instance->itemsPerCity>5) normal_time_limit = 200;

      if (time>normal_time_limit || params->normal_nonimprovement_counter>5000){
        params->normal_rounds = 0;
        params->ttp_rounds=1;
        params->probImprove = 0.1;
        params->quadratic_max = 100;
        params->quadratic_max_nonimprovement = 10;
      }
    }
    if (instance->numberOfNodes>5000){
      pgch_time_limit=150;
      params->normal_rounds = 1;
      params->ttp_rounds = 0;
      params->upper_improvement_prob = 0.1;
      // if (time>1000){
      if (time>300 || params->normal_nonimprovement_counter>5000){
        params->normal_rounds = 0;
        params->ttp_rounds=1;
        params->probImprove = 0.1;
        params->quadratic_max = 100;
        params->quadratic_max_nonimprovement = 10;
      }
    }
    if (instance->numberOfNodes>10000){
      pgch_time_limit=300;
      params->normal_rounds = 1;
      params->ttp_rounds = 0; 
      params->probImprove = 0;
      params->upper_quadratic_max = 10;
      params->upper_quadratic_max_nonimprovement = 3;
      params->upper_mbfs_max = 100;
      params->upper_improvement_prob = 0.01;
      // if (instance->itemsPerCity)
    }
    
    if (time<pgch_time_limit && params->pgch_nonimprovement_counter<500 && params->pgch_round_counter<5000){
      params->pgch_rounds=1;
      params->ttp_rounds=0;
      params->normal_rounds=0;
    }
}

int objective_copy(ttpObjective *src, ttpObjective *dest){
  dest->distance = src->distance;
  dest->finalweight = src->finalweight;
  dest->objective = src->objective;
  dest->profit = src->profit;
  dest->time = src->time;
  return 0;
}

void free_params(adaptive_params *params){
  CC_IFFREE(params->outIters, int);
}

void initInstance (ttpInstance *instance)
{
    instance->x = (double *) NULL;
    instance->y = (double *) NULL;
    instance->weight = (double *) NULL;
    instance->profit = (double *) NULL;
    instance->numberOfItems = 0;
    instance->numberOfNodes = 0;
    instance->maxSpeed = 0.0;
    instance->minSpeed = 0.0;
    instance->capacityOfKnapsack = 0.0;
    instance->rentingRatio=0.0;
    instance->profitRatios = (double *)NULL;
    instance->normalizedProfitRatios = (double *)NULL;
    instance->rDist = (double *)NULL;
    instance->normalizedRDist = (double *) NULL;
    instance->neighbors = (intvector *) NULL;
    instance->itemIndex = (int **) NULL;
    instance->packingIndex = (int **) NULL;
    instance->speed_factor = 0;
}

void freeInstance (ttpInstance *instance)
{
    CC_IFFREE (instance->x, double);
    CC_IFFREE (instance->y, double);
    CC_IFFREE (instance->weight, double);
    CC_IFFREE (instance->profit, double);
    CC_IFFREE (instance->rDist, double);
    CC_IFFREE (instance->normalizedRDist, double);
    CC_IFFREE (instance->profitRatios, double);
    CC_IFFREE (instance->normalizedProfitRatios, double);
    for (int i=0;i<instance->numberOfNodes;i++){
      CC_IFFREE(instance->itemIndex[i], int);
      CC_IFFREE(instance->packingIndex[i], int);
    }
    CC_IFFREE (instance->itemIndex, int *);
    CC_IFFREE (instance->packingIndex, int *);
}

void initSolution(ttpSolution *sol){
    sol->tour = (int *) NULL;
    sol->packing = (bool *) NULL;
    sol->isCityPacking=false;
}

void freeSolution(ttpSolution *sol){
    CC_IFFREE (sol->tour, int);
    CC_IFFREE (sol->packing, bool);
}

int generate_neighbors(CCdatagroup *dat, CCrandstate *rstate, ttpInstance *instance){
  int wantlist = 1;
  int   edgesCount;
  int   *edgesList;
  int local_quadtry = 200;
  int numCities=  instance->numberOfNodes;

  CCkdtree localkt;
  CCkdtree_build (&localkt, numCities, dat, (double *) NULL, rstate);

  // CCkdtree_k_nearest (&localkt, numCities, local_quadtry,
  //                         dat, (double *) NULL, 1, &edgesCount, &edgesList,
  //                         1, rstate);

  // CCkdtree_quadrant_k_nearest (&localkt, numCities, local_quadtry,
  //                         dat, (double *) NULL, 1, &edgesCount, &edgesList,
  //                         1, rstate);//CCkdtree_quadrant_k_nearest
  // printf("Num edges %d\n",edgesCount);
  CCedgegen_delaunay(numCities, dat, wantlist, &edgesCount, &edgesList);
  instance->neighbors = CC_SAFE_MALLOC(instance->numberOfNodes, intvector);
    for (int i=0;i<instance->numberOfNodes;i++) intvector_alloc(&instance->neighbors[i],8);

  for (int edgeIndex = 0; edgeIndex < edgesCount; edgeIndex++) {
        int cityId1 = edgesList[2*edgeIndex];
        int cityId2 = edgesList[2*edgeIndex+1];
        intvector_append(&instance->neighbors[cityId1],cityId2);
        intvector_append(&instance->neighbors[cityId2],cityId1);
  }
  // printf("Neighbors:\n");
  // for (int i=0;i<instance->numberOfNodes;i++){
  //   printf("%d ",instance->neighbors[i].len);
  //   intvector_print(&instance->neighbors[i], stdout);
  //   printf("\n");
  // }
  // exit(0);

  CC_IFFREE(edgesList, int);
  return 0;
}

int sol_from_file(ttpInstance *instance, ttpSolution *sol, char *fileName){
  sol->isCityPacking = true;
  FILE *file=fopen(fileName, "r");

  sol->packing = CC_SAFE_MALLOC(instance->numberOfItems,bool);
  for (int i=0;i<instance->numberOfItems;i++) sol->packing[i]=false;
  sol->tour = CC_SAFE_MALLOC(instance->numberOfNodes+1,int);

  int i=0;
  int num;
  while(fscanf(file, "%d", &num) > 0) {
      if (i<=instance->numberOfNodes) sol->tour[i] = num;
      else sol->packing[num] = true; 
      i++;
  }

  fclose(file);
  return 0;
}

int read_ttp_text (char *datname, ttpInstance *instance)
{
    int i;
    FILE *datin = (FILE *) NULL;
    
    instance->weight = (double *) NULL;
    instance->profit = (double *) NULL;
    instance->x = (double *) NULL;
    instance->y = (double *) NULL;

    datin = fopen (datname, "r");
    if (datin == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        goto FAILURE;
    }
    fscanf (datin, "%d%d%lf%lf%lf%lf", &(instance->numberOfNodes),  &(instance->numberOfItems), &(instance->maxSpeed), &(instance->minSpeed), &(instance->capacityOfKnapsack), &(instance->rentingRatio));

    instance->itemsPerCity = instance->numberOfItems/(instance->numberOfNodes-1);

    printf ("nNodes = %d, nItems = %d, maxSpeed = %lf, minSpeed = %lf, capacity = %lf, renting ratio = %lf\n", instance->numberOfNodes, instance->numberOfItems, instance->maxSpeed, instance->minSpeed, instance->capacityOfKnapsack, instance->rentingRatio); fflush (stdout);
    instance->x = CC_SAFE_MALLOC (instance->numberOfNodes, double);
    instance->y = CC_SAFE_MALLOC (instance->numberOfNodes, double);
    instance->weight = CC_SAFE_MALLOC (instance->numberOfItems, double);
    instance->profit = CC_SAFE_MALLOC (instance->numberOfItems, double);
    if (instance->x == (double *) NULL ||
        instance->y == (double *) NULL ||
        instance->weight ==(double *) NULL ||
        instance->profit == (double *) NULL
        ) {
        goto FAILURE;
    }
    for (i=0; i<instance->numberOfNodes; i++) {
        fscanf (datin, "%lf%lf", &(instance->x[i]), &(instance->y[i]));
    }
    for (i=0; i<instance->numberOfItems; i++) {
        fscanf (datin, "%lf%lf", &(instance->profit[i]), &(instance->weight[i]));
    }
    fclose (datin);
    datin = (FILE *) NULL;

    instance->profitRatios =  CC_SAFE_MALLOC(instance->numberOfItems, double);
    instance->normalizedProfitRatios =  CC_SAFE_MALLOC(instance->numberOfItems, double);
    for (int i = 0; i < instance->numberOfItems; i++)
        instance->profitRatios[i] = instance->profit[i] / instance->weight[i];
        instance->normalizedProfitRatios[i] = instance->profitRatios[i];
    normalizeMAD(instance->normalizedProfitRatios, instance->numberOfItems);
    instance->rDist = CC_SAFE_MALLOC(instance->numberOfNodes,double);
    instance->normalizedRDist = CC_SAFE_MALLOC(instance->numberOfNodes,double);
    instance->knapsackType = knapsackType(instance);
    instance->capacityFactor = capacityFactor(instance);

    return 0;
    
 FAILURE:
    CC_IFFREE (instance->x, double);
    CC_IFFREE (instance->y, double);
    CC_IFFREE (instance->weight, double);
    CC_IFFREE (instance->profit, double);
    if (datin != (FILE *) NULL) fclose (datin);
    return 1;
}

int read_ttp_tsplib (char *datname, CCdatagroup *dat, CCrandstate *rstate, ttpInstance *instance)
{
    int i;
    FILE *datin = (FILE *) NULL;
    
    instance->weight = (double *) NULL;
    instance->profit = (double *) NULL;
    instance->x = (double *) NULL;
    instance->y = (double *) NULL;

    datin = fopen (datname, "r");
    if (datin == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        goto FAILURE;
    }

    char mystr[256];
    fgets(mystr, 256, datin);
    fgets(mystr, 256, datin);

    fscanf(datin, "%*s %d", &(instance->numberOfNodes));
    fscanf(datin, "%*s %*s %*s %d", &(instance->numberOfItems));
    fscanf(datin, "%*s %*s %*s %lf", &(instance->capacityOfKnapsack));
    fscanf(datin, "%*s %*s %lf", &(instance->minSpeed));
    fscanf(datin, "%*s %*s %lf", &(instance->maxSpeed));
    fscanf(datin, "%*s %*s %lf", &(instance->rentingRatio));

    fgets(mystr, 256, datin);
    fgets(mystr, 256, datin);
    fgets(mystr, 256, datin);

    instance->itemsPerCity = instance->numberOfItems/(instance->numberOfNodes-1);

    printf ("nNodes = %d, nItems = %d, maxSpeed = %lf, minSpeed = %lf, capacity = %lf, renting ratio = %lf\n", instance->numberOfNodes, instance->numberOfItems, instance->maxSpeed, instance->minSpeed, instance->capacityOfKnapsack, instance->rentingRatio); fflush (stdout);
    instance->x = CC_SAFE_MALLOC (instance->numberOfNodes, double);
    instance->y = CC_SAFE_MALLOC (instance->numberOfNodes, double);
    instance->weight = CC_SAFE_MALLOC (instance->numberOfItems, double);
    instance->profit = CC_SAFE_MALLOC (instance->numberOfItems, double);
    if (instance->x == (double *) NULL ||
        instance->y == (double *) NULL ||
        instance->weight ==(double *) NULL ||
        instance->profit == (double *) NULL
        ) {
        goto FAILURE;
    }

    
    for (i=0; i<instance->numberOfNodes; i++) {
        fscanf (datin, "%*d%lf%lf", &(instance->x[i]), &(instance->y[i]));
    }

    fgets(mystr, 256, datin);
    fgets(mystr, 256, datin);

    for (i=0; i<instance->numberOfItems; i++) {
        fscanf (datin, "%*d%lf%lf%*d", &(instance->profit[i]), &(instance->weight[i]));
    }
    fclose (datin);
    datin = (FILE *) NULL;


    instance->profitRatios =  CC_SAFE_MALLOC(instance->numberOfItems, double);
    instance->normalizedProfitRatios =  CC_SAFE_MALLOC(instance->numberOfItems, double);
    for (int i = 0; i < instance->numberOfItems; i++)
    {
        instance->profitRatios[i] = instance->profit[i] / instance->weight[i];
        instance->normalizedProfitRatios[i] = instance->profitRatios[i];
    }
    normalizeMAD(instance->normalizedProfitRatios, instance->numberOfItems);
    instance->rDist = CC_SAFE_MALLOC(instance->numberOfNodes,double);
    instance->normalizedRDist = CC_SAFE_MALLOC(instance->numberOfNodes,double);
    instance->knapsackType = knapsackType(instance);
    instance->capacityFactor = capacityFactor(instance);
    generate_neighbors(dat, rstate, instance);

    instance->itemIndex = CC_SAFE_MALLOC(instance->numberOfNodes, int *);
    instance->packingIndex = CC_SAFE_MALLOC(instance->numberOfNodes, int *);
    instance->itemIndex[0]=NULL;
    instance->packingIndex[0]=NULL;
    for (int i=1;i<instance->numberOfNodes;i++){
      instance->itemIndex[i] = CC_SAFE_MALLOC(instance->itemsPerCity, int);
      instance->packingIndex[i] = CC_SAFE_MALLOC(instance->itemsPerCity, int);
      for (int j=0;j<instance->itemsPerCity;j++){
        instance->itemIndex[i][j]=getItemIndex(instance, i, j);
        instance->packingIndex[i][j]=getPackingIndex(instance, i, j);
      }
    }
    instance->speed_factor = (instance->maxSpeed-instance->minSpeed)/instance->capacityOfKnapsack;

    return 0;
    
 FAILURE:
    CC_IFFREE (instance->x, double);
    CC_IFFREE (instance->y, double);
    CC_IFFREE (instance->weight, double);
    CC_IFFREE (instance->profit, double);
    if (datin != (FILE *) NULL) fclose (datin);
    return 1;
}

int calcRDist(ttpInstance *instance, int *tour){
    double dist = 0;
    for (int i=instance->numberOfNodes-1;i>=0;i--){
        dist+=getDistance(instance,tour[i],tour[(i+1)%instance->numberOfNodes]);
        instance->rDist[tour[i]]=dist;
        instance->normalizedRDist[tour[i]] = dist;
    }
    normalizeMAD(&(instance->normalizedRDist[1]), instance->numberOfNodes-1);
    return 0;
}

double getDistance(ttpInstance *instance, int i, int j){
    double t1 = instance->x[i] - instance->x[j], t2 = instance->y[i] - instance->y[j];
    return (ceil (sqrt (t1 * t1 + t2 * t2)));
}

ttpObjective evaluateFromIndices(ttpInstance *instance, ttpSolution *sol, int startIndex, int stopIndex, double startWeight){
    double weight = startWeight;
    double time = 0;
    double distance = 0;
    double profit = 0;
    
    for (int i=startIndex;i<stopIndex+1;i++){
        if (i>0){
            for (int j=0;j<instance->itemsPerCity;j++){
                int itemIndex = instance->itemIndex[sol->tour[i]][j];
                int packingIdx = sol->isCityPacking ? itemIndex : instance->packingIndex[i][j];

                if (sol->packing[packingIdx]) {
                    weight+=instance->weight[itemIndex];
                    profit+=instance->profit[itemIndex];
                    // printf("tour idx %d, item number %d, packing index %d item idx %d weight %lf profit %lf packed\n",i,j,packingIdx, itemIdx, instance->weight[itemIdx], instance->profit[itemIdx]);
                }
                
            }
        }
        double dist = getDistance(instance, sol->tour[i], sol->tour[(i+1)%instance->numberOfNodes]);
        distance+=dist;
        double speed = getSpeed(instance, weight);
        // speed = weight>instance->capacityOfKnapsack ? instance->minSpeed/1E20 : speed;
        // time+=dist/speed;
        time += weight>instance->capacityOfKnapsack ? (1E20*dist)/instance->minSpeed : dist/speed;
        // double gotSpeed = instance->maxSpeed-weight*(instance->maxSpeed-instance->minSpeed)/instance->capacityOfKnapsack;
        // if (weight>instance->capacityOfKnapsack) printf("wc-diff %lf weight %lf capacity %lf maxSpeed %lf minSpeed %lf speed %lf gotSpeed %lf adjusted %lf time %lf\n", instance->capacityOfKnapsack-weight, weight, instance->capacityOfKnapsack, instance->maxSpeed, instance->minSpeed, getSpeed(instance, weight), gotSpeed, speed, dist/speed);
    }
    ttpObjective ob;
    ob.profit = profit;
    ob.time = time;
    ob.finalweight = weight;
    ob.distance=distance;
    ob.objective=profit-time*instance->rentingRatio;
    
    // if (weight>instance->capacityOfKnapsack) printf("wc-diff %lf objective %lf profit %lf time %lf\n", instance->capacityOfKnapsack-weight,ob.objective, profit, time);
    return ob;
}

ttpObjective evaluate_intermediate(ttpInstance *instance, ttpSolution *sol, ttpObjective *intermediate, int startIndex, bool update){
    double weight = 0;
    double time = 0;
    double distance = 0;
    double profit = 0;
    if (startIndex>0){
      weight = intermediate[startIndex-1].finalweight;
      time = intermediate[startIndex-1].time;
      distance = intermediate[startIndex-1].distance;
      profit = intermediate[startIndex-1].profit;
    }
    
    for (int i=startIndex;i<instance->numberOfNodes;i++){
        if (i>0){
            for (int j=0;j<instance->itemsPerCity;j++){
                // int packingIdx = sol->isCityPacking ? getItemIndex(instance, sol->tour[i],j) : getPackingIndex(instance, i,j);
                int itemIndex = instance->itemIndex[sol->tour[i]][j];
                int packingIdx = sol->isCityPacking ? itemIndex : instance->packingIndex[i][j];

                if (sol->packing[packingIdx]) {
                    // int itemIdx = getItemIndex(instance, sol->tour[i],j);
                    weight+=instance->weight[itemIndex];
                    profit+=instance->profit[itemIndex];
                    // printf("tour idx %d, item number %d, packing index %d item idx %d weight %lf profit %lf packed\n",i,j,packingIdx, itemIdx, instance->weight[itemIdx], instance->profit[itemIdx]);
                }
            }
        }
        double dist = getDistance(instance, sol->tour[i], sol->tour[(i+1)%instance->numberOfNodes]);
        distance+=dist;
        double speed = getSpeed(instance, weight);
        // speed = weight>instance->capacityOfKnapsack ? instance->minSpeed/1E20 : speed;
        // time+=dist/speed;
        time += weight>instance->capacityOfKnapsack ? (1E20*dist)/instance->minSpeed : dist/speed;
        if (update){
          intermediate[i].profit = profit;
          intermediate[i].time = time;
          intermediate[i].distance = distance;
          intermediate[i].finalweight = weight;
          intermediate[i].objective = profit-time*instance->rentingRatio;
        }
        // printf("eval_intermediate tour %d weight %lf distance %lf time %lf\n",i,weight,dist,time);
        
        // double gotSpeed = instance->maxSpeed-weight*(instance->maxSpeed-instance->minSpeed)/instance->capacityOfKnapsack;
        // if (weight>instance->capacityOfKnapsack) printf("wc-diff %lf weight %lf capacity %lf maxSpeed %lf minSpeed %lf speed %lf gotSpeed %lf adjusted %lf time %lf\n", instance->capacityOfKnapsack-weight, weight, instance->capacityOfKnapsack, instance->maxSpeed, instance->minSpeed, getSpeed(instance, weight), gotSpeed, speed, dist/speed);
    }
    ttpObjective ob;
    ob.profit = profit;
    ob.time = time;
    ob.finalweight = weight;
    ob.distance=distance;
    ob.objective=profit-time*instance->rentingRatio;
    // printf("eval_intermediate profit %lf\n",ob.profit);
    
    // if (weight>instance->capacityOfKnapsack) printf("wc-diff %lf objective %lf profit %lf time %lf\n", instance->capacityOfKnapsack-weight,ob.objective, profit, time);
    return ob;
}

ttpObjective evaluate(ttpInstance *instance, ttpSolution *sol){
    return evaluateFromIndices(instance, sol, 0, instance->numberOfNodes-1, 0);
}

ttpObjective evaluateFromIndicesVerbose(ttpInstance *instance, ttpSolution *sol, int startIndex, int stopIndex, double startWeight){
    double weight = startWeight;
    double time = 0;
    double distance = 0;
    double profit = 0;
    int *itemIndices = CC_SAFE_MALLOC(instance->numberOfItems,int);
    int ctr=0;
    
    for (int i=startIndex;i<stopIndex+1;i++){
        if (i>0){
            for (int j=0;j<instance->itemsPerCity;j++){
                // int packingIdx = sol->isCityPacking ? getItemIndex(instance, sol->tour[i],j) : getPackingIndex(instance, i,j);
                // int itemIdx = getItemIndex(instance, sol->tour[i],j);
                int itemIndex = instance->itemIndex[sol->tour[i]][j];
                int packingIdx = sol->isCityPacking ? itemIndex : instance->packingIndex[i][j];

                if (sol->packing[packingIdx]) {
                    
                    weight+=instance->weight[itemIndex];
                    profit+=instance->profit[itemIndex];
                    itemIndices[ctr++]=itemIndex;
                    // printf("eval city %d, item number %d, packing index %d item idx %d weight %lf profit %lf packed %d\n",sol->tour[i],j,packingIdx, itemIdx, weight, profit, sol->packing[packingIdx]);
                }
                
                
            }
        }
        double dist = getDistance(instance, sol->tour[i], sol->tour[(i+1)%instance->numberOfNodes]);
        distance+=dist;
        double speed = getSpeed(instance, weight);
        // speed = speed<instance->minSpeed ? instance->minSpeed/1E20 : speed;
        // time+=dist/speed;
        time += weight>instance->capacityOfKnapsack ? (1E20*dist)/instance->minSpeed : dist/speed;
        // if (weight>instance->capacityOfKnapsack) printf("wc-diff %lf\n", instance->capacityOfKnapsack-weight);
    }

    // qsort(itemIndices,ctr,sizeof(int),compareInt);
    // printf("evaluate verbose items:\n");
    // printIntArray(itemIndices,ctr);
    // printf("\n");
    ttpObjective ob;
    ob.profit = profit;
    ob.time = time;
    ob.finalweight = weight;
    ob.distance=distance;
    ob.objective=profit-time*instance->rentingRatio;
    return ob;
}

ttpObjective evaluateVerbose(ttpInstance *instance, ttpSolution *sol){
    return evaluateFromIndicesVerbose(instance, sol, 0, instance->numberOfNodes-1, 0);
}

/**
 * Calculate the distance from each city in the tour to the end of the tour
 */
double * reverseDistances(ttpInstance *instance, int *tour){
    double * rDist = CC_SAFE_MALLOC(instance->numberOfNodes, double);
    double dist = 0;
    for (int i=instance->numberOfNodes-1;i>=0;i--){
        dist+=getDistance(instance,tour[i],tour[(i+1)%instance->numberOfNodes]);
        rDist[tour[i]]=dist;
    }
    return rDist;
}


/**
 * Comparator for sorting doubleArrayEntries in descending order
 */
int comparingDouble(const void *a, const void *b){
    doubleArrayEntry *a_ = (doubleArrayEntry *) a;
    doubleArrayEntry *b_ = (doubleArrayEntry *) b;
    if (a_->value > b_->value) return -1;
    else if (a_->value < b_->value) return 1;
    else return 0;
}

/**
 * Comparator for sorting doubleArrayEntries in ascending order
 */
int comparingDoubleAsc(const void *a, const void *b){
    doubleArrayEntry *a_ = (doubleArrayEntry *) a;
    doubleArrayEntry *b_ = (doubleArrayEntry *) b;
    if (a_->value > b_->value) return 1;
    else if (a_->value < b_->value) return -1;
    else return 0;
}

void toCityPacking(ttpInstance *instance, ttpSolution *sol){
    bool *packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    for (int i=0;i<instance->numberOfNodes;i++){
        for (int j=0;j<instance->itemsPerCity;j++){
            int itemIndex = getItemIndex(instance, sol->tour[i], j);
            int packingIndex = getPackingIndex(instance,i,j);
            packing[itemIndex] = sol->packing[packingIndex];
        }
    }
    for (int i=0;i<instance->numberOfItems;i++) sol->packing[i] = packing[i];
    sol->isCityPacking = true;
    CC_IFFREE(packing, bool);
}

void toIndexPacking(ttpInstance *instance, ttpSolution *sol){
    bool *packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    for (int i=0;i<instance->numberOfNodes;i++){
        for (int j=0;j<instance->itemsPerCity;j++){
            int itemIndex = getItemIndex(instance, sol->tour[i], j);
            int packingIndex = getPackingIndex(instance,i,j);
            packing[packingIndex] = sol->packing[itemIndex];
        }
    }
    for (int i=0;i<instance->numberOfItems;i++) sol->packing[i] = packing[i];
    sol->isCityPacking = false;
    CC_IFFREE(packing, bool);
}

void updateWithIndexPacking(ttpInstance *instance, ttpSolution *sol, bool *packing){
    // sol->packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    for (int i=1;i<instance->numberOfNodes;i++){
        for (int j=0;j<instance->itemsPerCity;j++){
            int itemIndex = getItemIndex(instance, sol->tour[i],j);
            int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, i,j);
            sol->packing[packingIndex]= packing[itemIndex];
        }
    }
}

/**
 * Guess at whether this is a "bounded strongly corr", "uncorr similar weights", or "uncorr" problem
 */
int knapsackType(ttpInstance *instance){
    int minWeight = CCutil_MAXINT;
    int maxWeight = -(CCutil_MAXINT-1);
    bool corr = true;
    for (int i=0;i<instance->numberOfItems;i++){
        int weight = (int)instance->weight[i];
        minWeight = minWeight>weight ? weight : minWeight;
        maxWeight = maxWeight<weight ? weight : maxWeight;
        if ((int)(instance->profit[i]-instance->weight[i])%100 != 0) corr=false;
    }

    if (corr) {
        return 0;
    } else if (maxWeight-minWeight<=10){
        return 1;
    } else {
        return 2;
    }
}

/**
 * Get the capacity factor of the knapsack
 */
double capacityFactor(ttpInstance *instance){
    double tw = 0;
    for (int i=0;i<instance->numberOfItems;i++) tw+=instance->weight[i];
    return 11*instance->capacityOfKnapsack/tw;
    // return tw/instance->capacityOfKnapsack;
}

// 0 = linear     = 2 weights = x0 x1
// 1 = quadratic  = 3 weights = x0 x1 x0*x1
// 2 = 5-term 1   = 4 weights = x0 x1 x0*x1 x0**2
// 3 = 5-term 2   = 4 weights = x0 x1 x0*x1 x1**2
// 4 = 6-term     = 5 weights = x0 x1 x0*x1 x0**2 x1**2
#define HEURISTIC_ALG 1

#if HEURISTIC_ALG==0
double estimateParameter(ttpInstance *instance, int parameterIdx){
  printf("heuristic alg %d\n",HEURISTIC_ALG);
  double x0 = instance->capacityFactor;
  double ans;
  if (instance->knapsackType==0){
    if (parameterIdx==0){
      ans = x0/(0.48412*x0 + 0.9455312) + 0.603/(0.61676888*(x0*x0*x0) + 2.1574710976*(x0*x0) - 1.198778970912*x0 - 5.97609487371429) - 3.89285714285714/(0.48412*x0 + 0.9455312);
      return ans*(0.09789034114334358)+(0.4710828569975715);
    } else if (parameterIdx==1){
      ans = 2*(x0*x0*x0)/(0.723*(x0*x0*x0) + 2.26244343891403*(x0*x0)) - 8.91188687782805*(x0*x0)/(0.723*(x0*x0*x0) + 2.26244343891403*(x0*x0)) + 3.02356561085973*x0/(0.723*(x0*x0*x0) + 2.26244343891403*(x0*x0)) + 0.364554/(0.723*(x0*x0*x0*x0) + 2.26244343891403*(x0*x0*x0)) - 2.149869/(0.723*(x0*x0*x0) + 2.26244343891403*(x0*x0));
      return ans*(0.04692066043768883)+(-0.8753838441380654);
    } else{
      ans = -0.33951*x0 + 1.82307505 + 0.628/(7.69230769230769*x0 - 12.9076923076923) + 0.628/(7.69230769230769*x0 - 20.3153846153846);
      return ans*(0.06758875920266104)+(0.9131949540503761);
    }
  } else if (instance->knapsackType==1){
    if (parameterIdx==0){
      ans = 0.181*x0 - 0.423 - 1.419/((x0*x0*x0) - 5.4475138121547*(x0*x0)) - 2.255/x0;
      return ans*(0.06565900720164893)+(0.8431385263676947);
    } else if (parameterIdx==1){
      ans = -0.27*x0/(-0.959 + 0.404/x0) + 0.06426/(-0.959*x0 + 5.84713646153846 - 2.29304184615385/x0) + 1.38937846153846/(-0.959 + 0.404/x0);
      return ans*(0.0943657486178776)+(-0.5252631965173206);
    } else{
      ans = -0.267573696145125*x0 + 1.091 + 1.54969574036511/(x0 + 0.493);
      return ans*(0.08384775211462533)+(0.8421597539441523);
    }
  } else{
    if (parameterIdx==0){
      ans = x0/(0.519841*x0 + x0/(0.721*x0 + 0.177) - 0.397/(0.721*x0 + 0.177)) - 4.08598473282443/(0.519841*x0 + x0/(0.721*x0 + 0.177) - 0.397/(0.721*x0 + 0.177));
      return ans*(0.10604172693263764)+(0.8600932740123025);
    } else if (parameterIdx==1){
      ans = 0.183*x0 - 0.485189506 + 0.042761/(x0 - 4.3605652173913) - 1.90203715442357/x0;
      return ans*(0.15265547561062018)+(-0.47506948745845323);
    } else{
      ans = -0.205*x0 - 0.205*x0/(0.205*x0 + 0.894) + 1 + 1/(0.205*(x0*x0) + 1.20273*x0 + 1.346364) + 0.964091403699673/(0.205*x0 + 0.894);
      return ans*(0.0724182697392003)+(0.7791455900565144);
    }
  }
}
#elif HEURISTIC_ALG==1
/**
 * Estimate a parameter of the packing heuristic depending on the guessed knapsack type
 * @param parameterIdx represents the coefficient of x0, x1, or x0*x1, or the percent capacity used, if `parameterIdx=3`
 */
double estimateParameter(ttpInstance *instance, int parameterIdx){
  // printf("heuristic alg %d\n",HEURISTIC_ALG);
  double x0 = instance->capacityFactor;
  double ans;
  if (instance->knapsackType==0){
    if (parameterIdx==0){
      ans = 0.102*x0 + 0.0740359195402299*x0/(x0 - 0.696/(0.339 - 3.06607929515418/x0)) + 0.339 - 0.696/(0.227*x0 + 0.040179);
      return ans*(0.07794943796974876)+(0.4683691144892353);
    } else if (parameterIdx==1){
      ans = 0.17595*x0 - 0.315054 + 0.043020396/(0.414*x0 - 4.392) - 1.76328502415459/x0;
      return ans*(0.03793404146167837)+(-0.8774225327928984);
    } else if (parameterIdx==2){
      ans = 0.148/(0.028836*(x0*x0*x0*x0) - 0.351116*(x0*x0*x0) + 1.003793*(x0*x0)) + 0.258426966292135/(-0.178*(x0*x0*x0) + 1.351*(x0*x0)) + 0.239/(-0.178*(x0*x0) + 1.351*x0);
      return ans*(0.02550931072519169)+(0.050933684578813775);
    } else{
      ans = -0.233643044619423*(x0*x0)/(0.742976*x0 + 1.06047032474804/x0) + x0/(0.742976*x0 + 1.06047032474804/x0) + 1.584/(0.742976*x0 + 1.06047032474804/x0);
      return ans*(0.06703610589783689)+(0.9128793626171633);
    }
  } else if (instance->knapsackType==1){
    if (parameterIdx==0){
      ans = 0.365*x0 - 2.065 + 0.059/(-0.231775*(x0*x0) + 2.26619*x0) + 0.059/(x0 - 0.635);
      return ans*(0.08425172516917678)+(0.8204902478209194);
    } else if (parameterIdx==1){
      ans = 0.44*(x0*x0)/(0.787*x0 + 3.462) - 1.44002*x0/(0.787*x0 + 3.462) - 0.785 - 0.654481/(0.787*x0 + 3.462);
      return ans*(0.1435978222840473)+(-0.4916949349140941);
    } else if (parameterIdx==2){
      ans = -0.519509476031215 - 0.098/(1.68155497493115 + 2.82024574535697/(-0.098*(x0*x0) + 0.639*x0)) + 2.34207395332715/x0;
      return ans*(0.02172513230075182)+(0.23841387350236692);
    } else{
      ans = -0.243*x0 + 0.722454545454545 - 1.06461538461538/(-0.325*x0 - 0.267) - 0.616/(x0*x0*x0);
      return ans*(0.08759401224063153)+(0.8526498649711435);
    }
  } else{
    if (parameterIdx==0){
      ans = 1.46859903381643*x0/(x0 + 0.212 + 0.959/(0.696*x0 + 0.092416)) - 2.08695652173913/((x0*x0) + 0.212*x0 + 0.959*x0/(0.696*x0 + 0.092416)) - 5.27737198067633/(x0 + 0.212 + 0.959/(0.696*x0 + 0.092416));
      return ans*(0.12258165126028134)+(0.8589358651893645);
    } else if (parameterIdx==1){
      ans = 0.162052*x0 - 0.0639699 - 2.12612612612613/(x0 - 1 + 1.85810810810811/x0) - 0.944/x0;
      return ans*(0.20679174682604445)+(-0.3946503934450688);
    } else if (parameterIdx==2){
      ans = -0.497041420118343*(x0*x0)/(2*x0 - 0.489 - 0.569266589057043/(x0*x0)) + 0.994082840236686*x0/(2*x0 - 0.489 - 0.569266589057043/(x0*x0)) + 1.7636684303351 - 1.63844797178131/x0 - 2.5515794746564/(x0*x0);
      return ans*(0.02765484479999771)+(0.21892091773760267);
    } else{
      ans = -0.218026796589525*x0 + 0.544 + 3.10955135310566/x0 - 1.49413942516727/(x0*x0);
      return ans*(0.0835233858685778)+(0.7929762706832296);
    }
  }
}
#elif HEURISTIC_ALG==2
double estimateParameter(ttpInstance *instance, int parameterIdx){
  printf("heuristic alg %d\n",HEURISTIC_ALG);
  double x0 = instance->capacityFactor;
  double ans;
  if (instance->knapsackType==0){
    if (parameterIdx==0){
      ans = x0/(0.453*x0 + 0.889 + 0.732891832229581/x0) - 0.332/(0.293091*(x0*x0) + 0.477876788*x0 + 0.283220259452539 - 0.157428097130243/x0) - 4.02579591836735/(0.453*x0 + 0.889 + 0.732891832229581/x0);
      return ans*(0.2596415551307137)+(0.37479463296938265);
    } else if (parameterIdx==1){
      ans = 0.283*x0 - 1.633 + 0.06/(0.160178*x0 - 0.572872857142857 + 0.06/(0.283*x0 - 1.108));
      return ans*(0.08618514053284809)+(-0.8557279965342615);
    } else if (parameterIdx==2){
      ans = 0.1390544064*x0 - 0.9971829504 - 0.168137424/(4.06808510638298 - 0.725*x0) - 0.173680416/(0.831979310344828 - x0);
      return ans*(0.10722236993088358)+(-0.15562917162107853);
    } else if (parameterIdx==3){
      ans = -0.215*x0 + 1.93121 + 0.55598816/(-0.215*(x0*x0) + 1.12367*x0) - 0.07826/(1.12367 - 0.215*x0) - 2.86853232/x0;
      return ans*(0.05610011416635705)+(0.11641591370128762);
    } else{
      ans = -0.0479582089552239*(x0*x0)/(0.138*x0 + 0.468*x0/((x0*x0*x0) - 3.12975213675214*(x0*x0) + 2.12179487179487*x0)) + 0.0861009137436694*x0/(0.138*x0 + 0.468*x0/((x0*x0*x0) - 3.12975213675214*(x0*x0) + 2.12179487179487*x0)) + 1.278;
      return ans*(0.06814587458310797)+(0.9137136385807437);
    }
  } else if (instance->knapsackType==1){
    if (parameterIdx==0){
      ans = -0.02901767253*(x0*x0) + 0.757034508391437*x0 - 3.35464584185848 + 1.0891794448916/x0;
      return ans*(0.09674471127360099)+(0.8128563184477983);
    } else if (parameterIdx==1){
      ans = 0.4418*(x0*x0*x0)/((x0*x0) + 1.684*x0 + 0.706048) - 0.7001496*(x0*x0)/((x0*x0) + 1.684*x0 + 0.706048) - 1.441 + 0.854545454545454/((x0*x0) + 1.684*x0 + 0.706048);
      return ans*(0.1568855521867961)+(-0.506800298669995);
    } else if (parameterIdx==2){
      ans = -0.054*x0 - 0.787 + 0.046548/(-0.054*x0 + 0.417306 + 1.41080196399345/(x0 - 0.359)) + 2.42823057485964/x0;
      return ans*(0.047860009332280795)+(0.19949315805182893);
    } else if (parameterIdx==3){
      ans = -x0/(154.022771411071 - 14.7540402584409*x0) + 0.541 - 0.99/(x0*x0);
      return ans*(0.05476220126316566)+(0.0579521172448557);
    } else{
      ans = -0.261501423*x0 + x0/(x0 - 0.915 + 0.929/x0) + 0.929/(x0 - 0.915 + 0.929/x0);
      return ans*(0.08816427406641573)+(0.8536049974670996);
    }
  } else{
    if (parameterIdx==0){
      ans = 0.533*(x0*x0)/(0.472*(x0*x0) + 0.978) - 1.195611*x0/(0.472*(x0*x0) + 0.978) - 2.873/(0.472*(x0*x0) + 0.978);
      return ans*(0.16428982836165495)+(0.8426286885264234);
    } else if (parameterIdx==1){
      ans = 0.626*(x0*x0)/(0.215070643642072*(x0*x0) + 0.747252747252747*x0) - 3*x0/(0.215070643642072*(x0*x0) + 0.747252747252747*x0) + 0.391/(0.215070643642072*(x0*x0) + 0.747252747252747*x0);
      return ans*(0.23458009288248102)+(-0.4139663509604798);
    } else if (parameterIdx==2){
    //   ans = 0.03220184/(0.154*x0 - 1.005084 + 0.502866/x0) - 0.687843528/((x0*x0) - 0.546*x0) + 0.154/(x0 - 0.546);
    //   return ans*(0.03519900157559365)+(0.13498733399364055);
      return x0*(0.0030733071098402462)+(0.1180841448895192);
    } else if (parameterIdx==3){
      ans = 0.015876*(x0*x0*x0) - 0.331075584*(x0*x0) + 1.977111703392*x0 - 2.99108179633536 + 0.418499125622964/x0;
      return ans*(0.04569998202833745)+(0.12264072511539278);
    } else{
      ans = -0.185015290519878*x0 - x0/(-0.164179104477612*(x0*x0) - 0.241179104477612*x0 - 0.0833419104477612) + 0.0636452599388379/(x0 + 0.556);
      return ans*(0.08498937164732476)+(0.7935947690677406);
    }
  }
}
#elif HEURISTIC_ALG==3
double estimateParameter(ttpInstance *instance, int parameterIdx){
  printf("heuristic alg %d\n",HEURISTIC_ALG);
  double x0 = instance->capacityFactor;
  double ans;
  if (instance->knapsackType==0){
    if (parameterIdx==0){
      ans = 0.23*x0 - 0.924 + 0.06164/(x0 - 1.573 - 1.13859649122807/x0) - 0.866/x0;
      return ans*(0.10746598326463856)+(0.46156448817338547);
    } else if (parameterIdx==1){
      ans = 0.276*x0 - 1.524 + 0.02390022/(0.086595*(x0*x0) - 1.474546655*x0 + 6.195359398);
      return ans*(0.08125148032741815)+(-0.8430699222900591);
    } else if (parameterIdx==2){
      ans = 0.131588235294118*x0 - 0.625 + 0.028/(0.101*x0 - 0.952) + 0.028/(0.965699208443272 - 0.133245382585752*x0);
      return ans*(0.08433175199629035)+(0.1251727208564827);
    } else if (parameterIdx==3){
      ans = 0.24/(1.272 - x0) - 0.113/(0.902 + 0.24/(1.379 - 0.386*x0)) - 0.113/(0.902 + 0.24/(0.794 - 0.113*x0));
      return ans*(0.1476110838931197)+(-0.11606680077938529);
    } else{
      ans = -0.239*x0 + x0/(0.130494*(x0*x0) + 0.0703541750434782*x0 + 0.777578091365739 - 0.0609400434782609/x0) + 0.184986 + 0.239/(0.130494*(x0*x0) + 0.0703541750434782*x0 + 0.777578091365739 - 0.0609400434782609/x0);
      return ans*(0.06785822983792258)+(0.9133764525676066);
    }
  } else if (instance->knapsackType==1){
    if (parameterIdx==0){
      ans = 0.277727678571429*x0 + 1.4691787732042*x0/(0.923832923832924*x0 + 1.91703685503686) - 2.60969752495176;
      return ans*(0.10801354356780439)+(0.7937469385392443);
    } else if (parameterIdx==1){
      ans = 0.0429677419354839*x0/(0.0594516129032258 + 0.931/x0) - 1.15416129032258 - 0.009246024/(0.0594516129032258*x0 + 0.931) + 0.0188864477991419/(0.0594516129032258 + 0.931/x0) - 0.50731097282695/x0 + 0.248358478/(x0*x0);
      return ans*(0.16695615791852964)+(-0.5009050590720251);
    } else if (parameterIdx==2){
      ans = -0.105*x0 + 0.932 - 0.105/(0.512*x0 - 0.105*x0/(-0.046*(x0*x0) + 0.336*x0) + 0.512) - 0.426/x0;
      return ans*(0.04874444874179613)+(0.2718743396761879);
    } else if (parameterIdx==3){
      ans = 0.190552197596796*x0 - 0.948 + 0.132573/(0.177*(x0*x0) - 0.948*x0) + 0.132573/(x0 - 0.836) - 1.0/x0;
      return ans*(0.03334333982678032)+(-0.04599107442024427);
    } else{
      ans = -0.266888150609081*x0 + 1.06803875968992 + 1.664/(x0 + 0.954/(x0*x0*x0));
      return ans*(0.08868574090396958)+(0.8540593760272128);
    }
  } else{
    if (parameterIdx==0){
      ans = 0.469*(x0*x0)/(0.011236*(x0*x0*x0) + 0.154522984*(x0*x0) + 0.482730316*x0 - 0.052919016) - 1.892*x0/(0.011236*(x0*x0*x0) + 0.154522984*(x0*x0) + 0.482730316*x0 - 0.052919016);
      return ans*(0.1515889884521732)+(0.8259452772287567);
    } else if (parameterIdx==1){
      ans = -0.174*x0/(-0.004319849815424*(x0*x0) - 0.004319849815424*x0 - 0.326772) + 0.816/(-0.004319849815424*(x0*x0) - 0.004319849815424*x0 - 0.326772);
      return ans*(0.23606847923579524)+(-0.3933847198315476);
    } else if (parameterIdx==2){
      ans = -0.09911673*(x0*x0)/(0.030237384321*(x0*x0) - 0.012608989261857*x0 + 0.590889) + 1.0208145133*x0/(0.030237384321*(x0*x0) - 0.012608989261857*x0 + 0.590889) - 1.73621103117506/(0.030237384321*(x0*x0) - 0.012608989261857*x0 + 0.590889);
      return ans*(0.05840882545807337)+(0.27532713129446434);
    } else if (parameterIdx==3){
      ans = 0.22*x0 - 0.962 + 0.0484/(0.368*x0 - 1.213264) + 0.0484/(0.368*x0 - 1.968) - 0.984/x0;
      return ans*(0.04070329187034253)+(-0.05900792176140322);
    } else{
      ans = -0.237426*x0 + 0.934 + 0.706/(0.706*x0 - 1 + 0.876/x0);
      return ans*(0.08528958073697082)+(0.7939984526614366);
    }
  }
}
#elif HEURISTIC_ALG==4
double estimateParameter(ttpInstance *instance, int parameterIdx){
  printf("heuristic alg %d\n",HEURISTIC_ALG);
  double x0 = instance->capacityFactor;
  double ans;
  if (instance->knapsackType==0){
    if (parameterIdx==0){
      ans = -2*x0/(-0.763*x0 - 1.78 - 1.04575163398693/(0.193*x0 + 0.83)) + 8.91457142857143/(-0.763*x0 - 1.78 - 1.04575163398693/(0.193*x0 + 0.83));
      return ans*(0.29103383224801765)+(0.37867140994592735);
    } else if (parameterIdx==1){
      ans = 0.064*(x0*x0) - 0.443*x0 + 0.167/(-0.004096*(x0*x0*x0*x0*x0) - 0.103744*(x0*x0*x0*x0) + 0.758*(x0*x0*x0));
      return ans*(0.11631859636972328)+(-0.7757845539741349);
    } else if (parameterIdx==2){
      ans = 0.02768*(x0*x0) - 1.19724 - 0.181/(0.0256*(x0*x0*x0) - 1.179*x0 - 0.907);
      return ans*(0.2461797683188084)+(-0.12363813386747632);
    } else if (parameterIdx==3){
      ans = -0.236*x0 + x0/(x0 - 0.799*x0/(x0 - 0.236) + 2.70508898305085/((x0*x0) - 0.236*x0)) + 0.232224 - 0.044501104/(1.117 - 0.236*x0);
      return ans*(0.06253793260732617)+(0.12310469548874733);
    } else if (parameterIdx==4){
      ans = -0.292*x0 + 1.531 - 0.0452722063037249/(0.085264*(x0*x0) - 0.779348*x0 + 1.435734) + 0.079/(-0.235*(x0*x0) + 0.747*x0 + 0.961);
      return ans*(0.24782328616569343)+(-0.012750491067133957);
    } else{
      ans = -0.063*(x0*x0)/(0.21*x0 + 1.0/x0 + 89.4202886226221/(x0*x0*x0*x0)) + 1.329;
      return ans*(0.06777910989800047)+(0.9132848974001594);
    }
  } else if (instance->knapsackType==1){
    if (parameterIdx==0){
      ans = 0.22*x0 - 0.241*x0/(0.891/(0.22*x0 - 1.714) + 0.891/(-x0 - 0.089)) - 0.768 + 1.714/(0.891/(0.22*x0 - 1.714) + 0.891/(-x0 - 0.089));
      return ans*(0.09792129863021484)+(0.798516942911712);
    } else if (parameterIdx==1){
      ans = -0.484*x0/(0.081*x0 - 2.168) - 1.757 + 0.242/x0;
      return ans*(0.16756127683310157)+(-0.49844188298096753);
    } else if (parameterIdx==2){
      ans = 0.00145432*(x0*x0*x0) + 0.285/(-0.14*x0 + 0.577 - 0.061/(0.14*x0 - 0.447));
      return ans*(0.10269662315858692)+(0.23813638942263954);
    } else if (parameterIdx==3){
      ans = -0.00322161*(x0*x0*x0) + 0.001277535*(x0*x0*x0)/(0.161*x0 + 0.143*x0/(0.133*(x0*x0) + 0.035*x0 - 0.969) - 0.969);
      return ans*(0.07739120550999515)+(0.02056441159033432);
    } else if (parameterIdx==4){
      ans = -0.391 + 0.569/(0.218*x0 - 0.439 - 0.782/(-0.113590583*(x0*x0*x0) + 0.09541608972*(x0*x0) + x0 - 0.439) - 0.359/(x0 - 0.84));
      return ans*(0.04521965262554614)+(-0.023417341450388644);
    } else{
      ans = -0.238*x0 - 0.238*x0/(x0 - 0.696 + 2.03719951975684/x0) + 1 + 2.28114951975684/(x0 - 0.696 + 2.03719951975684/x0);
      return ans*(0.08878255613565315)+(0.8538683465716045);
    }
  } else{
    if (parameterIdx==0){
      ans = -1.08108108108108*x0/(-x0 + 1.17785630153121 - 2.22025912838634/x0 + 0.85962308598351/(x0*x0)) + 3.835/(-x0 + 1.17785630153121 - 2.22025912838634/x0 + 0.85962308598351/(x0*x0));
      return ans*(0.17908275744664975)+(0.8305256659917525);
    } else if (parameterIdx==1){
      ans = 0.630819672131148*x0/(0.261*x0 + 0.962/(0.171147540983607*x0 + 0.962)) - 2.91020125786163/(0.261*x0 + 0.962/(0.171147540983607*x0 + 0.962));
      return ans*(0.23136053304135645)+(-0.4100868749181984);
    } else if (parameterIdx==2){
      ans = 0.078*x0 + 0.291762533976986 + 0.078/(x0 - 1.87) + 0.078/(0.078*x0 - 0.827) - 0.935/x0;
      return ans*(0.10927830092593334)+(0.14408997410979613);
    } else if (parameterIdx==3){
    //   ans = -0.045*x0/(0.464*x0 - 2.78212290502793) + 0.01791*x0/(0.398*x0 - 3.78112290502793) - 0.125195530726257 - 0.06271209/(0.398*x0 - 3.78112290502793);
    //   return ans*(0.07986523902358482)+(0.11293204130726894);
      return x0*(-0.005761335678416693)+(0.14461938753856074);
    } else if (parameterIdx==4){
      ans = -0.076634628 + 0.018424/(0.031329*(x0*x0) - 0.257712*x0 + 0.502/(2.149 - x0)) + 0.330316/(1.43 - x0);
      return ans*(0.06791409471976426)+(0.0030738531127616);
    } else{
      ans = -0.176258992805755*x0 - 0.096589928057554 + 1.95128205128205/(0.309136*x0 + 0.423116 + 0.160441584/x0);
      return ans*(0.0860822236295349)+(0.7943625945655375);
    }
  }
}
#endif



double stdvsmean(ttpInstance *instance, int parameterIdx){
    if (instance->knapsackType==0){
        if (parameterIdx==0){
            return  0.10825761723848046/ 0.383766305969535;
        } else if (parameterIdx==1){
            return 0.121506613869339/0.7127071210041852;
        } else if (parameterIdx==2){
            return 0.01986299921736277 / 0.040154389497925554;
        } else {
            return 0.5545773627761429/0.09341010721905452;
        }
    //1: uncorr-similar-weights
    } else if (instance->knapsackType==1){
        if (parameterIdx==0){
            return 0.1657006227076693/ 0.7081549976053036;
        } else if (parameterIdx==1){
            return 0.08130687702799178/0.4027274749172905;
        } else if (parameterIdx==2){
            return 0.02792613718098652 /0.20128471208696877;
        } else {
            return 0.4135362857115528/0.2999198420207125;
        }
    //2: uncorr
    } else {
        if (parameterIdx==0){
            return 0.19431251387348641/0.7658850751111516;
        } else if (parameterIdx==1){
            return 0.11922766938497528/0.3224293423549689;
        } else if (parameterIdx==2){
            return 0.043320868882208395/0.19324114754374214;
        } else {
            return 0.4530550864034252/0.11405713750122952;
        }
    }
}


int lk_tour(CCdatagroup *dat, CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol, int *outIters, int **outtours, int numIters, int numOutputs, double timeBound){

  // CCrandstate rstate2;
  // int seed = (int) CCutil_real_zeit ();
  // CCutil_sprand (seed, &rstate2);

  int local_quadtry = 12;//2
  int run_silently = 1;//0
  int local_in_repeater = numIters;
  int local_number_runs = 0;
  double local_length_bound = -1.0;
  int local_kick_type = CC_LK_WALK_KICK;
  int *incycle = CC_SAFE_MALLOC (instance->numberOfNodes, int);
  // int *incycle_copy = CC_SAFE_MALLOC (instance->numberOfNodes, int);

  // local_time_bound = 30;
  // local_in_repeater = rand()%(instance->numberOfNodes + 1);
  // local_in_repeater = instance->numberOfNodes>MAXLKITERS ? MAXLKITERS : instance->numberOfNodes;

  double val;
  int tempcount, *templist;

  CCkdtree localkt;
  CCkdtree_build (&localkt, instance->numberOfNodes, dat, (double *) NULL, rstate);
  CCkdtree_qboruvka_tour (&localkt, instance->numberOfNodes, dat, incycle, &val, rstate);

  // copyInt(incycle,incycle_copy,instance->numberOfNodes);

  CCkdtree_quadrant_k_nearest (&localkt, instance->numberOfNodes, local_quadtry,
                          dat, (double *) NULL, 1, &tempcount, &templist,
                          run_silently, rstate);

  if (local_in_repeater == -1) local_in_repeater = instance->numberOfNodes;
  double objective = -CCutil_MAXDOUBLE;
  if (numOutputs>0 && outIters[0]==-1){
    for (int i=0;i<instance->numberOfNodes;i++) outtours[0][i]=incycle[i];
    outIters++;
    outtours++;
    numOutputs--;
  }
  CClinkern_tour (instance->numberOfNodes, dat, tempcount, templist, 100000000,
                   numIters, incycle, sol->tour, &val, run_silently,
                   timeBound, local_length_bound, (char *) NULL, local_kick_type,
                   rstate, outIters, outtours, numOutputs);

  // printIntArray(sol->tour,instance->numberOfNodes);
  // printf("\n");
  CCkdtree_free(&localkt);
  
  CC_IFFREE(templist, int);
  CC_IFFREE(incycle, int);
  return 0;
}