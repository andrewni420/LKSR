#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "util.h"
#include "ttp_util.h"
#include "kdtree.h"
#include "ttp_heuristic.h"
#include "ttp_optimisation.h"
#include "linkern.h"
#include "ttp.h"
#include "ttp_linkern.h"
#include "mbfs.h"
// #include "naive_optimisation.h"

char *global_filesuffix = (char *) NULL;
int global_outtime = 10;

int ttp_lk_tour(CCdatagroup *dat, CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol, mbfs_main *main, int num_iters, double time_bound, double *objective, bool silent){
  int local_quadtry = 12;//2
  int run_silently = silent ? 1 : 0;
  int local_in_repeater = -1;
  int local_number_runs = 0;
  double local_time_bound = -1.0;
  double local_length_bound = -1.0;
  int local_kick_type = CC_LK_WALK_KICK;
  int *incycle = CC_SAFE_MALLOC (instance->numberOfNodes, int);

  // local_time_bound = 1+(rand()%101)/100.0;
  // local_in_repeater = rand()%(instance->numberOfNodes + 1);


  double val;
  int tempcount, *templist;

  CCkdtree localkt;
  CCkdtree_build (&localkt, instance->numberOfNodes, dat, (double *) NULL, rstate);
  CCkdtree_qboruvka_tour (&localkt, instance->numberOfNodes, dat, incycle, &val, rstate);

  CCkdtree_quadrant_k_nearest (&localkt, instance->numberOfNodes, local_quadtry,
                          dat, (double *) NULL, 1, &tempcount, &templist,
                          run_silently, rstate);

  int tempcount2, *templist2;
  CCkdtree_quadrant_k_nearest (&localkt, instance->numberOfNodes, local_quadtry,
                          dat, (double *) NULL, 1, &tempcount2, &templist2,
                          run_silently, rstate);

  if (local_in_repeater == -1) local_in_repeater = instance->numberOfNodes;


  logisticRegressor classifier;
  initRegressor(&classifier);
  classifier.nWeights = 3;
  classifier.weight = CC_SAFE_MALLOC(3,double);
  classifier.weight[0] = estimateParameter(instance, 0);
  classifier.weight[1] = estimateParameter(instance, 1);
  classifier.weight[2] = estimateParameter(instance, 2);
  classifier.bias = estimateParameter(instance, 3);
  double capacity = estimatePercent(instance, 0)*instance->capacityOfKnapsack;

  CCttp_linkern_tour (instance->numberOfNodes, dat, tempcount2, templist2, num_iters,
                  num_iters, incycle, &val, run_silently,
                  time_bound, local_length_bound, (char *) NULL, local_kick_type,
                  rstate, instance, sol, main, &classifier, objective, &capacity, 0);

  // printIntArray(sol->tour,instance->numberOfNodes);
  // printf("\n");
  freeRegressor(&classifier);
  CCkdtree_free(&localkt);
  CC_IFFREE(templist, int);
  CC_IFFREE(templist2, int);
  CC_IFFREE(incycle, int);
  return 0;
}


int main(int argc, char** argv) {

  double initial_time = CCutil_zeit();
  global_outtime = (int) CCutil_real_zeit();
  global_filesuffix = argv[1];

  if (argc<3){
    printf("Need 2 arguments: ttp instance file, and random seed\n");
    return 1;
  }
  printf("problem %s\n", argv[1]);

  printf("zeit initial %lf",zeit_initial);
  init_params(&ttp_lk_params);
  // Initialization
  CCdatagroup dat;
  CCrandstate rstate;
  ttpInstance instance;
  initInstance(&instance);
  ttpSolution sol;
  initSolution(&sol);

  CCutil_init_datagroup(&dat);

  int seed = atoi(argv[2]);
  // int seed = (int) CCutil_real_zeit ();
  CCutil_sprand (seed, &rstate);
  srand(seed+100);
  srand(rand());
  srand(rand()*10);
  // srand(atoi(argv[2])+10);
  // srand(rand());
  // double time = CCutil_real_zeit();
  // srand(fmod(time,1)*1000000+rand()/2);
  int ncount; 
  CCutil_gettsplib (argv[1], &ncount, &dat);

  initInstance(&instance);
  read_ttp_tsplib(argv[1], &dat, &rstate, &instance);

  set_params(&instance, &ttp_lk_params, 0);


  initSolution(&sol);
  sol.isCityPacking = true;
  sol.packing = CC_SAFE_MALLOC(instance.numberOfItems, bool);
  sol.tour = CC_SAFE_MALLOC(instance.numberOfNodes+1,int);

  double objective = -CCutil_MAXDOUBLE;

  mbfs_main main;
  init_mbfs_main(&main);
  construct_mbfs_main(&main, &instance);
  read_mbfs_main(&main, &instance);

  lk_tour(&dat, &rstate, &instance, &sol, NULL, NULL, rand()%instance.numberOfNodes, 0, 1+(rand()%101)/100.);
  logisticRegressor classifier;
  initRegressor(&classifier);
  classifier.nWeights = 3;
  classifier.weight = CC_SAFE_MALLOC(3,double);
  classifier.weight[0] = estimateParameter(&instance, 0);
  classifier.weight[1] = estimateParameter(&instance, 1);
  classifier.weight[2] = estimateParameter(&instance, 2);
  classifier.bias = estimateParameter(&instance, 3);
  double capacity = estimatePercent(&instance, 0)*instance.capacityOfKnapsack;
  calcRDist(&instance, sol.tour);
  quadraticHillClimber(&instance, &sol, &classifier, ttp_lk_params.upper_quadratic_max_nonimprovement,ttp_lk_params.upper_quadratic_max,&capacity);

  ttpObjective obj = evaluate(&instance, &sol);
  objective = obj.objective;
  printf("Initial objective %lf\n",obj.objective);
  
  ttp_lk_tour(&dat, &rstate, &instance, &sol, &main, 1000000000, 590-(CCutil_zeit()-initial_time),&objective,false);//590-(CCutil_zeit()-initial_time)

  ttpOutput(&instance, &sol, argv[1], "LKSR", NULL);

  
  free_params(&ttp_lk_params);
  freeInstance(&instance);
  free_mbfs_main(&main);
  CCutil_freedatagroup(&dat);
  freeSolution(&sol);
  
  return 0;
}