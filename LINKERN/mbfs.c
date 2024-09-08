#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "util.h"
#include "math_util.h"
#include "ttp_util.h"
#include "mbfs.h"

#define LCIPR_DEFAULT_VALUE 1E6

int *city_to_tour_idx;

void init_mbfs_main(mbfs_main *main){
    main->iprSortedIdx = (mbfs_entry **) NULL;
    main->num_cities = 0;
    main->items_per_city = 0;
    main->ipr = (double *) NULL;
    city_to_tour_idx = (int *) NULL;
}

void init_mbfs(mbfs *mbfs){
    mbfs->packed = (intvector *) NULL;
    mbfs->unpacked = (intvector *) NULL;
    mbfs->huipr = (doubleArrayEntry *) NULL;
    mbfs->lcipr = (doubleArrayEntry *) NULL;
    mbfs->prefix_min_lcipr = (double *) NULL;
    mbfs->suffix_max_huipr = (double *) NULL;
    mbfs->intermediate = (ttpObjective *) NULL;
    intvector_init(&mbfs->to_pack);
    intvector_init(&mbfs->to_unpack);
    intvector_init(&mbfs->to_pack_loc);
    intvector_init(&mbfs->to_unpack_loc);
    mbfs->city_to_packing = (int *) NULL;
    mbfs->city_to_tour_idx = (int *) NULL;
    mbfs->mbfs_main = (mbfs_main *) NULL;
    mbfs->total_weight = 0;
    mbfs->total_profit=0;
}

void  free_mbfs_main(mbfs_main *main){
    for (int i=1;i<main->num_cities;i++){
        CC_IFFREE(main->iprSortedIdx[i], mbfs_entry);
    }
    CC_IFFREE(main->iprSortedIdx, mbfs_entry *);
    CC_IFFREE(city_to_tour_idx, int);
}
void free_mbfs(mbfs *mbfs){
    CC_IFFREE(mbfs->huipr, doubleArrayEntry);
    CC_IFFREE(mbfs->lcipr, doubleArrayEntry);
    CC_IFFREE(mbfs->prefix_min_lcipr, double);
    CC_IFFREE(mbfs->suffix_max_huipr, double);
    for (int i=1;i<mbfs->mbfs_main->instance->numberOfNodes;i++){
        intvector_free(&mbfs->packed[i]);
        intvector_free(&mbfs->unpacked[i]);
    }
    CC_IFFREE(mbfs->packed, intvector);
    CC_IFFREE(mbfs->unpacked, intvector);
    CC_IFFREE(mbfs->intermediate, ttpObjective);
    intvector_free(&mbfs->to_pack);
    intvector_free(&mbfs->to_unpack);
    intvector_free(&mbfs->to_pack_loc);
    intvector_free(&mbfs->to_unpack_loc);
    CC_IFFREE(mbfs->city_to_packing, int);
    CC_IFFREE(mbfs->city_to_tour_idx, int);
}

int construct_mbfs_main(mbfs_main *main, ttpInstance *instance){
    main->instance = instance;
    main->ipr = instance->profitRatios;
    main->num_cities = instance->numberOfNodes;
    main->iprSortedIdx = CC_SAFE_MALLOC(instance->numberOfNodes, mbfs_entry *);
    city_to_tour_idx = CC_SAFE_MALLOC(instance->numberOfNodes,int);
    if (main->iprSortedIdx== (mbfs_entry **) NULL || city_to_tour_idx == (int *) NULL) goto FAILURE1;
    for (int i=1;i<instance->numberOfNodes;i++){
        main->iprSortedIdx[i] = CC_SAFE_MALLOC(instance->itemsPerCity,mbfs_entry);
        if (main->iprSortedIdx[i]== (mbfs_entry *) NULL) goto FAILURE0;
    }
    return 0;
  FAILURE0:
    for (int i=1;i<instance->numberOfNodes;i++) CC_IFFREE(main->iprSortedIdx[i], mbfs_entry);
  FAILURE1:
    CC_IFFREE(main->iprSortedIdx, mbfs_entry *);
    CC_IFFREE(city_to_tour_idx, int);
    return 1;
}

int construct_mbfs(mbfs *mbfs, mbfs_main *main){
    ttpInstance *instance = main->instance;
    mbfs->huipr = CC_SAFE_MALLOC(instance->numberOfNodes,doubleArrayEntry);
    mbfs->lcipr = CC_SAFE_MALLOC(instance->numberOfNodes,doubleArrayEntry);
    mbfs->packed = CC_SAFE_MALLOC(instance->numberOfNodes,intvector);
    mbfs->unpacked = CC_SAFE_MALLOC(instance->numberOfNodes,intvector);
    mbfs->city_to_packing = CC_SAFE_MALLOC(instance->numberOfItems,int);
    mbfs->city_to_tour_idx = CC_SAFE_MALLOC(instance->numberOfNodes,int);
    mbfs->prefix_min_lcipr = CC_SAFE_MALLOC(instance->numberOfNodes,double);
    mbfs->suffix_max_huipr = CC_SAFE_MALLOC(instance->numberOfNodes,double);
    mbfs->intermediate = CC_SAFE_MALLOC(instance->numberOfNodes, ttpObjective);
    if (mbfs->huipr==(doubleArrayEntry *) NULL || 
        mbfs->lcipr==(doubleArrayEntry *) NULL || 
        mbfs->packed==(intvector *) NULL || 
        mbfs->unpacked==(intvector *) NULL || 
        mbfs->city_to_packing == (int *) NULL || 
        mbfs->suffix_max_huipr == (double *) NULL || 
        mbfs->prefix_min_lcipr == (double *) NULL || 
        mbfs->intermediate == (ttpObjective *) NULL) goto FAILURE1;

    int ret = 0;
    ret+=intvector_alloc(&mbfs->to_pack,instance->numberOfNodes);
    ret+=intvector_alloc(&mbfs->to_unpack,instance->numberOfNodes);
    ret+=intvector_alloc(&mbfs->to_pack_loc,instance->numberOfNodes);
    ret+=intvector_alloc(&mbfs->to_unpack_loc,instance->numberOfNodes);
    if (ret>0) goto FAILURE1;
 
    for (int i=0;i<instance->numberOfNodes;i++){
        intvector_init(&mbfs->packed[i]);
        intvector_init(&mbfs->unpacked[i]);
    }
    for (int i=1;i<instance->numberOfNodes;i++){
        ret+=intvector_alloc(&mbfs->packed[i],instance->itemsPerCity);
        ret+=intvector_alloc(&mbfs->unpacked[i],instance->itemsPerCity);
    }
    if (ret>0) goto FAILURE0;
    mbfs->mbfs_main = main;

    return 0;

  FAILURE0:
    for (int i=1;i<instance->numberOfNodes;i++){
        intvector_free(&mbfs->packed[i]);
        intvector_free(&mbfs->unpacked[i]);
    }
  FAILURE1:
    intvector_free(&mbfs->to_pack);
    intvector_free(&mbfs->to_unpack);
    intvector_free(&mbfs->to_pack_loc);
    intvector_free(&mbfs->to_unpack_loc);
    CC_IFFREE(mbfs->huipr, doubleArrayEntry);
    CC_IFFREE(mbfs->lcipr, doubleArrayEntry);
    CC_IFFREE(mbfs->packed, intvector);
    CC_IFFREE(mbfs->unpacked, intvector);
    CC_IFFREE(mbfs->city_to_packing, int);
    CC_IFFREE(mbfs->city_to_tour_idx, int);
    CC_IFFREE(mbfs->intermediate, ttpObjective);
    return 1;
}

int copy_mbfs_main(mbfs_main *src, mbfs_main *dest){
    dest->items_per_city = src->items_per_city;
    dest->num_cities = src->num_cities;
    for (int i=1;i<src->num_cities;i++){
        for (int j=0;j<src->items_per_city;j++){
            dest->iprSortedIdx[i][j] = src->iprSortedIdx[i][j];
        }
    }
    dest->ipr = src->ipr;
    return 0;
} 

int copy_mbfs(mbfs *src, mbfs *dest){
    int ret = 0;
    dest->mbfs_main = src->mbfs_main;
    int num_cities = src->mbfs_main->num_cities;
    for (int i=1;i<num_cities;i++){
        dest->huipr[i] = src->huipr[i];
        dest->lcipr[i] = src->lcipr[i];
        ret+=intvector_copy(&dest->packed[i],&src->packed[i]);
        ret+=intvector_copy(&dest->packed[i],&src->packed[i]);
    }
    int num_items = src->mbfs_main->items_per_city * num_cities;
    for (int i=0;i < num_items;i++){
        dest->city_to_packing[i] = src->city_to_packing[i];
    }
    return ret;
}


int lcipr(mbfs *mbfs, int city, doubleArrayEntry *ret){
    mbfs_main *main = mbfs->mbfs_main;
    ttpInstance *instance = main->instance;

    if (mbfs->packed[city].len==0) {
        ret->value=LCIPR_DEFAULT_VALUE;
        ret->index=instance->itemsPerCity;
    } else{
        int idx = intvector_peek(&mbfs->packed[city]);
        ret->value = main->iprSortedIdx[city][idx].ipr;
        ret->index = idx;
    }
    return 0;
}

int huipr(mbfs *mbfs, int city, doubleArrayEntry *ret){
    mbfs_main *main = mbfs->mbfs_main;
    ttpInstance *instance = main->instance;
    if (mbfs->unpacked[city].len==0) {
        ret->value=-1;
        ret->index=-1;
    } else{
        int idx = intvector_peek(&mbfs->unpacked[city]);
        ret->value = main->iprSortedIdx[city][idx].ipr;
        ret->index = idx;
    }
    return 0;
}

int prefix_min_lcipr(mbfs *mbfs, ttpSolution *sol, double init_ipr, int init_idx, bool early_stop){
    mbfs_main *main = mbfs->mbfs_main;
    ttpInstance *instance=  main->instance;
    double min_ipr = init_ipr;
    int ret = 0;
    for (int i=init_idx;i<instance->numberOfNodes;i++){
        int city = sol->tour[i];
        if (mbfs->lcipr[city].value<min_ipr){
            min_ipr = mbfs->lcipr[city].value;
            if (mbfs->packed[city].len>0){
                ret+=intvector_append(&mbfs->to_unpack,city);
                ret+=intvector_set(&mbfs->to_unpack_loc,city,mbfs->to_unpack.len-1);  
            }
            
        }
        if (early_stop && mbfs->prefix_min_lcipr[i]==min_ipr) break;
        mbfs->prefix_min_lcipr[i]=min_ipr;
    }
    if (ret>1) goto FAILURE;
    return 0;
  FAILURE:
    return 1;
}

int prefix_min_lcipr_default(mbfs *mbfs, ttpSolution *sol){
    return prefix_min_lcipr(mbfs, sol, LCIPR_DEFAULT_VALUE,1,false);
}

int suffix_max_huipr(mbfs *mbfs, ttpSolution *sol, double init_ipr, int init_idx, bool early_stop){
    mbfs_main *main = mbfs->mbfs_main;
    ttpInstance *instance=  main->instance;
    double max_ipr = init_ipr;
    int ret = 0;
    for (int i=init_idx;i>0;i--){
        int city = sol->tour[i];
        if (mbfs->huipr[city].value>max_ipr){
            max_ipr = mbfs->huipr[city].value;
            if (mbfs->unpacked[city].len>0){
                ret+=intvector_append(&mbfs->to_pack,city);
                ret+=intvector_set(&mbfs->to_pack_loc,city,mbfs->to_pack.len-1);
            }
        }
        if (early_stop && mbfs->suffix_max_huipr[i]==max_ipr) break;
        mbfs->suffix_max_huipr[i]=max_ipr;
    }
    if (ret>1) goto FAILURE;
    return 0;
  FAILURE:
    return 1;
}

int suffix_max_huipr_default(mbfs *mbfs, ttpSolution *sol){
    return suffix_max_huipr(mbfs, sol, -1,mbfs->mbfs_main->instance->numberOfNodes-1,false);
}

// Read the packed and unpacked indices at each city from a ttpSolution
int read_mbfs(mbfs *mbfs, ttpSolution *sol){
    // printf("READ MBFS\n");
    mbfs->total_weight=0;
    mbfs->total_profit=0;
    mbfs->sol = sol;
    mbfs_main *main = mbfs->mbfs_main;
    ttpInstance *instance = main->instance;

    // Read the city_to_packing map
    for (int i=1;i<instance->numberOfNodes;i++){
        for (int j=0;j<instance->itemsPerCity;j++){
            // int itemIndex = getItemIndex(instance, sol->tour[i],j);
            // int packingIndex = sol->isCityPacking ? itemIndex : getPackingIndex(instance, i,j);
            int itemIndex = instance->itemIndex[sol->tour[i]][j];
            int packingIndex = sol->isCityPacking ? itemIndex : instance->packingIndex[i][j];
            mbfs->city_to_packing[itemIndex] = packingIndex;
        }
    }

    // Read packed/unpacked and lcipr/huipr
    for (int i=1;i<instance->numberOfNodes;i++){
        intvector_clear(&mbfs->packed[i]);
        intvector_clear(&mbfs->unpacked[i]);
        // Go through item indices in increasing order of IPR
        for (int j=0;j<instance->itemsPerCity;j++){
            int itemIndex = main->iprSortedIdx[i][j].itemIndex;
            int packingIdx = mbfs->city_to_packing[itemIndex];
            if (sol->packing[packingIdx]){
                intvector_append(&mbfs->packed[i],j);
                mbfs->total_weight+=instance->weight[itemIndex];
                mbfs->total_profit+=instance->profit[itemIndex];
            } else {
                intvector_append(&mbfs->unpacked[i],j);
            }
        }
        // Reverse the packed indices for lcipr
        intvector_reverse(&mbfs->packed[i]);

        // Set lcipr and huipr. Equal to the top item in their respective packed/unpacked arrays if it exists, otherwise some default
        lcipr(mbfs, i, &mbfs->lcipr[i]);
        huipr(mbfs, i, &mbfs->huipr[i]);

        mbfs->city_to_tour_idx[sol->tour[i]]=i;
        city_to_tour_idx[sol->tour[i]]=i;
    }
    mbfs->city_to_tour_idx[sol->tour[0]]=0;
    // Compute the prefix min and suffix max of the lcipr and huipr sequences.
    prefix_min_lcipr(mbfs, sol, LCIPR_DEFAULT_VALUE, 1, false);
    suffix_max_huipr(mbfs, sol, -1, instance->numberOfNodes-1,false);

    mbfs->bestObj = evaluate_intermediate(instance, sol, mbfs->intermediate,0,true);

    return 0;
}


ttpObjective vir_eval_intermediate(ttpInstance *instance, ttpSolution *sol, ttpObjective *intermediate, int startIndex, int itemIndex, int packingIndex, bool update){
    double weight = 0;
    double time = 0;
    double distance = 0;
    if (startIndex>0){
      weight = intermediate[startIndex-1].finalweight;
      time = intermediate[startIndex-1].time;
      distance = intermediate[startIndex-1].distance;
    //   profit = intermediate[startIndex-1].profit;
    }
    double dprofit = 0;
    double dweight = 0;
    int itemCity=-1;
    if (itemIndex>=0){
        int itemCity = tourIndexFromPackingIndex(instance, packingIndex);
        dprofit = instance->profit[itemIndex];
        dweight = instance->weight[itemIndex];
        if (sol->packing[packingIndex]){
            dprofit*=-1;
            dweight*=-1;
        }
    }

    ttpObjective finalObj = intermediate[instance->numberOfNodes-1];
    double finalProfit = finalObj.profit;
    if (finalObj.finalweight+dweight>instance->capacityOfKnapsack && !update) {
        ttpObjective ret;
        ret.objective=-1E30*instance->rentingRatio;
        ret.finalweight=finalObj.finalweight+dweight;
        ret.time=1E30;
        return ret;
    }

    for (int i=startIndex;i<instance->numberOfNodes;i++){

        weight = intermediate[i].finalweight+(i>=itemCity ? dweight : 0);
        
        double ddist = intermediate[i].distance-distance;
        distance = intermediate[i].distance;
        double speed = getSpeed(instance, weight);
        
        time += weight>instance->capacityOfKnapsack ? (1E20*ddist)/instance->minSpeed : ddist/speed;
        // printf("mbfs tour %d weight %lf dist %lf time %lf\n",i,weight,ddist, ddist/speed);
    
        if (update && i>=itemCity){
          intermediate[i].profit += dprofit ;
        //   intermediate[i].time = time;
          intermediate[i].finalweight += dweight;
        //   intermediate[i].objective = intermediate[i].profit-time*instance->rentingRatio;
        }
        
        // double gotSpeed = instance->maxSpeed-weight*(instance->maxSpeed-instance->minSpeed)/instance->capacityOfKnapsack;
        // if (weight>instance->capacityOfKnapsack) printf("wc-diff %lf weight %lf capacity %lf maxSpeed %lf minSpeed %lf speed %lf gotSpeed %lf adjusted %lf time %lf\n", instance->capacityOfKnapsack-weight, weight, instance->capacityOfKnapsack, instance->maxSpeed, instance->minSpeed, getSpeed(instance, weight), gotSpeed, speed, dist/speed);
    }
    ttpObjective ob;
    ob.profit = finalProfit + dprofit;
    ob.time = time;
    ob.finalweight = weight;
    ob.distance=distance;
    ob.objective=ob.profit-time*instance->rentingRatio;
    // printf("mbfs profit %lf\n",ob.profit);
    
    // if (weight>instance->capacityOfKnapsack) printf("wc-diff %lf objective %lf profit %lf time %lf\n", instance->capacityOfKnapsack-weight,ob.objective, profit, time);
    return ob;
}

ttpObjective vir_eval(mbfs *mbfs, ttpSolution *sol, int city, int packingIndex, int itemIndex){
    ttpInstance *instance = mbfs->mbfs_main->instance;
    int tourIdx = mbfs->city_to_tour_idx[city];
    // sol->packing[packingIndex] = !sol->packing[packingIndex];
    // ttpObjective obj = evaluate_intermediate(instance, sol, mbfs->intermediate, tourIdx, false);
    ttpObjective obj = vir_eval_intermediate(instance, sol, mbfs->intermediate, tourIdx, itemIndex, packingIndex, false);
    // sol->packing[packingIndex] = !sol->packing[packingIndex];
    return obj;
}


ttpObjective vir_eval_test(mbfs *mbfs, ttpSolution *sol, int city, int packingIndex, int itemIndex){
    ttpInstance *instance = mbfs->mbfs_main->instance;
    int tourIdx = mbfs->city_to_tour_idx[city];
    sol->packing[packingIndex] = !sol->packing[packingIndex];
    ttpObjective obj = evaluate_intermediate(instance, sol, mbfs->intermediate, tourIdx, false);
    ttpObjective objtest = evaluate(instance, sol);
    if (abs(obj.objective-objtest.objective)>0.01){
        printf("Eval intermediate is wrong. Intermediate:\n");
        objective_print(&obj);
        printf("\nActual:\n");
        objective_print(&objtest);
        printf("\n");
        exit(0);
    }
    // ttpObjective obj = vir_eval_intermediate(instance, sol, mbfs->intermediate, tourIdx, itemIndex, packingIndex, false);
    sol->packing[packingIndex] = !sol->packing[packingIndex];
    return obj;
}

bool mbfs_rand_city(mbfs *mbfs, bool *to_pack, int* city_idx, bool sorted){
    int num_topack = mbfs->to_pack.len;
    int num_tounpack = mbfs->to_unpack.len;
    if (num_topack + num_tounpack==0) return false;
    if (sorted){
        if (num_tounpack==0 || (num_topack>0 && rand()%2==0)){
            *to_pack = true;
            *city_idx = mbfs->to_pack.len-1;
        } else {
            *to_pack = false;
            *city_idx = mbfs->to_unpack.len-1;
        }
        return true;
    } else {
        
        int rand_int = rand()%(num_topack+num_tounpack);
        if (rand_int>=num_topack){
            *to_pack = false;
            *city_idx = rand_int-num_topack;
        } else {
            *to_pack = true;
            *city_idx = rand_int;
        }
        return true;
    }
    
}

int reset_mbfs(mbfs *mbfs, ttpSolution *sol, int city){
    int ret = 0;
    ret+=intvector_clear(&mbfs->to_pack);
    ret+=intvector_clear(&mbfs->to_pack_loc);
    ret+=intvector_clear(&mbfs->to_unpack);
    ret+=intvector_clear(&mbfs->to_pack_loc);

    ret+=lcipr(mbfs, city, &mbfs->lcipr[city]);
    ret+=huipr(mbfs, city, &mbfs->huipr[city]);

    ret+= prefix_min_lcipr(mbfs, sol, LCIPR_DEFAULT_VALUE, 1, false);
    ret+=suffix_max_huipr(mbfs, sol, -1, mbfs->mbfs_main->instance->numberOfNodes-1,false);
    mbfs->bestObj = evaluate_intermediate(mbfs->mbfs_main->instance, sol, mbfs->intermediate,0,true);
    return ret;
}

int cmp_city(const void *a, const void *b)
{
  return city_to_tour_idx[* (int *)a] - city_to_tour_idx[* (int *)b];
}


ttpObjective run_mbfs(mbfs_main *main, ttpSolution *sol, int max_iter, double timebound, bool sorted){
    double starttime = CCutil_zeit();

    // sol_from_file(main->instance, sol, "test_sol.txt");

    mbfs mbfs;
    init_mbfs(&mbfs);
    construct_mbfs(&mbfs, main);
    read_mbfs(&mbfs, sol);
    if (sorted){
        intvector_qsort(&mbfs.to_pack,cmp_city);
        intvector_qsort(&mbfs.to_unpack,cmp_city);
    }
    ttpInstance *instance = main->instance;

    // for (int i=0;i<instance->numberOfNodes;i++) printf("(%lf, %lf)\n",mbfs.suffix_max_huipr[i], mbfs.prefix_min_lcipr[i]);

    

    // for (int i=0;i<mbfs.to_unpack.len;i++){
    //     int city = mbfs.to_unpack.arr[i];
    //     int idx = intvector_peek(&mbfs.packed[city]);
    //     int index = main->iprSortedIdx[city][idx].itemIndex;
    //     printf("insert item %d city %d lowpicked\n",index,city);
    // }

    // for (int i=0;i<mbfs.to_pack.len;i++){
    //     int city = mbfs.to_pack.arr[i];
    //     int idx = intvector_peek(&mbfs.unpacked[city]);
    //     int index = main->iprSortedIdx[city][idx].itemIndex;
    //     printf("insert item %d city %d highunpicked\n",index,city);
    // }
    // exit(0);

    bool to_pack;
    int city_idx;
    int iter=0;
    while (mbfs_rand_city(&mbfs, &to_pack, &city_idx, sorted)){
        if (timebound>=0 && (CCutil_zeit()-starttime>timebound)) break;
        // printf("MBFS ITERATION %d OBJECTIVE %lf\n",iter, mbfs.bestObj.objective);
        if (to_pack){
            int city = mbfs.to_pack.arr[city_idx];
            intvector_popidx(&mbfs.to_pack, &mbfs.to_pack_loc, city_idx);
            doubleArrayEntry huipr_entry = mbfs.huipr[city];
            if (huipr_entry.index>=instance->itemsPerCity) continue;
            mbfs_entry ipr_sorted_entry = main->iprSortedIdx[city][huipr_entry.index];
            if (mbfs.bestObj.finalweight+ipr_sorted_entry.weight>instance->capacityOfKnapsack) continue;
            int packingIdx = mbfs.city_to_packing[ipr_sorted_entry.itemIndex];

            // mbfs.bestObj = evaluate_intermediate(mbfs.mbfs_main->instance, sol, mbfs.intermediate,0,true);
            // ttpObjective testObj = vir_eval_test(&mbfs,sol, city, packingIdx, ipr_sorted_entry.itemIndex);
            ttpObjective newObj = vir_eval(&mbfs, sol, city, packingIdx, ipr_sorted_entry.itemIndex);
            // if (abs(testObj.objective-newObj.objective)>0.01){
            //     printf("tour %d city %d packing %d citypacking %d\n",mbfs.city_to_tour_idx[city],city,packingIdx,sol->isCityPacking);
            //     printf("vir_eval Differing objectives: Test\n");
            //     objective_print(&testObj);
            //     printf("\nActual:\n");
            //     objective_print(&newObj);
            //     printf("\n");
            //     exit(0);
            // }
            
            if (newObj.objective>mbfs.bestObj.objective){
                int tourIdx = mbfs.city_to_tour_idx[city];
                
                // mbfs.bestObj = evaluate_intermediate(mbfs.mbfs_main->instance, sol, mbfs.intermediate,0,true);
                mbfs.bestObj = vir_eval_intermediate(instance, sol, mbfs.intermediate, tourIdx, packingIdx, ipr_sorted_entry.itemIndex, true);
                sol->packing[packingIdx] = !sol->packing[packingIdx];
                // mbfs.bestObj = evaluate_intermediate(instance, sol, mbfs.intermediate, tourIdx, true);
                // testObj = evaluate_intermediate(mbfs.mbfs_main->instance, sol, mbfs.intermediate,0,true);
                // if (abs(testObj.objective-mbfs.bestObj.objective)>0.01){
                //     printf("tour %d city %d packing %d citypacking %d\n",tourIdx,city,packingIdx,sol->isCityPacking);
                //     printf("eval_intermediate Differing objectives: Test\n");
                //     objective_print(&testObj);
                //     printf("\nActual:\n");
                //     objective_print(&mbfs.bestObj);
                //     printf("\n");
                //     exit(0);
                // }
                //refresh
                intvector_append(&mbfs.packed[city],huipr_entry.index);
                intvector_bubble(&mbfs.packed[city],false);
                intvector_pop(&mbfs.unpacked[city]);
                
                reset_mbfs(&mbfs, sol, city);
                // reset_mbfs(&mbfs,sol);
                if (sorted){
                    intvector_qsort(&mbfs.to_pack,cmp_city);
                    intvector_qsort(&mbfs.to_unpack,cmp_city);
                }
            }
        } else {
            int city = mbfs.to_unpack.arr[city_idx];
            intvector_popidx(&mbfs.to_unpack, &mbfs.to_unpack_loc, city_idx);
            doubleArrayEntry lcipr_entry = mbfs.lcipr[city];
            if (lcipr_entry.index<0) continue;
            mbfs_entry ipr_sorted_entry = main->iprSortedIdx[city][lcipr_entry.index];
            int packingIdx = mbfs.city_to_packing[ipr_sorted_entry.itemIndex];

            // mbfs.bestObj = evaluate_intermediate(mbfs.mbfs_main->instance, sol, mbfs.intermediate,0,true);
            // ttpObjective testObj = vir_eval_test(&mbfs,sol,city, packingIdx, ipr_sorted_entry.itemIndex);
            ttpObjective newObj = vir_eval(&mbfs, sol, city, packingIdx, ipr_sorted_entry.itemIndex);
            // if (abs(testObj.objective-newObj.objective)>0.01){
            //     printf("tour %d city %d packing %d citypacking %d\n",mbfs.city_to_tour_idx[city],city,packingIdx,sol->isCityPacking);
            //     printf("vir_eval2 Differing objectives: Test\n");
            //     objective_print(&testObj);
            //     printf("\nActual:\n");
            //     objective_print(&newObj);
            //     printf("\n");
            //     exit(0);
            // }

            if (newObj.objective>mbfs.bestObj.objective){
                
                int tourIdx = mbfs.city_to_tour_idx[city];

                // mbfs.bestObj = evaluate_intermediate(mbfs.mbfs_main->instance, sol, mbfs.intermediate,0,true);
                mbfs.bestObj = vir_eval_intermediate(instance, sol, mbfs.intermediate, tourIdx, packingIdx, ipr_sorted_entry.itemIndex, true);
                sol->packing[packingIdx] = !sol->packing[packingIdx];
                // mbfs.bestObj = evaluate_intermediate(instance, sol, mbfs.intermediate, tourIdx, true);
                // testObj = evaluate_intermediate(mbfs.mbfs_main->instance, sol, mbfs.intermediate,0,true);
                // if (abs(testObj.objective-mbfs.bestObj.objective)>0.01){
                //     printf("tour %d city %d packing %d citypacking %d\n",tourIdx,city,packingIdx,sol->isCityPacking);
                //     printf("eval_intermediate2 Differing objectives: Test\n");
                //     objective_print(&testObj);
                //     printf("\nActual:\n");
                //     objective_print(&mbfs.bestObj);
                //     printf("\n");
                //     exit(0);
                // }
                intvector_append(&mbfs.unpacked[city],lcipr_entry.index);
                intvector_bubble(&mbfs.unpacked[city],true);
                intvector_pop(&mbfs.packed[city]);
                
                reset_mbfs(&mbfs, sol, city);
                // reset_mbfs(&mbfs,sol);
                if (sorted){
                    intvector_qsort(&mbfs.to_pack,cmp_city);
                    intvector_qsort(&mbfs.to_unpack,cmp_city);
                }
            }
        }
        iter++;
        if (max_iter>=0 && iter>max_iter) break;
    }

    // printf("BEST SOL %lf\n",mbfs.bestObj.objective);
    // exit(0);
    free_mbfs(&mbfs);
    return mbfs.bestObj;
}

bool pgch_update(mbfs *mbfs, ttpSolution *sol, city_info *info, int startIndex, int stopIndex, double *total_profit, double *total_weight){
    mbfs_main *main = mbfs->mbfs_main;
    ttpInstance *instance=  main->instance;
    double totalweight = mbfs->total_weight;
    double totalprofit = mbfs->total_profit;
    for (int i=startIndex;i<=stopIndex;i++){
        double lcipr = mbfs->prefix_min_lcipr[i];
        int city = sol->tour[i];
        intvector packed = mbfs->packed[city];
        for (int j=packed.len-1;j>=0;j--){
            mbfs_entry entry = main->iprSortedIdx[city][packed.arr[j]];
            if (entry.ipr<lcipr){
                int itemIndex = entry.itemIndex;
                int packingIndex = sol->isCityPacking ? itemIndex : instance->packingIndex[i][entry.itemNumber];
                sol->packing[packingIndex] = false;
                double weight = instance->weight[itemIndex];
                double profit = instance->profit[itemIndex];
                info[city].dweight-=weight;
                // info[city].dprofit-=profit;
                totalweight-=weight;
                totalprofit-=profit;
            } else break;
        }
    }

    for (int i=startIndex;i<=stopIndex;i++){
        double huipr = mbfs->suffix_max_huipr[i];
        int city = sol->tour[i];
        intvector unpacked = mbfs->unpacked[city];
        for (int j=unpacked.len-1;j>=0;j--){
            mbfs_entry entry = main->iprSortedIdx[city][unpacked.arr[j]];
            if (entry.ipr>huipr){
                int itemIndex = entry.itemIndex;
                int packingIndex = sol->isCityPacking ? itemIndex : instance->packingIndex[i][entry.itemNumber];
                double weight = instance->weight[itemIndex];
                double profit = instance->profit[itemIndex];
                if (totalweight+weight<=instance->capacityOfKnapsack){
                    sol->packing[packingIndex] = true;
                    info[city].dweight+=weight;
                    // info[city].dprofit+=profit;
                    totalweight+=weight;
                    totalprofit+=profit;
                }
            } else break;
        }
    }
    // if (totalweight>instance->capacityOfKnapsack){
    //     printf("WHY IS TOTAL WEIGHT LARGE IN PGCH UPDATE INIT %lf FINAL %lf CAP %lf\n",mbfs->total_weight, totalweight,instance->capacityOfKnapsack);
    //     exit(0);
    // }
    *total_profit=totalprofit;
    *total_weight=totalweight;
    if (abs(totalweight-mbfs->total_weight)<1 && (totalweight!=mbfs->total_weight)) {
        printf("HIII EQUAL? %d\n",totalweight==mbfs->total_weight);
    }
    return totalweight==mbfs->total_weight;
}

int reverse_segment(mbfs_main *main, int *src, int *dest, int posNumA, int posNumB){
    ttpInstance *instance=  main->instance;
    int numCities = instance->numberOfNodes;
    if (0<posNumA && posNumA<posNumB) {    
        for (int p=0;p<posNumA;p++) {
            dest[p]      = src[p];
        }
        for (int p = posNumA; p <= posNumB; p++) {
            int  q = posNumB-p+posNumA;
            dest[p]      = src[q];
        }
        dest[(posNumB+1)%numCities] = src[(posNumB+1)%numCities];
        for (int p=posNumB+1;p<instance->numberOfNodes;p++) {
            dest[p]      = src[p];
        }
    }
    else {
        int p = 1, last  ;
        dest[0]      = 0;
        if (posNumA == 0)
            last = numCities-1;
        else {  // posNumB < posNumA    
            for (int q = numCities-1; q >= posNumA; q--, p++) {
                dest[p]      = src[q];
            }
            last = posNumA-1; 
        }
        for (int q = posNumB+1; q <= last; q++, p++) {
            dest[p]      = src[q];
        }
        for (int q = posNumB; q > 0; q--, p++) {
            dest[p]      = src[q];
        }
    }  

    return 0;
}
// int real_reverse_segment(mbfs_main *main, int *tour,  int PosNumA, int PosNumB){
//     for (int p = 0; p < (PosNumB+1-PosNumA)/2; p++) {
//         SWAP(tour[PosNumA+p],tour[PosNumB-p], int);
//         SWAP(solPositions[PosNumA+p].addWeight,solPositions[PosNumB-p].addWeight);
//         SWAP(solCities[tour[PosNumA+p]].posNum,solCities[tour[PosNumB-p]].posNum);
//     }
// }


ttpObjective eval_via_info(mbfs *mbfs, ttpSolution *sol, ttpObjective *intermediate, city_info *info, int startIndex, int endIndex, double total_profit, bool sameweight){
    ttpInstance *instance = mbfs->mbfs_main->instance;
    double time = 0;
    double weight = 0;
    if (startIndex>0){
      weight = intermediate[startIndex-1].finalweight;
      time = intermediate[startIndex-1].time;
    }

    for (int i=startIndex;i<=endIndex;i++){
        int city = sol->tour[i];
        weight += info[city].dweight;
        double dist; 
        GET_DISTANCE(dist, instance, city, sol->tour[(i+1)%instance->numberOfNodes]);
        double speed = GET_SPEED(instance, weight);
        time += dist/speed;
    }
    if (sameweight){
        // double time2 = 0;
        // for (int i=endIndex+1;i<instance->numberOfNodes;i++){
        //     int city = sol->tour[i];
        //     weight += info[city].dweight;
        //     double dist; 
        //     GET_DISTANCE(dist, instance, city, sol->tour[(i+1)%instance->numberOfNodes]);
        //     double speed = GET_SPEED(instance, weight);
        //     time2 += dist/speed;
        // }
        // if (abs(time2-(intermediate[instance->numberOfNodes-1].time-intermediate[endIndex].time))>0.01){
        //     printf("UNEQUAL TIMES\n");
        //     printf("Actual time %lf\n",time2);
        //     printf("e-1 %lf e %lf e+1 %lf\n",intermediate[instance->numberOfNodes-1].time-intermediate[endIndex-1].time, intermediate[instance->numberOfNodes-1].time-intermediate[endIndex].time, intermediate[instance->numberOfNodes-1].time-intermediate[endIndex+1].time);
        //     exit(0);
        // }
        time+=intermediate[instance->numberOfNodes-1].time-intermediate[endIndex].time;
    } else {
        for (int i=endIndex+1;i<instance->numberOfNodes;i++){
            int city = sol->tour[i];
            weight += info[city].dweight;
            double dist; 
            GET_DISTANCE(dist, instance, city, sol->tour[(i+1)%instance->numberOfNodes]);
            double speed = GET_SPEED(instance, weight);
            time += dist/speed;
        }
    }

    ttpObjective ob;
    ob.profit = -1;
    ob.time = time;
    ob.finalweight = -1;
    ob.distance=-1;
    ob.objective=total_profit-time*instance->rentingRatio;
    
    return ob;
}


ttpObjective prep_eval(ttpInstance *instance, ttpSolution *sol, ttpObjective *intermediate, city_info *info){
    double weight = 0;
    double time = 0;
    double distance = 0;
    double profit = 0;
    
    for (int i=0;i<instance->numberOfNodes;i++){
        int city = sol->tour[i];
        double dweight = 0;
        double dprofit = 0;
        if (i>0){
            for (int j=0;j<instance->itemsPerCity;j++){
                // int packingIdx = sol->isCityPacking ? getItemIndex(instance, city,j) : getPackingIndex(instance, i,j);
                int itemIndex = instance->itemIndex[sol->tour[i]][j];
                int packingIdx = sol->isCityPacking ? itemIndex : instance->packingIndex[i][j];
                if (sol->packing[packingIdx]) {
                    // int itemIdx = getItemIndex(instance, city,j);
                    dweight+=instance->weight[itemIndex];
                    dprofit+=instance->profit[itemIndex];
                }
            }
        }
        weight += dweight;
        profit += dprofit;
        double dist = getDistance(instance, sol->tour[i], sol->tour[(i+1)%instance->numberOfNodes]);
        distance+=dist;
        double speed = getSpeed(instance, weight);
        time += weight>instance->capacityOfKnapsack ? (1E20*dist)/instance->minSpeed : dist/speed;
        intermediate[i].profit = profit;
        intermediate[i].time = time;
        intermediate[i].distance = distance;
        intermediate[i].finalweight = weight;
        intermediate[i].objective = profit-time*instance->rentingRatio;
        info[city].dweight = dweight;
        info[city].dprofit = dprofit;
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

ttpObjective pgch(mbfs_main *main, ttpSolution *sol, int maxIter){
    // printf("PGCH START\n");
    mbfs mbfs;
    init_mbfs(&mbfs);
    construct_mbfs(&mbfs, main);
    read_mbfs(&mbfs, sol);
    ttpInstance *instance = main->instance;
    ttpObjective initialObj = evaluate(instance, sol);
    // printf("WEIGHT %lf MBFS %lf\n",initialObj.finalweight,mbfs.total_weight);
    ttpObjective bestObj = initialObj;
    ttpObjective prevObj = bestObj;
    ttpObjective *intermediate = CC_SAFE_MALLOC(instance->numberOfNodes,ttpObjective);
    city_info *info = CC_SAFE_MALLOC(instance->numberOfNodes,city_info);
    city_info *info2 = CC_SAFE_MALLOC(instance->numberOfNodes,city_info);
    prep_eval(instance, sol, intermediate, info);
    // evaluate_intermediate(instance, sol, intermediate, 0, true);
    // printf("tour: (%d %d %d %d)\n",sol->tour[0],sol->tour[1],sol->tour[2],sol->tour[3]);
    // objective_print(&intermediate[0]);
    // printf("\n");
    
    
    int *city_to_tour_idx = CC_SAFE_MALLOC(instance->numberOfNodes, int);
    for (int i=0;i<instance->numberOfNodes;i++) city_to_tour_idx[sol->tour[i]]=i;

    ttpSolution newSol;
    newSol.isCityPacking = sol->isCityPacking;
    newSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    newSol.tour = CC_SAFE_MALLOC(instance->numberOfNodes+1, int);
    newSol.tour[instance->numberOfNodes]=0;
    int num_iter = 0;
    int a_todo = -1;
    int b_todo = -1;
    bool improved = false;
    ttpSolution cacheSol;
    cacheSol.isCityPacking = sol->isCityPacking;
    cacheSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    cacheSol.tour = CC_SAFE_MALLOC(instance->numberOfNodes+1, int);
    cacheSol.tour[instance->numberOfNodes]=0;
    double time = CCutil_zeit();
    double overalltime = 0;
    double copyBoolTime = 0;
    double updateTime = 0;
    double evalTime = 0;
    double revtime=0;
    int updatecities = 0;
    double prev_time = CCutil_zeit();
    int numCities = instance->numberOfNodes;
    int numItems = instance->numberOfItems;
    int iter = 0;
    do {
        a_todo = -1;
        b_todo = -1;
        for (int i=0;i<numCities;i++){
            info2[i].dweight=info[i].dweight;
            // info2[i].dprofit=info[i].dprofit;
        }
        // printf("PGCH ITERATION %d obj %lf\n",num_iter,bestObj.objective);
        iter++;
        prevObj = bestObj;
        for (int i=1;i<numCities;i++){
            // if (i%1000==999) printf("Iter %d PosNumA %d TIME %lf\n",iter, i,CCutil_zeit()-time);
            // printf("Average time %lf\n",(CCutil_zeit()-time)/num_iter);
            //PGCH ITERATION %d loc %d obj %lf\n / num_iter,i,bestObj.objective,

            // overalltime=0.99*overalltime+0.01*(CCutil_zeit()-prev_time);

            // printf("copyBool %le, revtime %le, adjtime %le per city or %le per iter, evaltime %le, total time %le ewm %le\n",copyBoolTime, revtime, updateTime, updateTime, evalTime, (CCutil_zeit()-time)/num_iter, overalltime);
            
            // prev_time = CCutil_zeit();
            
            int cityi = sol->tour[i];
            int posi = i;
            // if (maxIter>=0 && num_iter>maxIter) break;
            intvector neighbors = instance->neighbors[cityi];
            for (int j=0;j<neighbors.len;j++){
                
                int cityj = neighbors.arr[j];
                int posj = city_to_tour_idx[cityj];
                if (posj<posi) continue;

                // double prevTime = CCutil_zeit();
                reverse_segment(main, sol->tour, newSol.tour, posi, posj);
                // revtime=0.99*revtime + 0.01*(CCutil_zeit()-prevTime);

                // prevTime = CCutil_zeit();
                copyBool(sol->packing, newSol.packing, numItems);
                // copyBoolTime= 0.99*copyBoolTime+0.01*(CCutil_zeit()-prevTime);
                
                // prevTime = CCutil_zeit();
                
                double total_profit, total_weight;
                
                
                bool sameweight = pgch_update(&mbfs, &newSol, info2, posi, posj, &total_profit, &total_weight);
                
                
                
                
                // updateTime=0.99*updateTime+0.01*(CCutil_zeit()-prevTime);
                // updatecities+=posj-posi+1;

                // prevTime = CCutil_zeit();
                // ttpObjective obj;
                // if (total_weight>instance->capacityOfKnapsack){
                //     printf("WEIGHT SHOULD NOT BE THIS BIG %lf CAP %lf\n",total_weight,instance->capacityOfKnapsack);
                //     exit(1);
                //     obj.finalweight=total_weight;
                //     obj.objective=-1E30*instance->rentingRatio;
                // } else {
                // double inittime = CCutil_zeit();
                ttpObjective obj;
                // for (int kk=0;kk<1000;kk++){
                obj = eval_via_info(&mbfs, &newSol, intermediate, info2, posi-1, posj, total_profit, sameweight);
                // }
                // double endtime = CCutil_zeit();
                // printf("1000 iterations posA %d posB %d sameweight %d: %le\n",posi, posj, sameweight, endtime-inittime);
                // exit(0);
                // }
                 
                // evalTime=0.99*evalTime+0.01*(CCutil_zeit()-prevTime);

                // printf("tour: (%d %d %d %d)\n",newSol.tour[0],newSol.tour[1],newSol.tour[2],newSol.tour[3]);
                // ttpObjective obj = evaluate_intermediate(instance, &newSol, intermediate, posi-1, false);
                // ttpObjective obj2 = evaluate(instance, &newSol);
                // printf("posA %d cityA %d posB %d cityB %d gain: %lf best: %lf\n",posi,cityi,posj,cityj,obj.objective,bestObj.objective);
                // if (obj2.objective!=obj.objective){
                //     printf("UNEQUAL OBJECTIVES obj %lf obj2 %lf\n",obj.objective,obj2.objective);
                //     evaluate_intermediate(instance, &newSol, intermediate, posi-1, true);
                //     for (int kk=0;kk<instance->numberOfNodes;kk++){
                //         int city = sol->tour[kk];
                //         printf("(w %lf %lf)\n",intermediate[kk].finalweight - (kk==0 ? 0 : intermediate[kk-1].finalweight),info[city].dweight);
                //     }
                //     exit(0);
                // }
                

                if (obj.objective>=bestObj.objective){
                    improved = true;
                    a_todo = posi;
                    b_todo = posj;
                    // copyBool(newSol.packing, cacheSol.packing, instance->numberOfItems);
                    // copyInt(newSol.tour,cacheSol.tour,instance->numberOfNodes);
                    bestObj=obj;
                }

                for (int i=posi;i<=posj;i++){
                    int city = newSol.tour[i];
                    info2[city].dweight=info[city].dweight;
                    // info2[city].dprofit=info[city].dprofit;
                }

                if (maxIter>=0 && num_iter>maxIter) break;
                num_iter++;
            }
            // if (ttp_lk_params.fast_pgch && improved) {
            //     reverse_segment(main, sol->tour, cacheSol.tour, a_todo, b_todo);
            //     copyBool(sol->packing, cacheSol.packing, instance->numberOfItems);
            //     double total_profit, total_weight;
            //     pgch_update(&mbfs, &cacheSol, info2, a_todo, b_todo, &total_profit, &total_weight);

            //     copySolution(instance, sol, &cacheSol);
            //     for (int i=0;i<instance->numberOfNodes;i++) city_to_tour_idx[sol->tour[i]]=i;
            //     read_mbfs(&mbfs,sol);
            //     // ttpObjective tempObj = evaluate(instance, sol);
            //     // printf("WEIGHT %lf MBFS %lf\n",tempObj.finalweight,mbfs.total_weight);
            //     prep_eval(instance, sol, intermediate, info);
            //     // evaluate_intermediate(instance, sol, intermediate, 0, false);
            //     improved = false;
            // }
        }
        if (improved){
            reverse_segment(main, sol->tour, cacheSol.tour, a_todo, b_todo);
            copyBool(sol->packing, cacheSol.packing, instance->numberOfItems);
            double total_profit, total_weight;
            pgch_update(&mbfs, &cacheSol, info2, a_todo, b_todo, &total_profit, &total_weight);

            copySolution(instance, sol, &cacheSol);
            for (int i=0;i<instance->numberOfNodes;i++) city_to_tour_idx[sol->tour[i]]=i;
            read_mbfs(&mbfs,sol);
            // ttpObjective tempObj = evaluate(instance, sol);
            // printf("WEIGHT %lf MBFS %lf\n",tempObj.finalweight,mbfs.total_weight);
            prep_eval(instance, sol, intermediate, info);
            // evaluate_intermediate(instance, sol, intermediate, 0, false);
            improved = false;
        }
    } while (bestObj.objective>prevObj.objective*1.0001);

    
    free_mbfs(&mbfs);
    // printf("obj %lf\n",bestObj.objective);
    return bestObj;
}


/**
 * Comparator for sorting doubleArrayEntries in descending order
 */
int comparingMbfsEntry(const void *a, const void *b){
    mbfs_entry *a_ = (mbfs_entry *) a;
    mbfs_entry *b_ = (mbfs_entry *) b;
    if (a_->ipr < b_->ipr) return -1;
    else if (a_->ipr > b_->ipr) return 1;
    else return 0;
}

int read_mbfs_main(mbfs_main *main, ttpInstance *instance){
    for (int i=1;i<instance->numberOfNodes;i++){
        for (int j=0;j<instance->itemsPerCity;j++){
            main->iprSortedIdx[i][j].itemNumber = j;
            main->iprSortedIdx[i][j].city = i;
            int itemIndex = instance->itemIndex[i][j];
            // int itemIndex = getItemIndex(instance, i, j);
            main->iprSortedIdx[i][j].itemIndex = itemIndex;
            main->iprSortedIdx[i][j].weight = instance->weight[itemIndex];
            main->iprSortedIdx[i][j].profit = instance->profit[itemIndex];
            main->iprSortedIdx[i][j].ipr = instance->profitRatios[itemIndex];
        }
        qsort(main->iprSortedIdx[i], (size_t) instance->itemsPerCity, sizeof(main->iprSortedIdx[i][0]), comparingMbfsEntry);
    }
    return 0;
}


void mbfs_entry_print(mbfs_entry *entry){
    printf("((%d %d %d), w %lf p %lf ipr %lf)",entry->city, entry->itemNumber, entry->itemIndex, entry->weight, entry->profit, entry->ipr);
}