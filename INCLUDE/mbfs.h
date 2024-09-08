#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "util.h"
#include "ttp_util.h"

#ifndef MBFS_H

#define MBFS_H

/*
 * An entry for the iprSortedIdx
 * Fields:
 *  int itemNumber
 *  int city
 *  int itemIndex
 *  double ipr
 *  double profit
 *  double weight
 */
typedef struct mbfs_entry{
    int itemNumber;
    int city;
    int itemIndex;
    double ipr;
    double profit;
    double weight;
} mbfs_entry;


// Contains all the global information needed for mbfs. One mbfs_main per ttpInstance.
typedef struct mbfs_main{
    // For each city contains an mbfs_entry in order of increasing profit ratio
    mbfs_entry **iprSortedIdx;
    // Number of cities in the instance
    int num_cities;
    // Number of items per city
    int items_per_city;
    // Contains the profit ratio of the item at itemIndex. Points to the same place as instance->profitRatios
    double *ipr;
    // A pointer to the global ttpInstance
    ttpInstance *instance;
} mbfs_main;


/* Contains all the local information for mbfs. One mbfs per ttpSolution. 
 *  Fields:
 *  doubleArrayEntry *lcipr, *huipr;
 *  double *prefix_min_lcipr, *suffix_max_huipr;
 *  intvector *packed, *unpacked;
 *  int *city_to_packing;
 *  mbfs_main *mbfs_main;
 *  intvector topack, topack_loc, tounpack, tounpack_loc;
 *  ttpSolution *sol;
 */
typedef struct mbfs{
    /* For each city contains the lowest uncollected / highest collected [profit ratio, index of the entry in main->iprSortedIdx]
     * Otherwise has [instance->itemsPerCity, 1E6] for lcipr and [-1, -1] for huipr
     */
    doubleArrayEntry *lcipr, *huipr;

    // Prefix minimum sequence of the lcipr
    double *prefix_min_lcipr;
    // Suffix maximum sequence of the huipr
    double *suffix_max_huipr;

    // For each city contains the indices of the packed entries in main->iprSortedIdx from high to low
    intvector* packed;
    // For each city contains the indices of the unpacked entries in main->iprSortedIdx from low to high
    intvector* unpacked;

    // Maps the itemIdx to the packingIdx
    int *city_to_packing;
    // Maps the city to the index of its occurrence in the tour
    int *city_to_tour_idx;

    // Pointer to the global mbfs_main object
    mbfs_main *mbfs_main;

    // Contains a list of cities that have items that need to be packed/unpacked
    intvector to_unpack, to_pack;
    // Auxiliary arrays (for to_unpack and to_pack) to facilitate removing cities that have been checked
    intvector to_unpack_loc, to_pack_loc;

    // The ttpSolution we're working on
    ttpSolution *sol;

    // Cached intermediate objectives for faster evaluation.
    ttpObjective *intermediate;
    // Best seen objective
    ttpObjective bestObj;

    double total_weight;
    double total_profit;
} mbfs;

typedef struct city_info{
    double dweight;
    double dprofit;
} city_info;

void init_mbfs_main(mbfs_main *main), init_mbfs(mbfs *mbfs), free_mbfs_main(mbfs_main *main), free_mbfs(mbfs *mbfs), mbfs_entry_print(mbfs_entry *entry);

int construct_mbfs_main(mbfs_main *main, ttpInstance *instance), construct_mbfs(mbfs *mbfs, mbfs_main *main),
    read_mbfs_main(mbfs_main *main, ttpInstance *instance), 
    read_mbfs(mbfs *mbfs, ttpSolution *sol), reset_mbfs(mbfs *mbfs, ttpSolution *sol, int city),
    copy_mbfs_main(mbfs_main *src, mbfs_main *dest), copy_mbfs(mbfs *src, mbfs *dest);

bool mbfs_rand_city(mbfs *mbfs, bool *to_pack, int* city_idx, bool sorted);

ttpObjective pgch(mbfs_main *main, ttpSolution *sol, int maxIter);
ttpObjective run_mbfs(mbfs_main *main, ttpSolution *sol, int max_iter, double timebound, bool sorted);
ttpObjective vir_eval_intermediate(ttpInstance *instance, ttpSolution *sol, ttpObjective *intermediate, int startIndex, int itemIndex, int packingIndex, bool update);
// inline double eval_via_info(mbfs *mbfs, ttpSolution *sol, ttpObjective *intermediate, city_info *info, int startIndex, int endIndex, double total_profit, bool sameweight);
ttpObjective prep_eval(ttpInstance *instance, ttpSolution *sol, ttpObjective *intermediate, city_info *info);

#endif 
