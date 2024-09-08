/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                     CHAINED LIN-KERNIGHAN                                */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 22, 1995                                                    */
/*        May 1, 1998 (bico)                                                */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CClinkern_tour (int ncount, CCdatagroup *dat, int ecount,           */
/*      int *elist, int stallcount, int repeatcount, int *incycle,          */
/*      int *outcycle, double *val                                          */
/*      int silent, double time_bound, double length_bound,                 */
/*      char *saveit_name, int kicktype, CCrandstate *rstate)               */
/*    RUNS Chained Lin-Kernighan.                                           */
/*    -ncount (the number of nodes in the graph)                           */
/*    -dat (coordinate dat)                                                 */
/*    -ecount (the number of good edges - should not be 0)                  */
/*    -elist (the good edges in end1 end2 format)                           */
/*    -stallcount (the max number of 4-swaps without progress               */
/*    -repeatcount (the number of 4-swap kicks)                             */
/*    -incycle (a starting cycle, in node node node format - can be NULL)   */
/*    -outcycle (returns the cycle - can be NULL)                           */
/*    -run_slightly (if nonzero, then very little info will be printed)     */
/*    -time_bound (if postive, then the search will stop after the kick     */
/*       that puts the running time above this number of seconds)           */
/*    -length_bound (if postive, then the search will stop after the kick   */
/*       that puts the tour at or below this length)                        */
/*    -saveit_name (if non NULL then the tour will be saved to this file    */
/*       after every 10000 kicks - if it has improved)                      */
/*    -kicktype (specifies the type of kick used - should be one of         */
/*       CC_LK_RANDOM_KICK, CC_LK_GEOMETRIC_KICK, CC_LK_CLOSE_KICK, or      */
/*       CC_LK_WALK_KICK)                                                   */
/*                                                                          */
/*    NOTES: If incycle is NULL, then a random starting cycle is used. If   */
/*     outcycle is not NULL, then it should point to an array of length     */
/*     at least ncount.                                                     */
/*                                                                          */
/****************************************************************************/

#include <stdbool.h>
#include "machdefs.h"
#include "kdtree.h"
#include "util.h"
#include "macrorus.h"
#include "linkern.h"
#include "math_util.h"
#include "ttp_util.h"
#include "ttp_heuristic.h"
#include "ttp_optimisation.h"
#include "ttp_linkern.h"
#include "mbfs.h"
#include "ttp.h"

int numHillClimbers = 0, numZeroClimbers=0, numChecks=0, beneficialOpts=0,
    linearImprovements = 0, classifierImprovements = 0;
double linearTime=0, hillClimberTime=0, stepTime=0, checkTime=0;
adaptive_params ttp_lk_params;
//somehow print to a file to track the distribution of successes across gmaxes.

#define MAXDEPTH       25   /* Shouldn't really be less than 2.             */
#define KICK_MAXDEPTH  50
#define IMPROVE_SWITCH -1   /* When to start using IMPROVE_KICKS (-1 never) */
#define LONG_KICKER  
#define ACCEPT_TIES 
#undef  ACCEPT_BAD_TOURS 

#define GMAX_THRESHOLD ttp_lk_params.Gmax_threshold //-(CCutil_MAXDOUBLE-1)
#define CHECK_OBJECTIVE
// #undef CHECK_OBJECTIVE
#define NORMAL_LINKERN_ROUNDS 0
#define GA_POPSIZE 0
#define NORMAL_NONIMPROVEMENT 10
#define MBFS_MAX_ITER -1
#define QUADRATIC_MAX_NONIMPROVEMENT 50
#define QUADRATIC_MAX 500

#define USE_LESS_OR_EQUAL 
#define SUBTRACT_GSTAR 
#undef  SWITCH_LATE
#define LATE_DEPTH 10      /* Should be less than MAXDEPTH                 */

#define MAK_MORTON       
#undef  FULL_MAK_MORTON 
#undef  NODE_INSERTIONS

#undef  USE_HEAP
#undef  MARK_NEIGHBORS      /* Mark the good-edge neighbors after swaps     */
#define USE_LESS_MARKING    /* Do not mark the tour neighbors after swaps   */
#define MARK_LEVEL 10       /* Number of tour neighbors after 4-swap kick   */
#define BACKTRACK   4
#define MAX_BACK   12       /* Upper bound on the XXX_count entries         */
static const int backtrack_count[BACKTRACK] = {4, 3, 3, 2};
static const int weird_backtrack_count[3] = {4, 3, 3};

#define BIGINT 2000000000
#define Edgelen(n1, n2, D)  dist (n1, n2, D)
/*
#define Edgelen(n1, n2, D)  CCutil_dat_edgelen (n1, n2, D->dat) 
*/

#define FLIP(aprev, a, b, bnext, f, x) {                                   \
    CClinkern_flipper_flip ((x),(a), (b));                                 \
}
    // (f)->stack[(f)->counter].first = (a);                                  \
    // (f)->stack[(f)->counter++].last = (b);                                 \
// }

#define UNFLIP(aprev, a, b, bnext, f, x) {                                 \
    CClinkern_flipper_flip ((x), (b), (a));                                \
    (f)->counter--;                                                        \
}

// #define IMPROVE_PACKING(instance, sol, main, regressor, capacity){                    \
//     if (randUnif(0,1)< ttp_lk_params.probImprove){                                                     \
//         ttpObjective orig = evaluate(instance, sol);                            \
//         calcRDist(instance, sol->tour);\
//         double startTime = CCutil_zeit();                                       \
//         quadraticHillClimber(instance, sol, regressor, QUADRATIC_MAX_NONIMPROVEMENT, QUADRATIC_MAX, capacity);       \
//         ttpObjective postLinear = evaluate(instance, sol);                      \
//         linearTime+=CCutil_zeit()-startTime;                                    \
//         startTime = CCutil_zeit();                                              \
//         run_mbfs(main, sol, MBFS_MAX_ITER, true);         \
//         ttpObjective postClassifier = evaluate(instance, sol);                   \
//         if (postClassifier.objective>postLinear.objective) ++classifierImprovements;\
//         if (postLinear.objective>orig.objective) ++linearImprovements;\
//         hillClimberTime+=CCutil_zeit()-startTime;                               \
//         numHillClimbers++;                                                      \
//     }                                                                           \
// }      

#define IMPROVE_PACKING(instance, sol, main, regressor, capacity){                    \
    if (randUnif(0,1)< ttp_lk_params.probImprove){                                                     \
        calcRDist(instance, sol->tour);\
        initializedQuadraticClimber(instance, sol, regressor, ttp_lk_params.quadratic_max_nonimprovement, ttp_lk_params.quadratic_max, capacity);       \
        run_mbfs(main, sol, ttp_lk_params.mbfs_max, -1, false);         \
    }                                                                           \
}                                                                               

#ifdef USE_HEAP
#define MARK(xn, xQ, xF, xD, xG, xW)  turn ((xn), (xQ), (xF), (xD), (xG))
#else
#ifdef USE_LESS_MARKING
#define MARK(xn, xQ, xF, xD, xG, xW)  turn ((xn), (xQ), (xW))
#else
#define MARK(xn, xQ, xF, xD, xG, xW)  turn ((xn), (xQ), (xF), (xW))
#endif
#endif

#define markedge_add(n1, n2, E)    E->add_edges[n1 ^ n2] = 1
#define markedge_del(n1, n2, E)    E->del_edges[n1 ^ n2] = 1
#define unmarkedge_add(n1, n2, E)  E->add_edges[n1 ^ n2] = 0
#define unmarkedge_del(n1, n2, E)  E->del_edges[n1 ^ n2] = 0
#define is_it_added(n1, n2, E)     E->add_edges[n1 ^ n2]
#define is_it_deleted(n1, n2, E)   E->del_edges[n1 ^ n2]

static void
   ttp_lin_kernighan (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       double *val, int *win_cycle, flipstack *w, flipstack *fstack,
       CCptrworld *intptr_world, CCptrworld *edgelook_world, ttpInstance *instance, ttpSolution *sol, 
       mbfs_main *main, logisticRegressor *classifier, double *objective, double *capacity, int mode),
   look_ahead_noback (graph *G, distobj *D, adddel *E, CClk_flipper *F,
       int first, int last, int gain, edgelook *winner),
#ifdef USE_HEAP
   turn (int n, aqueue *Q, CClk_flipper *F, distobj *D, graph *G),
#else
#ifdef USE_LESS_MARKING
   turn (int n, aqueue *Q, CCptrworld *intptr_world),
#else
   turn (int n, aqueue *Q, CClk_flipper *F, CCptrworld *intptr_world),
#endif
#endif
   kickturn (int n, aqueue *Q, distobj *D, graph *G, CClk_flipper *F,
        CCptrworld *intptr_world),
   bigturn (graph *G, int n, int tonext, aqueue *Q, CClk_flipper *F,
        distobj *D, CCptrworld *intptr_world),
   first_kicker (graph *G, distobj *D, CClk_flipper *F, int *t1, int *t2),
   find_random_four (graph *G, distobj *D, CClk_flipper *F, int *t1, int *t2,
       int *t3, int *t4, int *t5, int *t6, int *t7, int *t8),
   find_close_four (graph *G, distobj *D, CClk_flipper *F, int *t1, int *t2,
       int *t3, int *t4, int *t5, int *t6, int *t7, int *t8),
   find_walk_four (graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8),
   randcycle (int ncount, int *cyc, CCrandstate *rstate),
   insertedge (graph *G, int n1, int n2, int w),
   initgraph (graph *G),
   freegraph (graph *G),
   init_adddel (adddel *E),
   free_adddel (adddel *E),
   init_aqueue (aqueue *Q),
   free_aqueue (aqueue *Q, CCptrworld *intptr_world),
#ifdef USE_HEAP
   add_to_active_queue (int n, aqueue *Q, distobj *D, graph *G, CClk_flipper *F),
#else
   add_to_active_queue (int n, aqueue *Q, CCptrworld *intptr_world),
#endif
   init_distobj (distobj *D),
   free_distobj (distobj *D),
   linkern_free_world (CCptrworld *intptr_world, CCptrworld *edgelook_world),
   free_flipstack (flipstack *f);

static int
   buildgraph (graph *G, int ncount, int ecount, int *elist, distobj *D),
   ttp_repeated_lin_kernighan (graph *G, distobj *D, int *cyc,
       int stallcount, int count, double *val, double time_bound,
       double length_bound,  char *saveit_name, int silent, int kicktype,
       CCptrworld *intptr_world, CCptrworld *edgelook_world,
        CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol, mbfs_main *main, CCdatagroup *dat, 
        logisticRegressor *classifier, double *objective, double *capacity, int mode), 
   ttp_weird_second_step (graph *G, distobj *D, adddel *E, aqueue *Q,
       CClk_flipper *F, int gain, int t1, int t2, flipstack *fstack,
       CCptrworld *intptr_world, CCptrworld *edgelook_world, 
       ttpInstance *instance, ttpSolution *solution, mbfs_main *main, logisticRegressor *classifier, double prevObj, double *capacity),
   ttp_step (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       int level, int gain, int *Gstar, int first, int last, flipstack *fstack,
       CCptrworld *intptr_world, CCptrworld *edgelook_world, ttpInstance *instance, 
       ttpSolution *solution, mbfs_main *main, logisticRegressor *classifier, double prevObj, double *capacity, int *Gmax),
   ttp_step_noback (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       int level, int gain, int *Gstar, int first, int last,
       flipstack *fstack, CCptrworld *intptr_world, ttpInstance *instance, ttpSolution *solution, double prevObj, int *Gmax),
   kick_step_noback (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       int level, int gain, int *Gstar, int first, int last, flipstack *win,
       flipstack *fstack, CCptrworld *intptr_world),
   find_geometric_four (graph *G, distobj *D, CClk_flipper *F, CCkdtree *kdt,
       int *t1, int *t2, int *t3, int *t4,
       int *t5, int *t6, int *t7, int *t8, CCrandstate *rstate),
   random_four_swap (graph *G, distobj *D, aqueue *Q, CClk_flipper *F,
       CCkdtree *kdt, int *delta, int kicktype, flipstack *win,
       flipstack *fstack, CCptrworld *intptr_world, CCrandstate *rstate),
   save_tour (int ncount, char *sname, CClk_flipper *F),
   build_adddel (adddel *E, int ncount),
   build_aqueue (aqueue *Q, int ncount, CCptrworld *intptr_world),
   pop_from_active_queue (aqueue *Q, CCptrworld *intptr_world),
   build_distobj (distobj *D, int ncount, CCdatagroup *dat),
   dist (int i, int j, distobj *D), 
   init_flipstack (flipstack *f, int total, int single);

static double
   ttp_improve_tour (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       int start, flipstack *fstack, CCptrworld *intptr_world,
       CCptrworld *edgelook_world, ttpInstance *instance, ttpSolution *sol, mbfs_main *main,
       logisticRegressor *classifier, double *capacity, double bestObj),
   kick_improve (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       flipstack *win, flipstack *fstack, CCptrworld *intptr_world),
   cycle_length (int ncount, int *cyc, distobj *D);

static edgelook
   *look_ahead (graph *G, distobj *D, adddel *E, CClk_flipper *F, int first,
       int last, int gain, int level, CCptrworld *edgelook_world),
   *weird_look_ahead  (graph *G, distobj *D, CClk_flipper *F, int gain, int t1,
       int t2, CCptrworld *edgelook_world),
   *weird_look_ahead2 (graph *G, distobj *D, CClk_flipper *F, int gain, int t2,
       int t3, int t4, CCptrworld *edgelook_world),
   *weird_look_ahead3 (graph *G, distobj *D, CClk_flipper *F, int gain, int t2,
       int t3, int t6, CCptrworld *edgelook_world);

ttpObjective evaluateFromFlipper(ttpInstance *instance, ttpSolution *sol, CClk_flipper *F), 
        evaluateFromFlipperIndices(ttpInstance *instance, ttpSolution *sol, CClk_flipper *F, CClk_childnode *start, CClk_childnode *stop, double startWeight);

CC_PTRWORLD_ROUTINES(intptr, intptralloc, intptr_bulkalloc, intptrfree)
CC_PTRWORLD_LISTFREE_ROUTINE(intptr, intptr_listfree, intptrfree)
CC_PTRWORLD_LEAKS_ROUTINE(intptr, intptr_check_leaks, this, int)

CC_PTRWORLD_ROUTINES(edgelook, edgelookalloc, edgelook_bulkalloc, edgelookfree)
CC_PTRWORLD_LISTFREE_ROUTINE(edgelook, edgelook_listfree, edgelookfree)
CC_PTRWORLD_LEAKS_ROUTINE(edgelook, edgelook_check_leaks, diff, int)




int CCttp_linkern_tour (int ncount, CCdatagroup *dat, int ecount,
        int *elist, int stallcount, int repeatcount, int *incycle, double *val,
        int silent, double time_bound, double length_bound,
        char *saveit_name, int kicktype, CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol, 
        mbfs_main *main, logisticRegressor *classifier, double *objective, double *capacity, int mode)
{
    int *outcycle = sol->tour;
    int rval = 0;
    int i;
    int *tcyc = (int *) NULL;
    double startzeit;
    graph G;
    distobj D;
    CCptrworld intptr_world;
    CCptrworld edgelook_world;

    if (silent == 0) {
        printf ("linkern ...\n"); fflush (stdout);
    }
    startzeit = CCutil_zeit ();

    initgraph (&G);
    init_distobj (&D);
    CCptrworld_init (&intptr_world);
    CCptrworld_init (&edgelook_world);
    G.rstate = rstate;

    if (ncount < 10 && repeatcount > 0) {
        printf ("Less than 10 nodes, setting repeatcount to 0\n");  
        fflush (stdout);
        repeatcount = 0;
    }

    if (((dat->norm) & CC_NORM_BITS) != CC_KD_NORM_TYPE) {
        if (kicktype == CC_LK_GEOMETRIC_KICK) {
            if (silent == 0) {
                printf ("Setting kick type to close\n"); fflush (stdout);
            }
            kicktype = CC_LK_CLOSE_KICK; 
        }
    }


    /* These bulkalloc's allocate sufficient objects that the individual
     * allocs will not fail, and thus do not need to be tested */
    rval = intptr_bulkalloc (&intptr_world, ncount);
    if (rval) {
        fprintf (stderr, "Unable to allocate initial intptrs\n");
        goto CLEANUP;
    }

    rval = edgelook_bulkalloc (&edgelook_world, MAX_BACK * (BACKTRACK + 3));
    if (rval) {
        fprintf (stderr, "Unable to allocate initial edgelooks\n");
        goto CLEANUP;
    }

    tcyc = CC_SAFE_MALLOC (ncount, int);
    if (tcyc == (int *) NULL) {
        fprintf (stderr, "out of memory in linkern\n");
        rval = 1; goto CLEANUP;
    }

    rval = build_distobj (&D, ncount, dat);
    if (rval) goto CLEANUP;
    
    rval = buildgraph (&G, ncount, ecount, elist, &D);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    if (incycle) {
        for (i = 0; i < ncount; i++) tcyc[i] = incycle[i];
    } else {
        randcycle (ncount, tcyc, G.rstate);
    }
    *val = cycle_length (ncount, sol->tour, &D);
    if (silent == 0) {
        printf ("Starting Cycle: %.0f\n", *val); fflush (stdout);
    }


    rval = ttp_repeated_lin_kernighan (&G, &D, tcyc, stallcount, repeatcount,
                 val, time_bound, length_bound, saveit_name, silent,
                 kicktype, &intptr_world, &edgelook_world, rstate, instance, sol, main, dat, classifier, 
                 objective, capacity, mode);

    
    if (rval) {
        fprintf (stderr, "repeated_lin_kernighan failed\n"); goto CLEANUP;
    }

    if (silent == 0) {
        printf ("Best cycle length: %.0f\n", *val);
        printf ("Lin-Kernighan Running Time: %.2f\n",
                  CCutil_zeit () - startzeit);
        fflush (stdout);
    }

    if (outcycle) {
        for (i = 0; i < ncount; i++) outcycle[i] = tcyc[i];
    }

CLEANUP:

    CC_IFFREE (tcyc, int);
    freegraph (&G);
    free_distobj (&D);
    linkern_free_world (&intptr_world, &edgelook_world);

    return rval;
}

#ifdef ACCEPT_BAD_TOURS
#define HEAT_FACTOR 0.999
#define HEAT_RESET 100000
#endif

ttpObjective do_pgch(mbfs_main *main, ttpInstance *instance, ttpSolution *sol, logisticRegressor *classifier, double timebound, double *capacity, int maxIters){
    double starttime = CCutil_zeit();
    calcRDist(instance, sol->tour);
    initializedQuadraticClimber(instance, sol, classifier, ttp_lk_params.upper_quadratic_max_nonimprovement, ttp_lk_params.upper_quadratic_max, capacity); 
    
    ttpObjective finalObj = pgch(main, sol, maxIters);
    double prevObj = finalObj.objective;
    
    finalObj = run_mbfs(main, sol, ttp_lk_params.upper_mbfs_max, timebound-(CCutil_zeit()-starttime), ttp_lk_params.upper_mbfs_sorted); 
    bool improve_tour = true;
    while (finalObj.objective>prevObj){
        if (timebound>=0 && (CCutil_zeit()-starttime>timebound)) break;
        prevObj = finalObj.objective;
        if (improve_tour) finalObj = pgch(main, sol, maxIters);
        else finalObj = run_mbfs(main, sol, ttp_lk_params.upper_mbfs_max, timebound-(CCutil_zeit()-starttime), ttp_lk_params.upper_mbfs_sorted); 
        improve_tour = !improve_tour;
    }
    return finalObj;
}

static int ttp_repeated_lin_kernighan (graph *G, distobj *D, int *cyc,
        int stallcount, int count, double *val, double time_bound,
        double length_bound, char *saveit_name, int silent, int kicktype,
        CCptrworld *intptr_world, CCptrworld *edgelook_world,
        CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol, mbfs_main *main, CCdatagroup *dat, 
        logisticRegressor *classifier, double *objective, double *capacity, int mode)
{
    double prevprinttime = CCutil_zeit();
    int rval    = 0;
    int round   = 0;
    int newtree = 0;
    int quitcount, hit, delta;
    int *win_cycle = (int *) NULL;
    CCkdtree kdt;
    flipstack winstack, fstack;
    double t, best = *val, oldbest = *val, bestObj = *objective, oldBestObj = *objective;
    double szeit = CCutil_zeit ();
#ifdef ACCEPT_BAD_TOURS
    double heat = *val / (20 * G->ncount), tdelta;
#endif
    int ncount = G->ncount;
    adddel E;
    CClk_flipper F;
    aqueue Q;

    init_aqueue (&Q);
    init_adddel (&E);
    rval = build_aqueue (&Q, ncount, intptr_world);
    if (rval) {
        fprintf (stderr, "build_aqueue failed\n"); goto CLEANUP;
    }
    rval = build_adddel (&E, ncount);
    if (rval) {
        fprintf (stderr, "build_adddel failed\n"); goto CLEANUP;
    }

    hit = 2 * (MAXDEPTH + 7 + KICK_MAXDEPTH);
    rval = init_flipstack (&fstack, hit, 0);
    if (rval) {
        fprintf (stderr, "init_flipstack failed\n"); goto CLEANUP;
    }
    rval = init_flipstack (&winstack, 500 + ncount / 50, hit);
    if (rval) {
        fprintf (stderr, "init_flipstack failed\n"); goto CLEANUP;
    }

    win_cycle = CC_SAFE_MALLOC (ncount, int);
    if (win_cycle == (int *) NULL) {
        fprintf (stderr, "out of memory in repeated_lin_kernighan\n");
        rval = 1; goto CLEANUP;
    }
    win_cycle[0] = -1;

    quitcount = stallcount;
    if (quitcount > count) quitcount = count;

    CClinkern_flipper_init (&F, ncount, cyc);
    fstack.counter = 0;
    winstack.counter = 0;
    win_cycle[0] = -1;

#ifdef USE_HEAP
    {
        int i;

        for (i = 0; i < ncount; i++) {
            add_to_active_queue (i, &Q, D, G, &F); 
        }
    }
#else
    {
        int *tcyc = (int *) NULL;
        int i;

        tcyc = CC_SAFE_MALLOC (ncount, int);
        if (tcyc == (int *) NULL) {
            fprintf (stderr, "out of memory in repeated_lin_kernighan\n");
            rval = 1; goto CLEANUP;
        }
        /* init active_queue with random order */
        randcycle (ncount, tcyc, G->rstate);
        for (i = 0; i < ncount; i++) {
            add_to_active_queue (tcyc[i], &Q, intptr_world);
        }
        CC_IFFREE (tcyc, int);
    }
#endif

    if (kicktype == CC_LK_GEOMETRIC_KICK) {
        rval = CCkdtree_build (&kdt, ncount, D->dat, (double *) NULL, rstate);
        if (rval) {
            fprintf (stderr, "CCkdtree_build failed\n"); goto CLEANUP;
        } else {
            newtree = 1;
        }
    }

    winstack.counter = 0;
    win_cycle[0] = -1;

    if (silent == 0) {
        if (quitcount > 0) {
            printf ("%4d Steps   Best: %.0f   %.2f seconds\n", round, best,
                                CCutil_zeit () - szeit);
        } else {
            printf ("LK Cycle: %.0f\n", best);
        }
        fflush (stdout);
    }
    ttpSolution newSol;
    newSol.packing = CC_SAFE_MALLOC(instance->numberOfItems, bool);
    newSol.tour = CC_SAFE_MALLOC(instance->numberOfNodes+1, int);
    newSol.isCityPacking = sol->isCityPacking;
    newSol.tour[instance->numberOfNodes]=0;

    int normal_nonimprovement_counter=0;
    bool improved = false;

    while (round < quitcount) {
        // Calculate the round
        int total_rounds = ttp_lk_params.normal_rounds+ttp_lk_params.ttp_rounds+ttp_lk_params.pgch_rounds;
        int round_type = -1;
        int pgch_rounds = ttp_lk_params.pgch_rounds;
        int ttp_rounds = ttp_lk_params.ttp_rounds;
        int normal_rounds = ttp_lk_params.normal_rounds;
        int round_index = round % total_rounds;
        if (round_index<pgch_rounds) round_type = 2;
        else if (round_index<pgch_rounds+ttp_rounds) round_type=1;
        else round_type=0;
        ////////////////
        // printf("ROUND TYPE %d\n",round_type);
        
// #ifdef CHECK_OBJECTIVE
        for (int i=0;i<instance->numberOfItems;i++) newSol.packing[i] = sol->packing[i];
        for (int i=0;i<instance->numberOfNodes;i++) newSol.tour[i] = sol->tour[i];
        CClinkern_flipper_cycle(&F, newSol.tour);
// #endif
        hit = 0;
        fstack.counter = 0;

        // Perturb and improve tour
        if (round_type!=2){
            if (IMPROVE_SWITCH == -1 || round < IMPROVE_SWITCH) {
                rval = random_four_swap (G, D, &Q, &F, &kdt, &delta, ttp_lk_params.kicktype,
                                        &winstack, &fstack, intptr_world, rstate);
                if (round_type==1){
                    calcRDist(instance, newSol.tour);
                    initializedQuadraticClimber(instance, &newSol, classifier, ttp_lk_params.quadratic_max_nonimprovement, ttp_lk_params.quadratic_max, capacity);                          
                    run_mbfs(main, &newSol, ttp_lk_params.mbfs_max, time_bound-(CCutil_zeit()-szeit), ttp_lk_params.mbfs_sorted); 
                }
                
                if (rval) {
                    fprintf (stderr, "random_four_swap failed\n"); goto CLEANUP;
                }
            } else {
                delta = kick_improve (G, D, &E, &Q, &F, &winstack, &fstack, intptr_world);
            }
        }
        fstack.counter = 0;
        t = best + delta;

#ifdef CHECK_OBJECTIVE
        ttpObjective initObj = evaluateFromFlipper(instance, &newSol, &F);
        double tObj = initObj.objective;
#else 
        double tObj=0;
#endif
        ttpObjective finalObj;
        if (round_type==2){
            int numIters = NUMLKITERS(instance->numberOfNodes);
            for (int i=0;i<instance->numberOfNodes+1;i++) newSol.tour[i]=0;
            for (int i=0;i<instance->numberOfItems;i++) newSol.packing[i]=false;
            lk_tour(dat, rstate, instance, &newSol, NULL, NULL, numIters, 0, 20);

            finalObj = do_pgch(main, instance, &newSol, classifier, time_bound-(CCutil_zeit()-szeit), capacity, -1);
            ttp_lk_params.pgch_round_counter++;
            if (finalObj.objective<=bestObj) ttp_lk_params.pgch_nonimprovement_counter++;
            // printf("pgch objective %lf\n",finalObj.objective);

        }
        else if (round_type==1){
            ttp_lin_kernighan (G, D, &E, &Q, &F, &t, win_cycle, &winstack, &fstack,
                       intptr_world, edgelook_world, instance, &newSol, main, classifier, 
                       capacity, &tObj, mode);
            calcRDist(instance, newSol.tour);
            initializedQuadraticClimber(instance, &newSol, classifier, ttp_lk_params.upper_quadratic_max_nonimprovement, ttp_lk_params.upper_quadratic_max, capacity); 
            finalObj = run_mbfs(main, &newSol, ttp_lk_params.upper_mbfs_max, time_bound-(CCutil_zeit()-szeit), ttp_lk_params.upper_mbfs_sorted); 
        } else {
            lin_kernighan (G, D, &E, &Q, &F, &t, win_cycle, &winstack, &fstack,
                    intptr_world, edgelook_world);

            CClinkern_flipper_cycle(&F, newSol.tour);
            newSol.tour[instance->numberOfNodes]=0;
            
            if (randUnif(0,1)<ttp_lk_params.upper_improvement_prob){
                if (finalObj.objective<bestObj){
                    copyInt(sol->tour,newSol.tour,instance->numberOfNodes+1);
                    copyBool(sol->packing,newSol.packing,instance->numberOfItems);
                }
                calcRDist(instance, newSol.tour);
                initializedQuadraticClimber(instance, &newSol, classifier, ttp_lk_params.upper_quadratic_max_nonimprovement, ttp_lk_params.upper_quadratic_max, capacity); 
                finalObj = run_mbfs(main, &newSol, ttp_lk_params.upper_mbfs_max, time_bound-(CCutil_zeit()-szeit), false); //ttp_lk_params.upper_mbfs_sorted
            } else {
                finalObj = evaluate(instance, &newSol);
                if (finalObj.objective>bestObj){
                    calcRDist(instance, newSol.tour);
                    initializedQuadraticClimber(instance, &newSol, classifier, ttp_lk_params.upper_quadratic_max_nonimprovement, ttp_lk_params.upper_quadratic_max, capacity); 
                    finalObj = run_mbfs(main, &newSol, ttp_lk_params.upper_mbfs_max, time_bound-(CCutil_zeit()-szeit), false); //ttp_lk_params.upper_mbfs_sorted
                }
            }
            if (finalObj.objective<=bestObj) ttp_lk_params.normal_nonimprovement_counter++;
        }
        
        tObj = finalObj.objective;

#ifdef ACCEPT_BAD_TOURS
        if (round % HEAT_RESET == HEAT_RESET - 1) {
            heat = oldbest / (20 * ncount);
            printf ("Reset Accept-Probablility\n");
            fflush (stdout);
        }
        tdelta = t - best;
        heat *= HEAT_FACTOR;
        if (t < best || (t > best && exp (-tdelta/heat) >
            (double) (CCutil_lprand (G->rstate) % ncount) / (double) ncount)) {
#else
#ifdef ACCEPT_TIES
    #ifdef CHECK_OBJECTIVE
        if (tObj >= bestObj) {
    #else
        // if (t <= best) {
            if (finalObj.objective>=bestObj){
    #endif
        // printf("round %d tObj %lf bestObj %lf\n", round, tObj, bestObjNormal.objective);
        
#else
        if (t < best) {
#endif /* ACCEPT_TIES */
#endif /* ACCEPT_BAD_TOURS */
            winstack.counter = 0;
            win_cycle[0] = -1;
            // Allow neutral mutations
            if (round_type==2){
                copyInt(newSol.tour,sol->tour,instance->numberOfNodes);
            }
            else if (round_type==1){
                CClinkern_flipper_cycle(&F,sol->tour);
            } else {
                copyInt(newSol.tour,sol->tour,instance->numberOfNodes+1);
            }
            copyBool(newSol.packing,sol->packing, instance->numberOfItems);
            best = t;
    #ifdef CHECK_OBJECTIVE
        if (tObj > bestObj) {
    #else 
        // if (t < best) {
            if (finalObj.objective > bestObj){
                
                printf("objective increased from %lf to %lf w/ dist %lf\n", bestObj, finalObj.objective, finalObj.distance);
                CClinkern_flipper_cycle(&F,sol->tour);
                bestObj=finalObj.objective;
    #endif
            best = t;
            improved = true;
                
#ifdef CHECK_OBJECTIVE
                printf("setting best to t: best = %lf t = %lf obj = %lf\n", best, t, tObj);
                bestObj = tObj;
#endif
                

                quitcount = round + stallcount;
                if (quitcount > count)
                    quitcount = count;
                hit++;
            }
#ifdef ACCEPT_BAD_TOURS
            else {
                if (silent == 0 && t > best) {
printf ("%4d Steps   Best: %.0f   %.2f seconds (Negative %.0f) (%.0f)\n",
                          round, t, CCutil_zeit () - szeit, t - best, oldbest);
                    fflush (stdout);
                }
                oldbest = best;
                best = t;
            }
#endif
        } else {
            win_cycle[0]=-1;
            winstack.counter=0;
            CClinkern_flipper_finish(&F);
            CClinkern_flipper_init(&F, ncount, sol->tour);
        }
        round++;

        double prevremtime = time_bound-prevprinttime;
        double remtime = time_bound-CCutil_zeit();
        if (improved && remtime<0.5*prevremtime){
            ttpOutput(instance, sol, global_filesuffix, "LKSR", NULL);
            prevprinttime = CCutil_zeit();
            improved = false;
        }
        
        if (silent == 0 && (hit || (round % 1000 == 999))) {
            printf ("%4d Steps   Best: %.3f   %.2f seconds\n",
                               round, bestObj, CCutil_zeit () - szeit);
            fflush(stdout);
        }
        
        // if (silent==0 && (improved && (round%1000 == 999))){
        //     ttpOutput(instance, sol, "", "", stdout);
        //     fflush (stdout);
        //     improved = false;
        // }

        // if (saveit_name && (round % 10000 == 9999) && best < oldbest) {
        //     rval = save_tour (ncount, saveit_name, &F);
        //     if (rval) {
        //         fprintf (stderr, "save_tour failed\n"); goto CLEANUP;
        //     }
        //     oldbest = best;
        // }
          
        double time = CCutil_zeit() - szeit;
        if (time_bound > 0.0 && time > time_bound) {
            printf ("STOP - timebound (%.2f seconds)\n", CCutil_zeit() - szeit);
            if (silent == 1) {
                printf ("STEPS: %d\n", round);
            }
            fflush (stdout);
            break;
        }
        set_params(instance, &ttp_lk_params, time);
        // if (length_bound > 0.0 && best <= length_bound) {
        //     printf ("STOP - length bound reached (%.0f)\n", length_bound);
        //     if (silent == 1) {
        //         printf ("STEPS: %d\n", round);
        //     }
        //     fflush (stdout);
        //     break;
        // }
        //only one round//
        // break;
    }
    if (silent == 0 && round > 0) {
        printf ("%4d Total Steps.\n", round); fflush (stdout);
    }

    CClinkern_flipper_cycle (&F, cyc);
    CClinkern_flipper_finish (&F);

//     if (saveit_name && best < oldbest) {
// /*
//         rval = CCutil_writecycle_edgelist (ncount, saveit_name, cyc,
//                                            D->dat, 0);
//         if (rval) {
//             fprintf (stderr, "could not write the cycle\n"); goto CLEANUP;
//         }
// */
//         rval = CCutil_writecycle (ncount, saveit_name, cyc, 0);
//         CCcheck_rval (rval, "CCutil_writecycle failed");

//         printf ("Wrote the last tour to %s\n", saveit_name); fflush (stdout);
//     }

    // t = cycle_length (ncount, cyc, D);
    // if (t != best) {
    //     printf ("WARNING: LK incremental counter was off by %.0f t: %.0f best: %.0f\n", t-best, t, best);
    //     fflush (stdout);
    //     best = t;
    // }
    
    *val = best;
    *objective = bestObj;

CLEANUP:

    free_aqueue (&Q, intptr_world);
    free_adddel (&E);
    free_flipstack (&fstack);
    free_flipstack (&winstack);
    CC_IFFREE (win_cycle, int);
    freeSolution(&newSol);
    if (newtree) CCkdtree_free (&kdt);
    return rval;
}


ttpObjective evaluateFromFlipperIndices(ttpInstance *instance, ttpSolution *sol, CClk_flipper *F, CClk_childnode *start, CClk_childnode *stop, double startWeight){
    double weight = startWeight;
    double time = 0;
    double distance = 0;
    double profit = 0;

    CClk_childnode *c = start;
    do {
        if (c!=start){
            for (int j=0;j<instance->itemsPerCity;j++){
                int itemIdx = getItemIndex(instance, c->name,j);
                if (sol->packing[itemIdx]) {
                    weight+=instance->weight[itemIdx];
                    profit+=instance->profit[itemIdx];
                    // if (weight>instance->capacityOfKnapsack){
                    //     printf("excessive weight\n");
                    //     exit(0);
                    // }
                }
                // printf("flipper eval city %d, item number %d, packing index %d item idx %d weight %lf profit %lf packed %d\n",c->name,j,itemIdx, itemIdx, weight, profit, sol->packing[itemIdx]);
            }
        }
        CClk_childnode *next = c->adj[!((F->reversed)^(c->parent->rev))];
        double dist = getDistance(instance, c->name, next->name);
        distance+=dist;
        double speed = getSpeed(instance, weight);
        speed = weight>instance->capacityOfKnapsack ? instance->minSpeed/1E20 : speed;
        time+=dist/speed;
        c = next;
    } while (c!=stop);
    
    ttpObjective ob;
    ob.profit = profit;
    ob.time = time;
    ob.finalweight = weight;
    ob.distance=distance;
    ob.objective=profit-time*instance->rentingRatio;
    return ob;
}

ttpObjective evaluateFromFlipper(ttpInstance *instance, ttpSolution *sol, CClk_flipper *F){
    return evaluateFromFlipperIndices(instance, sol, F, &(F->children[0]), &(F->children[0]), 0);
}

static void ttp_lin_kernighan (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, double *val, int *win_cycle, flipstack *win,
        flipstack *fstack, CCptrworld *intptr_world,
        CCptrworld *edgelook_world, ttpInstance *instance, ttpSolution *sol, mbfs_main *main, logisticRegressor *classifier, 
        double *capacity, double *objective, int mode)
{
    int start, i;
    double delta, totalwin = 0.0;
    double bestObj = *objective;

    while (1) {
        start = pop_from_active_queue (Q, intptr_world);
        if (start == -1) break;
        

        //based on delta's value, decide whether to check and enforce better objective.
        
        // printf("total win %lf delta %lf delta(distances) %lf objective %lf dist %lf totalwindist %lf best objective %lf\n", totalwin+delta, delta, prevobj.distance-obj.distance, obj.objective, obj.distance, *val-totalwin, bestObj);
        // printf("total win %lf delta %lf delta(distances) %lf dist %lf val %lf totalwindist %lf\n", totalwin+delta, delta, prevobj.distance-obj.distance, obj.distance,*val, *val-totalwin);
#ifdef CHECK_OBJECTIVE

        ttpObjective prevobj = evaluateFromFlipper(instance, sol, F);
        delta = ttp_improve_tour (G, D, E, Q, F, start, fstack, intptr_world,
                              edgelook_world, instance, sol, main, classifier, capacity, prevobj.objective);
        ttpObjective obj = evaluateFromFlipper(instance, sol, F);
        if (obj.objective>bestObj) {
            bestObj = obj.objective > bestObj ? obj.objective : bestObj;
#else
        delta = ttp_improve_tour (G, D, E, Q, F, start, fstack, intptr_world,
                              edgelook_world, instance, sol, classifier, linear, percent, 0);
        if (delta > 0){
#endif
            totalwin += delta;

            //
            // if (win->counter < win->max) {
            //     for (i = 0; i < fstack->counter; i++) {
            //         win->stack[win->counter].first = fstack->stack[i].first;
            //         win->stack[win->counter].last  = fstack->stack[i].last;
            //         win->stack[win->counter].firstprev =
            //                  fstack->stack[i].firstprev;
            //         win->stack[win->counter].lastnext =
            //                  fstack->stack[i].lastnext;
            //         win->counter++;
            //     }
            // } else if (win_cycle[0] == -1) {
            //     for (i = 0; i < fstack->counter; i++) {
            //         win->stack[win->counter].first = fstack->stack[i].first;
            //         win->stack[win->counter].last  = fstack->stack[i].last;
            //         win->counter++;
            //     }
            //     CClinkern_flipper_cycle (F, win_cycle);
            // }
            // fstack->counter = 0;
            //
        }
    }

    if (win_cycle[0] == -1) {
        for (i = 0; i < fstack->counter; i++) {
            win->stack[win->counter].first = fstack->stack[i].first;
            win->stack[win->counter].last  = fstack->stack[i].last;
            win->counter++;
        }
    }

    (*val) -= totalwin;
    (*objective) = bestObj;
}

ttpObjective test_sol(ttpInstance *instance, ttpSolution *sol, CClk_flipper *F, mbfs_main *main, logisticRegressor *classifier, double *capacity, double bestObj){
    ttpSolution copy;
    copy.isCityPacking = sol->isCityPacking;
    copy.tour = CC_SAFE_MALLOC(instance->numberOfNodes+1,int);
    copy.packing = CC_SAFE_MALLOC(instance->numberOfItems,bool);
    CClinkern_flipper_cycle (F, copy.tour);
    copy.tour[instance->numberOfNodes]=0;
    for (int i=0;i<instance->numberOfItems;i++) copy.packing[i]=sol->packing[i];
    IMPROVE_PACKING(instance, (&copy), main, classifier, capacity);
    ttpObjective finalObj = evaluate(instance, &copy);
    if (finalObj.objective>=bestObj){
        for (int i=0;i<instance->numberOfItems;i++) sol->packing[i]=copy.packing[i];
        for (int i=0;i<=instance->numberOfNodes;i++) sol->tour[i]=copy.tour[i];
    }
    CC_IFFREE(copy.packing, bool);
    CC_IFFREE(copy.tour,int);
    return finalObj;
}


static double ttp_improve_tour (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int t1, flipstack *fstack,
        CCptrworld *intptr_world, CCptrworld *edgelook_world, ttpInstance *instance, ttpSolution *sol, mbfs_main *main,
        logisticRegressor *classifier, double *capacity, double bestObj)
{
    // printf("improve tour\n");
    int t2 = CClinkern_flipper_next (F, t1);
    int gain, Gstar = 0, Gmax = -(CCutil_MAXINT-1);

    gain = Edgelen (t1, t2, D);
    markedge_del (t1, t2, E);

    // printf("--------STEP INIT--------\n");
    double st = CCutil_zeit();
    if (ttp_step (G, D, E, Q, F, 0, gain, &Gstar, t1, t2, fstack, intptr_world, edgelook_world, instance, sol, main, classifier, bestObj, capacity, &Gmax)
        == 0) {
        Gmax = ttp_weird_second_step (G, D, E, Q, F, gain, t1, t2, fstack,
                                   intptr_world, edgelook_world, instance, sol, main, classifier, bestObj, capacity); 
    }
    unmarkedge_del (t1, t2, E);
    stepTime+=CCutil_zeit()-st;
    
#ifdef CHECK_OBJECTIVE
    ttpObjective obj = evaluateFromFlipper(instance, sol, F);
    // ttpObjective obj = test_sol(instance, sol, F, main, classifier, capacity, bestObj);

    if (obj.objective > bestObj) {
        // printf("better obj %lf beneficial pct %lf\n", obj.objective, beneficialOpts/(double)numChecks);
    CClinkern_flipper_cycle (F, sol->tour);
#else
    if (Gstar>0) {
#endif
        
        IMPROVE_PACKING(instance, sol, main, classifier, capacity);
        
#ifdef CHECK_OBJECTIVE
        ttpObjective finalObj = evaluateFromFlipper(instance, sol, F);
        if (finalObj.objective<=obj.objective) numZeroClimbers++;
#endif
        MARK (t1, Q, F, D, G, intptr_world);
        MARK (t2, Q, F, D, G, intptr_world);
    }

    // return (double) Gstar;
    return (double) Gmax;
}

static int ttp_step (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
        int level, int gain, int *Gstar, int first, int last,
        flipstack *fstack, CCptrworld *intptr_world, 
        CCptrworld *edgelook_world, ttpInstance *instance, ttpSolution *solution,
        mbfs_main *main, logisticRegressor *classifier, double prevObj, double *capacity, int *Gmax)
{
    
    int val, this, newlast, hit = 0, oldG = gain, maxHit =0;
#if defined(MAK_MORTON) && defined(FULL_MAK_MORTON)
    int newfirst;
#endif
    edgelook *list, *e;

    if (level >= BACKTRACK) {
        return ttp_step_noback (G, D, E, Q, F, level, gain, Gstar, first, last,
                            fstack, intptr_world, instance, solution, prevObj, Gmax);
    }

    list = look_ahead (G, D, E, F, first, last, gain, level, edgelook_world);
    for (e = list; e; e = e->next) {
#if defined(MAK_MORTON) && defined(FULL_MAK_MORTON)
        if (e->mm) {
            this = e->other;
            newfirst = e->over;

            gain = oldG - e->diff;
            val = gain - Edgelen (newfirst, last, D);
            if (val > *Gstar) {
                *Gstar = val;
                hit++;
            }
            FLIP (this, newfirst, first, last, fstack, F);

            if (level < MAXDEPTH) {
                markedge_add (first, this, E);
                markedge_del (this, newfirst, E);
                hit += step (G, D, E, Q, F, level + 1, gain, Gstar, newfirst,
                             last, fstack, intptr_world, edgelook_world);
                unmarkedge_add (first, this, E);
                unmarkedge_del (this, newfirst, E);
            }

            if (!hit) {
                UNFLIP (this, newfirst, first, last, fstack, F);
            } else {
                MARK (this, Q, F, D, G, intptr_world);
                MARK (newfirst, Q, F, D, G, intptr_world);
                edgelook_listfree (edgelook_world, list);
                return 1;
            }
        } else 
#endif
        {
            this = e->other;
            newlast = e->over;

            gain = oldG - e->diff;
            val = gain - Edgelen (newlast, first, D);
            

            if (val > *Gstar) {
                *Gstar = val;
                hit++;
            }
            if (val>*Gmax){
                *Gmax = val;
                maxHit++;
            }
            // printf("updated Gstar to %d and Gmax to %d\n", *Gstar, *Gmax);
            int prevGmax = *Gmax;
            int prevGstar = *Gstar;

            FLIP (first, last, newlast, this, fstack, F);

            int nextHit = 0;

            if (level < MAXDEPTH) {
                markedge_add (last, this, E);
                markedge_del (this, newlast, E);
                nextHit += ttp_step (G, D, E, Q, F, level + 1, gain, Gstar, first,
                             newlast, fstack, intptr_world, edgelook_world, instance, solution, main, classifier, prevObj, capacity, Gmax);
                unmarkedge_add (last, this, E);
                unmarkedge_del (this, newlast, E);
            }

            // ttpObjective tempobj = evaluateFromFlipper(instance, solution, F);
            // if (tempobj.objective>=prevObj && nextHit==0)
            // printf("val %d Gmax %d objective difference %lf \n", val, prevGstar, tempobj.objective-prevObj);
            if (maxHit==0 && hit>0) {
                printf("MAX HIT ZERO HIT NONZERO Gstar %d Gmax %d val %d\n", *Gstar, *Gmax, val);
                exit(0);
            }
#ifdef CHECK_OBJECTIVE
            // printf("gmax: %d\n", prevGmax);
            if (maxHit>0 && nextHit==0 && prevGmax>GMAX_THRESHOLD){
            //     ttpObjective obj = evaluateFromFlipper(instance, solution, F);
            //     printf("val %d Gmax %d Gstar %d objective difference %lf \n", val, prevGmax,prevGstar,  obj.objective-prevObj);
            // }
            // if (hit>0 && nextHit==0){
                double initTime = CCutil_zeit();
                ttpObjective obj = evaluateFromFlipper(instance, solution, F);


                // ttpSolution sol2;
                // sol2.packing = CC_SAFE_MALLOC(instance->numberOfItems,bool);
                // sol2.tour = CC_SAFE_MALLOC(instance->numberOfNodes+1, int);
                // sol2.tour[instance->numberOfNodes]=0;
                // CClinkern_flipper_cycle(F, sol2.tour);
                // copyBool(solution->packing,sol2.packing,instance->numberOfItems);
                // sol2.isCityPacking = solution->isCityPacking;
                // double linear=59871127.663049, percent=1;
                // initializedLinearClimber(instance, &sol2, 10, 10, &linear, &percent);
                // classifierHillClimber(instance, &sol2, classifier, 10, 100, 100);
                // ttpObjective obj2 = evaluateFromFlipper(instance, &sol2, F);


                checkTime+=CCutil_zeit()-initTime;
                ++numChecks;
                // printf("step Gmax %d obj = %lf, prevObj = %lf\n", prevGmax, obj.objective, prevObj);
                hit = 0;
                if (obj.objective>prevObj){
                    ++beneficialOpts;
                    hit = 1;
                    *Gstar = val;
                    *Gmax = val;
                }
            }
#endif
            hit+=nextHit;

            if (hit==0) {
                UNFLIP (first, last, newlast, this, fstack, F);
                if (level==0) {
                    *Gstar=0;
                    *Gmax=0;
                }
            } else {
                // printf("add to queue: %d, %d\n", this, newlast);
                MARK (this, Q, F, D, G, intptr_world);
                MARK (newlast, Q, F, D, G, intptr_world);
                edgelook_listfree (edgelook_world, list);
                return 1;
            }
        }
    }
    edgelook_listfree (edgelook_world, list);
    return 0;
}

static int ttp_step_noback (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *fstack, CCptrworld *intptr_world, ttpInstance *instance, ttpSolution *solution, double prevObj, int *Gmax)
{
    edgelook e;

#ifdef    SUBTRACT_GSTAR
#ifdef    SWITCH_LATE
    if (level < LATE_DEPTH) {
        look_ahead_noback (G, D, E, F, first, last, gain - *Gstar, &e);
    } else {
        look_ahead_noback (G, D, E, F, first, last, gain - *Gstar - level, &e);
    }
#else
    look_ahead_noback (G, D, E, F, first, last, gain - *Gstar - level, &e);
#endif /* SWITCH_LATE */
#else
#ifdef    SWITCH_LATE
    if (level < LATE_DEPTH) {
        look_ahead_noback (G, D, E, F, first, last, gain, &e);
    } else {
        look_ahead_noback (G, D, E, F, first, last, gain - level, &e);
    }
#else
    look_ahead_noback (G, D, E, F, first, last, gain - level, &e);
#endif /* SWITCH_LATE */
#endif /* SUBTRACT_GSTAR */

    if (e.diff < BIGINT) {
#ifdef NODE_INSERTIONS
        if (e.ni) {
            int hit = 0;
            int newlast = e.other;
            int next = e.under;
            int prev = e.over;
            int val;

            gain -= e.diff;
            val = gain - Edgelen (newlast, first, D);

            if (val > *Gstar) {
                *Gstar = val;
                hit++;
            }

            FLIP (first, last, newlast, next, fstack, F);
            FLIP (newlast, prev, last, next, fstack, F);

            if (level < MAXDEPTH) {
                markedge_add (last, newlast, E);
                markedge_add (next, prev, E);
                markedge_del (newlast, prev, E);
                markedge_del (newlast, next, E);
                hit += step_noback (G, D, E, Q, F, level+1, gain, Gstar, first,
                                    newlast, fstack, intptr_world);
                unmarkedge_add (last, newlast, E);
                unmarkedge_add (next, prev, E);
                unmarkedge_del (newlast, prev, E);
                unmarkedge_del (newlast, next, E);
            }

            if (!hit) {
                UNFLIP (newlast, prev, last, next, fstack, F);
                UNFLIP (first, last, newlast, next, fstack, F);
                return 0;
            } else {
                MARK (newlast, Q, F, D, G, intptr_world);
                MARK (next, Q, F, D, G, intptr_world);
                MARK (prev, Q, F, D, G, intptr_world);
                return 1;
            }
        } else
#endif /* NODE_INSERTIONS */
        {
#ifdef MAK_MORTON
            if (e.mm) {
                int hit = 0;
                int maxHit=0;
                int this = e.other;
                int newfirst = e.over;
                int val;

                gain -= e.diff;
                val = gain - Edgelen (newfirst, last, D);
                if (val > *Gstar) {
                    *Gstar = val;
                    hit++;
                }
                if (val>*Gmax){
                    *Gmax = val;
                    maxHit++;
                }
                double prevGmax = *Gmax;
                // printf("updated Gstar to %d and Gmax to %d\n", *Gstar, *Gmax);
                FLIP (this, newfirst, first, last, fstack, F);
                int nextHit=0;

                if (level < MAXDEPTH) {
                    markedge_add (first, this, E);
                    markedge_del (this, newfirst, E);
                    nextHit = ttp_step_noback (G, D, E, Q, F, level + 1, gain, Gstar,
                                        newfirst, last, fstack, intptr_world, instance, solution, prevObj, Gmax);
                    unmarkedge_add (first, this, E);
                    unmarkedge_del (this, newfirst, E);
                }
    #ifdef CHECK_OBJECTIVE
                double objective = 0;
                if (maxHit>0 && nextHit==0 && prevGmax>GMAX_THRESHOLD){
                // if (hit>0 && nextHit==0){
                    double initTime = CCutil_zeit();
                    ttpObjective obj = evaluateFromFlipper(instance, solution, F);
                    checkTime+=CCutil_zeit()-initTime;
                    ++numChecks;
                    // printf("snb hit=%d, nextHit=%d, dist = %lf, obj = %lf, prevObj = %lf\n", hit, nextHit, obj.distance, obj.objective, prevObj);
                    hit = 0;
                    if (obj.objective>prevObj){
                        ++beneficialOpts;
                        hit = 1;
                        *Gstar = val;
                        *Gmax = val;
                    }
                }
    #endif
                hit+=nextHit;

                if (!hit) {
                    if (level==0) {
                        *Gstar=0;
                        *Gmax=0;
                    }
                    UNFLIP (this, newfirst, first, last, fstack, F);
                    return 0;
                } else {
                    // printf("add to queue: %d, %d\n", this, newfirst);
                    MARK (this, Q, F, D, G, intptr_world);
                    MARK (newfirst, Q, F, D, G, intptr_world);
                    return 1;
                }
            } else 
#endif  /* MAK_MORTON */
            {
                int hit = 0;
                int maxHit=0;
                int this = e.other;
                int newlast = e.over;
                int val;

                gain -= e.diff;
                val = gain - Edgelen (newlast, first, D);
                if (val > *Gstar) {
                    *Gstar = val;
                    hit++;
                }
                if (val>*Gmax){
                    *Gmax = val;
                    maxHit++;
                }
                double prevGmax = *Gmax;
                // printf("updated Gstar to %d and Gmax to %d\n", *Gstar, *Gmax);

                FLIP (first, last, newlast, this, fstack, F);

                int nextHit = 0;
                if (level < MAXDEPTH) {
                    markedge_add (last, this, E);
                    markedge_del (this, newlast, E);
                    nextHit = ttp_step_noback (G, D, E, Q, F, level + 1, gain, Gstar,
                                        first, newlast, fstack, intptr_world, instance, solution, prevObj, Gmax);
                    unmarkedge_add (last, this, E);
                    unmarkedge_del (this, newlast, E);
                }

                bool refresh = false;
#ifdef CHECK_OBJECTIVE
                if (maxHit>0 && nextHit==0 && prevGmax>GMAX_THRESHOLD){
                // if (hit>0 && nextHit==0){
                    double initTime = CCutil_zeit();
                    ttpObjective obj = evaluateFromFlipper(instance, solution, F);
                    checkTime+=CCutil_zeit()-initTime;
                    ++numChecks;
                    // printf("snb hit=%d, nextHit=%d, dist = %lf, obj = %lf, prevObj = %lf\n", hit, nextHit, obj.distance, obj.objective, prevObj);
                    hit = 0;
                    if (obj.objective>=prevObj){
                        ++beneficialOpts;
                        refresh = obj.objective>prevObj;
                        hit = 1;
                        *Gstar = val;
                        *Gmax = val;
                    }
                    
                }
#endif
                hit+=nextHit;

                if (!hit) {
                    if (level==0) {
                        *Gstar=0;
                        *Gmax=0;
                    }
                    UNFLIP (first, last, newlast, this, fstack, F);
                    return 0;
                } else {
                    // printf("add to queue: %d, %d\n", this, newlast);
                    MARK (this, Q, F, D, G, intptr_world);
                    MARK (newlast, Q, F, D, G, intptr_world);
                    return 1;
                }
            }
        }
    } else {
        return 0;
    }
}


static double kick_improve (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, flipstack *win, flipstack *fstack,
        CCptrworld *intptr_world)
{
    int t1, t2;
    int gain, Gstar = 0;
    int hit = 0;

    do {
        first_kicker (G, D, F, &t1, &t2);
        gain = Edgelen (t1, t2, D);
        markedge_del (t1, t2, E);
        hit = kick_step_noback (G, D, E, Q, F, 0, gain, &Gstar, t1, t2, win,
                                fstack, intptr_world);
        unmarkedge_del (t1, t2, E);
    } while (!hit);

    kickturn (t1, Q, D, G, F, intptr_world);
    kickturn (t2, Q, D, G, F, intptr_world);

    return (double) -Gstar;
}

#define G_MULT 1.5 

static int kick_step_noback (graph *G, distobj *D, adddel *E, aqueue *Q, 
        CClk_flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *win, flipstack *fstack, CCptrworld *intptr_world)
{
    edgelook winner;
    int val;
    int this, prev, newlast;
    int lastnext = CClinkern_flipper_next (F, last);
    int i;
    int cutoff = (int) (G_MULT * (double) gain);
    edge **goodlist = G->goodlist;

    winner.diff = BIGINT;
    for (i = 0; goodlist[last][i].weight < cutoff; i++) {
        this = goodlist[last][i].other;
        if (!is_it_deleted (last, this, E) && this != first &&
                                              this != lastnext) {
            prev = CClinkern_flipper_prev (F, this);
            if (!is_it_added (this, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (this, prev, D);
                if (val < winner.diff) {
                    winner.diff = val;
                    winner.other = this;
                    winner.over = prev;
                }
            }
        }
    }

    if (winner.diff < BIGINT) {
        this = winner.other;
        newlast = winner.over;
        gain -= winner.diff;
        *Gstar = gain - Edgelen (newlast, first, D);

        FLIP (first, last, newlast, this, fstack, F);
        kickturn (this, Q, D, G, F, intptr_world);
        kickturn (newlast, Q, D, G, F, intptr_world);
        if (win->counter < win->max) {
            win->stack[win->counter].first = last;
            win->stack[win->counter].last = newlast;
            win->counter++;
        }

        if (level < KICK_MAXDEPTH) {
            markedge_add (last, this, E);
            markedge_del (this, newlast, E);
            kick_step_noback (G, D, E, Q, F, level+1, gain, Gstar, first,
                              newlast, win, fstack, intptr_world);
            unmarkedge_add (last, this, E);
            unmarkedge_del (this, newlast, E);
        }
        return 1;
    } else {
        return 0;
    }
}

static int ttp_weird_second_step (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int len_t1_t2, int t1, int t2, flipstack *fstack,
        CCptrworld *intptr_world, CCptrworld *edgelook_world, ttpInstance *instance, 
        ttpSolution *solution, mbfs_main *main, logisticRegressor *classifier, double prevObj,double *capacity)
{
    int t3, t4, t5, t6, t7, t8;
    int oldG, gain, tG, Gstar = 0, val, hit, Gmax=-(CCutil_MAXINT-1), maxHit;
    int t3prev, t4next;
    edgelook *e, *f, *h, *list, *list2, *list3;

    list = weird_look_ahead (G, D, F, len_t1_t2, t1, t2, edgelook_world);
    for (h = list; h; h = h->next) {
        t3 = h->other;
        t4 = h->over;

        oldG = len_t1_t2 - h->diff;
  
        t3prev = CClinkern_flipper_prev (F, t3);
        t4next = CClinkern_flipper_next (F, t4);
  
        markedge_add (t2, t3, E);
        markedge_del (t3, t4, E);
        G->weirdmagic++;
        G->weirdmark[t1] = G->weirdmagic;
        G->weirdmark[t2] = G->weirdmagic;
        G->weirdmark[t3] = G->weirdmagic;
        G->weirdmark[t4next] = G->weirdmagic;

        list2 = weird_look_ahead2 (G, D, F, oldG, t2, t3, t4, edgelook_world);
        for (e = list2; e; e = e->next) {
            t5 = e->other;
            t6 = e->over;

            markedge_add (t4, t5, E);
            if (e->seq) {
                if (!e->side) {
                    gain = oldG - e->diff;
                    val = gain - Edgelen (t6, t1, D);
                    if (val > Gstar)
                        Gstar = val;
                    if (val>Gmax){
                        Gmax = val;
                        maxHit++;
                    }
                    FLIP (t1, t2, t6, t5, fstack, F);
                    FLIP (t2, t5, t3, t4, fstack, F);

                    markedge_del (t5, t6, E);
                    hit = ttp_step (G, D, E, Q, F, 2, gain, &Gstar, t1, t6, 
                                fstack, intptr_world, edgelook_world, 
                                instance, solution, main, classifier, prevObj, capacity, &Gmax);
                    unmarkedge_del (t5, t6, E);
#ifdef CHECK_OBJECTIVE
                    if (!hit && Gmax>GMAX_THRESHOLD) {
                        double initTime = CCutil_zeit();
                        ttpObjective obj = evaluateFromFlipper(instance, solution, F);
                        checkTime+=CCutil_zeit()-initTime;
                        ++numChecks;
                        // printf("wss hit=%d, nextHit=%d, dist = %lf, obj = %lf, prevObj = %lf\n", Gstar>0?1:0, hit, obj.distance, obj.objective, prevObj);
                        if (obj.objective>prevObj) {
                            ++beneficialOpts;
                            hit++;
                            Gstar = val;
                            Gmax = val;
                        }
                    }
#else 
                    if (!hit && Gstar) hit=1;
#endif

                    if (!hit) {
                        UNFLIP (t2, t5, t3, t4, fstack, F);
                        UNFLIP (t1, t2, t6, t5, fstack, F);
                    } else {
                        unmarkedge_add (t2, t3, E);
                        unmarkedge_del (t3, t4, E);
                        unmarkedge_add (t4, t5, E);
                        MARK (t3, Q, F, D, G, intptr_world);
                        MARK (t4, Q, F, D, G, intptr_world);
                        MARK (t5, Q, F, D, G, intptr_world);
                        MARK (t6, Q, F, D, G, intptr_world);
                        edgelook_listfree (edgelook_world, list);
                        edgelook_listfree (edgelook_world, list2);
                        return Gmax;
                    }
                } else {   
                    gain = oldG - e->diff;
                    val = gain - Edgelen (t6, t1, D);
                    if (val > Gstar)
                        Gstar = val;
                    if (val>Gmax){
                        Gmax = val;
                        maxHit++;
                    }
                    FLIP (t1, t2, t3, t4, fstack, F);
                    FLIP (t6, t5, t2, t4, fstack, F);
                    FLIP (t1, t3, t6, t2, fstack, F);

                    markedge_del (t5, t6, E);
                    hit = ttp_step (G, D, E, Q, F, 2, gain, &Gstar, t1, t6,
                                fstack, intptr_world, edgelook_world, 
                                instance, solution, main, classifier, prevObj, capacity, &Gmax);
                    unmarkedge_del (t5, t6, E);
#ifdef CHECK_OBJECTIVE
                    if (!hit && Gmax>GMAX_THRESHOLD) {
                        double initTime = CCutil_zeit();
                        ttpObjective obj = evaluateFromFlipper(instance, solution, F);
                        checkTime+=CCutil_zeit()-initTime;
                        ++numChecks;
                        // printf("wss hit=%d, nextHit=%d, dist = %lf, obj = %lf, prevObj = %lf\n", Gstar>0?1:0, hit, obj.distance, obj.objective, prevObj);
                        if (obj.objective>prevObj) {
                            ++beneficialOpts;
                            hit++;
                            Gstar = val;
                            Gmax = val;
                        }
                    }
#else 
                    if (!hit && Gstar) hit=1;
#endif

                    if (!hit) {
                        UNFLIP (t1, t3, t6, t2, fstack, F);
                        UNFLIP (t6, t5, t2, t4, fstack, F);
                        UNFLIP (t1, t2, t3, t4, fstack, F); 
                    } else {
                        unmarkedge_add (t2, t3, E);
                        unmarkedge_del (t3, t4, E);
                        unmarkedge_add (t4, t5, E);
                        MARK (t3, Q, F, D, G, intptr_world);
                        MARK (t4, Q, F, D, G, intptr_world);
                        MARK (t5, Q, F, D, G, intptr_world);
                        MARK (t6, Q, F, D, G, intptr_world);
                        edgelook_listfree (edgelook_world, list);
                        edgelook_listfree (edgelook_world, list2);
                        return Gmax;
                    }
                }
            } else {
                tG = oldG - e->diff;
                markedge_del (t5, t6, E);
                list3 = weird_look_ahead3 (G, D, F, tG, t2, t3, t6,
                                           edgelook_world);
                for (f = list3; f; f = f->next) {
                    t7 = f->other;
                    t8 = f->over;
                    gain = tG - f->diff;
                    if (!f->side) {
                        val = gain - Edgelen (t8, t1, D);
                        if (val > Gstar)
                            Gstar = val;
                        if (val>Gmax){
                            Gmax = val;
                            maxHit++;
                        }
                        FLIP (t1, t2, t8, t7, fstack, F);
                        FLIP (t2, t7, t3, t4, fstack, F);
                        FLIP (t7, t4, t6, t5, fstack, F);

                        markedge_add (t6, t7, E);
                        markedge_del (t7, t8, E);
                        hit = ttp_step (G, D, E, Q, F, 3, gain, &Gstar, t1, t8,
                                    fstack, intptr_world, edgelook_world, 
                                    instance, solution, main, classifier, prevObj, capacity, &Gmax);
                        unmarkedge_del (t6, t7, E);
                        unmarkedge_del (t7, t8, E);
#ifdef CHECK_OBJECTIVE
                        if (!hit && Gmax>GMAX_THRESHOLD) {
                            double initTime = CCutil_zeit();
                            ttpObjective obj = evaluateFromFlipper(instance, solution, F);
                            checkTime+=CCutil_zeit()-initTime;
                            ++numChecks;
                            // printf("wss hit=%d, nextHit=%d, dist = %lf, obj = %lf, prevObj = %lf\n", Gstar>0?1:0, hit, obj.distance, obj.objective, prevObj);
                            if (obj.objective>prevObj){
                                ++beneficialOpts;
                                hit++;
                                Gstar = val;
                                Gmax = val;
                            }
                            
                        }
#else 
                        if (!hit && Gstar) hit=1;
#endif
                        if (!hit) {
                            UNFLIP (t7, t4, t6, t5, fstack, F);
                            UNFLIP (t2, t7, t3, t4, fstack, F);
                            UNFLIP (t1, t2, t8, t7, fstack, F);
                        } else {
                            unmarkedge_add (t2, t3, E);
                            unmarkedge_del (t3, t4, E);
                            unmarkedge_add (t4, t5, E);
                            unmarkedge_del (t5, t6, E);
                            MARK (t3, Q, F, D, G, intptr_world);
                            MARK (t4, Q, F, D, G, intptr_world);
                            MARK (t5, Q, F, D, G, intptr_world);
                            MARK (t6, Q, F, D, G, intptr_world);
                            MARK (t7, Q, F, D, G, intptr_world);
                            MARK (t8, Q, F, D, G, intptr_world);
                            edgelook_listfree (edgelook_world, list);
                            edgelook_listfree (edgelook_world, list2);
                            edgelook_listfree (edgelook_world, list3);
                            return Gmax;
                        }
                    } else {
                        val = gain - Edgelen (t8, t1, D);
                        if (val > Gstar)
                            Gstar = val;
                        if (val>Gmax){
                            Gmax = val;
                            maxHit++;
                        }
                        FLIP (t1, t2, t6, t5, fstack, F);
                        FLIP (t1, t6, t8, t7, fstack, F);
                        FLIP (t3, t4, t2, t5, fstack, F);

                        markedge_add (t6, t7, E);
                        markedge_del (t7, t8, E);
                        hit = ttp_step (G, D, E, Q, F, 3, gain, &Gstar, t1, t8,
                                    fstack, intptr_world, edgelook_world, 
                                    instance, solution, main, classifier, prevObj, capacity, &Gmax);
                        unmarkedge_add (t6, t7, E);
                        unmarkedge_del (t7, t8, E);
#ifdef CHECK_OBJECTIVE
                        if (!hit && Gmax>GMAX_THRESHOLD) {
                            double initTime = CCutil_zeit();
                            ttpObjective obj = evaluateFromFlipper(instance, solution, F);
                            checkTime+=CCutil_zeit()-initTime;
                            ++numChecks;
                            // printf("wss hit=%d, nextHit=%d, dist = %lf, obj = %lf, prevObj = %lf\n", Gstar>0?1:0, hit, obj.distance, obj.objective, prevObj);
                            if (obj.objective>prevObj) {
                                ++beneficialOpts;
                                hit++;
                                Gstar = val;
                                Gmax = val;
                            }
                        }
#else 
                        if (!hit && Gstar) hit=1;
#endif
                        if (!hit) {
                            UNFLIP (t3, t4, t2, t5, fstack, F);
                            UNFLIP (t1, t6, t8, t7, fstack, F);
                            UNFLIP (t1, t2, t6, t5, fstack, F);
                        } else {
                            unmarkedge_add (t2, t3, E);
                            unmarkedge_del (t3, t4, E);
                            unmarkedge_add (t4, t5, E);
                            unmarkedge_del (t5, t6, E);
                            MARK (t3, Q, F, D, G, intptr_world);
                            MARK (t4, Q, F, D, G, intptr_world);
                            MARK (t5, Q, F, D, G, intptr_world);
                            MARK (t6, Q, F, D, G, intptr_world);
                            MARK (t7, Q, F, D, G, intptr_world);
                            MARK (t8, Q, F, D, G, intptr_world);
                            edgelook_listfree (edgelook_world, list);
                            edgelook_listfree (edgelook_world, list2);
                            edgelook_listfree (edgelook_world, list3);
                            return Gmax;
                        }
                    }
                }
                edgelook_listfree (edgelook_world, list3);
                unmarkedge_del (t5, t6, E);
            }
            unmarkedge_add (t4, t5, E);
        }
        edgelook_listfree (edgelook_world, list2);
        unmarkedge_add (t2, t3, E);
        unmarkedge_del (t3, t4, E);
    }
    edgelook_listfree (edgelook_world, list);
    return 0;
}

static edgelook *look_ahead (graph *G, distobj *D, adddel *E, CClk_flipper *F,
       int first, int last, int gain, int level, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL, *el;
    int i, val;
    int this, prev;
    int lastnext = CClinkern_flipper_next (F, last);
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
#if defined(MAK_MORTON) && defined(FULL_MAK_MORTON)
    int mm[MAX_BACK];
#endif
    int k, ahead = backtrack_count[level];
    edge **goodlist = G->goodlist;

    for (i = 0; i < ahead; i++) {
        value[i] = BIGINT;
#if defined(MAK_MORTON) && defined(FULL_MAK_MORTON)
        mm[i] = 0;
#endif
    }
    value[ahead] = -BIGINT;

#ifdef USE_LESS_OR_EQUAL
    // for (i = 0; goodlist[last][i].weight <= gain; i++) {
        for (i = 0; i<MAX_BACK; i++) {
#else
    for (i = 0; goodlist[last][i].weight < gain; i++) {
#endif
        this = goodlist[last][i].other;
        if (!is_it_deleted (last, this, E) && this != first &&
                                              this != lastnext) {
            prev = CClinkern_flipper_prev (F, this);
            if (!is_it_added (this, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (this, prev, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                    }
                    value[k] = val;
                    other[k] = this;
                    save[k] = prev;
                }
            }
        }
    }

#if defined(MAK_MORTON) && defined(FULL_MAK_MORTON)
    {
        int firstprev = CClinkern_flipper_prev (F, first);
        int next;
       
#ifdef USE_LESS_OR_EQUAL
        for (i = 0; goodlist[first][i].weight <= gain; i++) {
#else
        for (i = 0; goodlist[first][i].weight < gain; i++) {
#endif
            this = goodlist[first][i].other;
            if (!is_it_deleted (first, this, E) && this != last &&
                                                   this != firstprev) {
                next = CClinkern_flipper_next (F, this);
                if (!is_it_added (this, next, E)) {
                    val = goodlist[first][i].weight - Edgelen (this, next, D);
                    if (val < value[0]) {
                        for (k = 0; value[k+1] > val; k++) {
                            value[k] = value[k+1];
                            other[k] = other[k+1];
                            save[k] = save[k+1];
                            mm[k] = mm[k+1];
                        }
                        value[k] = val;
                        other[k] = this;
                        save[k] = next;
                        mm[k] = 1;
                    }
                }
            }
        }
    }
#endif

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc (edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->next = list;
#if defined(MAK_MORTON) && defined(FULL_MAK_MORTON)
            el->mm = mm[i];
#endif
            list = el;
        }
    }

    return list;
}

static void look_ahead_noback (graph *G, distobj *D, adddel *E, CClk_flipper *F,
        int first, int last, int gain, edgelook *winner)
{
    int val;
    int this, prev;
    int lastnext = CClinkern_flipper_next (F, last);
    int i;
#if defined(MAK_MORTON) || defined(NODE_INSERTIONS)
    int next;
#endif
    edge **goodlist = G->goodlist;

    winner->diff = BIGINT;
    for (i = 0; goodlist[last][i].weight < gain; i++) {
        this = goodlist[last][i].other;
        if (!is_it_deleted (last, this, E) && this != first &&
                                              this != lastnext) {
            prev = CClinkern_flipper_prev (F, this);
            if (!is_it_added (this, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (this, prev, D);
                if (val < winner->diff) {
                    winner->diff = val;
                    winner->other = this;
                    winner->over = prev;
#ifdef MAK_MORTON
                    winner->mm = 0;
#endif
#ifdef NODE_INSERTIONS
                    winner->ni = 0;
#endif
                }
#ifdef NODE_INSERTIONS
                next =  CClinkern_flipper_next (F, this);
                if (!is_it_added (this, next, E) &&
                    !is_it_deleted (prev, next, E)) {
                    val += (Edgelen (next, prev, D) - Edgelen (this, next, D));
                    if (val < winner->diff) {
                        winner->diff = val;
                        winner->other = this;
                        winner->over = prev;
                        winner->under = next;
                        winner->ni = 1;
                    }
                }
#endif
            }
        }
    }
#ifdef MAK_MORTON
    {
        int firstprev = CClinkern_flipper_prev (F, first);

        for (i = 0; goodlist[first][i].weight < gain; i++) {
            this = goodlist[first][i].other;
            if (!is_it_deleted (first, this, E) && this != last && 
                                                   this != firstprev) {
                next = CClinkern_flipper_next (F, this);
                if (!is_it_added (this, next, E)) {
                    val = goodlist[first][i].weight - Edgelen (this, next, D);
                    if (val < winner->diff) {
                        winner->diff = val;
                        winner->other = this;
                        winner->over = next;
                        winner->mm = 1;
#ifdef NODE_INSERTIONS
                        winner->ni = 0;
#endif
                    }
                }
            }
        }
    }
#endif
}

static edgelook *weird_look_ahead (graph *G, distobj *D, CClk_flipper *F,
        int gain, int t1, int t2, CCptrworld *edgelook_world)
{
    edgelook *list, *el;
    int i, this, next;
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val, ahead;
    edge **goodlist = G->goodlist;

    list = (edgelook *) NULL;
    ahead = weird_backtrack_count[0];
    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

#ifdef USE_LESS_OR_EQUAL
    for (i = 0; goodlist[t2][i].weight <= gain; i++) {
#else
    for (i = 0; goodlist[t2][i].weight < gain; i++) {
#endif
        this = goodlist[t2][i].other;
        if (this != t1) {
            next = CClinkern_flipper_next (F, this);
            val = goodlist[t2][i].weight - Edgelen (this, next, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                }
                value[k] = val;
                other[k] = this;
            save[k] = next;
            }
        }
    }
    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc (edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static edgelook *weird_look_ahead2 (graph *G, distobj *D, CClk_flipper *F,
       int gain, int t2, int t3, int t4, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL;
    edgelook *el;
    int i, t5, t6;
    int other[MAX_BACK], save[MAX_BACK], seq[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead = weird_backtrack_count[1];
    edge **goodlist = G->goodlist;
    int  *weirdmark = G->weirdmark;
    int  weirdmagic = G->weirdmagic;

    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

#ifdef USE_LESS_OR_EQUAL
    for (i = 0; goodlist[t4][i].weight <= gain; i++) {
#else
    for (i = 0; goodlist[t4][i].weight < gain; i++) {
#endif
        t5 = goodlist[t4][i].other;
        if (weirdmark[t5] != weirdmagic) {
            if (CClinkern_flipper_sequence (F, t2, t5, t3)) {
                t6 = CClinkern_flipper_prev (F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 0;
                }
                t6 = CClinkern_flipper_next (F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 1;
                }
            } else {
                t6 = CClinkern_flipper_prev (F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 0;
                    side[k] = 0;
                }
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc (edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->seq = seq[i];
            el->side = side[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static edgelook *weird_look_ahead3 (graph *G, distobj *D, CClk_flipper *F,
        int gain, int t2, int t3, int t6, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL;
    edgelook *el;
    int i, t7, t8;
    int other[MAX_BACK], save[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead = weird_backtrack_count[2];
    edge **goodlist = G->goodlist;
    int  *weirdmark = G->weirdmark;
    int  weirdmagic = G->weirdmagic;

    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

#ifdef USE_LESS_OR_EQUAL
    for (i = 0; goodlist[t6][i].weight <= gain; i++) {
#else 
    for (i = 0; goodlist[t6][i].weight < gain; i++) {
#endif
        t7 = goodlist[t6][i].other;   /* Need t7 != t2, t3, t2next, t3prev */
        if (weirdmark[t7] != weirdmagic &&
                   CClinkern_flipper_sequence (F, t2, t7, t3)) {
            t8 = CClinkern_flipper_prev (F, t7);
            val = goodlist[t6][i].weight - Edgelen (t7, t8, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                    side[k] = side[k+1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 0;
            }
            t8 = CClinkern_flipper_next (F, t7);
            val = goodlist[t6][i].weight - Edgelen (t7, t8, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                    side[k] = side[k+1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 1;
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc (edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->side = side[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static int save_tour (int ncount, char *sname, CClk_flipper *F)
{
    int rval = 0;
    int *ctemp = (int *) NULL;

    ctemp = CC_SAFE_MALLOC (ncount, int);
    if (ctemp == (int *) NULL) {
        fprintf (stderr, "out of memory in save_tour\n");
        rval = 1; goto CLEANUP;
    }
    CClinkern_flipper_cycle (F, ctemp);
    rval = CCutil_writecycle (ncount, sname, ctemp, 0);
    CCcheck_rval (rval, "CCutil_writecycle failed");

    printf ("Wrote the tour to %s\n", sname); fflush (stdout);

CLEANUP:

    CC_IFFREE (ctemp, int);
    return rval;
}

static double cycle_length (int ncount, int *cyc, distobj *D)
{
    int i;
    double val = 0.0;
    
    for (i = 1; i < ncount; i++) {
        val += (double) Edgelen (cyc[i - 1], cyc[i], D);
    }
    val += (double) Edgelen (cyc[0], cyc[ncount - 1], D);

    return val;
}

static int random_four_swap (graph *G, distobj *D, aqueue *Q, CClk_flipper *F,
       CCkdtree *kdt, int *delta, int kicktype, flipstack *win, 
       flipstack *fstack, CCptrworld *intptr_world, CCrandstate *rstate)
{
    int rval = 0;
    int t1, t2, t3, t4, t5, t6, t7, t8, temp;

    switch (kicktype) {
    case CC_LK_RANDOM_KICK:
        find_random_four (G, D, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);
        break;
    case CC_LK_WALK_KICK:
        find_walk_four (G, D, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);
        break;
    case CC_LK_CLOSE_KICK:
        find_close_four (G, D, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);
        break;
    case CC_LK_GEOMETRIC_KICK:
        rval = find_geometric_four (G, D, F, kdt, &t1, &t2, &t3, &t4,
                                                  &t5, &t6, &t7, &t8, rstate);
        if (rval) {
            fprintf (stderr, "find_geometric_four failed\n"); return 1;
        }
        break;
    default:
        fprintf (stderr, "unknown kick type %d\n", kicktype); return 1;
    }

    if (!CClinkern_flipper_sequence (F, t1, t3, t5)) {
        CC_SWAP (t3, t5, temp);
        CC_SWAP (t4, t6, temp);
    }
    if (!CClinkern_flipper_sequence (F, t1, t5, t7)) {
        CC_SWAP (t5, t7, temp);
        CC_SWAP (t6, t8, temp);
        if (!CClinkern_flipper_sequence (F, t1, t3, t5)) {
            CC_SWAP (t3, t5, temp);
            CC_SWAP (t4, t6, temp);
        }
    }
    FLIP (t1, t2, t5, t6, fstack, F);
    FLIP (t4, t3, t7, t8, fstack, F);
    FLIP (t1, t5, t6, t2, fstack, F);

    if (win->counter < win->max) {
        win->stack[win->counter].first = t2;
        win->stack[win->counter].last = t5;
        win->counter++;
    }
    if (win->counter < win->max) {
        win->stack[win->counter].first = t3;
        win->stack[win->counter].last = t7;
        win->counter++;
    }
    if (win->counter < win->max) {
        win->stack[win->counter].first = t5;
        win->stack[win->counter].last = t6;
        win->counter++;
    }

    bigturn (G ,t1, 0, Q, F, D, intptr_world);
    bigturn (G, t2, 1, Q, F, D, intptr_world);
    bigturn (G, t3, 0, Q, F, D, intptr_world);
    bigturn (G, t4, 1, Q, F, D, intptr_world);
    bigturn (G, t5, 0, Q, F, D, intptr_world);
    bigturn (G, t6, 1, Q, F, D, intptr_world);
    bigturn (G, t7, 0, Q, F, D, intptr_world);
    bigturn (G, t8, 1, Q, F, D, intptr_world);

    *delta = 
           Edgelen (t1, t6, D) + Edgelen (t2, t5, D) +
           Edgelen (t3, t8, D) + Edgelen (t4, t7, D) -
           Edgelen (t1, t2, D) - Edgelen (t3, t4, D) -
           Edgelen (t5, t6, D) - Edgelen (t7, t8, D);
    return 0;
}

#define HUNT_PORTION_LONG 0.001 

static void first_kicker (graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2)
{
#ifdef LONG_KICKER
    int longcount = (int) ((double) G->ncount * HUNT_PORTION_LONG) + 10;
    int i, best, try1, len, next, prev, nextl, prevl;
    int ncount = G->ncount;
    edge **goodlist = G->goodlist;

    try1 = CCutil_lprand (G->rstate) % ncount;
    next = CClinkern_flipper_next (F, try1);
    prev = CClinkern_flipper_prev (F, try1);
    nextl = Edgelen (try1, next, D);
    prevl = Edgelen (try1, prev, D);
    if (nextl >= prevl) {
        *t1 = try1;
        *t2 = next;
        best = nextl - goodlist[*t1][0].weight;
    } else {
        *t1 = prev;
        *t2 = try1;
        best = prevl - goodlist[*t1][0].weight;
    }

    for (i = 0; i < longcount; i++) {
        try1 = CCutil_lprand (G->rstate) % ncount;
        next = CClinkern_flipper_next (F, try1);
        prev = CClinkern_flipper_prev (F, try1);
        nextl = Edgelen (try1, next, D);
        prevl = Edgelen (try1, prev, D);
        if (nextl >= prevl) {
            len = nextl - goodlist[try1][0].weight;
            if (len > best) {
                *t1 = try1;
                *t2 = next;
            }
        } else {
            len = prevl - goodlist[try1][0].weight;
            if (len > best) {
                *t1 = prev;
                *t2 = try1;
            }
        }
    }
#else   /* LONG_KICKER */
    *t1 = CCutil_lprand (G->rstate) % G->ncount;
    *t2 = CClinkern_flipper_next (F, *t1);
#endif  /* LONG_KICKER */
}

static void find_random_four (graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8)
{
    int ncount = G->ncount;

    first_kicker (G, D, F, t1, t2);
    do {
        *t3 = CCutil_lprand (G->rstate) % ncount;
        *t4 = CClinkern_flipper_next (F, *t3);
    } while (*t3 == *t1 || *t3 == *t2 || *t4 == *t1);

    do {
        *t5 = CCutil_lprand (G->rstate) % ncount;
        *t6 = CClinkern_flipper_next (F, *t5);
    } while (*t5 == *t1 || *t5 == *t2 || *t5 == *t3 || *t5 == *t4 ||
             *t6 == *t1 || *t6 == *t3);

    do {
        *t7 = CCutil_lprand (G->rstate) % ncount;
        *t8 = CClinkern_flipper_next (F, *t7);
    } while (*t7 == *t1 || *t7 == *t2 ||
             *t7 == *t3 || *t7 == *t4 || *t7 == *t5 || *t7 == *t6 ||
             *t8 == *t1 || *t8 == *t3 || *t8 == *t5 );
}


#define HUNT_PORTION    0.03
#define RAND_TRYS       6    /* To find the 3 other edges */

static void find_close_four (graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8)
{
    int s1, s2, s3, s4, s5, s6, s7, s8;
    int i, k, try1, trydist;
    int count = (int) ((double) G->ncount * HUNT_PORTION) + 1 + RAND_TRYS;
    int trials[RAND_TRYS + 1];
    int tdist[RAND_TRYS + 1];

    first_kicker (G, D, F, &s1, &s2);


TRYAGAIN:

    for (k = 0; k < RAND_TRYS; k++) tdist[k] = BIGINT;
    tdist[RAND_TRYS] = -BIGINT;
    for (i = 0; i < count; i++) {
        try1 = CCutil_lprand (G->rstate) % G->ncount;
        trydist = Edgelen (try1, s1, D); 
        if (trydist < tdist[0]) {
            for (k = 0; tdist[k + 1] > trydist; k++) {
                tdist[k] = tdist[k + 1];
                trials[k] = trials[k + 1];
            }
            tdist[k] = trydist;
            trials[k] = try1;
        }
    }

    k = RAND_TRYS-1;
    do {
        if (k < 0) goto TRYAGAIN;
        s3 = trials[k--];
        s4 = CClinkern_flipper_next (F, s3);
    } while (s3 == s1 || s3 == s2 || s4 == s1);

    do {
        if (k < 0) goto TRYAGAIN;
        s5 = trials[k--];
        s6 = CClinkern_flipper_next (F, s5);
    } while (s5 == s1 || s5 == s2 || s5 == s3 ||
             s5 == s4 || s6 == s1 || s6 == s3);

    do {
        if (k < 0) goto TRYAGAIN;
        s7 = trials[k--];
        s8 = CClinkern_flipper_next (F, s7);
    } while (s7 == s1 || s7 == s2 || s7 == s3 || s7 == s4 ||
             s7 == s5 || s7 == s6 || s8 == s1 || s8 == s3 ||
             s8 == s5);

    *t1 = s1; *t2 = s2; *t3 = s3; *t4 = s4;
    *t5 = s5; *t6 = s6; *t7 = s7; *t8 = s8;
}

#define GEO_FACTOR  50
#define GEO_MAX    250

static int find_geometric_four (graph *G, distobj *D, CClk_flipper *F,
        CCkdtree *kdt, int *t1, int *t2, int *t3, int *t4, int *t5, int *t6,
        int *t7, int *t8, CCrandstate *rstate)
{
    int neigh[GEO_MAX];
    int temp, i, k, s1, s2, s3, s4, s5, s6, s7, s8;
    int trys, rval;

    first_kicker (G, D, F, &s1, &s2);
    trys = (G->ncount / GEO_FACTOR) + 25;
    if (trys > GEO_MAX) trys = GEO_MAX;
    if (trys > G->ncount - 1) trys = G->ncount - 1;

    rval = CCkdtree_node_k_nearest (kdt, G->ncount, s1, trys, D->dat,
                                    (double *) NULL, neigh, rstate);
    if (rval) {
        fprintf (stderr, "CCkdtree_node_k_nearest failed\n"); return rval;
    }

    for (i = trys; i > trys - 9; i--) {
        k = CCutil_lprand (G->rstate) % i;
        CC_SWAP (neigh[i - 1], neigh[k], temp);
    }

    k = trys - 1;
    do {
        s3 = neigh[k--];
        s4 = CClinkern_flipper_next (F, s3);
    } while (s3 == s2 || s4 == s1);

    do {
        s5 = neigh[k--];
        s6 = CClinkern_flipper_next (F, s5);
    } while (s5 == s2 || s5 == s4 || s6 == s1 || s6 == s3);

    do {
        s7 = neigh[k--];
        s8 = CClinkern_flipper_next (F, s7);
    } while (s7 == s2 || s7 == s4 || s7 == s6 || s8 == s1 || s8 == s3 ||
             s8 == s5);

    *t1 = s1; *t2 = s2; *t3 = s3; *t4 = s4;
    *t5 = s5; *t6 = s6; *t7 = s7; *t8 = s8;

    return 0;
}

#define WALK_STEPS 50

static void find_walk_four (graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8)
{
    int s1, s2, s3, s4, s5, s6, s7, s8;
    int old, n, i, j;

/*
    s1 = CCutil_lprand (G->rstate) % G->ncount;
    s2 = CClinkern_flipper_next (F, s1);
*/

    first_kicker (G, D, F, &s1, &s2);

    do {
        old = -1;
        n = s2;

        for (i = 0;  i < WALK_STEPS; i++) {
            j = CCutil_lprand (G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s3 = n;
        s4 = CClinkern_flipper_next (F, s3);

        n = s4;
        for (i = 0; i < WALK_STEPS; i++) {
            j = CCutil_lprand (G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s5 = n;
        s6 = CClinkern_flipper_next (F, s5);

        n = s6;
        for (i = 0; i < WALK_STEPS; i++) {
            j = CCutil_lprand (G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s7 = n;
        s8 = CClinkern_flipper_next (F, s7);
    } while (s1 == s3 || s1 == s4 || s1 == s5 || s1 == s6 || s1 == s7 ||
             s1 == s8 ||
             s2 == s3 || s2 == s4 || s2 == s5 || s2 == s6 || s2 == s7 ||
             s2 == s8 ||
             s3 == s5 || s3 == s6 || s3 == s7 || s3 == s8 ||
             s4 == s5 || s4 == s6 || s4 == s7 || s4 == s8 ||
             s5 == s7 || s5 == s8 ||
             s6 == s7 || s6 == s8);

    *t1 = s1;  *t2 = s2;  *t3 = s3;  *t4 = s4;
    *t5 = s5;  *t6 = s6;  *t7 = s7;  *t8 = s8;
}

#ifdef USE_HEAP

static void turn (int n, aqueue *Q, CClk_flipper *F, distobj *D, graph *G)

#else  /* USE_HEAP */

#ifdef USE_LESS_MARKING

static void turn (int n, aqueue *Q, CCptrworld *intptr_world)

#else /* USE_LESS_MARKING */

static void turn (int n, aqueue *Q, CClk_flipper *F, CCptrworld *intptr_world)

#endif /* USE_LESS_MARKING */
#endif /* USE_HEAP */
{
#ifdef USE_HEAP
    add_to_active_queue (n, Q, D, G, F);
#else
    add_to_active_queue (n, Q, intptr_world);
#endif

#ifdef MARK_NEIGHBORS
    {
       int i = 0;
       for (i = 0; i < bigG->degree[n]; i++) {
           if (CCutil_lprand (G->rstate) % 2) {
#ifdef USE_HEAP
               add_to_active_queue (G->goodlist[n][i].other, Q, D, G, F);
#else
               add_to_active_queue (bigG->goodlist[n][i].other, Q, intptr_world);

#endif
           }
       }
   }
#else 
#ifndef USE_LESS_MARKING
   {
        int k;
#ifdef USE_HEAP
        k = CClinkern_flipper_next (F, n);
        add_to_active_queue (k, Q, D, G, F);
        k = CClinkern_flipper_next (F, k);
        add_to_active_queue (k, Q, D, G, F);
        k = CClinkern_flipper_prev (F, n);
        add_to_active_queue (k, Q, D, G, F);
        k = CClinkern_flipper_prev (F, k);
        add_to_active_queue (k, Q, D, G, F);
#else
        k = CClinkern_flipper_next (F, n);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_next (F, k);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_prev (F, n);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_prev (F, k);
        add_to_active_queue (k, Q, intptr_world);
#endif /* USE_HEAP */
   }
#endif
#endif 
}

static void kickturn (int n, aqueue *Q, CC_UNUSED distobj *D,
        CC_UNUSED graph *G, CClk_flipper *F, CCptrworld *intptr_world)
{
#ifdef USE_HEAP
    add_to_active_queue (n, Q, D, G, F);
    {
        int k;
        k = CClinkern_flipper_next (F, n);
        add_to_active_queue (k, Q, D, G, F);
        k = CClinkern_flipper_next (F, k);
        add_to_active_queue (k, Q, D, G, F);
        k = CClinkern_flipper_prev (F, n);
        add_to_active_queue (k, Q, D, G, F);
        k = CClinkern_flipper_prev (F, k);
        add_to_active_queue (k, Q, D, G, F);
    }
#else
    add_to_active_queue (n, Q, intptr_world);
    {
        int k;
        k = CClinkern_flipper_next (F, n);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_next (F, k);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_prev (F, n);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_prev (F, k);
        add_to_active_queue (k, Q, intptr_world);
    }
#endif
}

static void bigturn (graph *G, int n, int tonext, aqueue *Q, CClk_flipper *F,
        CC_UNUSED distobj *D, CCptrworld *intptr_world)
{
    int i, k;

#ifdef USE_HEAP
    add_to_active_queue (n, Q, D, G, F);

    if (tonext) {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = CClinkern_flipper_next (F, k);
            add_to_active_queue (k, Q, D, G, F);
        }
    } else {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = CClinkern_flipper_prev (F, k);
            add_to_active_queue (k, Q, D, G, F);
        }
    }

    for (i = 0; i < G->degree[n]; i++) {
        add_to_active_queue (G->goodlist[n][i].other, Q, D, G, F);
    }
#else
    add_to_active_queue (n, Q, intptr_world);
    if (tonext) {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = CClinkern_flipper_next (F, k);
            add_to_active_queue (k, Q, intptr_world);
        }
    } else {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = CClinkern_flipper_prev (F, k);
            add_to_active_queue (k, Q, intptr_world);
        }
    }

    for (i = 0; i < G->degree[n]; i++) {
        add_to_active_queue (G->goodlist[n][i].other, Q, intptr_world);
    }
#endif
}

static void randcycle (int ncount, int *cyc, CCrandstate *rstate)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++) cyc[i] = i;
    for (i = ncount; i > 1; i--) {
        k = CCutil_lprand (rstate) % i;
        CC_SWAP (cyc[i - 1], cyc[k], temp);
    }
}

static void initgraph (graph *G)
{
    G->goodlist   = (edge **) NULL;
    G->edgespace  = (edge *) NULL;    
    G->degree     = (int *) NULL;
    G->weirdmark  = (int *) NULL;
    G->weirdmagic = 0;
    G->ncount     = 0;
}

static void freegraph (graph *G)
{
    if (G) {
        CC_IFFREE (G->goodlist, edge *);
        CC_IFFREE (G->edgespace, edge);
        CC_IFFREE (G->degree, int);
        CC_IFFREE (G->weirdmark, int);
        G->weirdmagic = 0;
        G->ncount     = 0;
    }
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist,
        distobj *D)
{
    int rval = 0;
    int n1, n2, w, i;
    edge *p;

    G->goodlist  = CC_SAFE_MALLOC (ncount, edge *);
    G->degree    = CC_SAFE_MALLOC (ncount, int);
    G->weirdmark = CC_SAFE_MALLOC (ncount, int);
    G->edgespace = CC_SAFE_MALLOC ((2 * ecount) + ncount, edge);
    if (G->goodlist == (edge **) NULL || G->degree == (int *) NULL ||
        G->edgespace == (edge *) NULL)  {
        fprintf (stderr, "out of memory in buildgraph\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        G->degree[i] = 1;
        G->weirdmark[i] = 0;
    }
    for (i = ecount - 1; i >= 0; i--) {
        G->degree[elist[2 * i]]++;
        G->degree[elist[(2 * i) + 1]]++;
    }

    for (i = 0, p = G->edgespace; i < ncount; i++) {
        G->goodlist[i] = p;
        p += (G->degree[i]);
        G->goodlist[i][G->degree[i] - 1].weight = BIGINT;
        G->degree[i] = 0;
    }

    for (i = ecount - 1; i >= 0; i--) {
        n1 = elist[2 * i];
        n2 = elist[(2 * i) + 1];
        w = Edgelen (n1, n2, D);
        insertedge (G, n1, n2, w);
        insertedge (G, n2, n1, w);
    }
    G->ncount     = ncount;
    G->weirdmagic = 0;

CLEANUP:

    if (rval) freegraph (G);
    return rval;
}

static void insertedge (graph *G, int n1, int n2, int w)
{
    int i;
    edge *e = G->goodlist[n1];

    for (i = G->degree[n1] - 1; i >= 0 && e[i].weight >= w; i--) {
        e[i + 1].weight = e[i].weight;
        e[i + 1].other  = e[i].other;
    } 
    e[i + 1].weight = w;
    e[i + 1].other  = n2;
    G->degree[n1]++;
}

static void linkern_free_world (CCptrworld *intptr_world,
        CCptrworld *edgelook_world)
{
    int total, onlist;

    if (intptr_check_leaks (intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs\n", total-onlist);
    }
    if (edgelook_check_leaks (edgelook_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding edgelooks\n", total-onlist);
    }
    CCptrworld_delete (intptr_world);
    CCptrworld_delete (edgelook_world);
}

static int init_flipstack (flipstack *f, int total, int single)
{
    f->counter = 0;
    f->max     = 0;
    f->stack   = (flippair *) NULL;

    f->stack = CC_SAFE_MALLOC (total + single, flippair);
    if (f->stack == (flippair *) NULL) {
        fprintf (stderr, "out of memory in init_flipstack\n"); return 1;
    }
    f->max = total;

    return 0;
}

static void free_flipstack (flipstack *f)
{
    f->counter = 0;
    f->max     = 0;
    CC_IFFREE (f->stack, flippair);
}

static void init_adddel (adddel *E)
{
    E->add_edges = (char *) NULL;
    E->del_edges = (char *) NULL;
}

static void free_adddel (adddel *E)
{
    if (E) {
        CC_IFFREE (E->add_edges, char);
        CC_IFFREE (E->del_edges, char);
    }
}

static int build_adddel (adddel *E, int ncount)
{
    int rval = 0;
    int i, M;

    i = 0;
    while ((1 << i) < ncount)
        i++;
    M = (1 << i);

    E->add_edges = CC_SAFE_MALLOC (M, char);
    E->del_edges = CC_SAFE_MALLOC (M, char);
    if (E->add_edges == (char *) NULL || E->del_edges == (char *) NULL) {
        fprintf (stderr, "out of memory in build_adddel\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < M; i++) {
        E->add_edges[i] = 0;
        E->del_edges[i] = 0;
    }

CLEANUP:

    if (rval) {
        free_adddel (E);
    }
    return rval;
}

static void init_aqueue (aqueue *Q)
{
    Q->active = (char *) NULL;
    Q->active_queue = (intptr *) NULL;
    Q->bottom_active_queue = (intptr *) NULL;
    Q->h = (CCdheap *) NULL;
}

static void free_aqueue (aqueue *Q, CCptrworld *intptr_world)
{
    if (Q) {
        CC_IFFREE (Q->active, char);
        intptr_listfree (intptr_world, Q->active_queue);
        Q->active_queue = (intptr *) NULL;
        Q->bottom_active_queue = (intptr *) NULL;
        if (Q->h) {
            CCutil_dheap_free (Q->h);
            Q->h = (CCdheap *) NULL;
        }
    }
}

static int build_aqueue (aqueue *Q, int ncount, CCptrworld *intptr_world)
{
    int rval = 0;
    int i;

    init_aqueue (Q);

    Q->active = CC_SAFE_MALLOC (ncount, char);
    if (Q->active == (char *) NULL) {
        fprintf (stderr, "out of memory in build_aqueue\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) Q->active[i] = 0;

#ifdef USE_HEAP
    Q->h = CC_SAFE_MALLOC (1, CCdheap);
    if (Q->h == (CCdheap *) NULL) {
        fprintf (stderr, "out of memory in build_aqueue\n");
        rval = 1; goto CLEANUP;
    }
    rval = CCutil_dheap_init (Q->h, ncount);
    if (rval) {
        fprintf (stderr, "CCutil_dheap_init failed\n"); goto CLEANUP;
    }
#endif

CLEANUP:

    if (rval) {
        free_aqueue (Q, intptr_world);
    }
    return rval;
}

#ifdef    USE_HEAP
static void add_to_active_queue (int n, aqueue *Q, distobj *D, graph *G,
        CClk_flipper *F)
{
    if (Q->active[n] == 0) { 
        int next = CClinkern_flipper_next (F, n);
        Q->active[n] = 1;
        Q->h->key[n] = G->goodlist[n][0].weight - Edgelen (n, next, D);
        CCutil_dheap_insert (Q->h, n);
    }
}

static int pop_from_active_queue (aqueue *Q, CCptrworld *intptr_world)
{
    int n;

    n = CCutil_dheap_deletemin (Q->h);
    if (n != -1) {
        Q->active[n] = 0;
    }
    return n;
}

#else  /* USE_HEAP */

static void add_to_active_queue (int n, aqueue *Q, CCptrworld *intptr_world)
{
    intptr *ip;

    /* intptralloc will not fail - the initial supply reserved with
     * intptr_bulkalloc is large enough */

    if (Q->active[n] == 0) { 
        Q->active[n] = 1;
        ip = intptralloc (intptr_world);
        ip->this = n;
        ip->next = (intptr *) NULL;
        if (Q->bottom_active_queue) {
            Q->bottom_active_queue->next = ip;
        } else {
            Q->active_queue = ip;
        }
        Q->bottom_active_queue = ip;
    }
}

static int pop_from_active_queue (aqueue *Q, CCptrworld *intptr_world)
{
    intptr *ip;
    int n = -1;

    if (Q->active_queue != (intptr *) NULL) {
        ip = Q->active_queue;
        n = ip->this;
        Q->active_queue = ip->next;
        if (ip == Q->bottom_active_queue) {
            Q->bottom_active_queue = (intptr *) NULL;
        }
        intptrfree (intptr_world, ip);
        Q->active[n] = 0;
    }
    return n;
}
#endif /* USE_HEAP */

static void init_distobj (distobj *D)
{
    D->dat = (CCdatagroup *) NULL;
    D->cacheind  = (int *) NULL;
    D->cacheval  = (int *) NULL;
    D->cacheM = 0;
}

static void free_distobj (distobj *D)
{
    if (D) {
         D->dat = (CCdatagroup *) NULL;
         CC_IFFREE (D->cacheind, int);
         CC_IFFREE (D->cacheval, int);
         D->cacheM = 0;
    }
}

static int build_distobj (distobj *D, int ncount, CCdatagroup *dat)
{
    int rval = 0;
    int i;

    init_distobj (D);
    D->dat = dat;

#ifndef BENTLEY_CACHE
    i = 0;
    while ((1 << i) < (ncount << 2))
        i++;
    D->cacheM = (1 << i);  
#else
    i = 0;
    while ((1 << i) < ncount)
        i++;
    D->cacheM = (1 << i);
#endif

    D->cacheind = CC_SAFE_MALLOC (D->cacheM, int);
    D->cacheval = CC_SAFE_MALLOC (D->cacheM, int);
    if (D->cacheind == (int *) NULL || D->cacheval == (int *) NULL) {
        fprintf (stderr, "out of memory in build_distobj\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < D->cacheM; i++) {
        D->cacheind[i] = -1;
    }

#ifndef BENTLEY_CACHE
    D->cacheM--;
#endif

CLEANUP:

    if (rval) {
        free_distobj (D);
    }
    return rval; 
}


static int dist (int i, int j, distobj *D)   /* As in Bentley's kdtree paper */
{
    int ind;

    if (i > j) {
        int temp;
        CC_SWAP (i, j, temp);
    }

#ifndef BENTLEY_CACHE
    ind = (((i << 8) + i + j) & (D->cacheM));
#else
    ind = i ^ j;
#endif

    if (D->cacheind[ind] != i) {
        D->cacheind[ind] = i;
        D->cacheval[ind] = CCutil_dat_edgelen (i, j, D->dat);
    }
    return D->cacheval[ind];
}
