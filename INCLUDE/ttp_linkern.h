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

#ifndef  __TTP_LINKERN_H
#define  __TTP_LINKERN_H

#include "util.h"
#include "mbfs.h"

#define CC_LK_RANDOM_KICK    (0)
#define CC_LK_GEOMETRIC_KICK (1)
#define CC_LK_CLOSE_KICK     (2)
#define CC_LK_WALK_KICK      (3)


extern int numHillClimbers, numZeroClimbers, numChecks, beneficialOpts;
extern double linearTime, hillClimberTime, stepTime, checkTime;
extern char *global_filesuffix;

int
    CCttp_linkern_tour (int ncount, CCdatagroup *dat, int ecount,
        int *elist, int stallcount, int repeatcount, int *incycle, double *val, int silent, double time_bound,
        double length_bound, char *saveit_name, int kicktype,
        CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol, mbfs_main *main,
        logisticRegressor *classifier, double *objective, double *capacity, int mode),
    CClinkern_path (int ncount, CCdatagroup *dat, int ecount,
        int *elist, int nkicks, int *inpath, int *outpath, double *val,
        int silent, CCrandstate *rstate),
    CClinkern_fixed (int ncount, CCdatagroup *dat, int ecount, int *elist,
        int nkicks, int *incycle, int *outcycle, double *val, int fcount,
        int *flist, int silent, CCrandstate *rstate);

ttpObjective do_pgch(mbfs_main *main, ttpInstance *instance, ttpSolution *sol, logisticRegressor *classifier, double timebound, double *capacity, int maxIters);

#endif  /* __LINKERN_H */



/****************************************************************************/
/*                                                                          */
/*                  FLIPPER HEADER (TWO-LIST)                               */
/*                                                                          */
/****************************************************************************/


#ifndef __FLIPPER_H
#define __FLIPPER_H

typedef struct CClk_parentnode {
    struct CClk_parentnode *adj[2];
    struct CClk_childnode  *ends[2];
    int                     size;
    int                     id;
    int                     rev;
} CClk_parentnode;

typedef struct CClk_childnode {
    struct CClk_parentnode *parent;
    struct CClk_childnode  *adj[2];
    int                     id;
    int                     name;
} CClk_childnode;

typedef struct CClk_flipper {
    CClk_parentnode        *parents;
    CClk_childnode         *children;
    int                     reversed;
    int                     nsegments;
    int                     groupsize;
    int                     split_cutoff;
} CClk_flipper;



int
    CClinkern_flipper_init (CClk_flipper *f, int ncount, int *cyc),
    CClinkern_flipper_next (CClk_flipper *f, int x),
    CClinkern_flipper_prev (CClk_flipper *f, int x),
    CClinkern_flipper_sequence (CClk_flipper *f, int x, int y, int z);
void
    CClinkern_flipper_flip (CClk_flipper *F, int x, int y),
    CClinkern_flipper_cycle (CClk_flipper *F, int *x),
    CClinkern_flipper_finish (CClk_flipper *F);

ttpObjective evaluateFromFlipper(ttpInstance *instance, ttpSolution *sol, CClk_flipper *F), 
    evaluateFromFlipperIndices(ttpInstance *instance, ttpSolution *sol, CClk_flipper *F, CClk_childnode *start, CClk_childnode *stop, double startWeight);

#endif  /* __FLIPPER_H */
