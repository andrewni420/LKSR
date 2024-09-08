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
/*                  FLIPPER HEADER (TWO-LIST)                               */
/*                                                                          */
/****************************************************************************/

#define MAK_MORTON   
#undef  NODE_INSERTIONS

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

#endif  /* __FLIPPER_H */


#ifndef  __LINKERN_H
#define  __LINKERN_H

#include "util.h"

#define CC_LK_RANDOM_KICK    (0)
#define CC_LK_GEOMETRIC_KICK (1)
#define CC_LK_CLOSE_KICK     (2)
#define CC_LK_WALK_KICK      (3)


typedef struct edge {
    int other;
    int weight;
} edge;

typedef struct edgelook {
    struct edgelook *next;
    int other;
    int diff;
    int over;
    int seq;
    int side;
#ifdef MAK_MORTON
    int mm;
#endif
#ifdef NODE_INSERTIONS
    int ni;
    int under;
#endif
} edgelook;

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct flippair {
    int firstprev;
    int first;
    int last;
    int lastnext;
} flippair;

typedef struct flipstack {
    flippair *stack;
    int counter;
    int max;
} flipstack;

typedef struct graph {
    edge **goodlist;
    edge *edgespace;
    int  *degree;
    int  *weirdmark;
    int   weirdmagic;
    int   ncount;
    CCrandstate *rstate;
} graph;

typedef struct distobj {
    CCdatagroup *dat;
    int       *cacheval;
    int       *cacheind;
    int        cacheM;
} distobj;

typedef struct adddel {
    char *add_edges;
    char *del_edges;
} adddel;

typedef struct aqueue {
    char *active;
    intptr *active_queue;
    intptr *bottom_active_queue;
    CCdheap *h;
} aqueue;


int
    CClinkern_tour (int ncount, CCdatagroup *dat, int ecount,
        int *elist, int stallcount, int repeatcount, int *incycle,
        int *outcycle, double *val, int silent, double time_bound,
        double length_bound, char *saveit_name, int kicktype,
        CCrandstate *rstate, int *outIters, int **outtours, int numOutputs),
    CClinkern_path (int ncount, CCdatagroup *dat, int ecount,
        int *elist, int nkicks, int *inpath, int *outpath, double *val,
        int silent, CCrandstate *rstate),
    CClinkern_fixed (int ncount, CCdatagroup *dat, int ecount, int *elist,
        int nkicks, int *incycle, int *outcycle, double *val, int fcount,
        int *flist, int silent, CCrandstate *rstate);
void 
    lin_kernighan (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       double *val, int *win_cycle, flipstack *w, flipstack *fstack,
       CCptrworld *intptr_world, CCptrworld *edgelook_world);

#endif  /* __LINKERN_H */

