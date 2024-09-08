#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "util.h"
#include "ttp_util.h"
#include "ttp_heuristic.h"
#include "kdtree.h"
#include "linkern.h"
#include "ttp_optimisation.h"
#include "ttp.h"

#define WEIGHT_IDX 1
#define PROFIT_IDX 2 

static int norm = CC_EUCLIDEAN;
static int seed = 0;
static int binary_in = 0;
static int binary_edges = 0;
static int tsplib_in = 1;
static int nnodes_want = 0;
static int gridsize = 0;
static int nearnum = 0;
static int quadtry = 2;
static int run_silently = 0;
static int in_repeater = -1;
static int number_runs = 0;
static double time_bound = -1.0;
static double length_bound = -1.0;
static int kick_type = CC_LK_WALK_KICK;

// int main(int argc, char** argv) {
//     ttpInstance instance;
//     initInstance(&instance);
//     read_ttp_text("a280_n279_bounded-strongly-corr_01.ttp",&instance);
//     printf("Hello World %d\n",10);
//     return 0;
// }

int ttpInit(CCdatagroup *dat, ttpInstance *instance, ttpSolution *sol, CCrandstate *rstate, char *fileName){
    zeit_initial = CCutil_real_zeit ();

    CCutil_init_datagroup(dat);
    int seed = (int) CCutil_real_zeit ();
    CCutil_sprand (seed, rstate);
    int ncount; 
    CCutil_gettsplib (fileName, &ncount, dat);
    CCutil_dat_getnorm (dat, &norm);
    int *incycle = CC_SAFE_MALLOC (ncount, int);
    sol->tour = CC_SAFE_MALLOC(ncount+1,int);

    double val;
    int tempcount, *templist;

    CCkdtree localkt;
    CCkdtree_build (&localkt, ncount, dat, (double *) NULL, rstate);
    CCkdtree_qboruvka_tour (&localkt, ncount, dat, incycle, &val, rstate);

    CCkdtree_quadrant_k_nearest (&localkt, ncount, quadtry,
                           dat, (double *) NULL, 1, &tempcount, &templist,
                           run_silently, rstate);

    if (in_repeater == -1) in_repeater = ncount;

    CClinkern_tour (ncount, dat, tempcount, templist, 100000000,
                   in_repeater, incycle, sol->tour, &val, run_silently,
                   time_bound, length_bound, (char *) NULL, kick_type,
                   rstate, NULL, NULL, 0);

    CC_IFFREE(incycle, int);

    sol->tour[ncount]=0;

    initInstance(instance);
    read_ttp_tsplib(fileName, dat, rstate, instance);
    return 0;
}

int cmp_int(const void *a, const void *b)
{
  return * (int *)a - * (int *)b;
}

int ttpOutput(ttpInstance *instance, ttpSolution *sol, char *ttpName, char *algName, FILE *file_override){
    int cur_time = (int) CCutil_real_zeit ();
    char *fileName = CC_SAFE_MALLOC(100, char);
    int ret;
    ret = sprintf(fileName, "%s.%s.%d", ttpName, algName, global_outtime);
    if (file_override == (FILE *) NULL) printf("%s:\n",fileName);
    if (ret<0){
        printf("Error formatting filename\n");
        return ret;
    } 

    FILE *fptr;
    if (file_override!=(FILE *) NULL){
        fptr = file_override;
    } else{
        fptr = fopen(fileName, "w");
    }
    

    fprintf(fptr,"[");
    bool first = true;
    for (int i=0;i<instance->numberOfNodes;i++){
        if (!first) fprintf(fptr,", ");
        fprintf(fptr,"%d",sol->tour[i]+1);
        first = false;
    }
    fprintf(fptr,"]\n");

    intvector collected;
    intvector_init(&collected);
    intvector_alloc(&collected,instance->numberOfItems);

    for (int i=1;i<instance->numberOfNodes;i++){
        for (int j=0;j<instance->itemsPerCity;j++){
            int itemIdx = getItemIndex(instance, sol->tour[i],j);
            int idx = sol->isCityPacking ? itemIdx : getPackingIndex(instance,i,j);
            if (sol->packing[idx]){
                intvector_append(&collected, itemIdx+1);
            }
            
        }
    }

    intvector_qsort(&collected, cmp_int);
    intvector_print(&collected, fptr);
    fprintf(fptr,"\n");

    if (file_override == (FILE *) NULL) fclose(fptr);

    CC_IFFREE(fileName, char);
    intvector_free(&collected);
    
}