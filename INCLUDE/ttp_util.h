

#ifndef TTP_UTIL_H

#define TTP_UTIL_H

#define CHECK_TIME() (zeit_initial>0.0 && zeit_bound > 0.0 && CCutil_real_zeit() - zeit_initial > zeit_bound)
#define GET_SPEED(instance, weight) ((instance)->maxSpeed-weight*(instance)->speed_factor)
#define GET_DISTANCE(res, instance, i, j) \
    do { \
        int __i__ = i; int __j__ = j; \
        double __t1__ = instance->x[__i__] - instance->x[__j__], __t2__ = instance->y[__i__] - instance->y[__j__]; \
        res = ceil (sqrt (__t1__ * __t1__ + __t2__ * __t2__)); \
    } while (0) 

#define SWAP(a,b, type) { \
    do { \
        type __temp__ = a; \
        a = b; \
        b = __temp__; \
    } while (0);}
#define MAXLKITERS 3000
#define NUMLKITERS(n) (n)>pow(MAXLKITERS,4./3.)?MAXLKITERS:pow((n),3./4.)

typedef struct intvector{
    int len;
    int cap;
    int *arr;
} intvector;

typedef struct ttpInstance {
    double *x;
    double *y;
    double *weight;
    double *profit;
    int numberOfNodes;
    int itemsPerCity;
    int numberOfItems;
    double maxSpeed;
    double minSpeed;
    double capacityOfKnapsack;
    double rentingRatio;
    double *profitRatios;
    double *normalizedProfitRatios;
    //Maps city to reverse distance
    double *rDist;
    //Maps city to normalized reverse distance
    double *normalizedRDist;
    int knapsackType;
    double capacityFactor;
    intvector *neighbors;

    int **itemIndex;
    int **packingIndex;
    double speed_factor;
} ttpInstance;

typedef struct ttpSolution{
    int *tour;
    bool *packing;
    bool isCityPacking;
} ttpSolution;

typedef struct ttpObjective
{
    double distance;
    double time;
    double profit;
    double finalweight;
    double objective;
} ttpObjective;

typedef struct doubleArrayEntry{
    double value;
    int index;
} doubleArrayEntry;




/* Fields
 * double probImprove;
 * double Gmax_threshold;
 * double check_objective_frequency;
 * int upper_quadratic_max_nonimprovement;
 * int upper_quadratic_max;
 * int upper_mbfs_max; 
 * int quadratic_max_nonimprovement;
 * int quadratic_max;
 * int mbfs_max; 
 * int numOutputs;
 * int *outIters;
 */
typedef struct adaptive_params{
    double probImprove;
    double Gmax_threshold;
    double check_objective_frequency;
    int upper_quadratic_max_nonimprovement;
    int upper_quadratic_max;
    int upper_mbfs_max; 
    bool upper_mbfs_sorted;
    int quadratic_max_nonimprovement;
    int quadratic_max;
    int mbfs_max; 
    bool mbfs_sorted;
    bool fast_pgch;
    int initial_pgch_max;
    int numOutputs;
    int normal_rounds;
    int ttp_rounds;
    int pgch_rounds;
    int kicktype;
    int *outIters;
    int pgch_round_counter;
    int pgch_nonimprovement_counter;
    int normal_nonimprovement_counter;
    double upper_improvement_prob;
} adaptive_params;

/****** Intvector *****/
int intvector_extend(intvector *vector, intvector *other), intvector_concat(intvector *vector, intvector *other), 
    intvector_append(intvector *vector, int val), intvector_conj(intvector *vector, int val), intvector_pop(intvector *vector),
    intvector_grow(intvector *vector), intvector_alloc_default(intvector *vector), intvector_alloc(intvector *vector, int cap),
    intvector_copy(intvector *dest, intvector *src), intvector_clear(intvector *vector), intvector_reverse(intvector *vector),
    intvector_peek(intvector *vector), intvector_popidx(intvector *vector, intvector *loc, int idx), intvector_remove(intvector *vector, intvector *loc, int item),
    intvector_bubble(intvector *vector, bool ascending), intvector_set(intvector *vector, int idx, int val);

void intvector_init(intvector *vector), intvector_free(intvector *vector), intvector_print(intvector *vector, FILE *outfile),
    intvector_qsort(intvector *vector, __compar_fn_t fn);
/*********************/

int objective_copy(ttpObjective *src, ttpObjective *dest);
int getItemIndex(ttpInstance *instance, int city, int itemNumber), tourIndexFromPackingIndex(ttpInstance *instance, int packingIndex),
    getPackingIndex(ttpInstance *instance, int tourIdx, int itemNumber), itemNumberFromPackingIndex(ttpInstance *instance, int packingIndex),
    getCityPackingIndex(ttpInstance *instance, int city, int itemNumber), getCityFromItemIndex(ttpInstance *instance, int itemIndex),
    comparingDouble(const void *a, const void *b), comparingDoubleAsc(const void *a, const void *b), itemNumberFromItemIndex(ttpInstance *instance, int itemIndex),
    calcRDist(ttpInstance *instance, int *tour);

bool isPacked(ttpInstance *instance, ttpSolution *sol, int tourIdx, int itemNumber);

double * reverseDistances(ttpInstance *instance, int * tour);

double getSpeed(ttpInstance *instance, double weight), getDistance(ttpInstance *instance, int i, int j);

void initInstance (ttpInstance *instance), freeInstance (ttpInstance *instance), initSolution(ttpSolution *sol), freeSolution(ttpSolution *sol),
     toCityPacking(ttpInstance *instance, ttpSolution *sol), toIndexPacking(ttpInstance *instance, ttpSolution *sol), init_params(adaptive_params *params),
     updateWithIndexPacking(ttpInstance *instance, ttpSolution *sol, bool *packing), objective_print(ttpObjective *obj), sol_print(ttpInstance *instance, ttpSolution *sol);

int read_ttp_text (char *datname, ttpInstance *instance), read_ttp_tsplib (char *datname, CCdatagroup *dat, CCrandstate *rstate, ttpInstance *instance), 
    copySolution(ttpInstance *instance, ttpSolution *dest, ttpSolution *src);

int knapsackType(ttpInstance *instance);
void set_params(ttpInstance *instance, adaptive_params *params, double time), free_params(adaptive_params *params);

ttpObjective evaluate(ttpInstance *instance, ttpSolution *sol), evaluateVerbose(ttpInstance *instance, ttpSolution *sol),
    evaluate_intermediate(ttpInstance *instance, ttpSolution *sol, ttpObjective *intermediate, int startIndex, bool update);

double capacityFactor(ttpInstance *instance);
double estimateParameter(ttpInstance *instance, int parameterIdx), stdvsmean(ttpInstance *instance, int parameterIdx);

int lk_tour(CCdatagroup *dat, CCrandstate *rstate, ttpInstance *instance, ttpSolution *sol, int *outIters, int **outtours, int numIters, int numOutputs, double timeBound);
int sol_from_file(ttpInstance *instance, ttpSolution *sol, char *fileName);

extern double zeit_initial;
extern double zeit_bound;
extern adaptive_params ttp_lk_params;

#endif 
