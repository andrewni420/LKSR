#include "util.h"
#include "ttp_util.h"


typedef struct improvement{
    double objective;
    int iteration;
    double time;
    improvement *next;
} improvement;

typedef struct opt_params{
    //Want to be able to interpolate between purely CClinkern + occasional checks vs check every step
    //Interpolate between no improvements ever (classifier/linear) to improvement checks at every opt (step?)

    //no check < check < improve

    // SUCCESSFUL PARADIGMS:

    //a280_n279 (the OG)
    //  - check at every step
    //  - improve (10 linear, 100/1000 classifier) after every improve_tour
    //  - improve after every kick

    //Linear only (fast)
    //  - no checks
    //  - linear 10 after every CClinkern

    //Check only (fast)
    //  - check at every step
    //  - no improvements

    //fastest to most accurate:
    //  no checking, basically CClinkern_tour
    //  Checking after some (all) lin_kernighans (improve_tour) (step)
    //  Linear (and classifier) after some (all) lin_kernighans (improve_tour) (step)
    //No point in taking something lower in the list without also having something higher in the list.
    //Mode = how far down in the list are we? Fractional part = probability in [some]? Should range over multiple orders of magnitude though

    //Maybe switch modes based on improvement duration? If it takes more than 10k then switch modes?

    improvement *improvementCurve;

    //How many checks have we done?
    int numChecks;

    //How many k-opt moves have been beneficial?
    int beneficialOpts;

    int numLinears;
    int linearImprovements;

    int numClassifiers;
    int classifierImprovements;

    double timeLR;
    double linearTime;
    double hillClimberTime;
    double stepTime;
    double checkTime;
    double CClinkernTime;
    double ttpLinkernTime;

    double probCheckOpt;
    double probImproveOpt;

    double probCheckStep;
    double probImproveStep;

    double probImproveLinkern;

    int classifierGens;
    int classifierNonImprovementDuration;

    int linearGens;

    double linear;
    double percent;

    ttpInstance *instance;

    int mode;
    //Mode switch activated by simple improvement duration checking
    //0 = no checking

    //1 = check after some/all lin kernighan
    //2 = linear after some/all lin kernighan
    //3 = 10/100 classifier after some/all lin kernighan
    //4 = 100/1000 classifier after some/all lin kernighan

    //5 = check after some/all improve_tour 
    //6 = linear after some/all improve_tour
    //7 = 10/100 classifier after some/all improve_tour
    //8 = 100/1000 classifier after some/all improve_tour

    //9 = check after some/all step
    //10 = linear after some/all step
    //11 = 10/100 classifier after some/all step
    //12 = 100/1000 classifier after some/all step
    
    //What governs the threshold for checking steps? - nope, just a flat out percentage. 
    //What governs improvment after each kick?
    //What governs the possibility of genetic algoritms? - I feel like this should be reserved for just finding the optimal solutions to the algorithms. using MPI
    int nonImprovementDuration;
    double nonImprovementTime;
    double startTime;

} opt_params;

int init_meta_params(opt_params *params, ttpInstance *instance), update_improvement(double objective, int iteration),
    json_improvement_curve(opt_params *params, FILE *outfile);

