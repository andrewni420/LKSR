

double sigmoid(double input), randNormal(double mu, double sigma), randUnif(double a, double b), 
    quickselect(double *arr, int n, int k), median(double *arr, int n), mad(double *arr, int n);

void fisherYates(int *arr, int n), range(int *arr, int n), normalizeMAD(double *arr, int n), 
    printIntArray(int *arr, int n), printDoubleArray(double *arr, int n), copyDouble(double *src, double *dest, int n),
    copyInt(int *src, int *dest, int n), copyBool(bool *src, bool *dest, int n), shuffle(int *array, size_t n),
    logrange(double *arr, int n, double loglow, double loghigh), reverseInt(int *src, int n), printBoolArray(bool *arr, int n),
    scaleQuartile(double *arr, int n, double q);

int compareInt (const void * a, const void * b);