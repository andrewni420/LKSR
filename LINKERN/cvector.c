#include <stdbool.h>
#include "util.h"
#include "ttp_util.h"

#define INITIAL_VECTOR_CAPACITY 16

void intvector_init(intvector *vector){
    vector->len = 0;
    vector->cap = 0;
    vector->arr = (int *) NULL;
}

int intvector_alloc(intvector *vector, int cap){
    vector->len = 0;
    vector->cap = cap;
    vector->arr = CC_SAFE_MALLOC(cap, int);
    if (vector->arr==(int *) NULL) goto FAILURE;
    return 0;

  FAILURE:
    return 1;
}

int intvector_alloc_default(intvector *vector){
    return intvector_alloc(vector, INITIAL_VECTOR_CAPACITY);
}

int intvector_grow(intvector *vector){
    int cap = 2*vector->cap;
    int *arr = (int *) realloc(vector->arr,cap * sizeof(int));
    if (arr == (int *) NULL) goto FAILURE;
    vector->arr = arr;
    vector->cap=cap;
    return 0;
  FAILURE:
    return 1;
}

int intvector_conj(intvector *vector, int val){
    vector->arr[vector->len++] = val;
    if (vector->len>=vector->cap) return intvector_grow(vector);
    return 0;
}

int intvector_append(intvector *vector, int val){return intvector_conj(vector,val);}

int intvector_pop(intvector *vector){
    int ret = vector->arr[vector->len-1];
    vector->len--;
    return ret;
}

int intvector_concat(intvector *vector, intvector *other){
    int total_len = vector->len+other->len;
    int *arr = vector->arr;
    if (total_len>vector->cap) {
        vector->arr = (int *) CCutil_reallocrus(vector->arr,total_len*sizeof(int));
        if (vector->arr == (int *) NULL) goto FAILURE;
    }
    for (int i=0;i<other->len;i++) vector->arr[i+vector->len] = other->arr[i];
    return 0;

  FAILURE:
    vector->arr = arr;
    return 1;
}

int intvector_extend(intvector *vector, intvector *other){return intvector_concat(vector,other);}

int intvector_copy(intvector *dest, intvector *src){
    int *arr = dest->arr;
    if (dest->cap<=src->len){
        dest->arr = CCutil_reallocrus(dest->arr, src->cap * sizeof(int));
        if (dest->arr == (int *) NULL) goto FAILURE;
        dest->cap = src->cap;
    }
    memcpy(dest->arr,src->arr,src->len*sizeof(int));
    dest->len = src->len;
    return 0;
  FAILURE:
    dest->arr = arr;
    return 1;
}

int intvector_clear(intvector *vector){
    vector->len=0;
    return 0;
}

int intvector_peek(intvector *vector){
    return vector->arr[vector->len-1];
}

int intvector_reverse(intvector *vector){
    int i=vector->len-1;
    int j=0;
    while (i>j){
        int temp = vector->arr[i];
        vector->arr[i]=vector->arr[j];
        vector->arr[j]=temp;
        i--;
        j++;
    }
    return 0;
}

int intvector_set(intvector *vector, int idx, int val){
    int *arr = vector->arr;
    if (idx>=vector->cap-1){
        vector->arr = CCutil_reallocrus(vector->arr, 2*(idx+1)*sizeof(int));
        if (vector->arr == (int *) NULL) goto FAILURE;
        vector->cap = 2*(idx+1);
    }
    vector->arr[idx]=val;
    if (idx>=vector->len) {
        for (int i=vector->len;i<idx;i++) vector->arr[i]=0;
        vector->len=idx+1;
    }
    return 0;
  FAILURE:
    vector->arr = arr;
    return 1;
}

int intvector_popidx(intvector *vector, intvector *loc, int idx){
    //arr = idx1:item1, idx2:item2
    //loc = item1:idx1, item2:idx2
    //idx:itemLast
    //itemLast:idx
    //item:-1
    int ret = 0;
    ret+=intvector_set(loc,vector->arr[idx],-1);
    vector->arr[idx] = intvector_peek(vector);
    ret+=intvector_set(loc,vector->arr[idx],-1);
    ret+=intvector_pop(vector);
    return ret;
}

int intvector_remove(intvector *vector, intvector *loc, int item){
    return intvector_popidx(vector,loc,loc->arr[item]);
}

// Is the array [item0, item1, ..., i, j] sorted in ascending/descending order?
bool _intvector_bubble_helper(int i, int j, bool ascending){
    if (ascending) return i<=j;
    else return i>=j;
}

int intvector_bubble(intvector *vector, bool ascending){
    int i = vector->len-2;
    int j = vector->len-1;
    int temp;
    while (i>=0 && !_intvector_bubble_helper(vector->arr[i], vector->arr[j], ascending)){
        temp = vector->arr[i];
        vector->arr[i--]=vector->arr[j];
        vector->arr[j--]=temp;
    }
}


void intvector_qsort(intvector *vector, __compar_fn_t fn){
    qsort(vector->arr,vector->len, sizeof(int), fn);
}
void intvector_print(intvector *vector, FILE *outfile){
    fprintf(outfile, "[");
    for (int i=0;i<vector->len-1;i++) fprintf(outfile, "%d, ",vector->arr[i]);
    if (vector->len>0) fprintf(outfile, "%d]",vector->arr[vector->len-1]);
    else fprintf(outfile, "]");
}

void intvector_free(intvector *vector){
    CC_IFFREE(vector->arr,int);
    vector->arr = (int *) NULL;
}


