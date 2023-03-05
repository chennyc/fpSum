#include "../include/util.h"

void print_binary(Lint n, uint size){
    uint temp = size-1;
    int i = size-1;
    uint b;
    while(i !=-1) {
        b = GET_BIT(n, temp);
        printf("%u", b);
        temp--;
        i -= 1;
    }
    printf("\n");
}