#include "../include/Rss_Op.h"
// extern "C" {
// #include "../aes_ni.h"
// }

void Rss_eqz3(Lint **res, Lint **a, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    uint i, j, k, index; // used for loops
    Lint a1 = 0, a2 = 0;
    switch (pid) {
    case 1:
        a1 = 1;
        a2 = 0;
        break;
    case 2:
        a1 = 0;
        a2 = 0;
        break;
    case 3:
        a1 = 0;
        a2 = 1;
        break;
    }

    Lint **edaBit_r = new Lint *[2];
    Lint **edaBit_b_2 = new Lint *[2];
    Lint **sum = new Lint *[2];
    Lint **a_2 = new Lint *[2];

    // will hold the k shares of a in (k_prime)
    Lint *c = new Lint[size];
   
    for (i = 0; i < 2; i++) {
        a_2[i] = new Lint[size];
        sum[i] = new Lint[size];

        edaBit_r[i] = new Lint[size];
        edaBit_b_2[i] = new Lint[size];
    }

    Rss_edaBit(edaBit_r, edaBit_b_2, size, ring_size, map, nodeNet);

    for (i = 0; i < size; i++) {
        sum[0][i] = a[0][i] + edaBit_r[0][i];
        sum[1][i] = a[1][i] + edaBit_r[1][i];
    }

    Rss_Open(c, sum, size, map, ring_size, nodeNet);

    for(i = 0; i < size; ++i) {
      sum[0][i] = c[i] ^ edaBit_b_2[0][i];
      sum[1][i] = c[i] ^ edaBit_b_2[1][i];
    }

    Rss_kor(a_2, sum, ring_size, size, ring_size, map, nodeNet);

    Rss_b2a3(res, a_2, ring_size, size, map, nodeNet);

    for(i = 0; i < size; ++i){
      res[0][i] = a1 - res[0][i];
      res[1][i] = a2 - res[1][i];
    }
    
    delete[] c;

    for (i = 0; i < 2; i++) {
        delete[] edaBit_r[i];
        delete[] edaBit_b_2[i];
        delete[] sum[i];
        delete[] a_2[i];
    }
    delete[] edaBit_r;
    delete[] edaBit_b_2;
    delete[] sum;
    delete[] a_2;
}

