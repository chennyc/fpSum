#include "../include/Rss_Op.h"

void Rss_b2a(Lint **res, Lint **a, uint ring_size, uint size, int *map, NodeNetwork *nodeNet) {
    int pid = nodeNet->getID();
    uint i;
    Lint **b = new Lint *[2];
    Lint **b_2 = new Lint *[2];
    Lint **sum = new Lint *[2];

    Lint *c = new Lint[size];

    for (i = 0; i < 2; i++) {
        b[i] = new Lint[size];
        b_2[i] = new Lint[size];
        sum[i] = new Lint[size];
    }

    Lint a1 = 0;
    Lint a2 = 0;
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
    Rss_RandBit(b, size, ring_size, map, nodeNet);

    for (i = 0; i < size; i++) {
        b_2[0][i] = GET_BIT(b[0][i], Lint(0));
        b_2[1][i] = GET_BIT(b[1][i], Lint(0));
    }

    for (i = 0; i < size; i++) {
        sum[0][i] = (a[0][i] ^ b_2[0][i]);
        sum[1][i] = (a[1][i] ^ b_2[1][i]);
    }

    Rss_Open_Bitwise(c, sum, size, map, ring_size, nodeNet);

    for (i = 0; i < size; i++) {
        res[0][i] = a1 * c[i] + b[0][i] - 2 * c[i] * b[0][i];
        res[1][i] = a2 * c[i] + b[1][i] - 2 * c[i] * b[1][i];
    }

    delete[] c;

    for (i = 0; i < 2; i++) {
        delete[] b[i];
        delete[] b_2[i];
        delete[] sum[i];
    }
    delete[] b;
    delete[] b_2;
    delete[] sum;
}