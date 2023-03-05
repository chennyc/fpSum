#include "../include/Rss_Op.h"

// REQUIREMENT: the MSB of x MUST BE ZERO
void Rss_truncPr(Lint **res, Lint **x, uint m, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    Lint i; // used for loops

    Lint **edaBit_r = new Lint *[2];
    Lint **sum = new Lint *[2];
    Lint **r_m_prime = new Lint *[2];
    Lint **b = new Lint *[2];
    Lint **r_km1 = new Lint *[2];

    Lint *c = new Lint[size];
    memset(c, 0, sizeof(Lint) * size);
    Lint *c_prime = new Lint[size];
    memset(c_prime, 0, sizeof(Lint) * size);

    for (i = 0; i < 2; i++) {
        edaBit_r[i] = new Lint[size];
        r_m_prime[i] = new Lint[size];
        r_km1[i] = new Lint[size];
        sum[i] = new Lint[size];

        b[i] = new Lint[size];
        memset(b[i], 0, sizeof(Lint) * size);
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

    // generating one edaBit
    Rss_edaBit_trunc(edaBit_r, r_m_prime, r_km1, size, ring_size, m, map, nodeNet);

    // computing the sum of x and edabit_r
    for (i = 0; i < size; i++) {
        sum[0][i] = (x[0][i] + edaBit_r[0][i]);
        sum[1][i] = (x[1][i] + edaBit_r[1][i]);
    }

    Rss_Open(c, sum, size, map, ring_size, nodeNet);

    for (i = 0; i < size; i++) {

        // (c / 2^m) mod 2^(k-m-1)
        c_prime[i] = (c[i] >> Lint(m)) & nodeNet->SHIFT[ring_size - m - 1];

        b[0][i] = r_km1[0][i] + ((c[i] * a1) >> Lint(ring_size - 1)) - 2 * ((c[i]) >> Lint(ring_size - 1)) * r_km1[0][i];
        b[1][i] = r_km1[1][i] + ((c[i] * a2) >> Lint(ring_size - 1)) - 2 * ((c[i]) >> Lint(ring_size - 1)) * r_km1[1][i];

        r_m_prime[0][i] = r_m_prime[0][i] - (r_km1[0][i] << Lint(ring_size - 1 - m));
        r_m_prime[1][i] = r_m_prime[1][i] - (r_km1[1][i] << Lint(ring_size - 1 - m));

        res[0][i] = (c_prime[i] * a1) - r_m_prime[0][i] + (b[0][i] << Lint(ring_size - m - 1));
        res[1][i] = (c_prime[i] * a2) - r_m_prime[1][i] + (b[1][i] << Lint(ring_size - m - 1));
    }

    delete[] c;
    delete[] c_prime;
    for (i = 0; i < 2; i++) {
        delete[] edaBit_r[i];
        delete[] sum[i];
        delete[] b[i];
        delete[] r_m_prime[i];
        delete[] r_km1[i];
    }
    delete[] edaBit_r;
    delete[] b;
    delete[] sum;
    delete[] r_m_prime;
    delete[] r_km1;
}

void Rss_truncPre(Lint **res, Lint **x, uint m, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    Lint i; // used for loops

    Lint **edaBit_r = new Lint *[2];
    Lint **r_m_prime = new Lint *[2];
    Lint **b = new Lint *[2];
    Lint **r_km1 = new Lint *[2];

    Lint *c = new Lint[size];
    memset(c, 0, sizeof(Lint) * size);
    Lint *c_prime = new Lint[size];
    Lint **bitLTres = new Lint *[2];
    
    memset(c_prime, 0, sizeof(Lint) * size);

    for (i = 0; i < 2; i++) {
        edaBit_r[i] = new Lint[size];
        r_m_prime[i] = new Lint[size];
        r_km1[i] = new Lint[size];

        b[i] = new Lint[size];
        memset(b[i], 0, sizeof(Lint) * size);

	bitLTres[i] = new Lint[size];
	memset(bitLTres[i], 0, sizeof(Lint) * size);
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

    // generating one edaBit
    Rss_edaBit_truncPre(x, c, bitLTres, edaBit_r, r_m_prime, r_km1, size, ring_size, m, map, nodeNet);

    // computing the sum of x and edabit_r

    for (i = 0; i < size; i++) {

        // (c / 2^m) mod 2^(k-m-1)
        c_prime[i] = (c[i] >> Lint(m)) & nodeNet->SHIFT[ring_size - m - 1];
	
        b[0][i] = r_km1[0][i] + ((c[i] * a1) >> Lint(ring_size - 1)) - 2 * ((c[i]) >> Lint(ring_size - 1)) * r_km1[0][i];
        b[1][i] = r_km1[1][i] + ((c[i] * a2) >> Lint(ring_size - 1)) - 2 * ((c[i]) >> Lint(ring_size - 1)) * r_km1[1][i];
    }

    for (i = 0; i < size; i++) {
        r_m_prime[0][i] = r_m_prime[0][i] - (r_km1[0][i] << Lint(ring_size - 1 - m));
        r_m_prime[1][i] = r_m_prime[1][i] - (r_km1[1][i] << Lint(ring_size - 1 - m));

        res[0][i] = (c_prime[i] * a1) - r_m_prime[0][i] + (b[0][i] << Lint(ring_size - m - 1)) - bitLTres[0][i];
        res[1][i] = (c_prime[i] * a2) - r_m_prime[1][i] + (b[1][i] << Lint(ring_size - m - 1)) - bitLTres[1][i];
    }

    delete[] c;
    delete[] c_prime;

    for (i = 0; i < 2; i++) {
        delete[] edaBit_r[i];
        delete[] b[i];
        delete[] r_m_prime[i];
        delete[] r_km1[i];
	delete[] bitLTres[i];
    }
    delete[] edaBit_r;
    delete[] b;
    delete[] bitLTres;
    delete[] r_m_prime;
    delete[] r_km1;
}



// computes [x/2^m] probabilistically
void truncPr_1(Lint **res, Lint **x, Lint m, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    Lint i, j, index, k; // used for loops

    uint n_rand_bits = size * (ring_size);

    Lint **edaBit_r = new Lint *[2];
    Lint **edaBit_b_2 = new Lint *[2];
    Lint **sum = new Lint *[2];
    Lint **r_m_prime = new Lint *[2];
    Lint **b = new Lint *[2];

    Lint *c = new Lint[size];
    memset(c, 0, sizeof(Lint) * size);
    Lint *c_prime = new Lint[size];
    memset(c_prime, 0, sizeof(Lint) * size);

    Lint *x_open = new Lint[size];
    Lint *final_open = new Lint[size];
    Lint *edaBit_r_open = new Lint[size];
    Lint *r_prime_open = new Lint[size];
    Lint *edaBit_b_open = new Lint[size];
    Lint *b_open = new Lint[size];

    for (i = 0; i < 2; i++) {
        edaBit_r[i] = new Lint[size];
        memset(edaBit_r[i], 0, sizeof(Lint) * size);
        edaBit_b_2[i] = new Lint[size];
        memset(edaBit_b_2[i], 0, sizeof(Lint) * size);
        sum[i] = new Lint[size];
        memset(sum[i], 0, sizeof(Lint) * size);

        r_m_prime[i] = new Lint[size];
        memset(r_m_prime[i], 0, sizeof(Lint) * size);
        b[i] = new Lint[size];
        memset(b[i], 0, sizeof(Lint) * size);
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
    Lint b1, b2;
    switch (pid) {
    case 1:
        b1 = -1;
        b2 = 0;
        break;
    case 2:
        b1 = 0;
        b2 = 0;
        break;
    case 3:
        b1 = 0;
        b2 = -1;
        break;
    }

    // generating one edaBit
    // Rss_edaBit(edaBit_r, edaBit_b_2, size, ring_size, ring_size, map, nodeNet);
    Rss_edaBit(edaBit_r, edaBit_b_2, size, ring_size, map, nodeNet);

    // for (i = 0; i < 2; i++) {
    //     memset(edaBit_r[i], 0, sizeof(Lint) * size);
    //     memset(edaBit_b_2[i], 0, sizeof(Lint) * size);
    // }

    // computing the sum of x and edabit_r
    for (i = 0; i < size; i++) {
        sum[0][i] = (x[0][i] + edaBit_r[0][i]);
        sum[1][i] = (x[1][i] + edaBit_r[1][i]);

        // sum[0][i] = edaBit_r[0][i] + x[0][i];
        // sum[1][i] = edaBit_r[1][i] + x[1][i];
    }

    Rss_Open(c, sum, size, map, ring_size, nodeNet);
    Rss_Open(x_open, x, size, map, ring_size, nodeNet);
    Rss_Open(edaBit_r_open, edaBit_r, size, map, ring_size, nodeNet);
    Rss_Open_Bitwise(edaBit_b_open, edaBit_b_2, size, map, ring_size, nodeNet);
    printf("\n");

    for (size_t i = 0; i < size; i++) {

        printf("c[%i]  : %llu\n", i, c[i]);
        // print_binary(c[i], ring_size);
        printf("x_open[%i] + edaBit_r  : %llu\n", i, x_open[i] + edaBit_r_open[i]);
        // print_binary(x_open[i], ring_size);

        printf("edaBit_r_open[%i]  : %llu\n", i, edaBit_r_open[i]);
        print_binary(edaBit_r_open[i], ring_size);

        printf("x_open[%i]  : %llu\n", i, x_open[i]);
        print_binary(x_open[i], ring_size);
    }

    for (i = 0; i < size; i++) {
        // (c / 2^m) mod 2^(k-m-1)
        // c_prime[i] = (c[i] >> Lint(m)) ;
        c_prime[i] = (c[i] >> Lint(m)) & nodeNet->SHIFT[ring_size - m - 1];

        // c_prime[i] = c[i]  & nodeNet->SHIFT[ring_size - m - 1];
        // c_prime[i] = c[i]  & nodeNet->SHIFT[m];
    }
    printf("\n");

    for (size_t i = 0; i < size; i++) {

        printf("c[%i]  : %llu\n", i, c[i]);
        print_binary(c[i], ring_size);
        printf("c_prime[%i]  : %llu\n", i, c_prime[i]);
        print_binary(c_prime[i], ring_size);

        printf("actual:[%i]  : %llu\n", i, c_prime[i]);

        printf("expected:  : %llu\n", (x_open[i] + edaBit_r_open[i]) >> Lint(m));
        // print_binary(x_open[i], ring_size);

        printf("edaBit_r_open[%i]  : %llu\n", i, edaBit_r_open[i]);
        print_binary(edaBit_r_open[i], ring_size);

        printf("x_open[%i]  : %llu\n", i, x_open[i]);
        print_binary(x_open[i], ring_size);
    }
    printf("\n");
    printf("\n");

    for (i = 0; i < size; i++) {

        // b[0][i] = GET_BIT(edaBit_b_2[0][i], Lint(ring_size - 1)) ^ ((c[i]*a1) >> Lint(ring_size - 1));
        // b[1][i] = GET_BIT(edaBit_b_2[1][i], Lint(ring_size - 1)) ^ ((c[i]*a2) >> Lint(ring_size - 1));
        printf("edaBit_b_2[0] : %llu\n", edaBit_b_2[0][i]);
        print_binary(edaBit_b_2[0][i], ring_size);

        printf("edaBit_b_2[1] : %llu\n", edaBit_b_2[1][i]);
        print_binary(edaBit_b_2[1][i], ring_size);

        printf("edaBit_b_2_(k-1)[0] : %llu\n", GET_BIT(edaBit_b_2[0][i], Lint(ring_size - 1)));
        printf("edaBit_b_2_(k-1)[1] : %llu\n", GET_BIT(edaBit_b_2[1][i], Lint(ring_size - 1)));

        printf("c[%i]  : %llu\n", i, c[i]);
        print_binary(c[i], ring_size);

        printf("c/2^k-1[%i]  : %llu\n", i, ((c[i]) >> Lint(ring_size - 1)));
        print_binary(((c[i]) >> Lint(ring_size - 1)), ring_size);
        printf("\n");

        // b[0][i] = GET_BIT(edaBit_b_2[0][i], Lint(ring_size - 1)) ^ ((c[i] ) >> Lint(ring_size - 1));
        // b[1][i] = GET_BIT(edaBit_b_2[1][i], Lint(ring_size - 1)) ^ ((c[i] ) >> Lint(ring_size - 1));
        b[0][i] = GET_BIT(edaBit_b_2[0][i], Lint(ring_size - 1)) ^ ((c[i] * a1) >> Lint(ring_size - 1));
        b[1][i] = GET_BIT(edaBit_b_2[1][i], Lint(ring_size - 1)) ^ ((c[i] * a2) >> Lint(ring_size - 1));

        // b[0][i] = GET_BIT(edaBit_b_2[0][i], Lint(ring_size - 1)) ^ ((c[i] & b1) >> Lint(ring_size - 1));
        // b[1][i] = GET_BIT(edaBit_b_2[1][i], Lint(ring_size - 1)) ^ ((c[i] & b2) >> Lint(ring_size - 1));
    }
    for (i = 0; i < size; i++) {
        printf("b[0] : %llu\n", b[0][i]);
        printf("b[1] : %llu\n", b[1][i]);
    }

    Rss_Open_Bitwise(b_open, b, size, map, ring_size, nodeNet);

    printf("\n");
    Lint temp;
    for (size_t i = 0; i < size; i++) {

        // this value is changing, which makes sense because edaBit_b is changing
        temp = GET_BIT(edaBit_b_open[i], Lint(ring_size - 1)) ^ ((c[i]) >> Lint(ring_size - 1));

        printf("edaBit_(k-1)[%i]  : %llu\n", i, GET_BIT(edaBit_b_open[i], Lint(ring_size - 1)));
        print_binary(edaBit_b_open[i], ring_size);

        // THIS SHOULD BE ZERO
        printf("c/2^k-1[%i]  : %llu\n", i, ((c[i]) >> Lint(ring_size - 1)));
        print_binary(((c[i]) >> Lint(ring_size - 1)), ring_size);
        printf("\n");

        printf("actual[%i]  : %llu\n", i, b_open[i]);
        print_binary(b_open[i], ring_size);
        printf("expected[%i]  : %llu\n", i, temp);
        print_binary(temp, ring_size);
        // print_binary(edaBit_r_open[i], ring_size);
    }
    printf("\n");

    for (i = 0; i < size; i++) {
        for (j = m; j < ring_size - 1; j++) {

            // computing the summation in the last step of the algorithm
            r_m_prime[0][i] = r_m_prime[0][i] + (GET_BIT(edaBit_b_2[0][i], Lint(j)) << Lint(j - m));
            r_m_prime[1][i] = r_m_prime[1][i] + (GET_BIT(edaBit_b_2[1][i], Lint(j)) << Lint(j - m));
        }
    }
    Rss_Open_Bitwise(r_prime_open, r_m_prime, size, map, ring_size, nodeNet);

    for (size_t i = 0; i < size; i++) {

        // printf("r_m_prime[0] : %llu\n", r_m_prime[0][i]);
        // print_binary(r_m_prime[0][i], ring_size);

        // printf("r_m_prime[1] : %llu\n", r_m_prime[1][i]);
        // print_binary(r_m_prime[1][i], ring_size);

        printf("\n");
        printf("r_m_prime[%i]  : %llu\n", i, r_prime_open[i]);
        print_binary(r_prime_open[i], ring_size);

        // printf("edaBit_b[%i]  : %llu\n", i, edaBit_b_open[i] );
        // print_binary(edaBit_b_open[i], ring_size);

        printf("expected[%i]  : %llu\n", i, ((edaBit_b_open[i] >> Lint(m) & nodeNet->SHIFT[ring_size - 1 - m])));
        print_binary((edaBit_b_open[i] >> Lint(m) & nodeNet->SHIFT[ring_size - 1 - m]), ring_size);
    }

    printf("\n");

    // for checking, clear lowest m bits of edabit_r_open, which should give the same as r_m_prime

    for (i = 0; i < size; ++i) {
        res[0][i] = (c_prime[i] * a1) - r_m_prime[0][i] + (b[0][i] << Lint(ring_size - m - 1));
        res[1][i] = (c_prime[i] * a2) - r_m_prime[1][i] + (b[1][i] << Lint(ring_size - m - 1));
        // a1/2 does not affect final output
        // res[0][i] = (c_prime[i] )  - r_m_prime[0][i] + (b[0][i] << Lint(ring_size - m - 1));
        // res[1][i] = (c_prime[i] )  - r_m_prime[1][i] + (b[1][i] << Lint(ring_size - m - 1));

        // res[0][i] =  r_m_prime[0][i] - ((c[i] >> m) * a1)  ;
        // res[1][i] =  r_m_prime[1][i] - ((c[i] >> m) * a2)  ;
        // res[0][i] = (-1) *((c[i] >> m) * a1)  + r_m_prime[0][i]  ;
        // res[1][i] = (-1) *((c[i] >> m) * a2)  + r_m_prime[1][i] ;
    }

    Rss_Open(final_open, res, size, map, ring_size, nodeNet);

    for (size_t i = 0; i < size; i++) {

        // printf("r_m_prime[0] : %llu\n", r_m_prime[0][i]);
        // print_binary(r_m_prime[0][i], ring_size);

        // printf("r_m_prime[1] : %llu\n", r_m_prime[1][i]);
        // print_binary(r_m_prime[1][i], ring_size);

        printf("\n");
        printf("actual[%i]  : %llu\n", i, final_open[i]);
        print_binary(final_open[i], ring_size);

        // printf("edaBit_b[%i]  : %llu\n", i, edaBit_b_open[i] );
        // print_binary(edaBit_b_open[i], ring_size);

        printf("expected[%i]  : %llu\n", i, c_prime[i] - r_prime_open[i] + (temp << Lint(ring_size - m - 1)));
        print_binary(c_prime[i] - r_prime_open[i] + (temp << Lint(ring_size - m - 1)), ring_size);
    }

    printf("\n");

    //cleanup
    delete[] x_open;
    delete[] b_open;
    delete[] final_open;
    delete[] edaBit_r_open;
    delete[] edaBit_b_open;
    delete[] r_prime_open;
    delete[] c;
    delete[] c_prime;
    for (i = 0; i < 2; i++) {
        delete[] edaBit_r[i];
        delete[] edaBit_b_2[i];
        delete[] sum[i];
        delete[] b[i];
        delete[] r_m_prime[i];
    }
    delete[] edaBit_r;
    delete[] edaBit_b_2;
    delete[] b;
    delete[] sum;
    delete[] r_m_prime;
}

// computes [x/2^m] probabilistically
void truncPr_2(Lint **res, Lint **x, Lint m, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    Lint i, j, index, k; // used for loops

    uint n_rand_bits = size * (ring_size);

    Lint **r_shares = new Lint *[2];
    for (i = 0; i < 2; i++) {
        r_shares[i] = new Lint[n_rand_bits];
    }
    Rss_RandBit(r_shares, n_rand_bits, ring_size, map, nodeNet);

    Lint **sum = new Lint *[2];
    Lint **r_m_prime = new Lint *[2];
    Lint **b = new Lint *[2];

    Lint **r = new Lint *[2];
    Lint **rprime = new Lint *[2];

    Lint *c = new Lint[size];
    memset(c, 0, sizeof(Lint) * size);
    Lint *c_prime = new Lint[size];
    memset(c_prime, 0, sizeof(Lint) * size);
    for (i = 0; i < 2; i++) {
        r[i] = new Lint[size];
        memset(r[i], 0, sizeof(Lint) * size);
        sum[i] = new Lint[size];

        rprime[i] = new Lint[size];
        memset(rprime[i], 0, sizeof(Lint) * size);

        r_m_prime[i] = new Lint[size];
        memset(r_m_prime[i], 0, sizeof(Lint) * size);
        b[i] = new Lint[size];
        memset(b[i], 0, sizeof(Lint) * size);
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

    for (j = 0; j < size; j++) {
        for (k = 0; k < ring_size; k++) {
            index = j * ring_size + k;
            r[0][j] = r[0][j] + (r_shares[0][index] << Lint(k));
            r[1][j] = r[1][j] + (r_shares[1][index] << Lint(k));
        }
    }

    for (j = 0; j < size; j++) {
        // should NOT be ring_size - 1
        for (k = m; k < ring_size; k++) {
            // this is for step 3
            index = j * ring_size + k;
            rprime[0][j] = rprime[0][j] + (r_shares[0][index] << Lint(k - m));
            rprime[1][j] = rprime[1][j] + (r_shares[1][index] << Lint(k - m));
        }
    }

    // computing the sum of x and edabit_r
    for (i = 0; i < size; i++) {

        sum[0][i] = r[0][i] - x[0][i];
        sum[1][i] = r[1][i] - x[1][i];

        // sum[0][i] = r[0][i] + x[0][i];
        // sum[1][i] = r[1][i] + x[1][i];
    }

    Rss_Open(c, sum, size, map, ring_size, nodeNet);

    for (i = 0; i < size; i++) {
        // (c / 2^m) mod 2^(k-m-1)
        // c_prime[i] = (c[i] >> Lint(m)) ;
        // c_prime[i] = (c[i] >> Lint(m)) & nodeNet->SHIFT[ring_size - m - 1];
        // c_prime[i] = c[i]  & nodeNet->SHIFT[ring_size - m - 1];
        c_prime[i] = c[i] & nodeNet->SHIFT[m];
    }

    for (i = 0; i < size; i++) {
        index = i * ring_size + ring_size - 1;
        printf("index: %u\n", index);
        b[0][i] = r_shares[0][index] ^ ((c[i] * a1) >> Lint(ring_size - 1));
        b[1][i] = r_shares[1][index] ^ ((c[i] * a2) >> Lint(ring_size - 1));
    }

    for (i = 0; i < size; ++i) {

        // res[0][i] = (c_prime[i] * a1)  - rprime[0][i] + (b[0][i] << Lint(ring_size - m - 1));
        // res[1][i] = (c_prime[i] * a2)  - rprime[1][i] + (b[1][i] << Lint(ring_size - m - 1));
        res[0][i] = (-1) * ((c[i] >> m) * a1) + rprime[0][i];
        res[1][i] = (-1) * ((c[i] >> m) * a2) + rprime[1][i];
    }

    //cleanup
    delete[] c;
    delete[] c_prime;
    for (i = 0; i < 2; i++) {
        delete[] sum[i];
        delete[] b[i];
        delete[] r_m_prime[i];
        delete[] r[i];
        delete[] rprime[i];
    }
    delete[] b;
    delete[] sum;
    delete[] r_m_prime;
    delete[] r;
    delete[] rprime;
}

// computes [x/2^m] probabilistically
void truncPr_2_test(Lint **res, Lint **x, Lint m, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    Lint i, j, index, k; // used for loops

    uint n_rand_bits = size * (ring_size);

    Lint **r_shares = new Lint *[2];
    for (i = 0; i < 2; i++) {
        r_shares[i] = new Lint[n_rand_bits];
    }
    Rss_RandBit(r_shares, n_rand_bits, ring_size, map, nodeNet);
    Lint **edaBit_r = new Lint *[2];
    Lint **edaBit_b_2 = new Lint *[2];

    Lint **sum = new Lint *[2];
    Lint **r_m_prime = new Lint *[2];
    Lint **b = new Lint *[2];

    Lint **r = new Lint *[2];
    Lint **rprime = new Lint *[2];

    Lint *c = new Lint[size];
    memset(c, 0, sizeof(Lint) * size);
    Lint *c_prime = new Lint[size];
    memset(c_prime, 0, sizeof(Lint) * size);
    for (i = 0; i < 2; i++) {

        edaBit_r[i] = new Lint[size];
        memset(edaBit_r[i], 0, sizeof(Lint) * size);
        edaBit_b_2[i] = new Lint[size];
        memset(edaBit_b_2[i], 0, sizeof(Lint) * size);
        r[i] = new Lint[size];
        memset(r[i], 0, sizeof(Lint) * size);
        sum[i] = new Lint[size];

        rprime[i] = new Lint[size];
        memset(rprime[i], 0, sizeof(Lint) * size);

        r_m_prime[i] = new Lint[size];
        memset(r_m_prime[i], 0, sizeof(Lint) * size);
        b[i] = new Lint[size];
        memset(b[i], 0, sizeof(Lint) * size);
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

    Rss_edaBit(edaBit_r, edaBit_b_2, size, ring_size, map, nodeNet);

    for (j = 0; j < size; j++) {
        for (k = 0; k < ring_size; k++) {
            index = j * ring_size + k;
            r[0][j] = r[0][j] + (r_shares[0][index] << Lint(k));
            r[1][j] = r[1][j] + (r_shares[1][index] << Lint(k));
        }
    }

    for (j = 0; j < size; j++) {
        // should NOT be ring_size - 1
        for (k = m; k < ring_size; k++) {
            // this is for step 3
            index = j * ring_size + k;
            rprime[0][j] = rprime[0][j] + (r_shares[0][index] << Lint(k - m));
            rprime[1][j] = rprime[1][j] + (r_shares[1][index] << Lint(k - m));
        }
    }

    // computing the sum of x and edabit_r
    for (i = 0; i < size; i++) {

        // sum[0][i] = r[0][i] - x[0][i];
        // sum[1][i] = r[1][i] - x[1][i];

        sum[0][i] = edaBit_r[0][i] - x[0][i];
        sum[1][i] = edaBit_r[1][i] - x[1][i];

        // sum[0][i] = r[0][i] + x[0][i];
        // sum[1][i] = r[1][i] + x[1][i];
    }

    Lint temp1, temp2;

    for (i = 0; i < size; i++) {
        temp1 = 0;
        temp2 = 0;
        // computing the lower half of the sum
        for (j = 0; j < m; j++) {

            temp1 = temp1 + GET_BIT(edaBit_b_2[0][i], Lint(j)) << Lint(j);
            temp2 = temp2 + GET_BIT(edaBit_b_2[1][i], Lint(j)) << Lint(j);
        }

        r_m_prime[0][i] = (edaBit_r[0][i] - temp1) >> Lint(m);
        r_m_prime[1][i] = (edaBit_r[1][i] - temp2) >> Lint(m);
    }

    // for (i = 0; i < size; i++) {
    //     for (j = m; j < ring_size ; j++) {

    //         // computing the summation in the last step of the algorithm
    //         r_m_prime[0][i] = r_m_prime[0][i] + (GET_BIT(edaBit_b_2[0][i], Lint(j)) << Lint(j - m));
    //         r_m_prime[1][i] = r_m_prime[1][i] + (GET_BIT(edaBit_b_2[1][i], Lint(j)) << Lint(j - m));
    //     }
    // }

    Rss_Open(c, sum, size, map, ring_size, nodeNet);

    for (i = 0; i < size; i++) {
        // (c / 2^m) mod 2^(k-m-1)
        // c_prime[i] = (c[i] >> Lint(m)) ;
        // c_prime[i] = (c[i] >> Lint(m)) & nodeNet->SHIFT[ring_size - m - 1];
        // c_prime[i] = c[i]  & nodeNet->SHIFT[ring_size - m - 1];
        c_prime[i] = c[i] & nodeNet->SHIFT[m];
    }

    // for (i = 0; i < size; i++) {
    //     index = i * ring_size + ring_size - 1;
    //     b[0][i] = r_shares[0][index] ^ ((c[i] * a1) >> Lint(ring_size - 1));
    //     b[1][i] = r_shares[1][index] ^ ((c[i] * a2) >> Lint(ring_size - 1));
    // }

    for (i = 0; i < size; ++i) {

        // res[0][i] = (c_prime[i] * a1)  - rprime[0][i] + (b[0][i] << Lint(ring_size - m - 1));
        // res[1][i] = (c_prime[i] * a2)  - rprime[1][i] + (b[1][i] << Lint(ring_size - m - 1));
        // res[0][i] = (-1) * ((c_prime[i] >> m) * a1) + rprime[0][i];
        // res[1][i] = (-1) * ((c_prime[i] >> m) * a2) + rprime[1][i];

        res[0][i] = (r_m_prime[0][i] - ((c[i] >> m) * a1));
        res[1][i] = (r_m_prime[1][i] - ((c[i] >> m) * a2));
    }

    //cleanup
    delete[] c;
    delete[] c_prime;
    for (i = 0; i < 2; i++) {
        delete[] edaBit_r[i];
        delete[] edaBit_b_2[i];

        delete[] sum[i];
        delete[] b[i];
        delete[] r_m_prime[i];
        delete[] r[i];
        delete[] rprime[i];
    }
    delete[] b;
    delete[] sum;
    delete[] r_m_prime;
    delete[] r;
    delete[] rprime;
    delete[] edaBit_r;
    delete[] edaBit_b_2;
}
