#include "../include/Rss_Op.h"
// extern "C" {
// #include "../aes_ni.h"
// }
// converts shares of a_k to shares of a_k'
// k --> ring_size
// k' --> ring_size_prime
// rewrite from convert
void Rss_shift(Lint **res, Lint **a, Lint **offset, uint bitLenOffset, uint size, uint ring_size, int *map, NodeNetwork *nodeNet)
{
    uint i, j, k, index;
    int pid = nodeNet->getID();
    Lint one[2];
    switch (pid)
    {
    case 1:
        one[0] = 0;
        one[1] = 0;
        break;
    case 2:
        one[0] = 0;
        one[1] = 1;
        break;
    case 3:
        one[0] = 1;
        one[1] = 0;
        break;
    }

    Lint ***offSetBits = new Lint **[2];
    Lint ***offSets = new Lint **[2];
    Lint **finalOffset = new Lint *[2];
    for (i = 0; i < 2; ++i)
    {
        offSetBits[i] = new Lint *[size];
        offSets[i] = new Lint *[bitLenOffset];
        finalOffset[i] = new Lint[size];
        memset(finalOffset[i], 0, sizeof(Lint) * size);
        for (j = 0; j < size; ++j)
        {
            offSetBits[i][j] = new Lint[bitLenOffset];
            memset(&offSetBits[i][j][0], 0, sizeof(Lint) * bitLenOffset);
            // offSets[i][j] = new Lint[bitLenOffset];
        }
        for (j = 0; j < bitLenOffset; ++j)
        {
            offSets[i][j] = new Lint[size];
        }
    }

    // printf("offset: before bitDec\n");
    // Rss_reveal(offset, ring_size, size,  map, nodeNet);

    Rss_bitDec3(offSetBits, offset, size, ring_size, bitLenOffset, map, nodeNet);

    // printf("offsetBits: after bitDec\n");
    // Rss_reveal(offSetBits, ring_size, size, bitLenOffset, map, nodeNet);
    // offset_i = 1 - bit_i + bit_i * 2^(i+1);
    for (i = 0; i < size; ++i)
    {
        Lint pow2 = 2;
        for (j = 0; j < bitLenOffset; ++j)
        {
            offSets[0][j][i] = one[0] - offSetBits[0][i][j] + offSetBits[0][i][j] * pow2;
            offSets[1][j][i] = one[1] - offSetBits[1][i][j] + offSetBits[1][i][j] * pow2;
            pow2 = pow2 * pow2;
        }
    }

    // carryMult
    Rss_CarryMult(finalOffset, offSets, ring_size, bitLenOffset, size, map, nodeNet);

    // printf("final offsets: after bitDec\n");
    // Rss_reveal(finalOffset, ring_size, size, map, nodeNet);

    Rss_Mult(res, finalOffset, a, size, ring_size, map, nodeNet);

    // printf("final shifted\n");
    // Rss_reveal(res, ring_size, size, map, nodeNet);

    // cleanup
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < size; ++j)
        {
            delete[] offSetBits[i][j];
        }
        for (j = 0; j < bitLenOffset; ++j)
        {
            delete[] offSets[i][j];
        }
        delete[] offSetBits[i];
        delete[] offSets[i];
    }
    delete[] offSetBits;
    delete[] offSets;
}

void Rss_newshift(Lint ***res, Lint ***a, Lint **offset, uint bitLenOffset, uint beta, uint size, uint ring_size, int *map, NodeNetwork *nodeNet)
{
    uint i, j, k, index;
    int pid = nodeNet->getID();
    Lint one[2];
    int w = pow(2, bitLenOffset);
    switch (pid)
    {
    case 1:
        one[0] = 0;
        one[1] = 0;
        break;
    case 2:
        one[0] = 0;
        one[1] = 1;
        break;
    case 3:
        one[0] = 1;
        one[1] = 0;
        break;
    }

    Lint ***offSetBits = new Lint **[2];
    Lint ***offSets = new Lint **[2];
    Lint **finalOffset = new Lint *[2];
    for (i = 0; i < 2; ++i)
    {
        offSetBits[i] = new Lint *[size];
        offSets[i] = new Lint *[bitLenOffset];
        finalOffset[i] = new Lint[size * (beta - 1)];
        memset(finalOffset[i], 0, sizeof(Lint) * size * (beta - 1));
        for (j = 0; j < size; ++j)
        {
            offSetBits[i][j] = new Lint[bitLenOffset];
            memset(&offSetBits[i][j][0], 0, sizeof(Lint) * bitLenOffset);
            // offSets[i][j] = new Lint[bitLenOffset];
        }
        for (j = 0; j < bitLenOffset; ++j)
        {
            offSets[i][j] = new Lint[size * (beta - 1)];
        }
    }

    Rss_bitDec3(offSetBits, offset, size, ring_size, bitLenOffset, map, nodeNet);

    for (i = 0; i < size; ++i)
    {
        for (k = 0; k < beta - 1; ++k)
        {
            for (j = 0; j < bitLenOffset; ++j)
            {
                Lint powofpow = pow(2, pow(2, j));
                offSets[0][j][i * (beta - 1) + k] = one[0] - offSetBits[0][i][j] + offSetBits[0][i][j] * powofpow;
                offSets[1][j][i * (beta - 1) + k] = one[1] - offSetBits[1][i][j] + offSetBits[1][i][j] * powofpow;
            }
        }
    }

    Rss_CarryMult(finalOffset, offSets, ring_size, bitLenOffset, size * (beta - 1), map, nodeNet);

    Lint **Mult_tmp1 = new Lint *[2];
    Mult_tmp1[0] = new Lint[size * (beta - 1)];
    Mult_tmp1[1] = new Lint[size * (beta - 1)];

    Lint **Mult_tmp2 = new Lint *[2];
    Mult_tmp2[0] = new Lint[size * (beta - 1)];
    Mult_tmp2[1] = new Lint[size * (beta - 1)];

    Lint **trunc_res = new Lint *[2];
    trunc_res[0] = new Lint[size * (beta - 1)];
    trunc_res[1] = new Lint[size * (beta - 1)];

    for (i = 0; i < size; ++i)
    {
        memcpy(&Mult_tmp1[0][i * (beta - 1)], a[0][i], sizeof(Lint) * (beta - 1));
        memcpy(&Mult_tmp1[1][i * (beta - 1)], a[1][i], sizeof(Lint) * (beta - 1));
    }

    Rss_Mult(Mult_tmp2, finalOffset, Mult_tmp1, size * (beta - 1), ring_size, map, nodeNet);
    Rss_truncPre(trunc_res, Mult_tmp2, bitLenOffset, size * (beta - 1), ring_size, map, nodeNet);

    for (i = 0; i < size; ++i)
    {
        for (j = 1; j < beta - 1; ++j)
        {
            res[0][i][j] = Mult_tmp2[0][i * (beta - 1) + j] - trunc_res[0][i * (beta - 1) + j] * pow(2, w) + trunc_res[0][i * (beta - 1) + j - 1];
            res[0][i][j] = Mult_tmp2[0][i * (beta - 1) + j] - trunc_res[0][i * (beta - 1) + j] * pow(2, w) + trunc_res[0][i * (beta - 1) + j - 1];
        }
    }

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < size; ++j)
        {
            delete[] offSetBits[i][j];
        }
        for (j = 0; j < bitLenOffset; ++j)
        {
            delete[] offSets[i][j];
        }
        delete[] offSetBits[i];
        delete[] offSets[i];
        delete[] Mult_tmp1[i];
        delete[] Mult_tmp2[i];
    }
    delete[] offSetBits;
    delete[] offSets;
    delete[] Mult_tmp1;
    delete[] Mult_tmp2;
}