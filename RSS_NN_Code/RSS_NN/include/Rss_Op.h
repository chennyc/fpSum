#ifndef RSS_OP_H_
#define RSS_OP_H_
#include <stdio.h>
#include <vector>

#include "../../connection/NodeNetwork.h"
#include "util.h"
#include "sys/time.h"
#include "cmath"
using namespace std;

template <typename T>
void rss_sqrt_inv(T *c, T *e, uint size, uint ring_size) {

  T c1, c2, temp, d_;
  uint i, j;

  for (i = 0 ; i < size; i++){
    c1 = T(1);
    c2 = T(1);
    d_ = T(4); // 100 - the first mask

    for (j = 2; j < ring_size - 1; j++) {
        temp = e[i] - (c1)*(c1);
        if (temp != T(0)) {
            //get the jth+1 bit of temp, place it in jth position, and add to c1
            c1 += (temp & (d_ << T(1))) >> T(1);
        }

        temp = T(1) - c1*c2;
        // get the jth bit of temp and add it to c2
        c2 += temp & d_;
        d_ = d_ << T(1);
    }
    // last round for the inv portion
    temp = T(1) - c1*c2;
    c[i] = c2 + (temp & d_);

    }
}

template <typename T>
void Rss_Open(T *res, T **a, uint size, int *map, uint ring_size,
              NodeNetwork *nodeNet)
// a[0] will be sent to next neighbor map[0], res will be filled by the received
// value from map[1] the r_size in the functions that including open
// functionality (open_s, multpub) refers to the ring size that we are doing
// computation over, this might be different from the basic ring_size. E.g., in
// randbits, multpub is working over ring_size+2.
{
    // communication
    int i;
    // uint bytes = (RING[ring_size] + 7)>>3;
    // for(i = 0; i<size; i++){
    // a[1] = a[1] & nodeNet->SHIFT[maskring_size];
    //}
    nodeNet->SendAndGetDataFromPeer(map[0], map[1], a[1], res, size, ring_size);
    for (i = 0; i < size; i++)
    {
        res[i] = a[0][i] + a[1][i] + res[i];
        res[i] = res[i] & nodeNet->SHIFT[ring_size];
        // res[i] = bitExtracted(res[i], nodeNet->RING[ring_size]);
    }
}

template <typename T>
void Rss_Open_Bitwise(T *res, T **a, uint size, int *map, uint ring_size, NodeNetwork *nodeNet)
// a[0] will be sent to next neighbor map[0], res will be filled by the received
// value from map[1] the r_size in the functions that including open
// functionality (open_s, multpub) refers to the ring size that we are doing
// computation over, this might be different from the basic ring_size. E.g., in
// randbits, multpub is working over ring_size+2.
{
    // communication
    int i;
    // uint bytes = (RING[ring_size] + 7)>>3;
    // for(i = 0; i<size; i++){
    // a[1] = a[1] & nodeNet->SHIFT[maskring_size];
    //}
    nodeNet->SendAndGetDataFromPeer(map[0], map[1], a[1], res, size, ring_size);
    for (i = 0; i < size; i++)
    {
        res[i] = a[0][i] ^ a[1][i] ^ res[i];
        res[i] = res[i] & nodeNet->SHIFT[ring_size];
        // res[i] = bitExtracted(res[i], nodeNet->RING[ring_size]);
    }
}
template <typename T>
void Rss_Open_s(T *res, T *a_0, T *a_1, uint size, int *map, uint ring_size, NodeNetwork *nodeNet)
// a[0] will be sent to next neighbor map[0], res will be filled by the received
// value from map[1]
{
    // uint bytes =(RING[ring_size] + 7)>>3;

    // communication
    nodeNet->SendAndGetDataFromPeer(map[0], map[1], a_1, res, size, ring_size);
    for (int i = 0; i < size; i++)
    {
        res[i] = a_0[i] + a_1[i] + res[i];
        res[i] = res[i] & nodeNet->SHIFT[ring_size];
    }
}
template <typename T>
void Rss_Open_Byte(uint8_t *res, uint8_t **a, uint size, int *map, uint ring_size, NodeNetwork *nodeNet)
// a[0] will be sent to next neighbor map[0], res will be filled by the received
// value from map[1] the r_size in the functions that including open
// functionality (open_s, multpub) refers to the ring size that we are doing
// computation over, this might be different from the basic ring_size. E.g., in
// randbits, multpub is working over ring_size+2.
{
    // communication
    int i;
    // uint bytes = (size+8-1)>>3;  //number of bytes need to be send/recv

    // uint bytes = (RING[ring_size] + 7)>>3;
    // for(i = 0; i<size; i++){
    // a[1] = a[1] & nodeNet->SHIFT[maskring_size];
    //}
    nodeNet->SendAndGetDataFromPeer_bit(map[0], map[1], a[1], res, size);
    for (i = 0; i < size; i++)
    {
        res[i] = a[0][i] ^ a[1][i] ^ res[i];
        // res[i] = res[i] & nodeNet->SHIFT[ring_size];
        // res[i] = bitExtracted(res[i], nodeNet->RING[ring_size]);
    }
}

template <typename T>
void Rss_MultPub(T * c, T** a, T** b, uint size, int *map, uint ring_size, NodeNetwork* nodeNet)
//  For party 1, a[0,1]=a_2,3; b[0,1]=b_2,3;  c[0,1] = c_2,3;
//  For party 2, a[0,1]=a_3,1; b[0,1]=b_3,1;  c[0,1] = c_3,1;
//  For party 3, a[0,1]=a_1,2; b[0,1]=b_1,2;  c[0,1] = c_1,2;
{
    int i, j, k; // used for loops

    // uint bytes = (nodeNet->RING[ring_size] +7) >> 3;
    uint bytes = (ring_size +7) >> 3;

    T **sendbuf = new T *[3];
    T **recvbuf = new T *[3];
    for(i = 0; i< 3; i++){
        sendbuf[i] = new T [size];
        memset(sendbuf[i],0,sizeof(T)*size);
        recvbuf[i] = new T [size];
        memset(recvbuf[i],0,sizeof(T)*size);
    }

    int pid = nodeNet->getID();
    T *v = new T [size];
    T *v_a = new T [size];

    T opa = 0;
    T opb = 0;
    switch(pid){
        case 1:
            opa = 1;
            opb = 1;
            break;
        case 2:
            opa = -1;
            opb = 1;
            break;
        case 3:
            opa = -1;
            opb = -1;
            break;
    }


    uint8_t *buffer = new uint8_t [bytes*size];
    nodeNet->prg_getrandom(0, bytes, size, buffer);
    for (i = 0 ; i<size; i++) {
        memcpy(v_a+i, buffer + i*bytes, bytes);
    }
    nodeNet->prg_getrandom(1, bytes, size, buffer);
    for (i = 0 ; i<size; i++) {
        memcpy(c+i, buffer + i*bytes, bytes);
    }

    for(i = 0; i<size; i++){
        v[i] = a[0][i] * b[0][i] + a[0][i] * b[1][i] + a[1][i] * b[0][i];
        c[i] = v[i] + opb*c[i] + opa*v_a[i];
    }

    //communication
    //move data into buf
    for(i = 1; i<= 3; i++){
        if(i == pid) continue;
        memcpy(sendbuf[i-1],c, sizeof(T)*size);
    }

    nodeNet->multicastToPeers(sendbuf, recvbuf, size, ring_size);

    memcpy(v_a, recvbuf[map[0]-1], sizeof(T)*size);
    memcpy(v, recvbuf[map[1]-1], sizeof(T)*size);

    for(i = 0; i<size; i++){
        //mask here
        c[i] = c[i] + v_a[i] + v[i];
        c[i] = c[i] & nodeNet->SHIFT[ring_size];
    }

        //free
    delete [] v;
    delete [] v_a;
    delete [] buffer;
    for(i = 0; i< 3; i++){
        delete [] sendbuf[i];
        delete [] recvbuf[i];
    }
    delete [] sendbuf;
    delete [] recvbuf;

}

template <typename T>
void Rss_Mult_Bitwise(T** c, T** a, T** b, uint size, uint ring_size, int *map, NodeNetwork *nodeNet)
//  For party 1, a[0,1]=a_2,3; b[0,1]=b_2,3;  c[0,1] = c_2,3;
//  For party 2, a[0,1]=a_3,1; b[0,1]=b_3,1;  c[0,1] = c_3,1;
//  For party 3, a[0,1]=a_1,2; b[0,1]=b_1,2;  c[0,1] = c_1,2;
{
    // uint bytes = (nodeNet->RING[ring_size] + 7) >> 3;
    uint bytes = (ring_size + 7) >> 3;
    int i;

    T *v = new T [size];

    uint8_t *buffer = new uint8_t [bytes*size];
    nodeNet->prg_getrandom(1, bytes, size, buffer);



    for (i = 0 ; i < size; i++) {
    //nodeNet->prg_getrandom(1, bytes, c[0]+i);
        memcpy(c[0]+i, buffer + i*bytes, bytes);
        v[i] =((a[0][i] & b[0][i]) ^ (a[0][i] & b[1][i]) ^ (a[1][i] & b[0][i])) ^ c[0][i];
    }
    //communication
    nodeNet->SendAndGetDataFromPeer(map[0], map[1], v, c[1], size, ring_size);
    nodeNet->prg_getrandom(0, bytes, size, buffer);

    for(i = 0; i < size; i++){
        c[1][i] = c[1][i] ^ c[0][i];
        //nodeNet->prg_getrandom(0, bytes, c[0]+i);
        memcpy(c[0]+i, buffer + i*bytes, bytes);
        c[0][i] = c[0][i] ^ v[i];
    }



    //free
    delete [] v;
    delete [] buffer;
}
template <typename T>
void Rss_Mult(T** c, T** a, T** b, uint size, uint ring_size, int *map, NodeNetwork *nodeNet)
//  For party 1, a[0,1]=a_2,3; b[0,1]=b_2,3;  c[0,1] = c_2,3;
//  For party 2, a[0,1]=a_3,1; b[0,1]=b_3,1;  c[0,1] = c_3,1;
//  For party 3, a[0,1]=a_1,2; b[0,1]=b_1,2;  c[0,1] = c_1,2;
{
    // uint bytes = (nodeNet->RING[ring_size] + 7) >> 3;
    uint bytes = (ring_size + 7) >> 3;
    int i;

    T *v = new T [size];

    uint8_t *buffer = new uint8_t [bytes*size];
    nodeNet->prg_getrandom(1, bytes, size, buffer);
    // memcpy(c[0], buffer, size*bytes);


    for (i = 0 ; i < size; i++) {
        memcpy(c[0]+i, buffer + i*bytes, bytes);
        v[i] = a[0][i] * b[0][i] + a[0][i] * b[1][i] + a[1][i] * b[0][i] - c[0][i];
    }
    //communication
    nodeNet->SendAndGetDataFromPeer(map[0], map[1], v, c[1], size, ring_size);
    nodeNet->prg_getrandom(0, bytes, size, buffer);


    for(i = 0; i < size; i++){
        c[1][i] = c[1][i] + c[0][i];
        //nodeNet->prg_getrandom(0, bytes, c[0]+i);
        memcpy(c[0]+i, buffer + i*bytes, bytes);
        c[0][i] = c[0][i] + v[i];
    }

    //free
    delete [] v;
    delete [] buffer;
}

void Rss_Mult_Byte(uint8_t** c, uint8_t** a, uint8_t** b, uint size, int *map, NodeNetwork *nodeNet) ;

void Rss_MatMult(Lint*** c, Lint*** a, Lint*** b, uint m, uint n, uint s, uint ring_size, int *map, NodeNetwork *nodeNet);

void Rss_MatMultArray(Lint** c, Lint** a, Lint** b, uint m, uint n, uint s, uint ring_size, int *map, NodeNetwork *nodeNet);

void Rss_MatMultArray_batch(Lint** c, Lint** a, Lint** b, uint m, uint n, uint s, uint ring_size, uint batch_size, uint weight_flag_a, uint weight_flag_b, int *map, NodeNetwork *nodeNet);

void Rss_RandBit(Lint** b, uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void Rss_RandBit3(Lint** b, uint size, uint ring_size, int *map, NodeNetwork* nodeNet);

template <typename T>
void Rss_RandBitT(T** b, uint size, uint ring_size, int *map, NodeNetwork* nodeNet){

    int pid = nodeNet->getID();
    uint i;
    uint bytes = (ring_size+9) >> 3;
    // printf("bytes : %llu\n", bytes );

    T **u = new T *[2];
    T **a = new T *[2];
    T **d = new T *[2];

    for(i = 0; i < 2; i++){
        u[i] = new T [size];
        a[i] = new T [size];
        d[i] = new T [size];
    }
    T *e = new T [size];
    T *c = new T [size];
    uint8_t *buffer = new uint8_t [bytes*size];

    // used to make a odd, we only add 1 to one share of a
    // All shares will be doubled
    T a1, a2;
    switch(pid){
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


    nodeNet->prg_getrandom(0, bytes, size, buffer);
    for (i = 0 ; i<size; i++) {
        memcpy(u[0]+i, buffer + i*bytes, bytes);
    }
    nodeNet->prg_getrandom(1, bytes, size, buffer);
    for (i = 0 ; i<size; i++) {
        memcpy(u[1]+i, buffer + i*bytes, bytes);
    }

    for(i = 0 ; i<size; i++){
        // ensuring [a] is odd
        a[0][i] = (u[0][i] << T(1)) + a1;
        a[1][i] = (u[1][i] << T(1)) + a2;

    }
    // squaring a
    Rss_MultPub(e, a, a, size, map, ring_size+2, nodeNet); //ringsize+2
    rss_sqrt_inv(c, e, size, ring_size+2);

    // effectively combines the two loops into one, eliminates d variable
    for (i = 0; i < size; i++) {
        b[0][i] = (c[i]*a[0][i] + a1) >> (1);
        b[1][i] = (c[i]*a[1][i] + a2) >> (1);

    }

    // freeing up
    delete [] c;
    delete [] buffer;
    delete [] e;
    for(i = 0; i < 2; i++){
        delete [] d[i];
        delete [] a[i];
        delete [] u[i];
    }
    delete [] d;
    delete [] a;
    delete [] u;
}
void Rss_edaBit(Lint** r, Lint** b_2, uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void Rss_edaBit(Lint** r, Lint** b_2, uint size, uint ring_size, uint bit_length, int *map, NodeNetwork* nodeNet);

template <typename T>
void Rss_GenerateRandomShares(T **res, T **res_bitwise, uint ring_size, uint size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    uint i, j;
    uint bytes = (ring_size + 7) >> 3;
    // printf("bytes : %u \n", bytes);
    uint p_index = pid - 1;
    uint numParties = nodeNet->getNumParties();
    // printf("numParties : %u \n", numParties);

    T temp0, temp1;

    // used since we have effectively double the number of values
    // since we need to represent both arithmetic and binary shares
    uint new_size = 2 * size;

    // generate a single random value, which happens to already be the sum of random bits *2^j
    // [shares (0,1)][party (0,1,2)][new_size (2*size)]
    T ***r_values = new T **[2];
    for (i = 0; i < 2; i++) {
        r_values[i] = new T *[numParties];
        for (j = 0; j < numParties; j++) {
            r_values[i][j] = new T[new_size];
            memset(r_values[i][j], 0, sizeof(T) * new_size);
        }
    }

    int gamma[2];
    switch (pid) {
    case 1:
        gamma[0] = 1;
        gamma[1] = 2;
        break;
    case 2:
        gamma[0] = 2;
        gamma[1] = 0;
        break;
    case 3:
        gamma[0] = 0;
        gamma[1] = 1;
        break;
    }

    T *r_bits = new T[size];
    memset(r_bits, 0, sizeof(T) * size);

    uint8_t *buffer = new uint8_t[bytes * new_size];
    // each party generating a unique random value
    nodeNet->prg_getrandom(bytes, size, buffer);

    memcpy_fast(r_bits, buffer, size * bytes);

    // for (i = 0; i < size; i++) {
    //     memcpy_fast(r_bits + i, buffer + i * bytes, bytes);
    //     // is this what we need to do to ensure we have a shorter value?
    //     // or do we need to do something at the end of the computation
    //     // r_bits[i] = r_bits[i] & nodeNet->SHIFT[ring_size];
    // }

    nodeNet->prg_getrandom(1, bytes, new_size, buffer);

    // store arithmetic and bitwise representation sequentially
    // calculating p_i's own individual shares

    memcpy_fast(r_values[1][p_index], buffer, size * bytes);
    memcpy_fast(r_values[1][p_index] + size, buffer + size * bytes, size * bytes);

    for (i = 0; i < size; i++) {
        // memcpy_fast(r_values[1][p_index] + 2 * i, buffer + (2 * i) * bytes, bytes);
        // memcpy_fast(r_values[1][p_index] + 2 * i + 1, buffer + (2 * i + 1) * bytes, bytes);
        // r_values[0][p_index][2 * i] = r_bits[i] ;
        r_values[0][p_index][1 * i] = r_bits[i] - r_values[1][p_index][1 * i];
        // r_values[0][p_index][2 * i + 1] = r_bits[i];
        r_values[0][p_index][size + i] = r_bits[i] ^ r_values[1][p_index][size + i];
    }

    // need to generate more random shares so that binary and arithmetic representations are different
    nodeNet->prg_getrandom(0, bytes, new_size, buffer);
    memcpy_fast(r_values[0][gamma[1]], buffer, size * bytes);
    memcpy_fast(r_values[0][gamma[1]] + size, buffer + size * bytes, size * bytes);

    // for (i = 0; i < size; i++) {
    //     memcpy_fast(r_values[0][gamma[1]] + (2 * i), buffer + (2 * i) * bytes, bytes);
    //     memcpy_fast(r_values[0][gamma[1]] + (2 * i + 1), buffer + (2 * i + 1) * bytes, bytes);
    // }

    //  sending r_values[0][p_index], receiving r_values[1][gamma[0]],
    nodeNet->SendAndGetDataFromPeer(map[0], map[1], r_values[0][p_index], r_values[1][gamma[0]], new_size, ring_size);

    for (i = 0; i < numParties - 1; i++) {
        // for (i = 0; i < numParties; i++) {
        memcpy_fast(res_bitwise[0] + i * size, r_values[0][i] + (size), size * sizeof(T));
        memcpy_fast(res_bitwise[1] + i * size, r_values[1][i] + (size), size * sizeof(T));
    }

    // memcpy_fast(res_bitwise[0], r_values[0][0] + (size), size * sizeof(T));
    // memcpy_fast(res_bitwise[1], r_values[1][0] + (size), size * sizeof(T));

    // memcpy_fast(res_bitwise[0] + size, r_values[0][1] + (size), size * sizeof(T));
    // memcpy_fast(res_bitwise[1] + size, r_values[1][1] + (size), size * sizeof(T));

    // memcpy_fast(res_bitwise[0] + 2 * size, r_values[0][2] + (size), size * sizeof(T));
    // memcpy_fast(res_bitwise[1] + 2 * size, r_values[1][2] + (size), size * sizeof(T));

    for (i = 0; i < size; i++) {
        // this is so we only have two parties generating shares
        for (j = 0; j < numParties - 1; j++) {
            // for (j = 0; j < numParties; j++) {

            // adding all the parties arithmetic shares together
            // memcpy_fast(res[0] + (3 * i + j), r_values[0][j] + (2 * i), sizeof(T));
            // memcpy_fast(res[1] + (3 * i + j), r_values[1][j] + (2 * i), sizeof(T));
            res[0][i] += r_values[0][j][1 * i];
            res[1][i] += r_values[1][j][1 * i];

            // memcpy_fast(res_bitwise[0] + (numParties * i + j), r_values[0][j] + (size + i), sizeof(T));
            // memcpy_fast(res_bitwise[1] + (numParties * i + j), r_values[1][j] + (size + i), sizeof(T));
        }
    }

    for (i = 0; i < 2; i++) {
        for (j = 0; j < numParties; j++) {
            delete[] r_values[i][j];
        }
        delete[] r_values[i];
    }
    delete[] r_values;
    delete[] buffer;
    delete[] r_bits;
}

template <typename T>
void Rss_edaBitT(T **r, T **b_2, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    uint numParties = nodeNet->getNumParties();
    // printf("numParties : %u\n",numParties);
    struct timeval start;
    struct timeval end;
    unsigned long timer;

    uint i;
    // need to multiply by the number of parties in the computation
    uint new_size = numParties * size;

    T **r_bitwise = new T *[2];
    for (i = 0; i < 2; i++) {
        r_bitwise[i] = new T[new_size];
        memset(r_bitwise[i], 0, sizeof(T) * new_size);

        // ensuring destinations are sanitized
        memset(r[i], 0, sizeof(T) * size);
        memset(b_2[i], 0, sizeof(T) * size);
    }
    // gettimeofday(&start, NULL); //start timer here

    Rss_GenerateRandomShares(r, r_bitwise, ring_size, size, map, nodeNet);
    // gettimeofday(&end, NULL); //stop timer here
    // timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;

    // printf("Runtime for Rss_GenerateRandomShares with data size %d = %.6lf ms\n", size, (double)(timer * 0.001));

    // gettimeofday(&start, NULL); //start timer here

    Rss_nBitAdd(b_2, r_bitwise, ring_size, size, map, nodeNet);

    // gettimeofday(&end, NULL); //stop timer here
    // timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;

    // printf("Runtime for Rss_nBitAdd with data size %d = %.6lf ms\n", size, (double)(timer * 0.001));

    for (i = 0; i < 2; i++) {
        delete[] r_bitwise[i];
    }
    delete[] r_bitwise;
}


void Rss_edaBit_trunc(Lint** r, Lint **r_prime, Lint **r_km1, uint size, uint ring_size, uint m, int *map, NodeNetwork* nodeNet);

void Rss_edaBit_truncPre(Lint **x, Lint *c, Lint **bitLTres, Lint **r, Lint **r_prime, Lint **r_km1, uint size, uint ring_size, uint m, int *map, NodeNetwork *nodeNet);



void test_rssop();


template <typename T>
void new_Rss_MSB(T **res, T **a, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {
    int pid = nodeNet->getID();
    uint i; // used for loops

    // only need to generate a single random bit per private value
    T **b = new T *[2];
    T **sum = new T *[2];
    T **u_2 = new T *[2];
    T **edaBit_r = new T *[2];
    T **edaBit_b_2 = new T *[2];
    T **rprime = new T *[2];

    T *c = new T[size];
    T *e = new T[size];

    // used for testing correctness
    T *res_check = new T[size];
    T *r_2_open = new T[size];
    T *u_2_open = new T[size];

    for (i = 0; i < 2; i++) {
        b[i] = new T[size];
        edaBit_r[i] = new T[size];
        edaBit_b_2[i] = new T[size];
        sum[i] = new T[size];
        u_2[i] = new T[size];

        rprime[i] = new T[size];
    }

    T a1 = 0;
    T a2 = 0;
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
    // stays the same
    Rss_RandBit3(b, size, ring_size, map, nodeNet);

    // Rss_edaBit(edaBit_r, edaBit_b_2, size, ring_size, ring_size - 1, map, nodeNet);
    // need to generate full edabit for final implementation
    Rss_edaBit(edaBit_r, edaBit_b_2, size, ring_size, map, nodeNet);

    for (i = 0; i < size; i++) {
        rprime[0][i] = edaBit_r[0][i] - (GET_BIT(edaBit_b_2[0][i], T(ring_size - 1)) << T(ring_size - 1));
        rprime[1][i] = edaBit_r[1][i] - (GET_BIT(edaBit_b_2[1][i], T(ring_size - 1)) << T(ring_size - 1));
        // combining w the next loop
        // combining w the previous loop
        sum[0][i] = (a[0][i] + edaBit_r[0][i]);
        sum[1][i] = (a[1][i] + edaBit_r[1][i]);
        // sum[0][i] = (a[0][i] + edaBit_r[0][i]) << T(1);
        // sum[1][i] = (a[1][i] + edaBit_r[1][i]) << T(1);
    }

    Rss_Open(c, sum, size, map, ring_size, nodeNet);

    for (i = 0; i < size; i++) {

        // used for alternate solution
        // prevents us from using it directly later on
        edaBit_b_2[0][i] = edaBit_b_2[0][i] & nodeNet->SHIFT[ring_size - 1];
        edaBit_b_2[1][i] = edaBit_b_2[1][i] & nodeNet->SHIFT[ring_size - 1];

        // definitely needed
        c[i] = c[i] & nodeNet->SHIFT[ring_size - 1];
        // c[i] = c[i] >> T(1);
    }

    // Rss_Open_Bitwise(r_2_open, edaBit_b_2, size, map, ring_size, nodeNet);
    // // this part is still correct
    // however, the edaBit_b_2 shares do get modified
    // which may not be desierable
    Rss_BitLT(u_2, c, edaBit_b_2, ring_size, size, map, nodeNet);

    // Rss_Open_Bitwise(u_2_open, u_2, size, map, ring_size, nodeNet);

    // for (int i = 0; i < size; i++) {

    //     res_check[i] = (c[i] < r_2_open[i]);
    //     if (!(u_2_open[i] == res_check[i])) {
    //         // printf("[%i] c < r_2 : %u   --- expected: %u\n", i, u_2_open[i], res_check[i]);
    //         // printf("c = %u --- edaBit_b_2 = %u\n", c[i], r_2_open[i]);
    //         // printf("BitLT ERROR at %d \n", i);
    //     }
    // }

    for (i = 0; i < size; ++i) {
        // cant do this because we modify edaBit_b_2 earlier

        sum[0][i] = a[0][i] - c[i] * a1 + rprime[0][i] - (u_2[0][i] << T(ring_size - 1)) + (b[0][i] << T(ring_size - 1));
        sum[1][i] = a[1][i] - c[i] * a2 + rprime[1][i] - (u_2[1][i] << T(ring_size - 1)) + (b[1][i] << T(ring_size - 1));

        // sum[0][i] = a[0][i] - c[i] * a1 + edaBit_r[0][i] - (u_2[0][i] << T(ring_size - 1)) + (b[0][i] << T(ring_size - 1));
        // sum[1][i] = a[1][i] - c[i] * a2 + edaBit_r[1][i] - (u_2[1][i] << T(ring_size - 1)) + (b[1][i] << T(ring_size - 1));
    }
    // opening sum
    Rss_Open(e, sum, size, map, ring_size, nodeNet);

    T e_bit;
    for (i = 0; i < size; ++i) {
        e_bit = GET_BIT(e[i], ring_size - 1); // getting the (k-1)th bit
        res[0][i] = e_bit * a1 + b[0][i] - (e_bit * b[0][i] << T(1));
        res[1][i] = e_bit * a2 + b[1][i] - (e_bit * b[1][i] << T(1));
    }

    // cleanup
    delete[] c;
    delete[] e;
    delete[] res_check;
    delete[] r_2_open;
    delete[] u_2_open;

    for (i = 0; i < 2; i++) {
        delete[] edaBit_r[i];
        delete[] edaBit_b_2[i];
        delete[] b[i];
        delete[] sum[i];
        delete[] u_2[i];
        delete[] rprime[i];
    }
    delete[] edaBit_r;
    delete[] edaBit_b_2;
    delete[] b;
    delete[] sum;
    delete[] u_2;
    delete[] rprime;
}

template <typename T>
void test_new_Rss_MSB(T **res, T **a, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {
    int pid = nodeNet->getID();
    uint i; // used for loops

    // only need to generate a single random bit per private value
    T **b = new T *[2];

    T **sum = new T *[2];
    T **sum2 = new T *[2];

    T **u_2 = new T *[2];
    T **u_2_2 = new T *[2];

    T **edaBit_r = new T *[2];
    T **edaBit_b_2 = new T *[2];
    T **edaBit_b_2_2 = new T *[2];

    T **rprime = new T *[2];
    T **r = new T *[2];

    T *c = new T[size];
    T *c2 = new T[size];
    T *e = new T[size];
    T *e2 = new T[size];

    // used for testing correctness
    T *res_check = new T[size];
    T *res_check_2 = new T[size];
    T *r_2_open = new T[size];
    T *u_2_open = new T[size];
    T *r_2_open_2 = new T[size];
    T *u_2_open_2 = new T[size];

    for (i = 0; i < 2; i++) {
        b[i] = new T[2*size];
        edaBit_r[i] = new T[size];

        edaBit_b_2[i] = new T[size];
        edaBit_b_2_2[i] = new T[size];

        sum[i] = new T[size];
        sum2[i] = new T[size];

        u_2[i] = new T[size];
        u_2_2[i] = new T[size];

        rprime[i] = new T[size];
        r[i] = new T[size];
    }

    T a1 = 0;
    T a2 = 0;
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
    // stays the same
    Rss_RandBit(b, 2*size, ring_size, map, nodeNet);

    Rss_edaBit(edaBit_r, edaBit_b_2, size, ring_size, ring_size - 1, map, nodeNet);
    // need to generate full edabit for final implementation
    // Rss_edaBit(edaBit_r, edaBit_b_2, size, ring_size, map, nodeNet);

    for (i = 0; i < size; i++) {
        r[0][i] = edaBit_r[0][i] + (GET_BIT(b[0][i], T(ring_size - 1)) << T(ring_size - 1));
        r[1][i] = edaBit_r[1][i] + (GET_BIT(b[1][i], T(ring_size - 1)) << T(ring_size - 1));
        // combining w the next loop
        // }

        // sum = 2(a _ )
        // for (i = 0; i < size; i++) {
        // combining w the previous loop
        // sum[0][i] = (a[0][i] + edaBit_r[0][i]);
        // sum[1][i] = (a[1][i] + edaBit_r[1][i]);
        sum[0][i] = (a[0][i] + edaBit_r[0][i]) ;
        sum[1][i] = (a[1][i] + edaBit_r[1][i]) ;

        sum2[0][i] = (a[0][i] + r[0][i]) ;
        sum2[1][i] = (a[1][i] + r[1][i]) ;
    }

    Rss_Open(c, sum, size, map, ring_size, nodeNet);
    Rss_Open(c2, sum2, size, map, ring_size, nodeNet);

    // c = c/2
    for (i = 0; i < size; i++) {

        // used for alternate solution
        // prevents us from using it directly later on
        // edaBit_b_2_2[0][i] = edaBit_b_2[0][i] & nodeNet->SHIFT[ring_size - 1];
        // edaBit_b_2_2[1][i] = edaBit_b_2[1][i] & nodeNet->SHIFT[ring_size - 1];

        edaBit_b_2_2[0][i] = edaBit_b_2[0][i] ;
        edaBit_b_2_2[1][i] = edaBit_b_2[1][i] ;
        // definitely needed
        c[i] = c[i] & nodeNet->SHIFT[ring_size - 1];
        c2[i] = c2[i] & nodeNet->SHIFT[ring_size - 1];

        // c2[i] = c2[i] ;
        // c[i] = c[i] ;
    }

    // printf("hi2\n");
    Rss_Open_Bitwise(r_2_open, edaBit_b_2, size, map, ring_size, nodeNet);
    Rss_Open_Bitwise(r_2_open_2, edaBit_b_2_2, size, map, ring_size, nodeNet);
    // // this part is still correct
    // however, the edaBit_b_2 shares do get modified
    // which may not be desierable
    Rss_BitLT(u_2, c, edaBit_b_2, ring_size, size, map, nodeNet);
    Rss_BitLT(u_2_2, c2, edaBit_b_2_2, ring_size, size, map, nodeNet);

    Rss_Open_Bitwise(u_2_open, u_2, size, map, ring_size, nodeNet);
    Rss_Open_Bitwise(u_2_open_2, u_2_2, size, map, ring_size, nodeNet);

    for (int i = 0; i < size; i++) {

        res_check[i] = (c[i] < r_2_open[i]);
        res_check_2[i] = (c2[i] < r_2_open_2[i]);
        if (!(u_2_open[i] == res_check[i]) && !(u_2_open_2[i] == res_check_2[i])) {
            printf("[%i] c < r_2 : %u   --- expected: %u\n", i, u_2_open[i], res_check[i]);
            printf("[%i] c2 < r_2_2 : %u   --- expected: %u\n", i, u_2_open_2[i], res_check_2[i]);
            printf("c = %u --- edaBit_b_2 = %u\n", c[i], r_2_open[i]);
            printf("c2 = %u --- edaBit_b_2_2 = %u\n", c2[i], r_2_open_2[i]);
            printf("BitLT ERROR at %u \n", i);
        }
    }

    // printf("hi3\n");
    for (i = 0; i < size; ++i) {
        // cant do this because we modify edaBit_b_2 earlier
        // sum[0][i] = a[0][i] - c[i] * a1 + (edaBit_r[0][i] - (GET_BIT(edaBit_b_2[0][i], T(ring_size - 1)) << T(ring_size - 1))) - (u_2[0][i] << T(ring_size - 1)) + (b[0][i] << T(ring_size - 1));
        // sum[1][i] = a[1][i] - c[i] * a2 + (edaBit_r[1][i] - (GET_BIT(edaBit_b_2[1][i], T(ring_size - 1)) << T(ring_size - 1))) - (u_2[1][i] << T(ring_size - 1)) + (b[1][i] << T(ring_size - 1));

        sum2[0][i] = a[0][i] - c2[i] * a1 + edaBit_r[0][i] - (u_2_2[0][i] << T(ring_size - 1)) + (b[0][size+i] << T(ring_size - 1));
        sum2[1][i] = a[1][i] - c2[i] * a2 + edaBit_r[1][i] - (u_2_2[1][i] << T(ring_size - 1)) + (b[1][size+i] << T(ring_size - 1));

        sum[0][i] = a[0][i] - c[i] * a1 + edaBit_r[0][i] - (u_2[0][i] << T(ring_size - 1)) + (b[0][size+i] << T(ring_size - 1));
        sum[1][i] = a[1][i] - c[i] * a2 + edaBit_r[1][i] - (u_2[1][i] << T(ring_size - 1)) + (b[1][size+i] << T(ring_size - 1));
    }
    // opening sum
    Rss_Open(e, sum, size, map, ring_size, nodeNet);
    Rss_Open(e2, sum2, size, map, ring_size, nodeNet);
    for (int i = 0; i < size; i++) {
        if (!(GET_BIT(e[i], ring_size - 1) == GET_BIT(e2[i], ring_size - 1))) {
            printf("[%i] e = %u --- e2 = %u\n", i, GET_BIT(e[i], ring_size - 1), GET_BIT(e2[i], ring_size - 1));
        }
    }

    T e_bit;
    T e_bit2;
    for (i = 0; i < size; ++i) {
        e_bit = GET_BIT(e[i], ring_size - 1);   // getting the (k-1)th bit
        res[0][i] = e_bit * a1 + b[0][size+i] - (e_bit * b[0][size+i] << T(1));
        res[1][i] = e_bit * a2 + b[1][size+i] - (e_bit * b[1][size+i] << T(1));

        // res[0][i] = e_bit2 * a1 + b[0][i] - (e_bit2 * b[0][i] << T(1));
        // res[1][i] = e_bit2 * a2 + b[1][i] - (e_bit2 * b[1][i] << T(1));
    }

    Rss_Open(res_check, res, size, map, ring_size, nodeNet);

    for (i = 0; i < size; ++i) {
        e_bit2 = GET_BIT(e2[i], ring_size - 1); // getting the (k-1)th bit
        // res[0][i] = e_bit * a1 + b[0][i] - (e_bit * b[0][i] << T(1));
        // res[1][i] = e_bit * a2 + b[1][i] - (e_bit * b[1][i] << T(1));

        res[0][i] = e_bit2 * a1 + b[0][size+i] - (e_bit2 * b[0][size+i] << T(1));
        res[1][i] = e_bit2 * a2 + b[1][size+i] - (e_bit2 * b[1][size+i] << T(1));
    }
    Rss_Open(res_check_2, res, size, map, ring_size, nodeNet);

    for (i = 0; i < size; i++)
    {
        if (!(res_check_2[i] == res_check[i]))
        {
            printf("actual:  : %llu\t", res_check[i]);
            printf("expected  : %llu\n", res_check_2[i]);
        }


    }


    // cleanup
    delete[] c2;
    delete[] e2;

    delete[] c;
    delete[] e;

    delete[] res_check;
    delete[] res_check_2;

    delete[] r_2_open;
    delete[] r_2_open_2;
    delete[] u_2_open;
    delete[] u_2_open_2;

    for (i = 0; i < 2; i++) {
        delete[] edaBit_r[i];
        delete[] edaBit_b_2[i];
        delete[] edaBit_b_2_2[i];
        delete[] b[i];
        delete[] sum[i];
        delete[] sum2[i];
        delete[] u_2[i];
        delete[] u_2_2[i];
        delete[] rprime[i];
        delete[] r[i];
    }
    delete[] edaBit_r;
    delete[] edaBit_b_2;
    delete[] edaBit_b_2_2;
    delete[] b;

    delete[] sum;
    delete[] sum2;

    delete[] u_2;
    delete[] u_2_2;

    delete[] rprime;
    delete[] r;
}

template <typename T>
void Rss_MSB(T **res, T **a, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    uint i, j, k, index; // used for loops
    // struct timeval start;
    // struct timeval end;
    // unsigned long timer;

    uint n_rand_bits = size * (ring_size + 1);

    T **r_shares = new T *[2];
    for (i = 0; i < 2; i++) {
        r_shares[i] = new T[n_rand_bits];
    }

    T **b = new T *[2];
    T **r = new T *[2];
    T **sum = new T *[2];
    T **rprime = new T *[2];
    T **r_2 = new T *[2];
    T **u_2 = new T *[2];

    T *c = new T[size];
    memset(c, 0, sizeof(T) * size);
    T *e = new T[size];
    memset(e, 0, sizeof(T) * size);

    // T *res_check = new T[size];
    // T *r_2_open = new T[size];
    // T *u_2_open = new T[size];

    for (i = 0; i < 2; i++) {
        b[i] = new T[size];
        r[i] = new T[size];
        sum[i] = new T[size];
        rprime[i] = new T[size];
        memset(rprime[i], 0, sizeof(T) * size);
        r_2[i] = new T[size];
        memset(r_2[i], 0, sizeof(T) * size);
        u_2[i] = new T[size];
        memset(u_2[i], 0, sizeof(T) * size);
    }

    T a1 = 0;
    T a2 = 0;
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

    // offline component start
    Rss_RandBit(r_shares, n_rand_bits, ring_size, map, nodeNet);
    for (i = 0; i < size; i++) {
        b[0][i] = r_shares[0][size * ring_size + i];
        b[1][i] = r_shares[1][size * ring_size + i];
    }

    for (j = 0; j < size; j++) {
        for (k = 0; k < ring_size - 1; k++) {
            // this is for step 3
            index = j * ring_size + k;
            rprime[0][j] = rprime[0][j] + (r_shares[0][index] << T(k));
            rprime[1][j] = rprime[1][j] + (r_shares[1][index] << T(k));
        }
        index = j * ring_size + k;
        r[0][j] = (rprime[0][j] + ((r_shares[0][index]) << T(k)));
        r[1][j] = (rprime[1][j] + ((r_shares[1][index]) << T(k)));
    }

    // for ( i = 0; i < size; i++)
    // {
    //     printf("rprime[0][%i]: %llu\n",i, rprime[0][i]);
    //     print_binary(rprime[0][i], 8*sizeof(T));
    //     printf("r[0][%i]: %llu\n",i, r[0][i]);
    //     print_binary(r[0][i], 8*sizeof(T));

    //     printf("rprime[1][%i]: %llu\n",i, rprime[1][i]);
    //     print_binary(rprime[1][i], 8*sizeof(T));
    //     printf("r[1][%i]: %llu\n",i, r[1][i]);
    //     print_binary(r[1][i], 8*sizeof(T));
    // }



    for (j = 0; j < size; j++) {
        for (k = 0; k < ring_size - 1; k++) {
            index = j * ring_size + k;

            r_2[0][j] = T(SET_BIT(r_2[0][j], T(k), GET_BIT(r_shares[0][index], T(0))));
            r_2[1][j] = T(SET_BIT(r_2[1][j], T(k), GET_BIT(r_shares[1][index], T(0))));
        }
    }
    // offline component ends
    // step 2
    for (i = 0; i < size; i++) {
        sum[0][i] = (a[0][i] + r[0][i]); // & nodeNet->SHIFT[1] ;
        sum[1][i] = (a[1][i] + r[1][i]); // & nodeNet->SHIFT[1] ;
    }

    // for (size_t i = 0; i < size; i++)
    // {
    //     printf("r[0][%i] : %llu\n", i, r[0][i]);
    //     print_binary(r[0][i], 32);
    //     printf("r[1][%i] : %llu\n", i, r[1][i]);
    //     print_binary(r[0][i], 32);
    // }

    Rss_Open(c, sum, size, map, ring_size, nodeNet);

    // step 3 -- getting the k-1 lowest bits of c[i]
    for (i = 0; i < size; i++) {
        c[i] = c[i] & nodeNet->SHIFT[ring_size - 1];
    }

    // Rss_Open_Bitwise(r_2_open, r_2, size, map, ring_size, nodeNet);
    // calculating c <? r_2, where r_2 is bitwise shared
    // gettimeofday(&start,NULL); //start timer here
    Rss_BitLT(u_2, c, r_2, ring_size, size, map, nodeNet);
    // gettimeofday(&end,NULL);//stop timer here
    // timer = 1000000 * (end.tv_sec-start.tv_sec)+ end.tv_usec-start.tv_usec;
    // printf("runtime for BitLT with data size %d = %ld us\n", size, timer);

    // Rss_Open_Bitwise(u_2_open, u_2, size, map, ring_size, nodeNet);

    // for (int i = 0; i < size; i++) {

    //     res_check[i] = (c[i] < r_2_open[i]);
    //     if (!(u_2_open[i] == res_check[i])) {
    //         printf("[%i] c < r_2 : %u   --- expected: %u\n", i, u_2_open[i], res_check[i]);
    //         printf("c = %u --- r_2_open = %u\n", c[i], r_2_open[i]);
    //         printf("BitLT ERROR at %d \n", i);
    //     }
    // }

    for (i = 0; i < size; ++i) {
        sum[0][i] = a[0][i] - c[i] * a1 + rprime[0][i] - (u_2[0][i] << T(ring_size - 1)) + (b[0][i] << T(ring_size - 1));
        sum[1][i] = a[1][i] - c[i] * a2 + rprime[1][i] - (u_2[1][i] << T(ring_size - 1)) + (b[1][i] << T(ring_size - 1));
    }
    // opening sum
    Rss_Open(e, sum, size, map, ring_size, nodeNet);

    T e_bit;
    for (i = 0; i < size; ++i) {
        e_bit = GET_BIT(e[i], ring_size - 1); // getting the (k-1)th bit
        res[0][i] = e_bit * a1 + b[0][i] - (e_bit * b[0][i] << T(1));
        res[1][i] = e_bit * a2 + b[1][i] - (e_bit * b[1][i] << T(1));
    }

    // cleanup
    delete[] c;
    delete[] e;
    // delete[] res_check;
    // delete[] r_2_open;
    // delete[] u_2_open;

    for (i = 0; i < 2; i++) {
        delete[] r_shares[i];
        delete[] b[i];
        delete[] r[i];
        delete[] sum[i];
        delete[] rprime[i];
        delete[] r_2[i];
        delete[] u_2[i];
    }
    delete[] r_shares;
    delete[] b;
    delete[] r;
    delete[] sum;
    delete[] rprime;
    delete[] r_2;
    delete[] u_2;
}

void Rss_GenerateRandomShares(Lint** res, Lint** r_i_values, uint ring_size, uint size, int *map, NodeNetwork* nodeNet);
void Rss_GenerateRandomShares(Lint **res, Lint **r_i_values, uint ring_size, uint bit_length, uint size, int *map, NodeNetwork *nodeNet);
void Rss_GenerateRandomShares_trunc(Lint **res, Lint **res_prime, Lint **r_i_values, uint ring_size, uint m, uint size, int *map, NodeNetwork *nodeNet);
void Rss_GenerateRandomShares_trunc_2(Lint **res, Lint **res_prime, Lint **r_i_values, uint ring_size, uint m, uint size, int *map, NodeNetwork *nodeNet);

void Rss_Convert(Lint **res, Lint** a, uint size, uint ring_size, uint ring_size_prime,int *map, NodeNetwork* nodeNet);

template <typename T>
void new_Rss_Convert(T **res, T **a, uint size, uint ring_size, uint ring_size_prime, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    uint i, j, k, index; // used for loops

    T **edaBit_r = new T *[2];
    T **edaBit_b_2 = new T *[2];
    T **r_2 = new T *[2];
    T **sum = new T *[2];
    T **a_2 = new T *[2];
    T **a_2_buffer = new T *[2];
    T **a_k_prime = new T *[2];
    T *c = new T[size];

    for (i = 0; i < 2; i++) {
        r_2[i] = new T[size];
        a_2[i] = new T[size];
        a_2_buffer[i] = new T[size*ring_size];
        sum[i] = new T[size];
        a_k_prime[i] = new T[size*ring_size];

        edaBit_r[i] = new T[size];
        edaBit_b_2[i] = new T[size];
    }

    // only need ring_size bit-length values for both parts of computation
    // first protection and b2a
    Rss_edaBitT(edaBit_r, edaBit_b_2, size, ring_size, map, nodeNet);

    for (i = 0; i < size; i++) {
        sum[0][i] = (a[0][i] - edaBit_r[0][i]);
        sum[1][i] = (a[1][i] - edaBit_r[1][i]);
    }

    // Rss_Open(c, sum, size, map, ring_size_prime, nodeNet);
    Rss_Open(c, sum, size, map, ring_size, nodeNet);

    Rss_BitAdd(a_2, c, edaBit_b_2, ring_size, size, map, nodeNet);
    // Rss_BitAdd(res, c, edaBit_b_2, ring_size, size, map, nodeNet);
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < size; ++j) {
            for(k = 0; k < ring_size; ++k) {
                a_2_buffer[i][j*ring_size + k] =GET_BIT(a_2[i][j], T(k));
            }
        }
    }
    Rss_b2a3(a_k_prime, a_2_buffer, ring_size_prime, size*ring_size, map, nodeNet);
    for (i = 0; i < 2; ++i) {
        for (j = 0; j < size; j++) {
            res[0][j] = 0, res[1][j] = 0;
            for (k = 0; k < ring_size; k++) {
                res[i][j] += (a_k_prime[i][j*ring_size + k] << T(k));
            }
        }
    }

    for (i = 0; i < 2; i++) {
        delete[] edaBit_r[i];
        delete[] edaBit_b_2[i];
        delete[] sum[i];
        delete[] r_2[i];
        delete[] a_2[i];
        delete[] a_k_prime[i];
    }
    delete[] edaBit_r;
    delete[] edaBit_b_2;
    delete[] sum;
    delete[] r_2;
    delete[] a_2;
    delete[] a_k_prime;
}

void Rss_nBitAdd_trunc(Lint **res, Lint **carry, Lint **r_bitwise, uint ring_size, uint m, uint size, int *map, NodeNetwork *nodeNet);
void Rss_BitAdd_trunc(Lint **res, Lint **carry, Lint **a, Lint **b, uint ring_size, uint m, uint size, int *map, NodeNetwork *nodeNet);
// overloading BitAdd
template <typename T>
void Rss_nBitAdd(T **res, T **r_bitwise, uint ring_size, uint size, int *map, NodeNetwork *nodeNet) {

    uint i, j;
    uint numParties = nodeNet->getNumParties();
    uint rounds = ceil(log2(numParties));

    T **a = new T *[2];

    T **b = new T *[2];
    for (i = 0; i < 2; i++) {
        a[i] = new T[size];
        b[i] = new T[size];
    }

    // always will be 2 rounds for 3pc
    for (j = 0; j < rounds; j++) {

        // if this is the first iteration, we copy r_bitwise into a and b
        if (j == 0) {
            // copy p_1 and p_2 into a, b respectively

            memcpy_fast(a[0], r_bitwise[0], size * sizeof(T));
            memcpy_fast(a[1], r_bitwise[1], size * sizeof(T));

            memcpy_fast(b[0], r_bitwise[0] + size, size * sizeof(T));
            memcpy_fast(b[1], r_bitwise[1] + size, size * sizeof(T));

            // for (i = 0; i < size; i++) {
            //     memcpy_fast(a[0] + i, r_bitwise[0] + 3 * i, sizeof(T));
            //     memcpy_fast(a[1] + i, r_bitwise[1] + 3 * i, sizeof(T));

            //     memcpy_fast(b[0] + i, r_bitwise[0] + 3 * i + 1, sizeof(T));
            //     memcpy_fast(b[1] + i, r_bitwise[1] + 3 * i + 1, sizeof(T));
            // }
            Rss_BitAdd(res, a, b, ring_size, size, map, nodeNet);
        } else {
            // we only need to copy r_bitwise into b

            // memcpy_fast(a[0], res[0], size * sizeof(T));
            // memcpy_fast(a[1], res[1], size * sizeof(T));

            // memcpy_fast(b[0], r_bitwise[0] + 2*size, size * sizeof(T));
            // memcpy_fast(b[1], r_bitwise[1] + 2*size, size * sizeof(T));




            // for (i = 0; i < size; i++) {
            //     memcpy_fast(a[0] + i, res[0] + i, sizeof(T));
            //     memcpy_fast(a[1] + i, res[1] + i, sizeof(T));

            //     memcpy_fast(b[0] + i, r_bitwise[0] + 3 * i + 2, sizeof(T));
            //     memcpy_fast(b[1] + i, r_bitwise[1] + 3 * i + 2, sizeof(T));
            // }


            // Rss_BitAdd(res, a, b, ring_size, size, map, nodeNet);
        }
    }

    for (i = 0; i < 2; i++) {
        delete[] a[i];
        delete[] b[i];
    }
    delete[] a;
    delete[] b;
}


template <typename T>
void CarryBuffer2(T **buffer, T **d, uint **index_array, uint size, uint k) {
    // prepares input u for multiplication
    // extracts p2i, p2i-1, and g2i
    // buffer and d are the same size (4 x size)

    // struct timeval start;
    // struct timeval end;
    // unsigned long timer;

    T i, j, mask1, mask2, mask1p, mask2p;

    // gettimeofday(&start, NULL); //start timer here

    for (i = 0; i < size; i++) {
        for (j = 0; j < k; j++) {

            // used to set the bits in the correct positions in buffer
            mask1 = 2 * j;
            mask2 = 2 * j + 1;

            // used to get the correct bits from d
            mask1p = index_array[0][j];
            mask2p = index_array[1][j];

            buffer[0][i] = SET_BIT(buffer[0][i], mask1, GET_BIT(d[0][i], mask1p));
            buffer[1][i] = SET_BIT(buffer[1][i], mask1, GET_BIT(d[1][i], mask1p));
            buffer[0][i] = SET_BIT(buffer[0][i], mask2, GET_BIT(d[2][i], mask1p));
            buffer[1][i] = SET_BIT(buffer[1][i], mask2, GET_BIT(d[3][i], mask1p));

            buffer[2][i] = SET_BIT(buffer[2][i], mask1, GET_BIT(d[0][i], mask2p));
            buffer[3][i] = SET_BIT(buffer[3][i], mask1, GET_BIT(d[1][i], mask2p));
            buffer[2][i] = SET_BIT(buffer[2][i], mask2, GET_BIT(d[0][i], mask2p));
            buffer[3][i] = SET_BIT(buffer[3][i], mask2, GET_BIT(d[1][i], mask2p));
        }
    }

    // gettimeofday(&end, NULL); //stop timer here
    // timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    // printf("Runtime for CarryBuffer2 with data size %d = %.6lf ms\n", size, (double)(timer * 0.001));
}

// alternative BitAdd implementation when both a and b are secret shared
// used in edaBit
template <typename T>
void Rss_BitAdd(T **res, T **a, T **b, uint ring_size, uint size, int *map, NodeNetwork *nodeNet) {

    T i, j;
    int pid = nodeNet->getID();

    T **d = new T *[4];
    for (i = 0; i < 4; i++) {
        d[i] = new T[size];
        memset(d[i], 0, sizeof(T) * size);
    }
    // struct timeval start;
    // struct timeval end;
    // unsigned long timer;

    // gettimeofday(&start, NULL); //start timer here
    Rss_Mult_Bitwise(res, a, b, size, ring_size, map, nodeNet);
    // gettimeofday(&end, NULL); //stop timer here
    // timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    // printf("Runtime for Rss_Mult_Bitwise with data size %d = %.6lf ms\n", size, (double)(timer * 0.001));

    for (i = 0; i < size; i++) {
        d[0][i] = a[0][i] ^ b[0][i];
        d[1][i] = a[1][i] ^ b[1][i];

        d[2][i] = res[0][i];
        d[3][i] = res[1][i];
    }

    // Rss_CircleOpL(d, ring_size, size, map, nodeNet);

    // gettimeofday(&start, NULL); //start timer here

    // Rss_CircleOpL_T(d, ring_size, size, map, nodeNet); // new version w Ts
    Rss_CircleOpL(d, ring_size, size, map, nodeNet); // original

    // gettimeofday(&end, NULL); //stop timer here
    // timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    // printf("Runtime for Rss_CircleOpL with data size %d = %.6lf ms\n", size, (double)(timer * 0.001));

    for (i = 0; i < size; i++) {

        res[0][i] = (a[0][i] ^ b[0][i]) ^ (d[2][i] << 1);
        res[1][i] = (a[1][i] ^ b[1][i]) ^ (d[3][i] << 1);
    }

    for (i = 0; i < 4; i++) {
        delete[] d[i];
    }
    delete[] d;
}

// a is public, b is bitwise-shared
// res will be a bitwise shared output
template <typename T>
void Rss_BitAdd(T **res, T *a, T **b, uint ring_size, uint size, int *map, NodeNetwork *nodeNet) {

    T i, j;
    int pid = nodeNet->getID();

    T **d = new T *[4];
    for (i = 0; i < 4; i++) {
        d[i] = new T[size];
        memset(d[i], 0, sizeof(T) * size);
    }

    struct timeval start;
    struct timeval end;
    unsigned long timer;

    T a1, a2;
    switch (pid) {
    case 1:
        a1 = -1;
        a2 = 0;
        break;
    case 2:
        a1 = 0;
        a2 = 0;
        break;
    case 3:
        a1 = 0;
        a2 = -1;
        break;
    }
    for (i = 0; i < size; i++) {
        d[0][i] = (a[i] & a1) ^ b[0][i];
        d[1][i] = (a[i] & a2) ^ b[1][i];

        d[2][i] = (a[i] & b[0][i]);
        d[3][i] = (a[i] & b[1][i]);
    }

    // gettimeofday(&start, NULL); //start timer here
    // Rss_CircleOpL_T(d, ring_size, size, map, nodeNet);
    Rss_CircleOpL(d, ring_size, size, map, nodeNet);
    // gettimeofday(&end, NULL); //stop timer here
    // timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    // printf("Runtime for Rss_CircleOpL with data size %d = %.6lf ms\n", size, (double)(timer * 0.001));

    // we only need the values in the g position (indices 2 and 3)

    for (i = 0; i < size; i++) {
        res[0][i] = ((a[i] & a1) ^ b[0][i]) ^ (d[2][i] << 1);
        res[1][i] = ((a[i] & a2) ^ b[1][i]) ^ (d[3][i] << 1);
    }
    for (i = 0; i < 4; i++) {
        delete[] d[i];
    }

    delete[] d;

}

template <typename T>
void Rss_CircleOpL(T **d, uint r_size, uint size, int *map, NodeNetwork *nodeNet) {

    T i, j, l, k, y, z, op_r; // used for loops

    // struct timeval start;
    // struct timeval end;
    // unsigned long timer;

    // struct timeval start2;
    // struct timeval end2;
    // unsigned long timer2;

    if (r_size > 1) {

        // just three nested for-loops
        // r_size <=> k in algorithm

        uint num = ((r_size + 7) >> 3) * size;
        uint n_uints = ((r_size + 7) >> 3);
        uint t_index;
        T mask2, mask1m8, mask2m8;
        T mask1p, mask2p;
        uint r_size_2 = pow(2, ceil(log2(r_size))); // rounding up to next power of two
        uint rounds = ceil(log2(r_size_2));

        uint **index_array = new uint *[2];
        T **buffer = new T *[4];
        // we need at most (r_size + 7)/8 bytes to store bits from the buffer
        uint8_t **a = new uint8_t *[2];
        uint8_t **b = new uint8_t *[2];
        uint8_t **u = new uint8_t *[2];

        for (i = 0; i < 2; i++) {
            index_array[i] = new uint[r_size_2];
            memset(index_array[i], 0, sizeof(uint) * r_size_2);
            buffer[2 * i] = new T[size];
            buffer[2 * i + 1] = new T[size];
            memset(buffer[2 * i], 0, sizeof(T) * size);
            memset(buffer[2 * i + 1], 0, sizeof(T) * size);
            // memsets are actually needed here since are ORing
            a[i] = new uint8_t[num];
            memset(a[i], 0, sizeof(uint8_t) * num);
            b[i] = new uint8_t[num];
            memset(b[i], 0, sizeof(uint8_t) * num);
            u[i] = new uint8_t[num];
            memset(u[i], 0, sizeof(uint8_t) * num);
        }

        for (i = 1; i <= rounds; i++) {
            // gettimeofday(&start, NULL); //start timer here

            op_r = 0; // number of operations in a round
            // equivalent to the new_ring_size in MSB

            // gettimeofday(&start2, NULL); //start timer here

            for (j = 1; j <= ceil(r_size_2 / pow(2, i)); j++) {

                y = uint(pow(2, i - 1) + j * pow(2, i)) % r_size_2;

                for (z = 1; z <= (pow(2, i - 1)); z++) {
                    // checking we have a valid set of indices to add to our set
                    if ((((y % r_size_2)) <= r_size) && (((y + z) % (r_size_2 + 1)) <= r_size)) {
                        // printf("y : %u\n", y);
                        // printf("y+z : %u\n", y+z);
                        index_array[0][op_r] = (y % r_size_2) - 1;
                        index_array[1][op_r] = ((y + z) % (r_size_2 + 1)) - 1;
                        op_r++;
                    }
                }
            }

            // gettimeofday(&end2, NULL); //stop timer here
            // timer2 = 1000000 * (end2.tv_sec - start2.tv_sec) + end2.tv_usec - start2.tv_usec;
            // printf("Runtime for index_array with data size %d = %.6lf ms\n", size, (double)(timer2 * 0.001));

            // updating parameters for optimization
            n_uints = ((2 * op_r + 7) >> 3);
            num = ((2 * op_r + 7) >> 3) * size;

            // extracting terms into buffer
            CarryBuffer2(buffer, d, index_array, size, op_r);

            // gettimeofday(&start2, NULL); //start timer here
            // Splitting the buffer into bytes

            // THIS DOESNT WORK
            // DO NOT TRY
            // memcpy_fast(a[0], buffer[0], size*n_uints);
            // memcpy_fast(a[1], buffer[1], size*n_uints);
            // memcpy_fast(b[0], buffer[2], size*n_uints);
            // memcpy_fast(b[1], buffer[3], size*n_uints);


            for (j = 0; j < size; ++j) {
                memcpy_fast(a[0] + j * n_uints, buffer[0] + j, n_uints);
                memcpy_fast(a[1] + j * n_uints, buffer[1] + j, n_uints);
                memcpy_fast(b[0] + j * n_uints, buffer[2] + j, n_uints);
                memcpy_fast(b[1] + j * n_uints, buffer[3] + j, n_uints);
            }

            // gettimeofday(&end2, NULL); //stop timer here
            // timer2 = 1000000 * (end2.tv_sec - start2.tv_sec) + end2.tv_usec - start2.tv_usec;
            // printf("Runtime for memcpy_fast with data size %d = %.6lf ms\n", size, (double)(timer2 * 0.001));

            // gettimeofday(&start2, NULL); //start timer here

            // bitwise multiplication
            Rss_Mult_Byte(u, a, b, num, map, nodeNet);
            // gettimeofday(&end2, NULL); //stop timer here
            // timer2 = 1000000 * (end2.tv_sec - start2.tv_sec) + end2.tv_usec - start2.tv_usec;
            // printf("Runtime for Rss_Mult_Byte with data size %d = %.6lf ms\n", size, (double)(timer2 * 0.001));

            // gettimeofday(&start2, NULL); //start timer here

            for (l = 0; l < size; ++l) {
                for (j = 0; j < op_r; ++j) {
                    // loop constants
                    t_index = (j >> 2) + (l * n_uints);
                    mask2 = index_array[1][j];
                    mask1m8 = (2 * j) & 7; // "&7" = %8, used for leftover bits
                    mask2m8 = (2 * j + 1) & 7;

                    d[0][l] = SET_BIT(d[0][l], mask2, GET_BIT(u[0][t_index], mask1m8));
                    d[1][l] = SET_BIT(d[1][l], mask2, GET_BIT(u[1][t_index], mask1m8));

                    // simplified from needing two separate loops
                    d[2][l] = SET_BIT(d[2][l], mask2, (GET_BIT(u[0][t_index], mask2m8) ^ GET_BIT(d[2][l], mask2)));
                    d[3][l] = SET_BIT(d[3][l], mask2, (GET_BIT(u[1][t_index], mask2m8) ^ GET_BIT(d[3][l], mask2)));
                }
            }

            // gettimeofday(&end2, NULL); //stop timer here
            // timer2 = 1000000 * (end2.tv_sec - start2.tv_sec) + end2.tv_usec - start2.tv_usec;
            // printf("Runtime for rearranging/circleOp with data size %d = %.6lf ms\n", size, (double)(timer2 * 0.001));

            // gettimeofday(&end, NULL); //stop timer here
            // timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
            // printf("Runtime for round %i with data size %d = %.6lf ms\n\n", i, size, (double)(timer * 0.001));
        }

        for (i = 0; i < 2; i++) {
            delete[] buffer[2 * i];
            delete[] buffer[2 * i + 1];
            delete[] a[i];
            delete[] b[i];
            delete[] u[i];
            delete[] index_array[i];
        }

        delete[] a;
        delete[] b;
        delete[] u;
        delete[] index_array;
        delete[] buffer;
    }
}

template <typename T>
void Rss_CircleOpL_T(T **d, uint r_size, uint size, int *map, NodeNetwork *nodeNet) {

    T i, j, l, y, z, op_r; // used for loops
    struct timeval start;
    struct timeval end;
    unsigned long timer;

    if (r_size > 1) {

        T mask2, mask1m8, mask2m8;
        uint r_size_2 = pow(2, ceil(log2(r_size))); // rounding up to next power of two
        uint rounds = ceil(log2(r_size_2));

        uint **index_array = new uint *[2];
        T **a_prime = new T *[2];
        T **b_prime = new T *[2];
        T **u_prime = new T *[2];

        for (i = 0; i < 2; i++) {
            index_array[i] = new uint[r_size_2];
            a_prime[i] = new T[size];
            b_prime[i] = new T[size];
            u_prime[i] = new T[size];
        }

        for (i = 1; i <= rounds; i++) {
            gettimeofday(&start, NULL); //start timer here

            op_r = 0; // number of operations in a round
            // equivalent to the new_ring_size in MSB

            for (j = 1; j <= ceil(r_size_2 / pow(2, i)); j++) {

                y = uint(pow(2, i - 1) + j * pow(2, i)) % r_size_2;

                for (z = 1; z <= (pow(2, i - 1)); z++) {
                    // checking we have a valid set of indices to add to our set
                    if ((((y % r_size_2)) <= r_size) && (((y + z) % (r_size_2 + 1)) <= r_size)) {
                        // printf("y : %u\n", y);
                        // printf("y+z : %u\n", y+z);
                        index_array[0][op_r] = (y % r_size_2) - 1;
                        index_array[1][op_r] = ((y + z) % (r_size_2 + 1)) - 1;
                        op_r++;
                    }
                }
            }

            // extracting terms into a_prime and b_prime directly
            CarryBuffer_T(a_prime, b_prime, d, index_array, size, op_r);

            // bitwise multiplication
            Rss_Mult_Bitwise(u_prime, a_prime, b_prime, size, r_size, map, nodeNet);

            // printf("adding g2j+1\n");
            for (l = 0; l < size; ++l) {
                for (j = 0; j < op_r; ++j) {

                    // loop constants
                    // DO NOT REMOVE mask2
                    // putting it directly in next  operations
                    // of SET BIT causes it to FAIL
                    // DONT ASK WHY

                    mask2 = index_array[1][j]; // CHECK
                    mask1m8 = (2 * j);         // "&7" = %8, used for leftover bits
                    mask2m8 = (2 * j + 1);

                    d[0][l] = SET_BIT(d[0][l], mask2, GET_BIT(u_prime[0][l], mask1m8));
                    d[1][l] = SET_BIT(d[1][l], mask2, GET_BIT(u_prime[1][l], mask1m8));

                    d[2][l] = SET_BIT(d[2][l], mask2, (GET_BIT(u_prime[0][l], mask2m8) ^ GET_BIT(d[2][l], mask2)));
                    d[3][l] = SET_BIT(d[3][l], mask2, (GET_BIT(u_prime[1][l], mask2m8) ^ GET_BIT(d[3][l], mask2)));
                }
            }
            gettimeofday(&end, NULL); //stop timer here
            timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
            // printf("Runtime for round with data size %d = %.6lf ms\n", size, (double)(timer * 0.001));
        }
        for (i = 0; i < 2; i++) {
            delete[] a_prime[i];
            delete[] b_prime[i];
            delete[] u_prime[i];
            delete[] index_array[i];
        }
        delete[] a_prime;
        delete[] b_prime;
        delete[] u_prime;
        delete[] index_array;
    }
}



template <typename T>
void CarryBuffer_T(T **a_prime, T **b_prime, T **d, uint **index_array, uint size, uint k) {

    T i, j, mask1, mask2, mask1p, mask2p;

    for (i = 0; i < size; i++) {
        for (j = 0; j < k; j++) {

            // used to set the bits in the correct positions in buffer
            mask1 = 2 * j;
            mask2 = 2 * j + 1;

            // used to get the correct bits from d
            mask1p = index_array[0][j];
            mask2p = index_array[1][j];

            a_prime[0][i] = SET_BIT(a_prime[0][i], mask1, GET_BIT(d[0][i], mask1p));
            a_prime[1][i] = SET_BIT(a_prime[1][i], mask1, GET_BIT(d[1][i], mask1p));
            a_prime[0][i] = SET_BIT(a_prime[0][i], mask2, GET_BIT(d[2][i], mask1p));
            a_prime[1][i] = SET_BIT(a_prime[1][i], mask2, GET_BIT(d[3][i], mask1p));

            b_prime[0][i] = SET_BIT(b_prime[0][i], mask1, GET_BIT(d[0][i], mask2p));
            b_prime[1][i] = SET_BIT(b_prime[1][i], mask1, GET_BIT(d[1][i], mask2p));
            b_prime[0][i] = SET_BIT(b_prime[0][i], mask2, GET_BIT(d[0][i], mask2p));
            b_prime[1][i] = SET_BIT(b_prime[1][i], mask2, GET_BIT(d[1][i], mask2p));
        }
    }
}

void OptimalBuffer_Lint(Lint **a_prime, Lint **b_prime, Lint **d, uint **index_array, uint size, uint k);

void Rss_b2kprime(Lint **res, Lint **a, uint ring_size, uint ring_size_prime, uint size, int *map, NodeNetwork *nodeNet);
void Rss_b2a(Lint **res, Lint **a, uint ring_size, uint size, int *map, NodeNetwork *nodeNet);

void Rss_BitLT(Lint **res, Lint* a, Lint** b, uint ring_size, uint size, int *map, NodeNetwork* nodeNet);
void Rss_LT(Lint **res, Lint** a, Lint** b, uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void new_Rss_LT(Lint **res, Lint** a, Lint** b, uint size, uint ring_size, int *map, NodeNetwork* nodeNet);

void Rss_CarryOut(Lint** res, Lint* a, Lint** b, uint ring_size, uint size, int *map, NodeNetwork* nodeNet);

void Rss_CarryOutAux(Lint** res,  Lint** d, uint r_size, uint size, int *map, NodeNetwork* nodeNet);
void new_Rss_CarryOutAux(Lint** res,  Lint** d, uint r_size, uint size, int *map, NodeNetwork* nodeNet);
void new_Rss_CarryOutAux_Lint(Lint** res,  Lint** d, uint r_size, uint size, int *map, NodeNetwork* nodeNet);

void CarryBuffer(Lint **buffer, Lint **u, uint size, uint r_size);
void OptimalBuffer(Lint **buffer, Lint **u, uint size, uint r_size, NodeNetwork* nodeNet);

void CircleOp(Lint **res, Lint *p, Lint *g, uint size, uint r_size);


void truncPr_1(Lint **res, Lint **x, uint m,  uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void Rss_truncPr(Lint **res, Lint **x, uint m,  uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void Rss_truncPre(Lint **res, Lint **x, uint m,  uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void truncPr_2(Lint **res, Lint **x, uint m,  uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void truncPr_2_test(Lint **res, Lint **x, uint m,  uint size, uint ring_size, int *map, NodeNetwork* nodeNet);



void rss_sqrt_inv(Lint *c, Lint *e, uint size, uint ring_size);

void invert(Lint *c, Lint *a, int size, int ring_size);

void rss_sqrt(Lint *c, Lint *e, int size, int ring_size);

void Rss_PrefixMult(Lint ***res, Lint ***input, uint ring_size, uint size, uint batch_size, int *map, NodeNetwork *nodeNet, NodeConfiguration* nodeConfig);
void Rss_PrefixAnd_batch(Lint ***res, Lint ***input, uint ring_size, uint size, uint batch_size, int *map, NodeNetwork *nodeNet, NodeConfiguration *nodeConfig);
void Rss_PrefixAnd(Lint **res, Lint **input, uint ring_size, uint size, int *map, NodeNetwork *nodeNet);
void Rss_PrefixOr(Lint **res, Lint **input, uint ring_size, uint size, int *map, NodeNetwork *nodeNet);


void Rss_JustLT(Lint **res, Lint** a, Lint** b, uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void Rss_int2mask(Lint ***res, uint k, Lint **input, uint ring_size, uint batch_size, int *map, NodeNetwork* nodeNet, NodeConfiguration* nodeConfig);
void Rss_fpSum(Lint *res_b, Lint *res_p, Lint **res_v, Lint **b, Lint **p, Lint ***v, int m, int e, int w, uint batch_size, uint ring_size, int *map, NodeNetwork* nodeNet, NodeConfiguration* nodeConfig, unsigned long* timer);
void Rss_allor(Lint ***res, Lint ***d, uint start, uint end, uint ring_size, uint batch_size, int *map, NodeNetwork *nodeNet);
void Rss_b2u(Lint ***res, Lint **a, uint ring_size, uint size, uint m, int *map, NodeNetwork *nodeNet);
void Rss_reveal(Lint **toBeReveal, uint ring_size, uint size, int *map, NodeNetwork *nodeNet);
void Rss_reveal(Lint ***toBeReveal, uint ring_size, uint size1, uint size2, int *map, NodeNetwork *nodeNet);
template <typename T>
void Rss_bitDec(T ***res, T **a, uint size, uint ring_size, uint bitLen, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    uint i, j, k, index; // used for loops
    uint n_rand_bits = size * ring_size;
    // first size*ring_size random bits used for first part of share conversion
    // 2nd half used for 2tok'
    T a1 = 0, a2 = 0;
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

    T **r_shares = new T *[2];
    for (i = 0; i < 2; i++) {
        r_shares[i] = new T[n_rand_bits];
    }
    // used for B2A component
    Rss_RandBit(r_shares, n_rand_bits, ring_size, map, nodeNet);

    T **edaBit_r = new T *[2];
    T **edaBit_b_2 = new T *[2];
    T **r_2 = new T *[2];
    T **sum = new T *[2];
    T **a_2 = new T *[2];

    // will hold the k shares of a in (k_prime)
    T **a_k_prime = new T *[2];
    T *c = new T[size];

    T *res_check = new T[size];

    for (i = 0; i < 2; i++) {
        r_2[i] = new T[size];
        a_2[i] = new T[size];
        sum[i] = new T[size];
        a_k_prime[i] = new T[n_rand_bits];

        edaBit_r[i] = new T[size];
        edaBit_b_2[i] = new T[size];
    }

    // only need ring_size bit-length values for both parts of computation
    // first protection and b2a
    Rss_edaBit(edaBit_r, edaBit_b_2, size, bitLen, map, nodeNet);

    for (i = 0; i < size; i++) {
        sum[0][i] = (a[0][i] - edaBit_r[0][i]);
        sum[1][i] = (a[1][i] - edaBit_r[1][i]);
    }

    // Rss_Open(c, sum, size, map, ring_size_prime, nodeNet);
    Rss_Open(c, sum, size, map, bitLen, nodeNet);

    Rss_BitAdd(a_2, c, edaBit_b_2, bitLen, size, map, nodeNet);
    // Rss_BitAdd(res, c, edaBit_b_2, ring_size, size, map, nodeNet);

    // 2tok' component (B2A)
    // a_2 is k-bits long

    for (j = 0; j < size; j++) {
        for (k = 0; k < bitLen; k++) {
            index = j * bitLen + k;

            r_2[0][j] = T(SET_BIT(r_2[0][j], T(k), GET_BIT(r_shares[0][index], T(0))));
            r_2[1][j] = T(SET_BIT(r_2[1][j], T(k), GET_BIT(r_shares[1][index], T(0))));
        }
    }

    // 2tok' component
    // a_2 is k-bits long

    for (i = 0; i < size; i++) {
        sum[0][i] = (a_2[0][i] ^ r_2[0][i]);
        sum[1][i] = (a_2[1][i] ^ r_2[1][i]);
    }

    // should this be bitwise-open? (like Open_Byte)
    // and what should the ring_size be?
    // memset(c, 0, sizeof(T) * size); // clearing c
    Rss_Open_Bitwise(c, sum, size, map, ring_size, nodeNet);

    for (i = 0; i < size; i++) {
        for (j = 0; j < bitLen; j++) {
            index = i * bitLen + j;

            a_k_prime[0][index] = a1 * GET_BIT(c[i], T(j)) + r_shares[0][index] - 2 * GET_BIT(c[i], T(j)) * r_shares[0][index];
            a_k_prime[1][index] = a2 * GET_BIT(c[i], T(j)) + r_shares[1][index] - 2 * GET_BIT(c[i], T(j)) * r_shares[1][index];
        }
    }
    //printf("a_k: /n");
    //Rss_reveal(a_k_prime, ring_size, n_rand_bits, map, nodeNet);
    // constructing the final shares of a

    //T resMask = (1 << ring_size) - 1;
    //printf("mask: %lx, resMask");
    for (j = 0; j < size; j++) {
        for (k = 0; k < bitLen; k++) {
            // this is for step 3
            index = j * bitLen + k;
            res[0][j][k] = a_k_prime[0][index];
            res[1][j][k] = a_k_prime[1][index];
        }
    }
    //printf("bitRes: /n");
    //Rss_reveal(res, ring_size, size, bitLen, map, nodeNet);

    delete[] c;
    delete[] res_check;

    for (i = 0; i < 2; i++) {
        delete[] r_shares[i];
        delete[] edaBit_r[i];
        delete[] edaBit_b_2[i];
        delete[] sum[i];
        delete[] r_2[i];
        delete[] a_2[i];
        delete[] a_k_prime[i];
    }
    delete[] r_shares;
    delete[] edaBit_r;
    delete[] edaBit_b_2;
    delete[] sum;
    delete[] r_2;
    delete[] a_2;
    delete[] a_k_prime;
}

template <typename T>
void Rss_bitDec3(T ***res, T **a, uint size, uint ring_size, uint bitLen, int *map, NodeNetwork *nodeNet) {

    int pid = nodeNet->getID();
    uint i, j, k, index; // used for loops
    //    uint n_rand_bits = size * ring_size;
    // first size*ring_size random bits used for first part of share conversion
    // 2nd half used for 2tok'
    T a1 = 0, a2 = 0;
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

    T **edaBit_r = new T *[2];
    T **edaBit_b_2 = new T *[2];
    T **sum = new T *[2];
    T **a_2 = new T *[2];

    // will hold the k shares of a in (k_prime)
    T *c = new T[size];

    for (i = 0; i < 2; i++) {
        a_2[i] = new T[size];
        sum[i] = new T[size];

        edaBit_r[i] = new T[size];
        edaBit_b_2[i] = new T[size];
    }

    // only need ring_size bit-length values for both parts of computation
    // first protection and b2a
    Rss_edaBit(edaBit_r, edaBit_b_2, size, bitLen, map, nodeNet);

    for (i = 0; i < size; i++) {
        sum[0][i] = (a[0][i] - edaBit_r[0][i]);
        sum[1][i] = (a[1][i] - edaBit_r[1][i]);
    }

    // Rss_Open(c, sum, size, map, ring_size_prime, nodeNet);
    Rss_Open(c, sum, size, map, bitLen, nodeNet);

    Rss_BitAdd(a_2, c, edaBit_b_2, bitLen, size, map, nodeNet);
    // Rss_BitAdd(res, c, edaBit_b_2, ring_size, size, map, nodeNet);


    T **b2ares = new T *[2];
    T **b2ainput = new T *[2];

    for (i = 0; i < 2; i++) {
      b2ares[i] = new T[size * bitLen]();
      b2ainput[i] = new T[size * bitLen]();
    }

    for(i = 0; i < 2; ++i) {
      for(j = 0; j < size; ++j) {
	for(k = 0; k < bitLen; ++k) {
	  //	  b2ainput[i][j*bitLen + k] = a_2[i][j] & (1 << k);
	  b2ainput[i][j*bitLen + k] = GET_BIT(a_2[i][j], (T)k);
	}
      }
    }

    //    Rss_reveal(b2ainput, ring_size, size * bitLen, map, nodeNet);
    Rss_b2a3(b2ares, b2ainput, ring_size, size * bitLen, map, nodeNet);
    //Rss_reveal(b2ares, ring_size, size * bitLen, map, nodeNet);

    for(i = 0; i < 2; ++i) {
      for(j = 0; j < size; ++j) {
	for(k = 0; k < bitLen; ++k) {
	  res[i][j][k] = b2ares[i][j*bitLen + k];
	}
      }
    }

    delete[] c;

    for (i = 0; i < 2; i++) {
      delete[] b2ares[i];
      delete[] b2ainput[i];
        delete[] edaBit_r[i];
        delete[] edaBit_b_2[i];
        delete[] sum[i];
        delete[] a_2[i];
    }
    delete[] b2ares;
    delete[] b2ainput;

    delete[] edaBit_r;
    delete[] edaBit_b_2;
    delete[] sum;
    delete[] a_2;
}
void Rss_CarryMult(Lint **res, Lint ***input, uint ring_size, uint size, uint batch_size, int *map, NodeNetwork *nodeNet);
void Rss_shift(Lint **res, Lint **a, Lint **offset, uint bitLenOffset, uint size, uint ring_size, int *map, NodeNetwork *nodeNet);
void Rss_newshift(Lint ***res, Lint ***a, Lint **offset, uint bitLenOffset, uint beta, uint size, uint ring_size, int *map, NodeNetwork *nodeNet);

template <typename T>
void Rss_b2a3(T **res, T **a, uint ring_size, uint size, int *map, NodeNetwork *nodeNet) {
    int pid = nodeNet->getID();
    uint i;
    T **b = new T *[2];
    T **t = new T *[2];
    uint bytes = (ring_size + 7) >> 3;
    T *c = new T[size]();

    uint8_t *buffer = new uint8_t [bytes*size];

    for (i = 0; i < 2; i++) {
      b[i] = new T[size]();
      t[i] = new T[size]();
    }

    switch (pid) {
    case 1:
      nodeNet->prg_getrandom(0, bytes, size, buffer);
      for(i = 0; i < size; ++i){
	      memcpy(&b[0][i], buffer + i * bytes, bytes);
      }
      break;
    case 2:
      nodeNet->getDataFromPeer(3, size, b[1], ring_size);
      break;
    case 3:
      nodeNet->prg_getrandom(1, bytes, size, buffer);
      for(i = 0; i < size; ++i){
	      memcpy(&b[1][i], buffer + i * bytes, bytes);
	      b[0][i] = a[0][i] * a[1][i] - b[1][i];
      }
      nodeNet->sendDataToPeer(2, size, b[0], ring_size);
      break;
    }
    //        printf("before lin 9\n");
    //        Rss_reveal(b, ring_size, size, map, nodeNet);
    //line 9
    switch (pid) {
    case 1:
      for(i = 0; i < size; ++i) {
	b[0][i] = a[0][i] - 2 * b[0][i];
	b[1][i] = 0 - 2 * b[1][i];
      }
      break;
    case 2:
      for(i = 0; i < size; ++i) {
	b[0][i] = 0 - 2 * b[0][i];
	b[1][i] = a[1][i] - 2 * b[1][i];
      }
      break;
    case 3:
      for(i = 0; i < size; ++i) {
	b[0][i] = a[0][i] - 2 * b[0][i];
	b[1][i] = a[1][i] - 2 * b[1][i];
      }
      break;
    }

    //line 10
    switch (pid) {
    case 1:
      nodeNet->prg_getrandom(1, bytes, size, buffer);
      for(i = 0; i < size; ++i){
	memcpy(&c[i], buffer + i * bytes, bytes);
	t[0][i] = b[0][i] * a[1][i] - c[i];
      }
      nodeNet->sendDataToPeer(3, size, t[0], ring_size);
      nodeNet->getDataFromPeer(2, size, t[1], ring_size);
      for(i = 0; i < size; ++i){
	t[1][i] += c[i];
      }
      break;
    case 2:
      nodeNet->prg_getrandom(1, bytes, size, buffer);
      for(i = 0; i < size; ++i){
	memcpy(&t[1][i], buffer + i * bytes, bytes);
	t[0][i] = b[1][i] * a[0][i] - t[1][i];
      }
      nodeNet->sendDataToPeer(1, size, t[0], ring_size);

      nodeNet->prg_getrandom(0, bytes, size, buffer);
      for(i = 0; i < size; ++i){
	memcpy(&c[i], buffer + i * bytes, bytes);
	t[0][i] += c[i];
      }
      break;
    case 3:
      nodeNet->prg_getrandom(0, bytes, size, buffer);
      for(i = 0; i < size; ++i){
	memcpy(&t[0][i], buffer + i * bytes, bytes);
      }
      nodeNet->getDataFromPeer(1, size, t[1], ring_size);
      break;
    }

    //    Rss_reveal(b, ring_size, size, map, nodeNet);
    //    Rss_reveal(t, ring_size, size, map, nodeNet);
    //line 14
    switch (pid) {
    case 1:
      for(i = 0; i < size; ++i) {
	res[0][i] = b[0][i] - 2 * t[0][i];
	res[1][i] = b[1][i] - 2 * t[1][i] + a[1][i];
      }
      break;
    case 2:
      for(i = 0; i < size; ++i) {
	res[0][i] = b[0][i] - 2 * t[0][i] + a[0][i];
	res[1][i] = b[1][i] - 2 * t[1][i];
      }
      break;
    case 3:
      for(i = 0; i < size; ++i) {
	res[0][i] = b[0][i] - 2 * t[0][i];
	res[1][i] = b[1][i] - 2 * t[1][i];
      }
      break;
    }

    delete[] c;

    for (i = 0; i < 2; i++) {
        delete[] b[i];
        delete[] t[i];
    }
    delete[] b;
    delete[] t;
    delete[] buffer;
}
void Rss_kor(Lint **res, Lint **input, uint ring_size, uint size, uint productRangeSize, int *map, NodeNetwork *nodeNet);
void Rss_reveal(bool flag, Lint **toBeReveal, uint ring_size, uint size, int *map, NodeNetwork *nodeNet);
void Rss_eqz3(Lint **res, Lint **a, uint ring_size, uint size, int *map, NodeNetwork *nodeNet);
void Rss_MassDotProduct(Lint*** res, Lint***superA, Lint*** bitArray, int batch_size, int alpha, int beta, uint ring_size, int *map, NodeNetwork *nodeNet);
void Rss_SASum(Lint **res, Lint ***superA, int w, uint batch_size, int alpha, uint ring_size, int *map, NodeNetwork *nodeNet);
void Rss_superSum_batch(Lint ***res, Lint ***input, int w, uint batch_size, int alpha, uint ring_size, int *map, NodeNetwork *nodeNet);
void Rss_superSum(Lint **res, Lint **input, int w, int alpha, uint ring_size, int *map, NodeNetwork *nodeNet);

void Rss_sa2fl(Lint *b, Lint *p, Lint **v, Lint **sa, int m, int e, int w, uint ring_size, int *map, NodeNetwork* nodeNet, int pid);
void Rss_ORead(Lint** res, Lint**superA, Lint** bitArray,  int alpha, int beta, uint ring_size, int *map, NodeNetwork *nodeNet);
void Rss_normal(Lint *b, Lint **v, Lint *p, Lint **vp, int w, int beta, int m, uint ring_size, int pid, int *map, NodeNetwork *nodeNet);
void Rss_ORead2(Lint** res, Lint**bits, Lint** bitArray,  int k, int kp, uint ring_size, int *map, NodeNetwork *nodeNet);



Lint bitExtracted(Lint number, uint k);
uint8_t smallBitExtracted(uint8_t a, uint k);
uint FeedBufferBits(Lint *buffer, Lint source, Lint bitLen, Lint ring_size, uint buffer_point);
Lint TakeBufferBits(Lint *buffer, Lint bitLen, Lint ring_size, uint buffer_point);
Lint bitExtractedRange(Lint number, uint length, uint start);

#endif
