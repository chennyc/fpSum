#include "../include/Rss_Op.h"

void Rss_b2u(Lint ***res, Lint **a, uint ring_size, uint size, uint m, int *map, NodeNetwork *nodeNet) {
  //size = batch_size
  int pid = nodeNet->getID();
  uint i, j;
  int q = ceil(log2(m));
  int allorResSize = pow(2, q);
  Lint one[2]; //one: 0, 0, 1
  switch(pid){
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
  //generate randbits
  Lint **r = new Lint *[2];
  Lint ***bits = new Lint **[2];
  Lint ***allorRes = new Lint**[2];
  for (i = 0; i < 2; i++) {
    r[i] = new Lint[size*q];
    bits[i] = new Lint*[q];
    allorRes[i] = new Lint*[allorResSize];
    for(j = 0; j < q; ++j) {
      bits[i][j] = new Lint[size];
    }
    for(j = 0; j < allorResSize; ++j) {
      allorRes[i][j] = new Lint[size];
    }
  }
  Rss_RandBit3(r, size*q, ring_size, map, nodeNet);
  //reform bits
  for(i = 0; i < 2; ++i) {
    for(j = 0; j < q; ++j) {
      memcpy(&bits[i][j][0], &r[i][j*size], sizeof(Lint)*size);
    }
  }
  //allor
  Rss_allor(allorRes, bits, q-1, 0, ring_size, size, map, nodeNet);
  for(i = 0; i < allorResSize; ++i) {
    for(j = 0; j < size; ++j) {
      allorRes[0][i][j] = one[0] - allorRes[0][i][j];
      allorRes[1][i][j] = one[1] - allorRes[1][i][j];
    }
  }
  //printf("allor:\n");
  //Rss_reveal(allorRes, ring_size, allorResSize, size, map, nodeNet);

  //printf("randbits:\n");
  //Rss_reveal(bits, ring_size, q, size, map, nodeNet);
  //open
  Lint **open_tmp1 = new Lint *[2];
  Lint *c = new Lint[size];
  open_tmp1[0] = new Lint [size];
  open_tmp1[1] = new Lint [size];
  memset(open_tmp1[0], 0, sizeof(Lint)*size);
  memset(open_tmp1[1], 0, sizeof(Lint)*size);

  Lint power2 = 1;
  for(i = 0; i < q; ++i) {
    for(j = 0; j < size; ++j) {
      open_tmp1[0][j] += bits[0][i][j] * power2;
      open_tmp1[1][j] += bits[1][i][j] * power2;
    }
    power2 *= 2;
  }

  //printf("sum of bits:\n");
 //Rss_reveal(open_tmp1, ring_size, size, map, nodeNet);

  //printf("a:\n");
  //Rss_reveal(a, ring_size, size, map, nodeNet);

  Lint mask = (1 << q) - 1;
  //printf("mask:%lu\n", mask);

  for(i = 0; i < size; ++i) {
    open_tmp1[0][i] += a[0][i];
    open_tmp1[1][i] += a[1][i];

    open_tmp1[0][i] &= mask;
    open_tmp1[1][i] &= mask;
  }
  Rss_Open(c, open_tmp1, size, map, ring_size, nodeNet);
  /*
  for(i = 0; i < m; ++i) {
    for(j = 0; j < size; ++j) {
      res[0][i][j] = allorRes[0][(c[j]-i) % allorResSize][j];
      res[1][i][j] = allorRes[1][(c[j]-i) % allorResSize][j];
    }
  }
  */
  for(i = 0; i < size; ++i) {
    for(j = 0; j < m; ++j) {
      res[0][i][j] = allorRes[0][(c[i]-j) % allorResSize][i];
      res[1][i][j] = allorRes[1][(c[i]-j) % allorResSize][i];
    }
  }


  // for(int i = 0; i < size; ++i){
  //    printf("c[%d] is %lu \n", i, c[i]);
  //}

  //    printf("b2u res:\n");
  //Rss_reveal(res, ring_size, m, size, map, nodeNet);


    //cleanup
  for(i = 0; i < 2; i++){
    for(j = 0; j < q; ++j) {
      delete [] bits[i][j];
    }
    for(j = 0; j < allorResSize; ++j) {
      delete [] allorRes[i][j];
    }
    delete [] open_tmp1[i];
    delete [] r[i];
    delete [] bits[i];
    delete [] allorRes[i];
  }
  delete [] c;
  delete [] open_tmp1;
  delete [] r;
  delete [] bits;
  delete [] allorRes;
}
