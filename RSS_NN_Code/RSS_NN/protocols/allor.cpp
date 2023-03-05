#include "../include/Rss_Op.h"

void Rss_allor(Lint ***res, Lint ***d, uint start, uint end, uint ring_size, uint batch_size, int *map, NodeNetwork *nodeNet) {
  //start = k-2, end = 0
  uint i, j;
  int pid = nodeNet->getID();
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
  if(start == end) {
    for(i = 0; i < batch_size; ++i) {
      res[0][0][i] = d[0][start][i];
      res[1][0][i] = d[1][start][i];
      res[0][1][i] = one[0] - d[0][start][i];
      res[1][1][i] = one[1] - d[1][start][i];
    }
    return;
  }
  int k = start - end + 1;
  int ill = floor(k/2);
  int size1 = pow(2, ill);
  int size2 = pow(2, k-ill);
  int totalSize = pow(2, k);
  Lint ***res1 = new Lint**[2];
  Lint ***res2 = new Lint**[2];
  for(i = 0; i < 2; ++i) {
    res1[i] = new Lint*[size1];
    res2[i] = new Lint*[size2];
    for(j = 0; j < size1; ++j) {
      res1[i][j] = new Lint[batch_size];
    }
    for(j = 0; j < size2; ++j) {
      res2[i][j] = new Lint[batch_size];
    }
  }
  //recusive
  Rss_allor(res1, d, end+ill-1, end, ring_size, batch_size, map, nodeNet);
  Rss_allor(res2, d, start, end+ill, ring_size, batch_size, map, nodeNet);  
  //printf("res1: %d to %d\n", end+ill-1, end);
  //Rss_reveal(res1, ring_size, size1, batch_size, map, nodeNet);
  //printf("res2: %d to %d\n", start, end+ill);
  //Rss_reveal(res2, ring_size, size2, batch_size, map, nodeNet);   

  Lint **Mult_tmp1 = new Lint *[2];
  Mult_tmp1[0] = new Lint [batch_size*totalSize];
  Mult_tmp1[1] = new Lint [batch_size*totalSize];
  
  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint [batch_size*totalSize];
  Mult_tmp2[1] = new Lint [batch_size*totalSize];

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint [batch_size*totalSize];
  Mult_tmp3[1] = new Lint [batch_size*totalSize];

  for(i = 0; i < size2; ++i) {
    for(j = 0; j < size1; ++j) {
      memcpy(&Mult_tmp1[0][(i*size1+j)*batch_size], &res2[0][i][0], sizeof(Lint)*batch_size);
      memcpy(&Mult_tmp1[1][(i*size1+j)*batch_size], &res2[1][i][0], sizeof(Lint)*batch_size);

      memcpy(&Mult_tmp2[0][(i*size1+j)*batch_size], &res1[0][j][0], sizeof(Lint)*batch_size);
      memcpy(&Mult_tmp2[1][(i*size1+j)*batch_size], &res1[1][j][0], sizeof(Lint)*batch_size);      
    }
  }
  Rss_Mult(Mult_tmp3, Mult_tmp2, Mult_tmp1, batch_size*totalSize, ring_size, map, nodeNet);
  for(i = 0; i < size2; ++i) {
    for(j = 0; j < size1; ++j) {
      memcpy(&res[0][size1*i+j][0], &Mult_tmp3[0][(i*size1+j)*batch_size], sizeof(Lint)*batch_size);
      memcpy(&res[1][size1*i+j][0], &Mult_tmp3[1][(i*size1+j)*batch_size], sizeof(Lint)*batch_size);
      for(int i2 = 0; i2 < batch_size; ++i2) {
	res[0][size1*i+j][i2] = res2[0][i][i2] + res1[0][j][i2] - res[0][size1*i+j][i2];
	res[1][size1*i+j][i2] = res2[1][i][i2] + res1[1][j][i2] - res[1][size1*i+j][i2];
      }
    }
  }
  //cleanup
  for(i = 0; i < 2; i++){
    for(j = 0; j < size1; ++j) {
      delete [] res1[i][j];
    }
    for(j = 0; j < size2; ++j) {
      delete [] res2[i][j];
    }
    delete [] res1[i];
    delete [] res2[i];
    delete [] Mult_tmp1[i];
    delete [] Mult_tmp2[i];
    delete [] Mult_tmp3[i];

  }
  delete [] Mult_tmp1;
  delete [] Mult_tmp2;
  delete [] Mult_tmp3;

  delete [] res1;
  delete [] res2;
}
