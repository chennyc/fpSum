#include "../include/Rss_Op.h"

void Rss_reveal(Lint **toBeReveal, uint ring_size, uint size, int *map, NodeNetwork *nodeNet) {
  Lint* dtmp = new Lint[size];
  Rss_Open(dtmp, toBeReveal, size, map, ring_size, nodeNet);
  printf("\n");
  for(int i = 0; i < size; ++i){
    printf(" %llu ", dtmp[i]);
  }
  printf("\n");
  delete [] dtmp;
}

void Rss_reveal(bool flag, Lint **toBeReveal, uint ring_size, uint size, int *map, NodeNetwork *nodeNet) {
  Lint* dtmp = new Lint[size];
  Rss_Open_Bitwise(dtmp, toBeReveal, size, map, ring_size, nodeNet);
  printf("\n");
  for(int i = 0; i < size; ++i){
    printf(" %llu ", dtmp[i]);
  }
  printf("\n");
  delete [] dtmp;
}

void Rss_reveal(Lint ***toBeReveal, uint ring_size, uint size1, uint size2, int *map, NodeNetwork *nodeNet) {

  Lint* dtmp = new Lint[size1 * size2];
  Lint** toBeOpen = new Lint*[2];
  for(int i = 0; i < 2; ++i) {
    toBeOpen[i] = new Lint[size1 * size2];
  }
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < size1; ++j) {
      memcpy(&toBeOpen[i][j*size2], &toBeReveal[i][j][0], sizeof(Lint)*size2);
    }
  }
  Rss_Open(dtmp, toBeOpen, size1*size2, map, ring_size, nodeNet);
  for(int i = 0; i < size1; ++i){
    printf("\n");
    for(int j = 0; j < size2; ++j) {
      printf(" %llu ", dtmp[i*size2 + j]);
    }
  }
  printf("\n");
  delete [] dtmp;
  delete [] toBeOpen[0];
  delete [] toBeOpen[1];
  delete [] toBeOpen;
}
