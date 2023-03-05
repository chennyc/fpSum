#ifndef NEURAL_OP_H_
#define NEURAL_OP_H_
#include <stdio.h>
#include <math.h>

#include "../../connection/NodeNetwork.h"
#include "Rss_Op.h"
#include "sys/time.h"


// neural network operations

// original
void im2col(Lint** res, Lint** a, uint channels, uint height, uint width, uint ksize, uint stride, uint batch_size, int *map, NodeNetwork* nodeNet);

// tensorflow version
void extract_patch(Lint** res, Lint** data_im, int channels, int height, int width, int ksize, int stride, uint batch_size, int pad_flag);

void pad_original(Lint** res, Lint** data_im, int channels, int height, int width, int ksize, int stride, uint batch_size);




void ReLU(Lint** res, Lint** a,  uint size, uint ring_size, int *map, NodeNetwork* nodeNet, uint flag);
void old_ReLU(Lint** res, Lint** a,  uint size, uint ring_size, int *map, NodeNetwork* nodeNet);
void eda_ReLU(Lint** res, Lint** a,  uint size, uint ring_size, int *map, NodeNetwork* nodeNet);

void MaxPool(Lint** res, Lint** a, uint c, uint m, uint n, uint batch_size, uint ring_sze, int *map, NodeNetwork* nodeNet,uint flag);
void old_MaxPool(Lint** res, Lint** a, uint c, uint m, uint n, uint batch_size, uint ring_sze, int *map, NodeNetwork* nodeNet);
void eda_MaxPool(Lint** res, Lint** a, uint c, uint m, uint n, uint batch_size, uint ring_sze, int *map, NodeNetwork* nodeNet);

void add_biases(Lint** res, Lint **a, Lint **b, uint m, uint n, uint batch_size,  int *map, NodeNetwork* nodeNet);



void conv(Lint** res, Lint **a, Lint **b, uint K, uint M, uint N, uint F, uint stride, uint batch_size, uint ring_size, int *map, NodeNetwork* nodeNet);

void conv_dw(Lint** res, Lint **a, Lint **b, uint K, uint M, uint F, uint stride, uint batch_size, uint ring_size, int *map, NodeNetwork* nodeNet);
void max_pool_bench(Lint** res, Lint **a, Lint **b, uint M, uint F, uint win, uint stride, uint batch_size, uint ring_size, int *map, NodeNetwork* nodeNet);

void fire_module(Lint** res, Lint **a, Lint **b, uint M,  uint F, uint squeeze, uint expand, uint stride, uint batch_size, uint ring_size, int *map, NodeNetwork* nodeNet);


#endif
