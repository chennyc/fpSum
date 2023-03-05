#ifndef BENCHMARK_NN_H_
#define BENCHMARK_NN_H_
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "sys/time.h"
#include "Rss_Op.h"
#include "neural_ops.h"
#include "svm_ops.h"
#include "init.h"
extern "C"{
#include "aes_ni.h"
}



void bench_main(NodeNetwork* nodeNet, NodeConfiguration* nodeConfig, uint batch_size);

void bench_mobilenet(NodeNetwork* nodeNet, NodeConfiguration* nodeConfig,int *map, uint batch_size);
void bench_squeeze(NodeNetwork* nodeNet, NodeConfiguration* nodeConfig,int *map, uint batch_size);
void bench_resnet(NodeNetwork* nodeNet, NodeConfiguration* nodeConfig,int *map, uint batch_size);



#endif



