#ifndef BENCHMARK_FPSUM_H_
#define BENCHMARK_FPSUM_H_
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include "sys/time.h"
#include "Rss_Op.h"
extern "C"{
#include "aes_ni.h"
}

void benchmark_fpSum(NodeNetwork* nodeNet, NodeConfiguration* nodeConfig, int m, int e, int w, uint size, uint batch_size);

#endif



