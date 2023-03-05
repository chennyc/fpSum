#ifndef Program_H_
#define Program_H_
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



// #include "benchmark_main.h"

#define GET_BIT_TEST(X, N) ( ( (X) >> (N) ) & 1 )


void test(NodeNetwork*, NodeConfiguration*, uint, uint);

void test2(NodeNetwork *nNet, NodeConfiguration *nodeConfig, int size);

#endif
