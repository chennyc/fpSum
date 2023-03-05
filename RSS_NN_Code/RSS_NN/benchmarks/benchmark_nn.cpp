#include "../include/benchmark_nn.h"

void bench_main(NodeNetwork *nodeNet, NodeConfiguration *nodeConfig, uint batch_size) {

    int bits;
    int pid = nodeConfig->getID();
    int flag = 0;
    int ring_size = nodeNet->RING;

    int total;

    printf("ring_size = %i\n", ring_size);
    printf("8*sizeof(Lint) = %i\n", 8 * sizeof(Lint));
    printf("sizeof(Lint) = %i\n", sizeof(Lint));
    printf("batch_size = %i\n", batch_size);

    // checking if we can access this built-in function
    // uint aaa = _pext_u32((unsigned) ring_size, 0x55555555);
    printf("hello, I am %d\n", pid);

    int bytes = (ring_size + 2 + 8 - 1) / 8;
    int map[2];
    switch (pid) {
    case 1:
        map[0] = 3;
        map[1] = 2;
        break;
    case 2:
        map[0] = 1;
        map[1] = 3;
        break;
    case 3:
        map[0] = 2;
        map[1] = 1;
        break;
    }

    bench_mobilenet(nodeNet, nodeConfig, map, batch_size);
}

void bench_mobilenet(NodeNetwork *nodeNet, NodeConfiguration *nodeConfig, int *map, uint batch_size) {

    int ring_size = nodeNet->RING;
    uint size = 51380224 * batch_size;

    uint total = 10;

    struct timeval start;
    struct timeval end;
    unsigned long timer;

    //make up data
    Lint **a = new Lint *[2];
    Lint **b = new Lint *[2];
    Lint **c = new Lint *[2];
    Lint **d = new Lint *[2];

    for (int i = 0; i < 2; i++) {
        c[i] = new Lint[size];
        memset(c[i], 0, sizeof(Lint) * size);
        d[i] = new Lint[size];
        memset(d[i], 0, sizeof(Lint) * size);
    }

    Lint *res = new Lint[size];
    // memset(res, 0, sizeof(Lint) * size);

    generateData(a, size, ring_size, nodeNet);
    generateData(b, size, ring_size, nodeNet);

    gettimeofday(&start, NULL); // start timer here
    for (int i = 0; i < total; i++) {

        // (K) Filter shape - (K x K x M x N)
        // (F) Input size - (F x F x M)
        // (G) Output size - (F x F x N)
        conv(c, a, b, 3, 3, 32, 224, 2, batch_size, ring_size, map, nodeNet);
        // (K) Filter shape - (K x K x M)
        // (F) Input size - (F x F x M)
        // (G) Output size - (F x F x M)
        conv_dw(c, a, b, 3, 32, 112, 1, batch_size, ring_size, map, nodeNet);

        conv(c, a, b, 1, 32, 64, 112, 1, batch_size, ring_size, map, nodeNet);
        conv_dw(c, a, b, 3, 64, 112, 2, batch_size, ring_size, map, nodeNet);

        conv(c, a, b, 1, 64, 128, 56, 1, batch_size, ring_size, map, nodeNet);
        conv_dw(c, a, b, 3, 128, 56, 1, batch_size, ring_size, map, nodeNet);

        conv(c, a, b, 1, 128, 128, 56, 1, batch_size, ring_size, map, nodeNet);
        conv_dw(c, a, b, 3, 128, 56, 2, batch_size, ring_size, map, nodeNet);

        conv(c, a, b, 1, 128, 256, 28, 1, batch_size, ring_size, map, nodeNet);
        conv_dw(c, a, b, 3, 256, 28, 1, batch_size, ring_size, map, nodeNet);

        conv(c, a, b, 1, 256, 256, 28, 1, batch_size, ring_size, map, nodeNet);
        conv_dw(c, a, b, 3, 256, 28, 2, batch_size, ring_size, map, nodeNet);

        conv(c, a, b, 1, 256, 512, 14, 1, batch_size, ring_size, map, nodeNet);

        for (int i = 0; i < 5; i++) {
            conv_dw(c, a, b, 3, 512, 14, 1, batch_size, ring_size, map, nodeNet);
            conv(c, a, b, 1, 512, 512, 14, 1, batch_size, ring_size, map, nodeNet);
        }

        conv_dw(c, a, b, 3, 512, 14, 2, batch_size, ring_size, map, nodeNet);
        conv(c, a, b, 1, 512, 1024, 7, 1, batch_size, ring_size, map, nodeNet);

        conv_dw(c, a, b, 3, 1024, 7, 2, batch_size, ring_size, map, nodeNet);
        conv(c, a, b, 1, 1024, 1024, 7, 1, batch_size, ring_size, map, nodeNet);

        // avg pool

        // FC layer
        Rss_MatMultArray_batch(c, a, b, 1000, 1024, 1, ring_size, batch_size, 0, 0, map, nodeNet);

        // softmax
        // replace with arg_max
        eda_arg_max(c, a, b, ring_size, 1000, batch_size, map, nodeNet);
    }
    gettimeofday(&end, NULL); // stop timer here
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;

    printf("total runtime for MobileNet exeuction = %lf ms\n", (double)(timer * 0.001));
    printf("average prediction time = %lf ms\n", (double)(timer * 0.001) / total);

    delete[] res;

    for (int i = 0; i < 2; i++) {
        delete[] a[i];
        delete[] b[i];
        delete[] c[i];
        delete[] d[i];
    }
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
}


void bench_squeeze(NodeNetwork *nodeNet, NodeConfiguration *nodeConfig, int *map, uint batch_size) {

    int ring_size = nodeNet->RING;
    uint size = 51380224 * batch_size; // need to recompute

    uint total = 10;

    struct timeval start;
    struct timeval end;
    unsigned long timer;

    //make up data
    Lint **a = new Lint *[2];
    Lint **b = new Lint *[2];
    Lint **c = new Lint *[2];
    Lint **d = new Lint *[2];

    for (int i = 0; i < 2; i++) {
        c[i] = new Lint[size];
        memset(c[i], 0, sizeof(Lint) * size);
        d[i] = new Lint[size];
        memset(d[i], 0, sizeof(Lint) * size);
    }

    Lint *res = new Lint[size];
    // memset(res, 0, sizeof(Lint) * size);

    generateData(a, size, ring_size, nodeNet);
    generateData(b, size, ring_size, nodeNet);

    gettimeofday(&start, NULL); // start timer here
    for (int i = 0; i < total; i++) {

        // (K) Filter shape - (K x K x M x N)
        // (F) Input size - (F x F x M)
        // (G) Output size - (F x F x N)
        conv(c, a, b, 3, 3, 64, 32, 1, batch_size, ring_size, map, nodeNet);


        // max_pool <4,2,2,1>*(32/2)^2 * 64
        // max_pool_bench(c, a, b, 64, 32, 3, 2, batch_size, ring_size, map, nodeNet);

        fire_module(c, a, b, 64, 16, 32, 64, 1, batch_size, ring_size, map, nodeNet);

        fire_module(c, a, b, 128, 16, 32, 64, 1, batch_size, ring_size, map, nodeNet);




        max_pool_bench(c, a, b, 128, 16, 3, 2, batch_size, ring_size, map, nodeNet);
        // max_pool <4,2,2,1>*(16/2)^2 * 128

        fire_module(c, a, b, 128, 8, 32, 128, 1, batch_size, ring_size, map, nodeNet);

        fire_module(c, a, b, 256, 8, 32, 128, 1, batch_size, ring_size, map, nodeNet);

        // size = 1, filters = 10, stride = 1, input = (8x8x256)
        conv(c, a, b, 1, 256, 10, 8, 1, batch_size, ring_size, map, nodeNet);

        // avg pool

    }
    gettimeofday(&end, NULL); // stop timer here
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;

    printf("total runtime for squeezenet exeuction = %lf ms\n", (double)(timer * 0.001));
    printf("average prediction time = %lf ms\n", (double)(timer * 0.001) / total);

    delete[] res;

    for (int i = 0; i < 2; i++) {
        delete[] a[i];
        delete[] b[i];
        delete[] c[i];
        delete[] d[i];
    }
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

}





void bench_resnet(NodeNetwork *nodeNet, NodeConfiguration *nodeConfig, int *map, uint batch_size) {

    int ring_size = nodeNet->RING;
    uint size = 51380224 * batch_size; // need to recompute

    uint total = 10;

    struct timeval start;
    struct timeval end;
    unsigned long timer;

    //make up data
    Lint **a = new Lint *[2];
    Lint **b = new Lint *[2];
    Lint **c = new Lint *[2];
    Lint **d = new Lint *[2];

    for (int i = 0; i < 2; i++) {
        c[i] = new Lint[size];
        memset(c[i], 0, sizeof(Lint) * size);
        d[i] = new Lint[size];
        memset(d[i], 0, sizeof(Lint) * size);
    }

    Lint *res = new Lint[size];
    // memset(res, 0, sizeof(Lint) * size);

    generateData(a, size, ring_size, nodeNet);
    generateData(b, size, ring_size, nodeNet);

    gettimeofday(&start, NULL); // start timer here
    for (int i = 0; i < total; i++) {



    }
    gettimeofday(&end, NULL); // stop timer here
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;

    printf("total runtime for NN exeuction = %lf ms\n", (double)(timer * 0.001));
    printf("average prediction time = %lf ms\n", (double)(timer * 0.001) / total);

    delete[] res;

    for (int i = 0; i < 2; i++) {
        delete[] a[i];
        delete[] b[i];
        delete[] c[i];
        delete[] d[i];
    }
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

}
