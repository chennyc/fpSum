#include "../include/Rss_Op.h"

void Rss_int2mask(Lint ***res, uint k, Lint **input, uint ring_size, uint batch_size, int *map, NodeNetwork* nodeNet, NodeConfiguration* nodeConfig){
    //res:[2][batch_size][k]
    int i, j;
    int pid = nodeConfig->getID();
    int compSize = k * batch_size;
    Lint **diff = new Lint*[2];
    Lint **compRes = new Lint*[2];
    for (i = 0; i < 2; i++) {
        diff[i] = new Lint [compSize];
        compRes[i] = new Lint [compSize];
    }

    Lint **constK = new Lint*[2];
    constK[0] = new Lint[k];
    constK[1] = new Lint[k];

    for(i = 0; i < k; ++i) {
        switch(pid){
            case 1:
                constK[0][i] = 0;
                constK[1][i] = 0;
                break;
            case 2:
                constK[0][i] = 0;
                constK[1][i] = i;
                break;
            case 3:
                constK[0][i] = i;
                constK[1][i] = 0;
                break;
        }
    }
    for(i = 0; i < batch_size; ++i) {
        for(j = 0; j < k; ++j) {
            diff[0][i * k + j] = constK[0][j] - input[0][i];
            diff[1][i * k + j] = constK[1][j] - input[1][i];
        }
    }
    Rss_MSB(compRes, diff, compSize, ring_size, map, nodeNet);
    for(i = 0; i < batch_size; ++i) {
        memcpy(&res[0][i][0], &compRes[0][i*k], sizeof(Lint) * k);
        memcpy(&res[1][i][0], &compRes[1][i*k], sizeof(Lint) * k);
    }
    for(i = 0; i < batch_size; ++i) {
        for(j = k-1; j >0; --j) {
            res[0][i][j] = res[0][i][j] - res[0][i][j-1];
            res[1][i][j] = res[1][i][j] - res[1][i][j-1];
        }
    }
    for(i = 0; i < 2; i++){
        delete [] diff[i];
        delete [] compRes[i];
        delete [] constK[i];
    }
    delete [] constK;
    delete [] compRes;
    delete [] diff;
}
