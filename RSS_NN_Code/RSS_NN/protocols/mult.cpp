#include "../include/Rss_Op.h"
// extern "C"{
// #include "../aes_ni.h"
// }
void Rss_Mult_Byte(uint8_t** c, uint8_t** a, uint8_t** b, uint size, int *map, NodeNetwork *nodeNet) 
{
    // size == how many bytes we're multiplying
    // uint bytes = (size+8-1)>>3;  //number of bytes need to be send/recv
    // uint bytes = size;  //number of bytes need to be send/recv

    int i = 0;
    uint8_t *v = new uint8_t [size];
    nodeNet->prg_getrandom(1, 1, size, c[0]);  //<-- COMMENT OUT FOR TESTING
    // for(i = 0; i < 2; i++) memset(c[i], 0, sizeof(uint8_t)*size); //<-- COMMENT IN FOR TESTING
    
    for(i = 0 ; i < size; i++){  //do operations byte by byte
        v[i] =((a[0][i] & b[0][i]) ^ (a[0][i] & b[1][i]) ^ (a[1][i] & b[0][i])) ^ c[0][i];
    }

    //communication
    nodeNet->SendAndGetDataFromPeer_bit(map[0], map[1], v, c[1], size);
    for(i = 0; i < size; i++){
        c[1][i] = c[1][i] ^ c[0][i];
    }
    nodeNet->prg_getrandom(0, 1, size, c[0]); //<-- COMMENT OUT FOR TESTING
    for(i = 0; i < size; i++){
        c[0][i] = c[0][i] ^ v[i];
    }

    //free
    delete [] v;
}




