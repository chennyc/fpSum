#ifndef AES_NI_H
#define AES_NI_H
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>     //for int8_t
#include <string.h>     //for memcmp
#include <tmmintrin.h>
#include <wmmintrin.h>  //for intrinsics for AES-NI
#include <inttypes.h>
#include <unistd.h>
#include <stdint.h>     //for int8_t
#include <tmmintrin.h>
#include <wmmintrin.h>  //for intrinsics for AES-NI
#include <inttypes.h>

#include "../../setringsize.h"

void offline_prg(uint8_t * dest, uint8_t * src, void * ri);
__m128i * offline_prg_keyschedule(uint8_t * src);
void prg_aes_ni(Lint* destination, uint8_t * seed, __m128i * key);
void test_aes();
int print_u128_u(__uint128_t);
int print_u128_u2(__uint128_t);
void print_128(__uint128_t *A, int size);
void print_1283(__uint128_t A);
#endif
