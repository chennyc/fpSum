#include "../include/benchmark_fpSum.h"

void benchmark_fpSum(NodeNetwork *nodeNet, NodeConfiguration *nodeConfig, int m, int e, int w, uint size, uint batch_size)
{
  int bits;
  int i, j, k;
  int pid = nodeConfig->getID();
  int flag = 0;
  int total;
  int ring_size = nodeNet->RING;

  int bytes = (ring_size + 8 - 1) / 8;

  printf("hello, I am %d\n", pid);
  printf("ring_size = %i\n", ring_size);
  printf("8*sizeof(Lint) = %i\n", 8 * sizeof(Lint));
  printf("sizeof(Lint) = %i\n", sizeof(Lint));

  int map[2];
  Lint one[2]; // one: 0, 0, 1
  switch (pid)
  {
  case 1:
    map[0] = 3;
    map[1] = 2;
    one[0] = 0;
    one[1] = 0;
    break;
  case 2:
    map[0] = 1;
    map[1] = 3;

    one[0] = 0;
    one[1] = 1;
    break;
  case 3:
    map[0] = 2;
    map[1] = 1;

    one[0] = 1;
    one[1] = 0;
    break;
  }
  // setup prg key (will be used by all parties, only for data makeup)
  __m128i *key_prg;
  uint8_t key_raw[] = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
  key_prg = offline_prg_keyschedule(key_raw);
  // setup prg seed(k1, k2, k3)
  uint8_t k1[] = {0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34};
  uint8_t k2[] = {0xa2, 0x34, 0x6f, 0x67, 0x10, 0x1b, 0x13, 0xa3, 0x56, 0x45, 0x90, 0xb2, 0x13, 0xe3, 0x23, 0x24};

  // make up data
  Lint **a;
  Lint **b;
  Lint **c; // c = a op b
  Lint **d;
  a = new Lint *[2];
  b = new Lint *[2];
  c = new Lint *[2];
  d = new Lint *[2];

  c[0] = new Lint[batch_size];
  memset(c[0], 0, sizeof(Lint) * batch_size);
  c[1] = new Lint[batch_size];
  memset(c[1], 0, sizeof(Lint) * batch_size);
  d[0] = new Lint[batch_size];
  memset(d[0], 0, sizeof(Lint) * batch_size);
  d[1] = new Lint[batch_size];
  memset(d[1], 0, sizeof(Lint) * batch_size);

  Lint **Data1;
  Lint **Data2;

  Data1 = new Lint *[3];
  for (int i = 0; i < 3; i++)
  {
    Data1[i] = new Lint[batch_size];
    memset(Data1[i], 0, sizeof(Lint) * batch_size);
  }

  Data2 = new Lint *[3];
  for (int i = 0; i < 3; i++)
  {
    Data2[i] = new Lint[batch_size];
    memset(Data2[i], 0, sizeof(Lint) * batch_size);
  }

  for (Lint i = 0; i < batch_size; i++)
  {
    prg_aes_ni(Data1[0] + i, k1, key_prg);
    prg_aes_ni(Data1[1] + i, k1, key_prg);
    prg_aes_ni(Data2[0] + i, k2, key_prg);
    prg_aes_ni(Data2[1] + i, k2, key_prg);
    Data1[2][i] = i - Data1[0][i] - Data1[1][i];
    Data2[2][i] = i - Data2[0][i] - Data2[1][i];
    // printf("Data[2][%i] : %i\n", i, Data2[2][i]);
  }
  free(key_prg);

  Lint *res = new Lint[batch_size];
  memset(res, 0, sizeof(Lint) * batch_size);
  Lint *res2 = new Lint[batch_size];
  memset(res2, 0, sizeof(Lint) * batch_size);
  Lint *res3 = new Lint[batch_size];
  memset(res3, 0, sizeof(Lint) * batch_size);

  Lint *res_check = new Lint[batch_size];
  memset(res_check, 0, sizeof(Lint) * batch_size);

  // assign data
  switch (pid)
  {
  case 1:
    a[0] = Data1[1];
    a[1] = Data1[2];
    b[0] = Data2[1];
    b[1] = Data2[2];
    break;
  case 2:
    a[0] = Data1[2];
    a[1] = Data1[0];
    b[0] = Data2[2];
    b[1] = Data2[0];
    break;
  case 3:
    a[0] = Data1[0];
    a[1] = Data1[1];
    b[0] = Data2[0];
    b[1] = Data2[1];
    break;
  }
  /*
  //make up data for prefixmult
  Lint ***x = new Lint **[2];
  Lint **xx = new Lint *[2];
  Lint **xxx = new Lint *[2];

  int len = 10;
  for(i = 0; i < 2; ++i) {
      x[i] = new Lint *[len];
xx[i] = new Lint[batch_size];
xxx[i] = new Lint[batch_size];

      for(j = 0; j < len; ++j) {
          x[i][j] = new Lint [batch_size];
          for(k = 0; k < batch_size; ++k) {
              x[i][j][k] = a[i][j+1];
          }
      }
  }


  {//for carrymult test
  printf("x before prefix\n");
  Rss_reveal(x, ring_size, len, batch_size, map, nodeNet);
  //delete [] dtmp;

  Rss_CarryMult(xx, x, ring_size, len, batch_size, map, nodeNet);

  printf("xx after carryMult\n");
  Rss_reveal(xx, ring_size, batch_size, map, nodeNet);

  Rss_PrefixMult(x, x, ring_size, len, batch_size, map, nodeNet, nodeConfig);
  printf("x after prefix \n");

  for(i = 0; i < 2; ++i) {
    for(j = 0; j <batch_size; ++j) {
xx[i][j] = a[i][j+1];
    }
  }
  Rss_shift(xxx, xx, xx, 8, batch_size, ring_size, map, nodeNet);
  printf("xxx after shift\n");
  Rss_reveal(xxx, ring_size, batch_size, map, nodeNet);

  Rss_reveal(x, ring_size, len, batch_size, map, nodeNet);
  for(i = 0; i < 2; ++i) {
    for(j = 0; j < len; ++j){
delete[] x[i][j];
    }
    delete [] x[i];
    delete [] xx[i];
    delete [] xxx[i];
  }
  delete[] x;
  delete[] xx;
  delete[] xxx;
  } //for carryMult test
  */


  /*
  int beta = ceil((double)(m + 1) / w) + 1;
  Lint **fpRes = new Lint *[2];
  Lint **fb = new Lint *[2];
  Lint **fp = new Lint *[2];
  Lint ***fv = new Lint **[2];

  Lint *res_b = new Lint[2];
  Lint *res_p = new Lint[2];
  Lint **res_v = new Lint*[2];
  for (i = 0; i < 2; ++i)
  {
    fpRes[i] = new Lint[batch_size];
    fb[i] = new Lint[batch_size];
    fp[i] = new Lint[batch_size];
    memset(&fp[i][0], 0, sizeof(Lint)*batch_size);
    fv[i] = new Lint*[batch_size];
    for(j = 0; j < batch_size; ++j) {
      fv[i][j] = new Lint[beta];
      memset(&fv[i][j][0], 0, sizeof(Lint)*beta);
    }
    res_v[i] = new Lint[batch_size];
  }

  struct timeval start;
  struct timeval end;
  unsigned long timers[6] = {0, 0, 0, 0, 0, 0};
  unsigned long timer;
  int rep = 1000;
  if (batch_size > pow(2, 14)) {
    rep = 100;
  }
  gettimeofday(&start, NULL);
  for (i = 0; i < rep; ++i)
  {
    //printf("start  %d th run \n", i);
    Rss_fpSum(res_b, res_p, res_v, fb, fp, fv, m, e, w, batch_size, ring_size, map, nodeNet, nodeConfig, timers);
  }
  gettimeofday(&end, NULL);
  timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  printf("runtime for fpSum with data size %d (m=%d, e=%d, w=%d) = %.6lf ms\n", batch_size, m, e, w, (double)(timer * 0.001) / rep);
  printf("runtime for shift timer = %.6lf ms\n", (double)(timers[0] * 0.001) / rep);
  printf("runtime for split timer = %.6lf ms\n", (double)(timers[1] * 0.001) / rep);
  printf("runtime for b2u timer = %.6lf ms\n", (double)(timers[2] * 0.001) / rep);
  printf("runtime for massDotProduct timer = %.6lf ms\n", (double)(timers[3] * 0.001) / rep);
  printf("runtime for superSum = %.6lf ms\n", (double)(timers[4] * 0.001) / rep);
  printf("runtime for sa2fl timer = %.6lf ms\n", (double)(timers[5] * 0.001) / rep);

  for (i = 0; i < 2; ++i)
  {
    delete[] fpRes[i];
    delete[] fb[i];
    delete[] fp[i];
    delete[] fv[i];
  }
  delete[] fpRes;
  delete[] fb;
  delete[] fp;
  delete[] fv;
*/
  /*
  { // check prefix-mult vs EQZ
    Lint ***x = new Lint **[2];
    int len = 132; // length for prefix
    int preifx_batch_size = 1;
    for (i = 0; i < 2; ++i)
    {
      x[i] = new Lint *[len];
      for (j = 0; j < len; ++j)
      {
        x[i][j] = new Lint[preifx_batch_size];
        for (k = 0; k < preifx_batch_size; ++k)
        {
          x[i][j][k] = a[i][j + 1];
        }
      }
    }

    Lint **eqzRes = new Lint *[2];
    for (i = 0; i < 2; ++i)
    {
      eqzRes[i] = new Lint[len]();
    }

    len = 15;
    gettimeofday(&start, NULL);
    Rss_PrefixMult(x, x, ring_size, len, preifx_batch_size, map, nodeNet, nodeConfig);
    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("runtime for prefixmult with data size len %d = %.6lf ms\n", len, (double)(timer * 0.001));

    gettimeofday(&start, NULL);
    Rss_eqz3(eqzRes, a, len, ring_size, map, nodeNet);
    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("runtime for eqz with data size len %d = %.6lf ms\n", len, (double)(timer * 0.001));

    len = 7;
    gettimeofday(&start, NULL);
    Rss_PrefixMult(x, x, ring_size, len, preifx_batch_size, map, nodeNet, nodeConfig);
    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("runtime for prefixmult with data size len %d = %.6lf ms\n", len, (double)(timer * 0.001));

    gettimeofday(&start, NULL);
    Rss_eqz3(eqzRes, a, len, ring_size, map, nodeNet);
    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("runtime for eqz with data size len %d = %.6lf ms\n", len, (double)(timer * 0.001));

    len = 127;
    gettimeofday(&start, NULL);
    Rss_PrefixMult(x, x, ring_size, len, preifx_batch_size, map, nodeNet, nodeConfig);
    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("runtime for prefixmult with data size len %d = %.6lf ms\n", len, (double)(timer * 0.001));

    gettimeofday(&start, NULL);
    Rss_eqz3(eqzRes, a, len, ring_size, map, nodeNet);
    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("runtime for eqz with data size len %d = %.6lf ms\n", len, (double)(timer * 0.001));

    len = 63;
    gettimeofday(&start, NULL);
    Rss_PrefixMult(x, x, ring_size, len, preifx_batch_size, map, nodeNet, nodeConfig);
    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("runtime for prefixmult with data size len %d = %.6lf ms\n", len, (double)(timer * 0.001));

    gettimeofday(&start, NULL);
    Rss_eqz3(eqzRes, a, len, ring_size, map, nodeNet);
    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("runtime for eqz with data size len %d = %.6lf ms\n", len, (double)(timer * 0.001));

    len = 132;
    for (i = 0; i < 2; ++i)
    {
      delete[] eqzRes[i];
    }
    delete[] eqzRes;
    for (i = 0; i < 2; ++i)
    {
      for (j = 0; j < len; ++j)
      {
        delete[] x[i][j];
      }
      delete[] x[i];
    }
    delete[] x;
  } // check prefix-mult vs EQZ END
  */

/*
  {//test Trunc
    Lint **Trucres = new Lint *[2];
    for(i = 0; i < 2; ++i) {
      Trucres[i] = new Lint [batch_size];
      memset(Trucres[i], 0, sizeof(Lint) * batch_size);
    }

    Lint **buf = new Lint *[2]; //-1
    for(i = 0; i < 2; ++i) {
      buf[i] = new Lint [batch_size];
      memset(buf[i], 0, sizeof(Lint) * batch_size);
      for(j = 0; j < batch_size; ++j) {
        buf[i][j] = 0 - one[i];
      }
    }

    Lint **buf2 = new Lint *[2]; //res
    Lint **sign = new Lint *[2]; //res
    for(i = 0; i < 2; ++i) {
      buf2[i] = new Lint [batch_size];
      memset(buf2[i], 0, sizeof(Lint) * batch_size);
      sign[i] = new Lint [batch_size];
      memset(sign[i], 0, sizeof(Lint) * batch_size);
    }

    for(i = 0; i < batch_size; ++i) {
      a[0][i] = one[0] << 18;
      a[1][i] = one[1] << 18;
    }

    printf("a: \n");
    Rss_reveal(a, ring_size, batch_size,  map, nodeNet);

    for(i = 0; i < batch_size; ++i) {
      //a[0][i] = a[0][i] - one[0];
      //a[1][i] = a[1][i] - one[1];
    }

    printf("a after -1: \n");
    Rss_reveal(a, ring_size, batch_size,  map, nodeNet);

    printf("-1: \n");
    Rss_reveal(buf, ring_size, batch_size,  map, nodeNet);

    Rss_Mult(buf2, a, buf, batch_size, ring_size, map, nodeNet);

    printf("(a)*(-1): \n");
    Rss_reveal(buf2, ring_size, batch_size,  map, nodeNet);

    new_Rss_MSB(sign, buf2, batch_size, ring_size, map, nodeNet);
    printf("sign: \n");
    Rss_reveal(sign, ring_size, batch_size,  map, nodeNet);
    for(i = 0; i < batch_size; ++i) {
      sign[0][i] = 0 - sign[0][i];
      sign[1][i] = 0 - sign[1][i];
    }

    Rss_truncPre(Trucres, buf, 16, batch_size, ring_size, map, nodeNet);

    printf("Truc res: \n");
    Rss_reveal(Trucres, ring_size, batch_size,  map, nodeNet);

    for(i = 0; i < 2; ++i) {
      delete [] Trucres[i];
      delete [] buf[i];
      delete [] buf2[i];
      delete [] sign[i];
    }
    delete[] Trucres;
    delete [] buf;
    delete [] buf2;
    delete [] sign;
  }//END OF Trunc
*/

/*
  {
  //test b2u
  int b2uLen = 120;
  Lint ***b2ures = new Lint **[2];
  for(i = 0; i < 2; ++i) {
      b2ures[i] = new Lint *[batch_size];
      for(j = 0; j < batch_size; ++j) {
          b2ures[i][j] = new Lint [b2uLen];
      }
  }
  Rss_b2u(b2ures, a, ring_size, batch_size, b2uLen, map, nodeNet);
  printf("test b2u  \n");
  Rss_reveal(a, ring_size, batch_size,  map, nodeNet);
  Rss_reveal(b2ures, ring_size, batch_size, b2uLen, map, nodeNet);
  for(i = 0; i < 2; ++i) {
    for(j = 0; j < batch_size; ++j){
      delete[] b2ures[i][j];
    }
    delete [] b2ures[i];
  }
  delete[] b2ures;
} end of b2u
*/


  //test b2a3
  {
  struct timeval start;
  struct timeval end;
  unsigned long timer;
  Lint **b2ares = new Lint *[2];
  Lint **b2ainput = new Lint *[2];
  for(i = 0; i < 2; ++i) {
    b2ares[i] = new Lint [batch_size]();
    b2ainput[i] = new Lint [batch_size]();
  }

  for(i = 0; i < batch_size; ++i) {
    b2ainput[0][i] = one[0];
    //      printf("%lld, %lld \n", a[0][i], a[1][i]);
    b2ainput[1][i] = one[1];
    //      printf("%lld, %lld \n", a[0][i], a[1][i]);
  }
  //    printf("%lld, %lld \n", one[0], one[1]);
  //    printf("%lld, %lld \n", b2ainput[0][0], b2ainput[1][0]);
  printf("b2a test\n");
  //    Rss_reveal(b2ainput, ring_size, batch_size,  map, nodeNet);
  //    Rss_reveal(b2ares, ring_size, batch_size,  map, nodeNet);
  gettimeofday(&start, NULL);
  int rep = 1000;
  for(i = 0; i < rep; ++i){
    Rss_b2a(b2ares, b2ainput, ring_size, batch_size, map, nodeNet);
  }
  gettimeofday(&end, NULL);
  timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  printf("runtime for b2a with data size %d = %.6lf ms\n", batch_size, (double)(timer * 0.001) / rep);

  gettimeofday(&start, NULL);
  for(i = 0; i < rep; ++i){
    Rss_b2a3(b2ares, b2ainput, ring_size, batch_size, map, nodeNet);
  }
  gettimeofday(&end, NULL);
  timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  printf("runtime for b2a3 with data size %d = %.6lf ms\n", batch_size, (double)(timer * 0.001) / rep);


  //Rss_reveal(b2ares, ring_size, batch_size,  map, nodeNet);

  for(i = 0; i < 2; ++i) {
    delete [] b2ares[i];
    delete [] b2ainput[i];
  }
  delete[] b2ares;
  delete[] b2ainput;
} // b2a
/*
{
  //test bitdec
  int tmpRingsize = ring_size;
  Lint ***bitDecRes = new Lint**[2];
  for(i = 0; i < 2; ++i) {
      bitDecRes[i] = new Lint* [batch_size];
for(j = 0; j < batch_size; ++j) {
  bitDecRes[i][j] = new Lint [ring_size]();
  //memset(bitDecRes[i][j], 0, sizeof(Lint) * ring_size);
      }
  }
  gettimeofday(&start, NULL);
  rep = 1000;

  printf("a: \n");
  Rss_reveal(a, ring_size, batch_size,  map, nodeNet);
  for(i = 0; i < rep; ++i){
    Rss_bitDec(bitDecRes, a, batch_size, ring_size, 5, map, nodeNet);
  }
  gettimeofday(&end, NULL);
  timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  printf("runtime for bitDec with data size %d (m=%d, e=%d, w=%d) = %.6lf ms\n", batch_size, m, e, w, (double)(timer * 0.001) / rep);
    Rss_reveal(bitDecRes, ring_size, batch_size, ring_size, map, nodeNet);
  for(i = 0; i < 2; ++i) {
    for(j = 0; j < batch_size; ++j){
delete[] bitDecRes[i][j];
    }
    delete [] bitDecRes[i];
  }
  delete[] bitDecRes;
} //bitDec
*/

  /*
  //test kor
  {
    Lint **korRes = new Lint *[2];
    for(i = 0; i < 2; ++i) {
korRes[i] = new Lint [batch_size];
    }

    Lint **kora = new Lint *[2];
    for(i = 0; i < 2; ++i) {
kora[i] = new Lint [batch_size];
    }

    Lint **kordata = new Lint *[3];
    for(i = 0; i < 3; ++i) {
kordata[i] = new Lint [batch_size]();
    }
    for(j = 0; j < batch_size; ++j) {
kordata[0][j] = j;
    }


    switch(pid){
      case 1:
          kora[0] = kordata[1];
          kora[1] = kordata[2];
          break;
      case 2:
          kora[0] = kordata[2];
          kora[1] = kordata[0];
          break;
      case 3:
          kora[0] = kordata[0];
          kora[1] = kordata[1];
          break;
    }


    printf("input a: \n");
    Rss_reveal(kora, ring_size, batch_size, map, nodeNet);

    Rss_kor(korRes, kora, ring_size, batch_size, ring_size, map, nodeNet);

    printf("kor: \n");
    Rss_reveal(korRes, ring_size, batch_size, map, nodeNet);

    memset(korRes[0], 0, sizeof(Lint) * batch_size);
    memset(korRes[1], 0, sizeof(Lint) * batch_size);

    Lint bitGet;
    Lint bitTmp = 1;
    bitGet = bitExtractedRange(bitTmp, ring_size/2, 0);
    printf("bitGet: %d \n", bitGet);
    FeedBufferBits(korRes[0], bitGet, ring_size/2, ring_size, 0);

    bitGet = bitExtractedRange(bitTmp, ring_size/2, ring_size/2);
    printf("bitGet: %d \n", bitGet);
    FeedBufferBits(korRes[1], bitGet, ring_size/2, ring_size, 0);
    printf("\n check buffer, M1_0, after move everything in: \n");
    for(int i = 0; i < batch_size; ++i) {
      printf(" %d ", korRes[0][i]);
    }
    printf("\n check buffer, M1_1,  after move everything in: \n");

    for(int i = 0; i < batch_size; ++i) {
      printf(" %d ", korRes[1][i]);
    }


    for(i = 0; i < 2; ++i) {
delete [] korRes[i];
    }
    delete[] korRes;
    for(i = 0; i < 3; ++i) {
delete [] kordata[i];
    }
    delete[] kordata;
    delete[] kora;
  }
  */

  /*
  {//test eqz
    Lint **eqzRes = new Lint *[2];
    for(i = 0; i < 2; ++i) {
      eqzRes[i] = new Lint [batch_size]();
    }

    printf("input a: \n");
    Rss_reveal(a, ring_size, batch_size, map, nodeNet);

    Rss_eqz3(eqzRes, a, batch_size, ring_size, map, nodeNet);

    printf("eqz: \n");
    Rss_reveal(eqzRes, ring_size, batch_size, map, nodeNet);

    for(i = 0; i < 2; ++i) {
      delete [] eqzRes[i];
    }
    delete[] eqzRes;

  }//test eqz
  */

  /*
  { //test edabit
    Lint **edabitR = new Lint *[2];
    Lint **edabitB = new Lint *[2];
    for(i = 0; i < 2; ++i) {
      edabitR[i] = new Lint [1]();
      edabitB[i] = new Lint [1]();
    }

    Lint **edabitB2 = new Lint *[2];
    for(i = 0; i < 2; ++i) {
      edabitB2[i] = new Lint [1 * ring_size]();
    }

    Rss_edaBit(edabitR, edabitB, 1, ring_size, map, nodeNet);

    for(i = 0; i < 2; ++i) {
      for(j = 0; j < 1; ++j) {
        for(k = 0; k < ring_size; ++k) {
          edabitB2[i][j * ring_size + k] = edabitB[i][j] >> k & 1;
        }
      }
    }

    Rss_b2a3(edabitB2, edabitB2, ring_size, 1 * ring_size, map, nodeNet);

    Lint* dtmp = new Lint[ring_size];
    Rss_Open(dtmp, edabitB2, ring_size, map, ring_size, nodeNet);

    Lint edabitRes = 0;
    for(i = 0; i < ring_size; ++i) {
      edabitRes += dtmp[i] << i;
    }
    printf(" %llu ", edabitRes);

    Rss_reveal(edabitR, ring_size, 1, map, nodeNet);

    // test the case k < l

    uint bit_length = ring_size / 2;
    for(i = 0; i < 2; ++i) {
      memset(edabitR[i], 0, sizeof(Lint) * 1);
      memset(edabitB[i], 0, sizeof(Lint) * 1);
      memset(edabitB2[i], 0, sizeof(Lint) * ring_size);
    }
    Rss_edaBit(edabitR, edabitB, 1, ring_size, bit_length, map, nodeNet);
    for(i = 0; i < 2; ++i) {
      for(j = 0; j < 1; ++j) {
        for(k = 0; k < bit_length; ++k) {
          edabitB2[i][j * bit_length + k] = edabitB[i][j] >> k & 1;
        }
      }
    }
    Rss_b2a3(edabitB2, edabitB2, ring_size, 1 * bit_length, map, nodeNet);

    memset(dtmp, 0, sizeof(Lint) * bit_length);
    Rss_Open(dtmp, edabitB2, bit_length, map, ring_size, nodeNet);

    edabitRes = 0;
    for(i = 0; i < bit_length; ++i) {
      edabitRes += dtmp[i] << i;
    }
    printf(" %llu ", edabitRes);

    Rss_reveal(edabitR, ring_size, 1, map, nodeNet);

    for(i = 0; i < 2; ++i) {
      delete [] edabitR[i];
      delete [] edabitB[i];
      delete [] edabitB2[i];
    }
    delete[] edabitR;
    delete[] edabitB;
    delete[] edabitB2;
    delete [] dtmp;
  }
  */

  // Free
  delete[] res;
  delete[] res2;
  delete[] res3;
  delete[] res_check;
  for (int i = 0; i < 3; i++)
  {

    delete[] Data1[i];
    delete[] Data2[i];
  }

  delete[] Data1;
  delete[] Data2;
  delete[] a;
  delete[] b;

  delete[] c[0];
  delete[] c[1];
  delete[] c;

  delete[] d[0];
  delete[] d[1];
  delete[] d;
}
