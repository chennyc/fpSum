#include "../include/Rss_Op.h"
void Rss_fpSum(Lint *res_b, Lint *res_p, Lint **res_v, Lint **b, Lint **p, Lint ***v, int m, int e, int w, uint batch_size, uint ring_size, int *map, NodeNetwork *nodeNet, NodeConfiguration *nodeConfig, unsigned long *timer)
{
  uint i, j, k;
  int alpha = ceil((double)(pow(2, e) + m) / w);
  int beta = ceil((double)(m + 1) / w) + 1;
  int logw = log2(w);
  int pid = nodeConfig->getID();
  Lint one[2]; // one: 0, 0, 1

  struct timeval start;
  struct timeval end;

  switch (pid)
  {
  case 1:
    one[0] = 0;
    one[1] = 0;
    break;
  case 2:
    one[0] = 0;
    one[1] = 1;
    break;
  case 3:
    one[0] = 1;
    one[1] = 0;
    break;
  }

  Lint **Mult_tmp1 = new Lint *[2];
  Mult_tmp1[0] = new Lint[batch_size * beta];
  Mult_tmp1[1] = new Lint[batch_size * beta];

  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint[batch_size * beta];
  Mult_tmp2[1] = new Lint[batch_size * beta];

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint[batch_size * beta];
  Mult_tmp3[1] = new Lint[batch_size * beta];

  // trunc p into pHigh and pLow
  Lint **pHigh = new Lint *[2];
  Lint **pLow = new Lint *[2];
  Lint **z = new Lint *[2];
  for (i = 0; i < 2; i++)
  {
    pHigh[i] = new Lint[batch_size];
    pLow[i] = new Lint[batch_size];
    z[i] = new Lint[batch_size];
  }

  Rss_truncPre(pHigh, p, logw, batch_size, ring_size, map, nodeNet);
  for (i = 0; i < batch_size; ++i)
  {
    pLow[0][i] = p[0][i] - pHigh[0][i] * w;
    pLow[1][i] = p[1][i] - pHigh[1][i] * w;
  }

  Rss_eqz3(z, p, batch_size, ring_size, map, nodeNet);

  for (i = 0; i < batch_size; ++i)
  {
    v[0][i][beta - 2] = v[0][i][beta - 2] + ((one[0] - z[0][i]) << (m - w * (beta - 2)));
    v[1][i][beta - 2] = v[1][i][beta - 2] + ((one[1] - z[1][i]) << (m - w * (beta - 2)));
  }

  Lint ***vblock = new Lint **[2];
  for (i = 0; i < 2; ++i)
  {
    vblock[i] = new Lint *[batch_size];
    for (j = 0; j < batch_size; ++j)
    {
      vblock[i][j] = new Lint[beta];
      memset(&vblock[i][j][0], 0, sizeof(Lint)*beta);
    }
  }
  // p1.2
  // shift v by [pLow]
  gettimeofday(&start, NULL);
  /*
  Lint ***intmask = new Lint**[2];
  for(i = 0; i < 2; ++i) {
      intmask[i] = new Lint*[batch_size];
      for(j = 0; j < batch_size; ++j) {
          intmask[i][j] = new Lint[w];
      }
  }
  //Rss_reveal(intmask, ring_size, batch_size, w, map, nodeNet);
  //Rss_shift(Mult_tmp2, v, pLow, logw, batch_size, ring_size, map, nodeNet);
  { //using b2u for shift
    Rss_b2u(intmask, pLow, ring_size, batch_size, w, map, nodeNet);

    for(i = 0; i < batch_size; ++i) {
      for(j = 0; j < w; ++j) {
  Mult_tmp1[0][i] += intmask[0][i][j] << j;
      }
    }
    Rss_Mult(Mult_tmp2, Mult_tmp1, v, batch_size, ring_size, map, nodeNet);
    for(i = 0; i < 2; ++i) {
      memcpy(&v[i][0], &Mult_tmp2[i][0], sizeof(Lint) * batch_size);
    }
  }
  */
  { // using bitDec for shift
    // Rss_shift(v, v, pLow, logw, batch_size, ring_size, map, nodeNet);
    Rss_newshift(vblock, v, pLow, logw, beta, batch_size, ring_size, map, nodeNet);
  }

  gettimeofday(&end, NULL);
  timer[0] += 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  // printf("runtime for shift with data size %d = %.6lf ms\n", batch_size, (double)(timer * 0.001));
  //     printf("%.6lf\n", (double)(timer * 0.001));

  // split

  //  gettimeofday(&start, NULL);
  //  for(i = 0; i < beta - 1; ++i) {
  //      Rss_truncPre(truncRes, v, w, batch_size, ring_size, map, nodeNet);
  //      for(j = 0; j < batch_size; ++j) {
  //          vblock[0][j][i] = v[0][j] - truncRes[0][j] << w;
  //          vblock[1][j][i] = v[1][j] - truncRes[1][j] << w;
  //      }
  //      memcpy(&v[0][0], &truncRes[0][0], sizeof(Lint) * batch_size);
  //      memcpy(&v[1][0], &truncRes[1][0], sizeof(Lint) * batch_size);
  //  }
  //  gettimeofday(&end, NULL);
  //  timer[1] += 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  // printf("runtime for split  with data size %d = %.6lf ms\n", batch_size, (double)(timer * 0.001));
  // printf("%.6lf\n", (double)(timer * 0.001));

  //for (i = 0; i < batch_size; ++i)
 // {
  //  vblock[0][i][beta - 1] = v[0][i];
   // vblock[1][i][beta - 1] = v[1][i];
  //}
  // b * vblock
  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < batch_size; ++j)
    {
      for (k = 0; k < beta; ++k)
      {
        Mult_tmp1[i][j * beta + k] = b[i][j];
      }
      memcpy(&Mult_tmp2[i][j * beta], &vblock[i][j][0], sizeof(Lint) * beta);
    }
  }

  Rss_Mult(Mult_tmp3, Mult_tmp1, Mult_tmp2, batch_size * beta, ring_size, map, nodeNet);
  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < batch_size; ++j)
    {
      memcpy(&vblock[i][j][0], &Mult_tmp2[i][j * beta], sizeof(Lint) * beta);
    }
  }

  // p2.1
  Lint ***superA = new Lint **[2];
  for (i = 0; i < 2; i++)
  {
    superA[i] = new Lint *[batch_size];
    for (j = 0; j < batch_size; j++)
    {
      superA[i][j] = new Lint[alpha];
    }
  }

  Lint ***intmask2 = new Lint **[2];
  for (i = 0; i < 2; ++i)
  {
    intmask2[i] = new Lint *[batch_size];
    for (j = 0; j < batch_size; ++j)
    {
      intmask2[i][j] = new Lint[alpha];
    }
  }

  gettimeofday(&start, NULL);
  // Rss_int2mask(intmask2, alpha, pHigh, ring_size, batch_size, map, nodeNet, nodeConfig);
  Rss_b2u(intmask2, pHigh, ring_size, batch_size, alpha, map, nodeNet);
  gettimeofday(&end, NULL);
  timer[2] += 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  // printf("runtime for 2nd int2mask with data size %d = %.6lf ms\n", batch_size, (double)(timer * 0.001));
  //  printf("%.6lf\n", (double)(timer * 0.001));

  gettimeofday(&start, NULL);
  Rss_MassDotProduct(superA, vblock, intmask2, batch_size, alpha, beta, ring_size, map, nodeNet);
  gettimeofday(&end, NULL);
  timer[3] += 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;

  Lint **sumAccum = new Lint *[2];
  for (i = 0; i < 2; ++i)
  {
    sumAccum[i] = new Lint[alpha];
  }
  // superSum
  gettimeofday(&start, NULL);
  Rss_SASum(sumAccum, superA, w, batch_size, alpha, ring_size, map, nodeNet);
  gettimeofday(&end, NULL);
  timer[4] += 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  // Normalization
  gettimeofday(&start, NULL);
  Rss_sa2fl(res_b, res_p, res_v, sumAccum, m, e, w, ring_size, map, nodeNet, pid);
  gettimeofday(&end, NULL);
  timer[5] += 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
  // clearup
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < batch_size; ++j)
    {
      delete[] vblock[i][j];
    }
    delete[] vblock[i];
    delete[] sumAccum[i];
  }
  delete[] vblock;
  delete[] sumAccum;
  // clearup

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < batch_size; ++j)
    {
      delete[] intmask2[i][j];
      delete[] superA[i][j];
    }
    delete[] superA[i];
    delete[] intmask2[i];
    delete[] pLow[i];
    delete[] pHigh[i];
    delete[] z[i];
    delete[] Mult_tmp1[i];
    delete[] Mult_tmp2[i];
    delete[] Mult_tmp3[i];
  }
  delete[] superA;
  delete[] intmask2;
  delete[] pLow;
  delete[] pHigh;
  delete[] z;
  delete[] Mult_tmp1;
  delete[] Mult_tmp2;
  delete[] Mult_tmp3;
}

void Rss_MassDotProduct(Lint ***res, Lint ***superA, Lint ***bitArray, int batch_size, int alpha, int beta, uint ring_size, int *map, NodeNetwork *nodeNet)
{
  uint bytes = (ring_size + 7) >> 3;
  int i, j, k;

  Lint *v = new Lint[alpha * batch_size];
  memset(v, 0, sizeof(Lint) * alpha * batch_size);
  Lint *recv_buf = new Lint[alpha * batch_size];
  memset(recv_buf, 0, sizeof(Lint) * alpha * batch_size);

  uint8_t *buffer = new uint8_t[bytes * alpha * batch_size];
  nodeNet->prg_getrandom(1, bytes, alpha * batch_size, buffer);

  for (i = 0; i < batch_size; ++i)
  {
    for (j = 0; j < alpha; ++j)
    {
      if (j <= beta - 2)
      {
        for (k = 0; k <= j; ++k)
        {
          v[i * alpha + j] += bitArray[0][i][j - k] * superA[0][i][k] + bitArray[0][i][j - k] * superA[1][i][k] + bitArray[1][i][j - k] * superA[0][i][k];
        }
      }

      if (j >= beta - 1 && j <= alpha - beta)
      {
        for (k = 0; k <= beta - 1; ++k)
        {
          v[i * alpha + j] += bitArray[0][i][j - k] * superA[0][i][k] + bitArray[0][i][j - k] * superA[1][i][k] + bitArray[1][i][j - k] * superA[0][i][k];
        }
      }

      if (j >= alpha - beta + 1)
      {
        for (k = 0; k <= alpha - 1 - j; ++k)
        {
          v[i * alpha + j] += bitArray[0][i][j - beta + 1 + k] * superA[0][i][beta - 1 - k] + bitArray[0][i][j - beta + 1 + k] * superA[1][i][beta - 1 - k] + bitArray[1][i][j - beta + 1 + k] * superA[0][i][beta - 1 - k];
        }
      }
      memcpy(&res[0][i][j], buffer + (i * alpha + j) * bytes, bytes);
      v[i * alpha + j] = v[i * alpha + j] - res[0][i][j];
    }
  }

  nodeNet->SendAndGetDataFromPeer(map[0], map[1], v, recv_buf, batch_size * alpha, ring_size);
  for (i = 0; i < batch_size; i++)
  {
    memcpy(res[1][i], recv_buf + i * alpha, sizeof(Lint) * alpha);
  }
  nodeNet->prg_getrandom(0, bytes, batch_size * alpha, buffer);
  for (i = 0; i < batch_size; i++)
  {
    for (j = 0; j < alpha; j++)
    {
      res[1][i][j] = res[1][i][j] + res[0][i][j];
      memcpy(res[0][i] + j, buffer + (i * alpha + j) * bytes, bytes);
      res[0][i][j] = res[0][i][j] + v[i * alpha + j];
    }
  }

  // free
  delete[] v;
  delete[] recv_buf;
  delete[] buffer;
}
