#include "../include/Rss_Op.h"
void Rss_SASum(Lint **res, Lint ***superA, int w, uint batch_size, int alpha, uint ring_size, int *map, NodeNetwork *nodeNet)
{
  uint i, j, k;
  if (batch_size <= pow(2, w - 2))
  {
    for (k = 0; k < 2; ++k)
    {
      for (i = 0; i < batch_size; ++i)
      {
        for (j = 0; j < alpha; ++j)
        {
          res[k][j] += superA[k][i][j];
        }
      }
    }
    Rss_superSum(res, res, w, alpha, ring_size, map, nodeNet);
    return;
  }
  //printf("batch size is too large to be processed in one batch \n");
  int chunk_size = pow(2, w - 2);
  int chunk_no = batch_size / chunk_size;

  Lint ***allChunks = new Lint **[2];
  for (i = 0; i < 2; ++i)
  {
    allChunks[i] = new Lint *[chunk_no];
    for (j = 0; j < chunk_no; ++j)
    {
      allChunks[i][j] = new Lint [alpha];
    }
  }

  for (k = 0; k < 2; ++k)
  {
    for (i = 0; i < batch_size; ++i)
    {
      for (j = 0; j < alpha; ++j)
      {
        allChunks[k][i / chunk_size][j] += superA[k][i][j];
      }
    }
  }

  Rss_superSum_batch(allChunks, allChunks, w, chunk_no, alpha, ring_size, map, nodeNet);
  Rss_SASum(res, allChunks, w, chunk_no, alpha, ring_size, map, nodeNet);
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < chunk_no; ++j)
    {
      delete[] allChunks[i][j];
    }
    delete[] allChunks[i];
  }
  delete[] allChunks;
}

void Rss_superSum_batch(Lint ***res, Lint ***input, int w, uint batch_size, int alpha, uint ring_size, int *map, NodeNetwork *nodeNet)
{
  uint i, j, k;
  Lint one[2]; // one: 0, 0, 1
  int pid = nodeNet->getID();
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
  // truc
  Lint **truncRes = new Lint *[2];
  Lint **input_tmp = new Lint *[2];
  Lint **r = new Lint *[2];
  Lint **b = new Lint *[2];
  Lint **y = new Lint *[2];
  for (i = 0; i < 2; ++i)
  {
    truncRes[i] = new Lint[alpha*batch_size];
    input_tmp[i] = new Lint[alpha*batch_size];
    r[i] = new Lint[alpha*batch_size];
    b[i] = new Lint[alpha*batch_size];
    y[i] = new Lint[alpha*batch_size];
  }
  for (i = 0; i < 2; ++i) {
    for(j = 0; j < batch_size; ++j) {
      memcpy(&input_tmp[i][j*alpha], input[i][j], sizeof(Lint)*alpha);
    }
  }

  new_Rss_MSB(b, input_tmp, alpha*batch_size, ring_size, map, nodeNet);
  for (i = 0; i < alpha*batch_size; ++i)
  {
    b[0][i] = one[0] - 2 * b[0][i];
    b[1][i] = one[1] - 2 * b[1][i];
  }

  Rss_Mult(y, b, input_tmp, alpha*batch_size, ring_size, map, nodeNet);

  Rss_truncPre(truncRes, y, w, alpha*batch_size, ring_size, map, nodeNet);

  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < alpha*batch_size; ++j)
    {
      r[i][j] = y[i][j] - (truncRes[i][j] << w);
    }
  }

  for (i = 0; i < 2; ++i)
  {
    for (j = 1; j < alpha; ++j)
    {
      for(k = 0; k < batch_size; ++k) {
        res[i][k][j] = r[i][k*alpha + j] + truncRes[i][k*alpha + j - 1];
      }
    }
  }

  for (i = 0; i < 2; i++)
  {
    delete[] truncRes[i];
    delete[] input_tmp[i];
    delete[] r[i];
    delete[] b[i];
    delete[] y[i];
  }
  delete[] input_tmp;
  delete[] r;
  delete[] b;
  delete[] y;
  delete[] truncRes;
}

void Rss_superSum(Lint **res, Lint **input, int w, int alpha, uint ring_size, int *map, NodeNetwork *nodeNet)
{
  uint i, j, k;
  Lint one[2]; // one: 0, 0, 1
  int pid = nodeNet->getID();
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
  // truc
  Lint **truncRes = new Lint *[2];
  Lint **r = new Lint *[2];
  Lint **b = new Lint *[2];
  Lint **y = new Lint *[2];
  for (i = 0; i < 2; ++i)
  {
    truncRes[i] = new Lint[alpha];
    r[i] = new Lint[alpha];
    b[i] = new Lint[alpha];
    y[i] = new Lint[alpha];
  }

  new_Rss_MSB(b, input, alpha, ring_size, map, nodeNet);
  for (i = 0; i < alpha; ++i)
  {
    b[0][i] = one[0] - 2 * b[0][i];
    b[1][i] = one[1] - 2 * b[1][i];
  }

  Rss_Mult(y, b, input, alpha, ring_size, map, nodeNet);

  Rss_truncPre(truncRes, y, w, alpha, ring_size, map, nodeNet);

  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < alpha; ++j)
    {
      r[i][j] = y[i][j] - (truncRes[i][j] << w);
    }
  }

  for (i = 0; i < 2; ++i)
  {
    for (j = 1; j < alpha; ++j)
    {
      res[i][j] = r[i][j] + truncRes[i][j - 1];
    }
  }

  for (i = 0; i < 2; i++)
  {
    delete[] truncRes[i];
    delete[] r[i];
    delete[] b[i];
    delete[] y[i];
  }
  delete[] r;
  delete[] b;
  delete[] y;
  delete[] truncRes;
}
