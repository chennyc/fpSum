#include "../include/Rss_Op.h"

void Rss_PrefixMult(Lint ***res, Lint ***input, uint ring_size, uint size, uint batch_size, int *map, NodeNetwork *nodeNet, NodeConfiguration *nodeConfig)
{
  // res : [2] * [size] * [batch_size]
  uint i, j;
  uint numParties = nodeNet->getNumParties();
  uint rounds = ceil(log2(size));
  uint halfSize = pow(2, rounds - 1);
  int pid = nodeConfig->getID();
  Lint one[2]; // one: 0, 0, 1
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
  Mult_tmp1[0] = new Lint[halfSize * batch_size];
  Mult_tmp1[1] = new Lint[halfSize * batch_size];

  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint[halfSize * batch_size];
  Mult_tmp2[1] = new Lint[halfSize * batch_size];

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint[halfSize * batch_size];
  Mult_tmp3[1] = new Lint[halfSize * batch_size];

  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < size; ++j)
    {
      memcpy(&res[i][j][0], &input[i][j][0], sizeof(Lint) * batch_size);
    }
  }

  int index = 0;
  int offset = 2;
  for (i = 0; i < rounds; ++i)
  {
    int pow2i = pow(2, i);
    index = pow2i - 1;
    int multIndex = 0;
    while (index < size)
    {
      for (j = 1; j <= pow2i && index + j < size; ++j)
      {
        // Mult_tmp1[0][multIndex] = res[0][][];
        memcpy(&Mult_tmp1[0][multIndex], &res[0][index][0], sizeof(Lint) * batch_size);
        memcpy(&Mult_tmp1[1][multIndex], &res[1][index][0], sizeof(Lint) * batch_size);

        memcpy(&Mult_tmp2[0][multIndex], &res[0][index + j][0], sizeof(Lint) * batch_size);
        memcpy(&Mult_tmp2[1][multIndex], &res[1][index + j][0], sizeof(Lint) * batch_size);

        multIndex += batch_size;
      }
      index = index + offset;
    }
    // do mult
    Rss_Mult(Mult_tmp3, Mult_tmp1, Mult_tmp2, multIndex, ring_size, map, nodeNet);

    // take out mult result
    index = pow2i - 1;
    multIndex = 0;
    while (index < size)
    {
      for (j = 1; j <= pow2i && index + j < size; ++j)
      {
        memcpy(&res[0][index + j][0], &Mult_tmp3[0][multIndex], sizeof(Lint) * batch_size);
        memcpy(&res[1][index + j][0], &Mult_tmp3[1][multIndex], sizeof(Lint) * batch_size);
        multIndex += batch_size;
      }
      index = index + offset;
    }
    offset *= 2;
  }

  for (i = 0; i < 2; i++)
  {
    delete[] Mult_tmp1[i];
    delete[] Mult_tmp2[i];
    delete[] Mult_tmp3[i];
  }
  delete[] Mult_tmp1;
  delete[] Mult_tmp2;
  delete[] Mult_tmp3;
}

void Rss_PrefixOr_batch(Lint ***res, Lint ***input, uint ring_size, uint size, uint batch_size, int *map, NodeNetwork *nodeNet, NodeConfiguration *nodeConfig)
{
  // res : [2] * [size] * [batch_size]
  uint i, j;
  uint numParties = nodeNet->getNumParties();
  uint rounds = ceil(log2(size));
  uint halfSize = pow(2, rounds - 1);
  int pid = nodeConfig->getID();
  Lint one[2]; // one: 0, 0, 1
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
  Mult_tmp1[0] = new Lint[halfSize * batch_size];
  Mult_tmp1[1] = new Lint[halfSize * batch_size];

  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint[halfSize * batch_size];
  Mult_tmp2[1] = new Lint[halfSize * batch_size];

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint[halfSize * batch_size];
  Mult_tmp3[1] = new Lint[halfSize * batch_size];

  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < size; ++j)
    {
      memcpy(&res[i][j][0], &input[i][j][0], sizeof(Lint) * batch_size);
    }
  }

  int index = 0;
  int offset = 2;
  for (i = 0; i < rounds; ++i)
  {
    int pow2i = pow(2, i);
    index = pow2i - 1;
    int multIndex = 0;
    while (index < size)
    {
      for (j = 1; j <= pow2i && index + j < size; ++j)
      {
        // Mult_tmp1[0][multIndex] = res[0][][];
        memcpy(&Mult_tmp1[0][multIndex], &res[0][index][0], sizeof(Lint) * batch_size);
        memcpy(&Mult_tmp1[1][multIndex], &res[1][index][0], sizeof(Lint) * batch_size);

        memcpy(&Mult_tmp2[0][multIndex], &res[0][index + j][0], sizeof(Lint) * batch_size);
        memcpy(&Mult_tmp2[1][multIndex], &res[1][index + j][0], sizeof(Lint) * batch_size);

        multIndex += batch_size;
      }
      index = index + offset;
    }
    // do mult
    Rss_Mult(Mult_tmp3, Mult_tmp1, Mult_tmp2, multIndex, ring_size, map, nodeNet);

    // take out mult result
    index = pow2i - 1;
    multIndex = 0;
    while (index < size)
    {
      for (j = 1; j <= pow2i && index + j < size; ++j)
      {
        for (int k = 0; k < batch_size; ++k)
        {
          res[0][index + j][k] = Mult_tmp1[0][multIndex + k] + Mult_tmp1[0][multIndex + k] - Mult_tmp3[0][multIndex + k];
          res[1][index + j][k] = Mult_tmp1[1][multIndex + k] + Mult_tmp1[1][multIndex + k] - Mult_tmp3[1][multIndex + k];
        }
        multIndex += batch_size;
      }
      index = index + offset;
    }
    offset *= 2;
  }

  for (i = 0; i < 2; i++)
  {
    delete[] Mult_tmp1[i];
    delete[] Mult_tmp2[i];
    delete[] Mult_tmp3[i];
  }
  delete[] Mult_tmp1;
  delete[] Mult_tmp2;
  delete[] Mult_tmp3;
}

void Rss_PrefixOr(Lint **res, Lint **input, uint ring_size, uint size, int *map, NodeNetwork *nodeNet)
{
  // res : [2] * [size] *
  uint i, j;
  uint rounds = ceil(log2(size));
  uint halfSize = pow(2, rounds - 1);
  int pid = nodeNet->getID();
  Lint one[2]; // one: 0, 0, 1
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
  Mult_tmp1[0] = new Lint[halfSize];
  Mult_tmp1[1] = new Lint[halfSize];

  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint[halfSize];
  Mult_tmp2[1] = new Lint[halfSize];

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint[halfSize];
  Mult_tmp3[1] = new Lint[halfSize];

  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < size; ++j)
    {
      //memcpy(&res[i][j], &input[i][j], sizeof(Lint) );
      res[i][j] = input[i][j];
    }
  }

  int index = 0;
  int offset = 2;
  for (i = 0; i < rounds; ++i)
  {
    int pow2i = pow(2, i);
    index = pow2i - 1;
    int multIndex = 0;
    while (index < size)
    {
      for (j = 1; j <= pow2i && index + j < size; ++j)
      {
        // Mult_tmp1[0][multIndex] = res[0][][];
        Mult_tmp1[0][multIndex] = res[0][index];
        Mult_tmp1[1][multIndex] = res[1][index];

        Mult_tmp2[0][multIndex] = res[0][index + j];
        Mult_tmp2[1][multIndex] = res[1][index + j];

        multIndex ++;
      }
      index = index + offset;
    }
    // do mult
    Rss_Mult(Mult_tmp3, Mult_tmp1, Mult_tmp2, multIndex, ring_size, map, nodeNet);

    // take out mult result
    index = pow2i - 1;
    multIndex = 0;
    while (index < size)
    {
      for (j = 1; j <= pow2i && index + j < size; ++j)
      {
        res[0][index + j] = Mult_tmp1[0][multIndex] + Mult_tmp1[0][multIndex] - Mult_tmp3[0][multIndex];
        res[1][index + j] = Mult_tmp1[1][multIndex] + Mult_tmp1[1][multIndex] - Mult_tmp3[1][multIndex];
        multIndex ++;
      }
      index = index + offset;
    }
    offset *= 2;
  }

  for (i = 0; i < 2; i++)
  {
    delete[] Mult_tmp1[i];
    delete[] Mult_tmp2[i];
    delete[] Mult_tmp3[i];
  }
  delete[] Mult_tmp1;
  delete[] Mult_tmp2;
  delete[] Mult_tmp3;
}

void Rss_PrefixAnd(Lint **res, Lint **input, uint ring_size, uint size, int *map, NodeNetwork *nodeNet)
{
  // res : [2] * [size] *
  uint i, j;
  uint rounds = ceil(log2(size));
  uint halfSize = pow(2, rounds - 1);
  int pid = nodeNet->getID();
  Lint one[2]; // one: 0, 0, 1
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
  Mult_tmp1[0] = new Lint[halfSize];
  Mult_tmp1[1] = new Lint[halfSize];

  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint[halfSize];
  Mult_tmp2[1] = new Lint[halfSize];

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint[halfSize];
  Mult_tmp3[1] = new Lint[halfSize];

  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < size; ++j)
    {
      //memcpy(&res[i][j], &input[i][j], sizeof(Lint) );
      res[i][j] = input[i][j];
    }
  }

  int index = 0;
  int offset = 2;
  for (i = 0; i < rounds; ++i)
  {
    int pow2i = pow(2, i);
    index = pow2i - 1;
    int multIndex = 0;
    while (index < size)
    {
      for (j = 1; j <= pow2i && index + j < size; ++j)
      {
        // Mult_tmp1[0][multIndex] = res[0][][];
        Mult_tmp1[0][multIndex] = res[0][index];
        Mult_tmp1[1][multIndex] = res[1][index];

        Mult_tmp2[0][multIndex] = res[0][index + j];
        Mult_tmp2[1][multIndex] = res[1][index + j];

        multIndex ++;
      }
      index = index + offset;
    }
    // do mult
    Rss_Mult(Mult_tmp3, Mult_tmp1, Mult_tmp2, multIndex, ring_size, map, nodeNet);

    // take out mult result
    index = pow2i - 1;
    multIndex = 0;
    while (index < size)
    {
      for (j = 1; j <= pow2i && index + j < size; ++j)
      {
        res[0][index + j] = Mult_tmp3[0][multIndex];
        res[1][index + j] = Mult_tmp3[1][multIndex];
        multIndex ++;
      }
      index = index + offset;
    }
    offset *= 2;
  }

  for (i = 0; i < 2; i++)
  {
    delete[] Mult_tmp1[i];
    delete[] Mult_tmp2[i];
    delete[] Mult_tmp3[i];
  }
  delete[] Mult_tmp1;
  delete[] Mult_tmp2;
  delete[] Mult_tmp3;
}

// Actually, this is Product of size numbers
void Rss_CarryMult(Lint **res, Lint ***input, uint ring_size, uint size, uint batch_size, int *map, NodeNetwork *nodeNet)
{
  // res : [2] * [size] * [batch_size]
  uint i, j;
  uint halfSize = size / 2 + size % 2;
  int pid = nodeNet->getID();
  Lint one[2]; // one: 0, 0, 1
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
  Mult_tmp1[0] = new Lint[halfSize * batch_size];
  Mult_tmp1[1] = new Lint[halfSize * batch_size];

  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint[halfSize * batch_size];
  Mult_tmp2[1] = new Lint[halfSize * batch_size];

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint[size * batch_size];
  Mult_tmp3[1] = new Lint[size * batch_size];

  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < size; ++j)
    {
      memcpy(&Mult_tmp3[i][j * batch_size], &input[i][j][0], sizeof(Lint) * batch_size);
    }
  }

  int curSize = size;
  while (curSize > 1)
  {
    int halfCurSize = curSize / 2;
    for (i = 0; i < 2; ++i)
    {
      for (j = 0; j < halfCurSize; ++j)
      {
        memcpy(&Mult_tmp1[i][j * batch_size], &Mult_tmp3[i][2 * j * batch_size], sizeof(Lint) * batch_size);
        memcpy(&Mult_tmp2[i][j * batch_size], &Mult_tmp3[i][(2 * j + 1) * batch_size], sizeof(Lint) * batch_size);
      }
    }
    Rss_Mult(Mult_tmp3, Mult_tmp1, Mult_tmp2, halfCurSize * batch_size, ring_size, map, nodeNet);
    if (curSize % 2 == 1)
    {
      memcpy(&Mult_tmp3[0][halfCurSize * batch_size], &Mult_tmp3[0][(curSize - 1) * batch_size], sizeof(Lint) * batch_size);
      memcpy(&Mult_tmp3[1][halfCurSize * batch_size], &Mult_tmp3[1][(curSize - 1) * batch_size], sizeof(Lint) * batch_size);
    }
    curSize = halfCurSize + curSize % 2;
  }

  memcpy(&res[0][0], &Mult_tmp3[0][0], sizeof(Lint) * batch_size);
  memcpy(&res[1][0], &Mult_tmp3[1][0], sizeof(Lint) * batch_size);

  for (i = 0; i < 2; i++)
  {
    delete[] Mult_tmp1[i];
    delete[] Mult_tmp2[i];
    delete[] Mult_tmp3[i];
  }
  delete[] Mult_tmp1;
  delete[] Mult_tmp2;
  delete[] Mult_tmp3;
}

void Rss_CarryMult_batch(Lint **res, Lint ***input, uint ring_size, uint size, uint batch_size, int *map, NodeNetwork *nodeNet)
{
  // res : [2] * [size] * [batch_size]
  uint i, j;
  uint halfSize = size / 2 + size % 2;
  int pid = nodeNet->getID();
  Lint one[2]; // one: 0, 0, 1
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
  Mult_tmp1[0] = new Lint[halfSize * batch_size];
  Mult_tmp1[1] = new Lint[halfSize * batch_size];

  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint[halfSize * batch_size];
  Mult_tmp2[1] = new Lint[halfSize * batch_size];

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint[size * batch_size];
  Mult_tmp3[1] = new Lint[size * batch_size];

  for (i = 0; i < 2; ++i)
  {
    for (j = 0; j < size; ++j)
    {
      memcpy(&Mult_tmp3[i][j * batch_size], &input[i][j][0], sizeof(Lint) * batch_size);
    }
  }

  int curSize = size;
  while (curSize > 1)
  {
    int halfCurSize = curSize / 2;
    for (i = 0; i < 2; ++i)
    {
      for (j = 0; j < halfCurSize; ++j)
      {
        memcpy(&Mult_tmp1[i][j * batch_size], &Mult_tmp3[i][2 * j * batch_size], sizeof(Lint) * batch_size);
        memcpy(&Mult_tmp2[i][j * batch_size], &Mult_tmp3[i][(2 * j + 1) * batch_size], sizeof(Lint) * batch_size);
      }
    }
    Rss_Mult(Mult_tmp3, Mult_tmp1, Mult_tmp2, halfCurSize * batch_size, ring_size, map, nodeNet);
    if (curSize % 2 == 1)
    {
      memcpy(&Mult_tmp3[0][halfCurSize * batch_size], &Mult_tmp3[0][(curSize - 1) * batch_size], sizeof(Lint) * batch_size);
      memcpy(&Mult_tmp3[1][halfCurSize * batch_size], &Mult_tmp3[1][(curSize - 1) * batch_size], sizeof(Lint) * batch_size);
    }
    curSize = halfCurSize + curSize % 2;
  }

  memcpy(&res[0][0], &Mult_tmp3[0][0], sizeof(Lint) * batch_size);
  memcpy(&res[1][0], &Mult_tmp3[1][0], sizeof(Lint) * batch_size);

  for (i = 0; i < 2; i++)
  {
    delete[] Mult_tmp1[i];
    delete[] Mult_tmp2[i];
    delete[] Mult_tmp3[i];
  }
  delete[] Mult_tmp1;
  delete[] Mult_tmp2;
  delete[] Mult_tmp3;
}

void Rss_kor(Lint **res, Lint **input, uint ring_size, uint size, uint productRangeSize, int *map, NodeNetwork *nodeNet)
{
  // productRangeSize: how many bits we are going to do the product, i.e., prod(input) = input_(input_(productRangeSize)-1) * ... * input_(1) * input_(0)
  // res : [2] * [size]
  uint i, j;
  uint halfSize = productRangeSize * size / 2;
  uint halfSize_Lint = (halfSize + ring_size - 1) / ring_size;
  int pid = nodeNet->getID();

  Lint **Mult_tmp1 = new Lint *[2];
  Mult_tmp1[0] = new Lint[halfSize_Lint]();
  Mult_tmp1[1] = new Lint[halfSize_Lint]();

  Lint **Mult_tmp2 = new Lint *[2];
  Mult_tmp2[0] = new Lint[halfSize_Lint]();
  Mult_tmp2[1] = new Lint[halfSize_Lint]();

  Lint **Mult_tmp3 = new Lint *[2];
  Mult_tmp3[0] = new Lint[halfSize_Lint]();
  Mult_tmp3[1] = new Lint[halfSize_Lint]();

  for (i = 0; i < 2; ++i)
  {
    memcpy(&res[i][0], &input[i][0], sizeof(Lint) * size);
  }

  uint curSize = productRangeSize;
  Lint tmp = 0;
  while (curSize > 1)
  {

    // Rss_reveal(res, ring_size, size, map, nodeNet);
    uint bufferP0 = 0;
    uint bufferP1 = 0;
    int halfCurSize = curSize / 2;
    for (i = 0; i < size; ++i)
    {
      tmp = bitExtractedRange(res[0][i], halfCurSize, 0);
      FeedBufferBits(Mult_tmp1[0], tmp, halfCurSize, ring_size, bufferP0);

      tmp = bitExtractedRange(res[1][i], halfCurSize, 0);
      FeedBufferBits(Mult_tmp1[1], tmp, halfCurSize, ring_size, bufferP0);

      tmp = bitExtractedRange(res[0][i], halfCurSize, halfCurSize);
      FeedBufferBits(Mult_tmp2[0], tmp, halfCurSize, ring_size, bufferP0);

      tmp = bitExtractedRange(res[1][i], halfCurSize, halfCurSize);
      FeedBufferBits(Mult_tmp2[1], tmp, halfCurSize, ring_size, bufferP0);

      bufferP0 += halfCurSize;
    }
    uint mult_size = bufferP0 / ring_size;
    mult_size += bufferP0 % ring_size == 0 ? 0 : 1;

    // printf("Xor input after reaseemble: \n");
    // Rss_reveal(Mult_tmp1, ring_size, halfSize_Lint, map, nodeNet);
    // Rss_reveal(Mult_tmp2, ring_size, halfSize_Lint, map, nodeNet);

    Rss_Mult_Bitwise(Mult_tmp3, Mult_tmp1, Mult_tmp2, mult_size, ring_size, map, nodeNet);

    // Mult_tmp1 + Mult_tmp2 in bitwise
    for (int i = 0; i < mult_size; ++i)
    {
      Mult_tmp1[0][i] ^= Mult_tmp2[0][i];
      Mult_tmp1[1][i] ^= Mult_tmp2[1][i];
    }
    // Mult_tmp1 + Mult_tmp2 - (Mult_tmp1 x Mult_tmp2) in bitwise
    for (int i = 0; i < mult_size; ++i)
    {
      Mult_tmp3[0][i] ^= Mult_tmp1[0][i];
      Mult_tmp3[1][i] ^= Mult_tmp1[1][i];
    }

    //      memset(&res[0][0], 0, sizeof(Lint) * size);
    // memset(&res[1][0], 0, sizeof(Lint) * size);

    // printf("before take extra bit, curSize is %d,  halfCurSize is: %d\n", curSize, halfCurSize);
    // Rss_reveal(res, ring_size, size, map, nodeNet);

    if (curSize % 2 == 1)
    {
      for (i = 0; i < size; ++i)
      {
        res[0][i] = (res[0][i] & (1 << 2 * halfCurSize)) >> halfCurSize;
        res[1][i] = (res[1][i] & (1 << 2 * halfCurSize)) >> halfCurSize;
      }
    }
    else
    {
      memset(&res[0][0], 0, sizeof(Lint) * size);
      memset(&res[1][0], 0, sizeof(Lint) * size);
    }

    // printf("After taking extra size: \n");
    // Rss_reveal(res, ring_size, size, map, nodeNet);

    // printf("Xor result before reaseemble: \n");
    // Rss_reveal(Mult_tmp3, ring_size, halfSize_Lint, map, nodeNet);

    for (i = 0; i < size; ++i)
    {
      res[0][i] += TakeBufferBits(Mult_tmp3[0], halfCurSize, ring_size, bufferP1);
      res[1][i] += TakeBufferBits(Mult_tmp3[1], halfCurSize, ring_size, bufferP1);
      bufferP1 += halfCurSize;
    }
    // printf("cur round final res: \n");
    // Rss_reveal(res, ring_size, size, map, nodeNet);

    curSize = halfCurSize + curSize % 2;
    memset(&Mult_tmp1[0][0], 0, sizeof(Lint) * mult_size);
    memset(&Mult_tmp1[1][0], 0, sizeof(Lint) * mult_size);
    memset(&Mult_tmp2[0][0], 0, sizeof(Lint) * mult_size);
    memset(&Mult_tmp2[1][0], 0, sizeof(Lint) * mult_size);
    memset(&Mult_tmp3[0][0], 0, sizeof(Lint) * mult_size);
    memset(&Mult_tmp3[1][0], 0, sizeof(Lint) * mult_size);
  }
  for (i = 0; i < 2; i++)
  {
    delete[] Mult_tmp1[i];
    delete[] Mult_tmp2[i];
    delete[] Mult_tmp3[i];
  }
  delete[] Mult_tmp1;
  delete[] Mult_tmp2;
  delete[] Mult_tmp3;
}

// take out number_{start + length -1}, ..., number_{start}
Lint bitExtractedRange(Lint number, uint length, uint start)
{
  return (((1 << length) - 1) & (number >> (start)));
}

uint FeedBufferBits(Lint *buffer, Lint source, Lint bitLen, Lint ring_size, uint buffer_point)
{
  // bitLen: length of bits that going to be added into buffer
  // buffer_point: the start posistion(in bits) of free space in buffer
  uint index = buffer_point / ring_size;
  uint sub_index = buffer_point % ring_size;

  buffer[index] += source << sub_index;
  if (sub_index + bitLen > ring_size)
  {
    uint extra_len = sub_index + bitLen - ring_size;
    buffer[index + 1] += source >> (bitLen - extra_len);
  }
  return buffer_point + bitLen;
}

Lint TakeBufferBits(Lint *buffer, Lint bitLen, Lint ring_size, uint buffer_point)
{
  uint index = buffer_point / ring_size;
  uint sub_index = buffer_point % ring_size;
  Lint res = 0;
  if (sub_index + bitLen <= ring_size)
  {
    res = bitExtractedRange(buffer[index], bitLen, sub_index);
  }
  else
  {
    uint extra_len = sub_index + bitLen - ring_size;
    uint cur_len = bitLen - extra_len;
    res = bitExtractedRange(buffer[index], cur_len, sub_index);
    res += bitExtractedRange(buffer[index + 1], extra_len, 0) << cur_len;
    //    printf("one op: %lld vs two op: %lld \n", (bitExtractedRange(buffer[index], cur_len, sub_index) + bitExtractedRange(buffer[index+1], extra_len, 0)) << cur_len, res);
  }
  return res;
}
