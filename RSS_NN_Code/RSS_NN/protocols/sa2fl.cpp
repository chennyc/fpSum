#include "../include/Rss_Op.h"
void Rss_sa2fl(Lint *b, Lint *p, Lint **v, Lint **sa, int m, int e, int w, uint ring_size, int *map, NodeNetwork *nodeNet, int pid)
{
  int i, j, k;
  int alpha = ceil((double)(pow(2, e) + m) / w);
  int beta = ceil((double)(m + 1) / w) + 1;

  Lint one[2];
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

  Lint **c = new Lint *[2];
  Lint **d = new Lint *[2];
  Lint **cp = new Lint *[2];
  Lint **dp = new Lint *[2];
  Lint **tmp_sa = new Lint *[2];
  Lint **vp = new Lint *[2];

  for (i = 0; i < 2; ++i)
  {
    c[i] = new Lint[alpha - beta];
    d[i] = new Lint[alpha - beta];
    tmp_sa[i] = new Lint[alpha];
    cp[i] = new Lint[alpha - beta];
    dp[i] = new Lint[alpha];
    vp[i] = new Lint[beta];
  }

  // move sa[\beta] - sa[\alpha-1] to tmp_sa
  for (i = beta; i < alpha; ++i)
  {
    tmp_sa[0][i - beta] = sa[0][i];
    tmp_sa[1][i - beta] = sa[1][i];
  }

  Rss_eqz3(c, tmp_sa, alpha - beta, ring_size, map, nodeNet);

  for (i = 0; i < alpha - beta; ++i)
  {
    c[0][i] = one[0] - c[0][i];
    c[1][i] = one[1] - c[1][i];
  }

  Rss_PrefixAnd(d, c, ring_size, alpha - beta, map, nodeNet);

  for (i = beta; i <= alpha - 2; ++i)
  {
    dp[0][i] = d[0][i - beta] - d[0][i - beta + 1];
    dp[1][i] = d[1][i - beta] - d[1][i - beta + 1];
  }

  dp[0][alpha - 1] = one[0] - d[0][alpha - beta - 1];
  dp[1][alpha - 1] = one[0] - d[1][alpha - beta - 1];

  dp[0][beta - 1] = d[0][0];
  dp[1][beta - 1] = d[1][0];

  Rss_ORead(vp, sa, dp, alpha, beta, ring_size, map, nodeNet);

  Lint *pp = new Lint[2];

  Rss_normal(b, v, pp, vp, w, beta, m, ring_size, pid, map, nodeNet);

  for (i = 0; i < alpha - beta; ++i)
  {
    p[0] += dp[0][i + beta - 1] * i * w;
    p[1] += dp[1][i + beta - 1] * i * w;
  }

  p[0] += pp[0];
  p[1] += pp[1];

  // clean
  for (i = 0; i < 2; i++)
  {
    delete[] c[i];
    delete[] d[i];
    delete[] cp[i];
    delete[] dp[i];
    delete[] tmp_sa[i];
    delete[] vp[i];
  }
  delete[] c;
  delete[] d;
  delete[] cp;
  delete[] dp;
  delete[] tmp_sa;
  delete[] vp;
}

void Rss_ORead(Lint **res, Lint **superA, Lint **bitArray, int alpha, int beta, uint ring_size, int *map, NodeNetwork *nodeNet)
{
  int i, j;
  uint bytes = (ring_size + 7) >> 3;

  Lint *v = new Lint[beta];
  memset(v, 0, sizeof(Lint) * beta);

  uint8_t *buffer = new uint8_t[bytes * beta];
  nodeNet->prg_getrandom(1, bytes, beta, buffer);

  for (i = 0; i < beta; ++i)
  {
    for (j = i; j <= alpha - beta + i; ++j)
    {
      v[i] += bitArray[0][j + beta - 1 - i] * superA[0][j] + bitArray[0][j + beta - 1 - i] * superA[1][j] + bitArray[1][j + beta - 1 - i] * superA[0][j];
    }
    memcpy(&res[0][i], buffer + i * bytes, bytes);
    v[i] = v[i] - res[0][i];
  }

  nodeNet->SendAndGetDataFromPeer(map[0], map[1], v, res[1], beta, ring_size);
  nodeNet->prg_getrandom(0, bytes, beta, buffer);

  for (i = 0; i < beta; ++i)
  {
    res[1][i] = res[1][i] + res[0][i];
    memcpy(&res[0][i], buffer + i * bytes, bytes);
    res[0][i] = res[0][i] + v[i];
  }

  delete[] v;
}

void Rss_normal(Lint *b, Lint **v, Lint *p, Lint **vp, int w, int beta, int m, uint ring_size, int pid, int *map, NodeNetwork *nodeNet)
{
  int i, j;
  Lint **tmp = new Lint *[2];
  Lint **tmp_b = new Lint *[2];
  Lint **tmp_c = new Lint *[2];
  Lint ***bits = new Lint **[2];
  Lint **bitsp = new Lint *[2];
  Lint **dp = new Lint *[2];
  int ell = w * beta;
  LLint **vell = new LLint *[2];
  LLint *sell = new LLint [2];
  LLint **vpell = new LLint *[2];
  Lint **bits_tmp = new Lint *[2];
  for (i = 0; i < 2; ++i)
  {
    tmp[i] = new Lint[1];
    tmp_b[i] = new Lint[1];
    tmp_c[i] = new Lint[1];
    bits[i] = new Lint *[1];
    bits[i][0] = new Lint[w*beta];
    bitsp[i] = new Lint[w*beta];
    dp[i] = new Lint[w*beta];

    vell[i] = new LLint[beta];
    bits_tmp[i] = new Lint[w*beta];
    vpell[i] = new LLint[beta];
  }

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

  for(i = 0; i < 2; ++i) {
    for(j = 0; j < beta; ++j) {
      vpell[i][j] = vp[i][j];
    }
  }
  //new_Rss_Convert(vell, vp, beta, ring_size,  ell, map, nodeNet);
  new_Rss_Convert(vell, vpell, beta, ring_size,  ell, map, nodeNet);

  sell[0] = 0;
  sell[1] = 0;

  for(i = 0; i < beta; ++i) {
    sell[0] += vell[0][i] << (w*i);
    sell[1] += vell[1][i] << (w*i);
  }

  // msb
  tmp[0][0] = sell[0];
  tmp[1][0] = sell[1];

  new_Rss_MSB(tmp_b, tmp, 1, ring_size, map, nodeNet); // new

  tmp_b[0][0] = one[0] - 2 * tmp_b[0][0];
  tmp_b[1][0] = one[1] - 2 * tmp_b[1][0];

  b[0] = tmp_b[0][0];
  b[1] = tmp_b[1][0];

  // mult
  Rss_Mult(tmp_c, tmp, tmp_b, 1, ring_size, map, nodeNet);


  // BitDec
  Rss_bitDec3(bits, tmp_c, 1, ring_size, w*beta, map, nodeNet); //TODO check ring size

  for (i = m; i < ell - 1; ++i) {
    bits_tmp[0][i-m] = bits[0][0][i];
    bits_tmp[1][i-m] = bits[1][0][i];
  }

  Rss_PrefixOr(bits_tmp, bits_tmp, ring_size, ell - 1 - m, map, nodeNet);

  dp[0][beta*w - 1] =  bits[0][0][beta*w - 1];
  dp[1][beta*w - 1] =  bits[1][0][beta*w - 1];

  //for(i = m; i < ell - 2; ++i) {
  //  dp[0][i] = bits_tmp[0][i - m+1] - bits_tmp[0][i - m]; //TODO: Double check here
  //  dp[1][i] = bits_tmp[1][i - m+1] - bits_tmp[1][i - m];
  //}

  for (i = 0; i < 2; ++i)
  {
    memcpy(bitsp[i], bits[i][0], sizeof(Lint) * ell);
  }
  Rss_ORead2(dp, bitsp, dp, ell, m, ring_size, map, nodeNet);

  // clear
  for (i = 0; i < 2; ++i)
  {
    delete[] dp[i];
    delete[] tmp[i];
    delete[] tmp_b[i];
    delete[] tmp_c[i];
    delete[] bits[i][0];
    delete[] bits[i];
    delete[] bitsp[i];
  }
  delete[] dp;
  delete[] tmp_c;
  delete[] bits;
  delete[] bitsp;
  delete[] tmp;
  delete[] tmp_b;
}

void Rss_ORead2(Lint **res, Lint **bits, Lint **bitArray, int k, int kp, uint ring_size, int *map, NodeNetwork *nodeNet)
{
  int i, j;
  uint bytes = (ring_size + 7) >> 3;

  Lint *v = new Lint[kp + 1];
  memset(v, 0, sizeof(Lint) * (kp + 1));
  Lint *recv_buf = new Lint[kp + 1];
  memset(recv_buf, 0, sizeof(Lint) * (kp + 1));

  uint8_t *buffer = new uint8_t[bytes * (kp + 1)];
  nodeNet->prg_getrandom(1, bytes, kp + 1, buffer);

  for (i = 0; i < kp; ++i)
  {
    for (j = i; j <= k - kp - 1 + i; ++j)
    {
      v[i] += bitArray[0][j + kp - i] * bits[0][j] + bitArray[0][j + kp - i] * bits[1][j] + bitArray[1][j + kp - i] * bits[0][j];
    }
    memcpy(&res[0][i], buffer + i * bytes, bytes);
    v[i] = v[i] - res[0][i];
  }
  v[kp] = bitArray[0][kp] * bits[0][kp] + bitArray[0][kp] * bits[1][kp] + bitArray[1][kp] * bits[0][kp];
  memcpy(&bitArray[0][kp], buffer + kp * bytes, bytes);
  v[kp] = v[kp] - bitArray[0][kp];

  // memcpy(send_buf, v, sizeof(Lint)*kp);

  nodeNet->SendAndGetDataFromPeer(map[0], map[1], v, recv_buf, kp + 1, ring_size);
  nodeNet->prg_getrandom(0, bytes, kp + 1, buffer);
  memcpy(&res[1][0], recv_buf, sizeof(Lint) * kp);
  bitArray[1][kp] = recv_buf[kp];

  for (i = 0; i < kp; ++i)
  {
    res[1][i] = res[1][i] + res[0][i];
    memcpy(&res[0][i], buffer + i * bytes, bytes);
    res[0][i] = res[0][i] + v[i];
  }
  bitArray[1][kp] = bitArray[1][kp] + bitArray[0][kp];
  memcpy(&bitArray[0][kp], buffer + kp * bytes, bytes);
  bitArray[0][kp] = bitArray[0][kp] + v[kp];

  delete[] recv_buf;
  delete[] v;
}
