#include "../include/neural_ops.h"

// does the actual "convolution" step required prior to matrix multiplication
// ksize = kernelsize
// we assume the size of res is sufficient to handle
void im2col(Lint **res, Lint **a, uint channels, uint height, uint width, uint ksize, uint stride, uint batch_size, int *map, NodeNetwork *nodeNet) {

    for (uint i = 0; i < batch_size; i++) {
        uint pad = 0;

        uint c, h, w;
        uint height_col = (height + 2 * pad - ksize) / stride + 1;
        uint width_col = (width + 2 * pad - ksize) / stride + 1;
        uint channels_col = channels * ksize * ksize;

        uint w_offset, h_offset, c_im, im_row, im_col, col_index;

        // used for batch optimization
        uint in_batch_size = channels * height * width;
        uint out_batch_size = channels_col * height_col * width_col;

        for (c = 0; c < channels_col; ++c) {
            w_offset = c % ksize;
            h_offset = (c / ksize) % ksize;
            c_im = c / ksize / ksize;
            for (h = 0; h < height_col; ++h) {
                for (w = 0; w < width_col; ++w) {

                    im_row = h_offset + h * stride;
                    im_col = w_offset + w * stride;
                    col_index = (c * height_col + h) * width_col + w;

                    if (im_row < 0 || im_col < 0 || (im_row >= height) || (im_col >= width)) {
                        res[0][col_index] = 0;
                        res[1][col_index] = 0;

                        // check if correct
                        // res[0][col_index + i * out_batch_size] = 0;
                        // res[1][col_index + i * out_batch_size] = 0;

                    } else {
                        res[0][col_index + i * out_batch_size] = a[0][i * in_batch_size + im_col + width * (im_row + height * c_im)];
                        res[1][col_index + i * out_batch_size] = a[1][i * in_batch_size + im_col + width * (im_row + height * c_im)];
                    }
                }
            }
        }
    }
}

void extract_patch(Lint **res, Lint **data_im, int channels, int height, int width, int ksize, int stride, uint batch_size, int pad_flag) {


    int c, h, w, pad;
    int pad_along_height, pad_along_width, pad_top, pad_bottom, pad_left, pad_right;

    int height_col, width_col, channels_col;

    if (pad_flag == 0) {
        // we have VALID padding (so none at all)
        pad_top = 0;
        pad_left = 0;
        pad = 0;

        height_col = floor((height + 2 * pad - ksize) / stride + 1);
        width_col = floor((width + 2 * pad - ksize) / stride + 1);

    } else {
        // SAME PADDING

        height_col = ceil(float(height) / stride);
        width_col = ceil(float(width) / stride);

        if ((height % stride) == 0) {
            pad_along_height = max(ksize - stride, 0);

        } else {
            pad_along_height = max(ksize - (height % stride), 0);
        }
        if ((width % stride) == 0) {
            pad_along_width = max(ksize - stride, 0);

        } else {
            pad_along_width = max(ksize - (width % stride), 0);
        }

        if (!(stride > 1)) {

            pad_top = floor(float(pad_along_height) / 2.0);
            pad_bottom = pad_along_height - pad_top;
            pad_left = floor((float(pad_along_width) / 2.0));
            pad_right = pad_along_width - pad_left;

        } else {
            pad_left = channels;
            // pad_left = 8;
            pad_top = 1;
        }
    }

    channels_col = channels * width_col * height_col;

    uint in_batch_size = channels * height * width;

    uint out_batch_size = ksize * ksize * channels * width_col * height_col;


    for (uint i = 0; i < batch_size; i++) {

        for (c = 0; c < width_col * height_col; ++c) {
            int w_offset = c % width_col;
            int h_offset = (c / height_col) % height_col;
            int c_im = c / width_col / width_col;

            for (h = 0; h < ksize; ++h) {
                for (w = 0; w < ksize * channels; ++w) {
                    int im_row = h_offset * stride + h;
                    int im_col = w_offset * stride * channels + w;
                    int col_index = (c * ksize + h) * ksize * channels + w;

                    im_row -= pad_top;
                    im_col -= pad_left;

                    if (im_row < 0 || im_col < 0 || (im_row >= height) || (im_col >= width * channels)) {
                        res[0][col_index + i * out_batch_size] = 0;
                        res[1][col_index + i * out_batch_size] = 0;

                    } else {
                        res[0][col_index + i * out_batch_size] = data_im[0][i * in_batch_size + im_col + width * channels * (im_row + height * c_im)];
                        res[1][col_index + i * out_batch_size] = data_im[1][i * in_batch_size + im_col + width * channels * (im_row + height * c_im)];
                    }
                }
            }
        }
    }
}

void pad_original(Lint **res, Lint **data_im, int channels, int height, int width, int ksize, int stride, uint batch_size) {

    int padw = (ksize / 2);

    int c, h, w;

    int height_col = height + 2;
    int width_col = width + 2;

    int pad_left = channels;
    int pad_top = 1;
    int channels_col = channels;

    vector<double> data_col((height_col) * (width_col)*channels, 0);

    uint in_batch_size = channels * height * width;
    uint out_batch_size = channels_col * height_col * width_col;

    for (uint i = 0; i < batch_size; i++) {

        for (c = 0; c < 1; ++c) {
            int w_offset = c % width_col;
            int h_offset = (c / height_col) % height_col;
            int c_im = c / width_col / width_col;

            // cout << "w_offset : " << w_offset << endl;
            // cout << "h_offset : " << h_offset << endl;
            // cout << "c_im : " << c_im << endl;

            for (h = 0; h < height_col; ++h) {
                for (w = 0; w < width_col * channels; ++w) {
                    int im_row = h_offset * stride + h;
                    int im_col = w_offset * stride * channels + w;
                    // int col_index = (c * height_col + h) * width_col + w;
                    int col_index = (c * height_col + h) * width_col * channels + w;
                    // cout << "c_im : " << c_im << "\t";

                    // cout << "im_row : " << im_row << "\t";
                    // cout << "im_col : " << im_col << "\t";
                    // cout << "col_index : " << col_index << "\t";
                    // cout << "result : " << pad_get_pixel(data_im, height, width, channels, im_row, im_col, c_im, pad_top, pad_left) << endl;

                    im_row -= pad_top;
                    im_col -= pad_left;

                    if (im_row < 0 || im_col < 0 || (im_row >= height) || (im_col >= width * channels)) {
                        res[0][col_index + i * out_batch_size] = 0;
                        res[1][col_index + i * out_batch_size] = 0;

                    } else {
                        res[0][col_index + i * out_batch_size] = data_im[0][i * in_batch_size + im_col + width * channels * (im_row + height * c_im)];
                        res[1][col_index + i * out_batch_size] = data_im[1][i * in_batch_size + im_col + width * channels * (im_row + height * c_im)];
                    }

                    // data_col[col_index] = pad_get_pixel(data_im, height, width, channels, im_row, im_col, c_im, pad_top, pad_left);
                }
            }
        }
    }

    // return data_col;
}

void MaxPool(Lint **res, Lint **a, uint c, uint m, uint n, uint batch_size, uint ring_size, int *map, NodeNetwork *nodeNet, uint flag) {
    if (flag == 0) {
        old_MaxPool(res, a, c, m, n, batch_size, ring_size, map, nodeNet);
    } else {
        eda_MaxPool(res, a, c, m, n, batch_size, ring_size, map, nodeNet);
    }
}
void old_MaxPool(Lint **res, Lint **a, uint c, uint m, uint n, uint batch_size, uint ring_size, int *map, NodeNetwork *nodeNet) {
    // takes in matrix of dimension m x n and pools according to the window
    // output matrix dimensions are calculated based off of inputs
    // we exclude the z direction in function since it is one in all cases
    uint i, j, k, l, index; // used for loops
    // uint bytes = ((ring_size) + 7) >> 3;

    // uint size = (c*m*n)/2;
    uint size = (c * m * n) / 2;
    // printf("size : %u\n", size);

    // Lint *res2 = new Lint [100000];
    // memset(res2,0,sizeof(Lint)*100000);

    Lint **x0 = new Lint *[2];
    Lint **x1 = new Lint *[2];
    Lint **temp = new Lint *[2];
    for (i = 0; i < 2; i++) {
        x0[i] = new Lint[size * batch_size];
        memset(x0[i], 0, sizeof(Lint) * (size * batch_size));
        x1[i] = new Lint[size * batch_size];
        memset(x1[i], 0, sizeof(Lint) * (size * batch_size));
        temp[i] = new Lint[size * batch_size];
        memset(temp[i], 0, sizeof(Lint) * (size * batch_size));
    }

    for (j = 0; j < batch_size; j++) {

        for (i = 0; i < size; i++) {
            x0[0][j * size + i] = a[0][j * size * 2 + 2 * i];
            x0[1][j * size + i] = a[1][j * size * 2 + 2 * i];
            x1[0][j * size + i] = a[0][j * size * 2 + 2 * i + 1];
            x1[1][j * size + i] = a[1][j * size * 2 + 2 * i + 1];
        }
    }

    Rss_LT(temp, x0, x1, size * batch_size, ring_size, map, nodeNet); // c is the problem

    for (i = 0; i < 2; i++) {
        memset(x0[i], 0, sizeof(Lint) * (size * batch_size));
        memset(x1[i], 0, sizeof(Lint) * (size * batch_size));
        // memset(res[i],0,sizeof(Lint)*(size));
    }
    size = size / 2;
    // printf("size : %u\n", size);

    // Rss_Open(res2, temp, 16*12*12, map, ring_size, nodeNet);
    // for(int i =0; i< 30; i++) {
    //     printf("res[%i]  : %llu\n", i, res2[i]);
    // }
    for (l = 0; l < batch_size; l++) {
        k = 0;
        for (i = 0; i < c * m / 2; i++) {
            for (j = 0; j < n / 2; j++) {

                // printf("j + i*n : %u\n", j + i*n);
                // printf("j + m/2 + i*n : %u\n", j + m/2 +i*n);
                x0[0][l * size + k] = temp[0][2 * l * size + j + i * n];
                x0[1][l * size + k] = temp[1][2 * l * size + j + i * n];
                x1[0][l * size + k] = temp[0][2 * l * size + j + m / 2 + i * n];
                x1[1][l * size + k] = temp[1][2 * l * size + j + m / 2 + i * n];

                k++;
            }
        }
    }

    Rss_LT(res, x0, x1, size * batch_size, ring_size, map, nodeNet);

    // delete [] res2;
    // cleanup
    for (i = 0; i < 2; i++) {
        delete[] x0[i];
        delete[] temp[i];
        delete[] x1[i];
    }
    delete[] x0;
    delete[] x1;
    delete[] temp;
}

void eda_MaxPool(Lint **res, Lint **a, uint c, uint m, uint n, uint batch_size, uint ring_size, int *map, NodeNetwork *nodeNet) {
    // takes in matrix of dimension m x n and pools according to the window
    // output matrix dimensions are calculated based off of inputs
    // we exclude the z direction in function since it is one in all cases
    uint i, j, k, l, index; // used for loops
    // uint bytes = ((ring_size) + 7) >> 3;

    // uint size = (c*m*n)/2;
    uint size = (c * m * n) / 2;
    // printf("size : %u\n", size);

    // Lint *res2 = new Lint [100000];
    // memset(res2,0,sizeof(Lint)*100000);

    Lint **x0 = new Lint *[2];
    Lint **x1 = new Lint *[2];
    Lint **temp = new Lint *[2];
    for (i = 0; i < 2; i++) {
        x0[i] = new Lint[size * batch_size];
        memset(x0[i], 0, sizeof(Lint) * (size * batch_size));
        x1[i] = new Lint[size * batch_size];
        memset(x1[i], 0, sizeof(Lint) * (size * batch_size));
        temp[i] = new Lint[size * batch_size];
        memset(temp[i], 0, sizeof(Lint) * (size * batch_size));
    }

    for (j = 0; j < batch_size; j++) {

        for (i = 0; i < size; i++) {
            x0[0][j * size + i] = a[0][j * size * 2 + 2 * i];
            x0[1][j * size + i] = a[1][j * size * 2 + 2 * i];
            x1[0][j * size + i] = a[0][j * size * 2 + 2 * i + 1];
            x1[1][j * size + i] = a[1][j * size * 2 + 2 * i + 1];
        }
    }

    new_Rss_LT(temp, x0, x1, size * batch_size, ring_size, map, nodeNet); // c is the problem

    for (i = 0; i < 2; i++) {
        memset(x0[i], 0, sizeof(Lint) * (size * batch_size));
        memset(x1[i], 0, sizeof(Lint) * (size * batch_size));
        // memset(res[i],0,sizeof(Lint)*(size));
    }
    size = size / 2;
    // printf("size : %u\n", size);

    // Rss_Open(res2, temp, 16*12*12, map, ring_size, nodeNet);
    // for(int i =0; i< 30; i++) {
    //     printf("res[%i]  : %llu\n", i, res2[i]);
    // }
    for (l = 0; l < batch_size; l++) {
        k = 0;
        for (i = 0; i < c * m / 2; i++) {
            for (j = 0; j < n / 2; j++) {

                // printf("j + i*n : %u\n", j + i*n);
                // printf("j + m/2 + i*n : %u\n", j + m/2 +i*n);
                x0[0][l * size + k] = temp[0][2 * l * size + j + i * n];
                x0[1][l * size + k] = temp[1][2 * l * size + j + i * n];
                x1[0][l * size + k] = temp[0][2 * l * size + j + m / 2 + i * n];
                x1[1][l * size + k] = temp[1][2 * l * size + j + m / 2 + i * n];

                k++;
            }
        }
    }

    new_Rss_LT(res, x0, x1, size * batch_size, ring_size, map, nodeNet);

    // delete [] res2;
    // cleanup
    for (i = 0; i < 2; i++) {
        delete[] x0[i];
        delete[] temp[i];
        delete[] x1[i];
    }
    delete[] x0;
    delete[] x1;
    delete[] temp;
}

void ReLU(Lint **res, Lint **a, uint size, uint ring_size, int *map, NodeNetwork *nodeNet, uint flag) {
    if (flag == 0) {
        old_ReLU(res, a, size, ring_size, map, nodeNet);
    } else {
        eda_ReLU(res, a, size, ring_size, map, nodeNet);
    }
}

// does not need to be changed for batch optimization
void old_ReLU(Lint **res, Lint **a, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {
    // takes in matrix of dimension m x n and calculates the ReLU for each element
    // output dimensions same as input
    uint i, j, k, index; // used for loops
    uint bytes = ((ring_size) + 7) >> 3;

    Lint **zero = new Lint *[2];
    for (i = 0; i < 2; i++) {
        zero[i] = new Lint[size];
        memset(zero[i], 0, sizeof(Lint) * (size));
    }
    Rss_LT(res, a, zero, size, ring_size, map, nodeNet);
    // cleanup
    for (i = 0; i < 2; i++) {
        delete[] zero[i];
    }
    delete[] zero;
}

// does not need to be changed for batch optimization
void eda_ReLU(Lint **res, Lint **a, uint size, uint ring_size, int *map, NodeNetwork *nodeNet) {
    // takes in matrix of dimension m x n and calculates the ReLU for each element
    // output dimensions same as input
    uint i, j, k, index; // used for loops
    uint bytes = ((ring_size) + 7) >> 3;

    Lint **zero = new Lint *[2];
    for (i = 0; i < 2; i++) {
        zero[i] = new Lint[size];
        memset(zero[i], 0, sizeof(Lint) * (size));
    }
    new_Rss_LT(res, a, zero, size, ring_size, map, nodeNet);
    // cleanup
    for (i = 0; i < 2; i++) {
        delete[] zero[i];
    }
    delete[] zero;
}

void add_biases(Lint **res, Lint **a, Lint **b, uint m, uint n, uint batch_size, int *map, NodeNetwork *nodeNet) {

    uint i, j, k;

    for (k = 0; k < batch_size; k++) {
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                res[0][k * m * n + i * n + j] = a[0][k * m * n + i * n + j] + b[0][i];
                res[1][k * m * n + i * n + j] = a[1][k * m * n + i * n + j] + b[1][i];
            }
        }
    }
}

// (K) Filter shape - (K x K x M x N)
// (F) Input size - (F x F x M)
// (G) Output size - (F x F x N
void conv(Lint **res, Lint **a, Lint **b, uint K, uint M, uint N, uint F, uint stride, uint batch_size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    // printf("F*F*K*K*M*N:\t %llu\n", F*F*K*K*M*N*batch_size);
    // need to divide by stride^2
    Rss_Mult(res, a, b, (F * F * K * K * M * N * batch_size) / (stride * stride), ring_size, map, nodeNet);

    ReLU(res, a, F * F * N * batch_size / (stride * stride), ring_size, map, nodeNet, 1);
}

void conv_dw(Lint **res, Lint **a, Lint **b, uint K, uint M, uint F, uint stride, uint batch_size, uint ring_size, int *map, NodeNetwork *nodeNet) {
    //  need to ensure that when generating data we have enough generated to satisfy a large batch size

    // printf("F*F*K*K*M:\t %llu\n", F*F*K*K*M*batch_size);
    // printf("F*F*M*M:\t %llu\n", F*F*M*M*batch_size);

    // 3 x 3 dw conv
    Rss_Mult(res, a, b, (F * F * K * K * M * batch_size) / (stride * stride), ring_size, map, nodeNet);

    ReLU(res, a, F * F * M * batch_size / (stride * stride), ring_size, map, nodeNet, 1);

    // 1 x 1 conv, N is the num_filters
    Rss_Mult(res, a, b, (F * F * M * M * batch_size) / (stride * stride), ring_size, map, nodeNet);

    ReLU(res, a, F * F * M * batch_size / (stride * stride), ring_size, map, nodeNet, 1);
}

void fire_module(Lint **res, Lint **a, Lint **b, uint M, uint F, uint squeeze, uint expand, uint stride, uint batch_size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    conv(res, a, b, 1, M, squeeze, F, stride, batch_size, ring_size, map, nodeNet);

    conv(res, a, b, 1, squeeze, expand, F, stride, batch_size, ring_size, map, nodeNet);

    conv(res, a, b, 3, squeeze, expand, F, stride, batch_size, ring_size, map, nodeNet);
}

void max_pool_bench(Lint **res, Lint **a, Lint **b, uint M, uint F, uint win, uint stride, uint batch_size, uint ring_size, int *map, NodeNetwork *nodeNet) {

    vector<uint> rounds;

    if (win == 2) {
        rounds.push_back(2);
        rounds.push_back(1);
    } else if (win == 3) {

        rounds.push_back(4);
        rounds.push_back(2);
        rounds.push_back(2);
        rounds.push_back(1);
    }
    for (size_t i = 0; i < rounds.size(); i++) {
        new_Rss_LT(res, a, b, (rounds.at(i) * M * M * F * batch_size) / (stride * stride), ring_size, map, nodeNet);
    }

    // this is assuming a window size of 3 x 3

    // new_Rss_LT(res, a, b, (4*M*M*F*batch_size)/(stride*stride), ring_size, map, nodeNet);

    // new_Rss_LT(res, a, b, (2*M*M*F*batch_size)/(stride*stride), ring_size, map, nodeNet);

    // new_Rss_LT(res, a, b, (2*M*M*F*batch_size)/(stride*stride), ring_size, map, nodeNet);

    // new_Rss_LT(res, a, b, (1*M*M*F*batch_size)/(stride*stride), ring_size, map, nodeNet);
}