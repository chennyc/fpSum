#include "../include/Rss_Op.h"
// #include "../Rss_Op.h"
// extern "C"{
// #include "../aes_ni.h"
// }

void Rss_RandBit(Lint** b, uint size, uint ring_size, int *map, NodeNetwork* nodeNet){

    int pid = nodeNet->getID();
    uint i;
    uint bytes = (ring_size+9) >> 3;
    // printf("bytes : %llu\n", bytes );

    Lint **u = new Lint *[2];
    Lint **a = new Lint *[2];
    Lint **d = new Lint *[2];

    for(i = 0; i < 2; i++){
        u[i] = new Lint [size];
        a[i] = new Lint [size];
        d[i] = new Lint [size];
    }
    Lint *e = new Lint [size];
    Lint *c = new Lint [size];
    uint8_t *buffer = new uint8_t [bytes*size];

    // used to make a odd, we only add 1 to one share of a
    // All shares will be doubled
    Lint a1, a2;
    switch(pid){
        case 1:
            a1 = 1;
            a2 = 0;
            break;
        case 2:
            a1 = 0;
            a2 = 0;
            break;
        case 3:
            a1 = 0;
            a2 = 1;
            break;
    }


    nodeNet->prg_getrandom(0, bytes, size, buffer);
    for (i = 0 ; i<size; i++) {
        memcpy(u[0]+i, buffer + i*bytes, bytes);
    }
    nodeNet->prg_getrandom(1, bytes, size, buffer);
    for (i = 0 ; i<size; i++) {
        memcpy(u[1]+i, buffer + i*bytes, bytes);
    }

    for(i = 0 ; i<size; i++){
        // ensuring [a] is odd
        a[0][i] = (u[0][i] << Lint(1)) + a1;
        a[1][i] = (u[1][i] << Lint(1)) + a2;

    }
    // squaring a
    Rss_MultPub(e, a, a, size, map, ring_size+2, nodeNet); //ringsize+2
    rss_sqrt_inv(c, e, size, ring_size+2);

    // effectively combines the two loops into one, eliminates d variable
    for (i = 0; i < size; i++) {
        b[0][i] = (c[i]*a[0][i] + a1) >> (1);
        b[1][i] = (c[i]*a[1][i] + a2) >> (1);

    }

    // freeing up
    delete [] c;
    delete [] buffer;
    delete [] e;
    for(i = 0; i < 2; i++){
        delete [] d[i];
        delete [] a[i];
        delete [] u[i];
    }
    delete [] d;
    delete [] a;
    delete [] u;
}

void Rss_RandBit3(Lint** b, uint size, uint ring_size, int *map, NodeNetwork* nodeNet){
// This randbit is implemented by using local randomness to generate shares of r_1 over Z_2 and converting them to the larger ring using B2A.
    int pid = nodeNet->getID();
    uint i;
    uint bytes = (size+7) >> 3;
    // printf("bytes : %llu\n", bytes );

    Lint **a = new Lint *[2];

    for(i = 0; i < 2; i++){
        a[i] = new Lint [size];
    }
    uint8_t *buffer0 = new uint8_t [bytes];
    uint8_t *buffer1 = new uint8_t [bytes];

    // used to make a odd, we only add 1 to one share of a
    // All shares will be doubled
    Lint a1, a2;
    switch(pid){
        case 1:
            a1 = 1;
            a2 = 0;
            break;
        case 2:
            a1 = 0;
            a2 = 0;
            break;
        case 3:
            a1 = 0;
            a2 = 1;
            break;
    }


    nodeNet->prg_getrandom(0, bytes, 1, buffer0);
    nodeNet->prg_getrandom(1, bytes, 1, buffer1);
    for (i = 0 ; i< size; i++) {
        a[0][i] = (buffer0[i/8] >> (i % 8)) & 1;
        a[1][i] = (buffer1[i/8] >> (i % 8)) & 1;
    }

    // b2a
    Rss_b2a3(b, a, ring_size, size, map, nodeNet);

    // freeing up
    delete [] buffer0;
    delete [] buffer1;
    for(i = 0; i < 2; i++){
        delete [] a[i];
    }
    delete [] a;
}

void rss_sqrt_inv(Lint *c, Lint *e, uint size, uint ring_size) {

  Lint c1, c2, temp, d_;
  uint i, j;

  for (i = 0 ; i < size; i++){
    c1 = Lint(1);
    c2 = Lint(1);
    d_ = Lint(4); // 100 - the first mask

    for (j = 2; j < ring_size - 1; j++) {
        temp = e[i] - (c1)*(c1);
        if (temp != Lint(0)) {
            //get the jth+1 bit of temp, place it in jth position, and add to c1
            c1 += (temp & (d_ << Lint(1))) >> Lint(1);
        }

        temp = Lint(1) - c1*c2;
        // get the jth bit of temp and add it to c2
        c2 += temp & d_;
        d_ = d_ << Lint(1);
    }
    // last round for the inv portion
    temp = Lint(1) - c1*c2;
    c[i] = c2 + (temp & d_);

    }
}

// not used

void invert(Lint *c, Lint *a, int size, int ring_size){

    // note, should really check if a is odd
    // if even abort

    Lint res;
    Lint temp;

    for (int i = 0; i < size; i++) {
        res = 1;
        // can start with j  = 1, since a is odd
        for (int j = 1; j < ring_size; j++) {
            // res += {Z2<K>((Z2<K>(1) - Z2<K>::Mul(*this, res)).get_bit(j)) << j};

            temp = Lint(1) - a[i]*res;
            // printf("%d th temp: %u \n", j, temp);


            // https://codeforwin.org/2016/01/c-program-to-get-value-of-nth-bit-of-number.html
            // b = temp.get_bit(j)
            temp = (temp >> j) & Lint(1); // getting the jth bit

            res += (temp << j);
        }
        c[i] = res;
        // printf("%d th c: %u \n", i, c[i]);
    }

}



// not used
void rss_sqrt(Lint *c, Lint *e, int size, int ring_size){

    // calculates the smallest root of some

    Lint d_;
    Lint temp;

    for(int i = 0 ; i <size; i++){


        // make sure e[i] is odd square
        // assert(e[i] % Lint(8) == Lint(1));
        //printf("e[%i] %u\n", i, e[i]);

        // the least significant bit of the root must be 1
        c[i] = 1;
        d_ = 0;

        // starting from 2 since the first two bits will be the same for all square roots
        for (int j = 2; j < ring_size-1; j++) {


            // printf("%i: c[%i]: %llu \n", i, j, c[i]);


            temp = e[i] - (c[i])*(c[i]);

            // printf("%d th temp: %u \n", j, temp);
            if(temp == Lint(0))
            {
                //printf("stop");
                break;
            }
            // d will be the j-th (j>0) least significant bit of the root
            d_ = (temp >> (j+1) ) & 1; //getting the jth+1 bit
            //printf("%d th d: %u \n", j, d_);
            // left shift d for j bits and add it to res;
            c[i] += (d_ << j);
            //printf("%d th c: %u \n", j, c[i]);
        }
        //printf("\nc[%i] %u\n", i, c[i]);
        //printf("c[%i]square is  %u\n", i, c[i]*c[i]); // c^2 should equal to e
    //  printf("bits of c[%i] %i\n", i, countBits(c[i]));
    //  res[i] = c[i];
        // printf("%d th c: %u \n", i, c[i]);

    }


}
