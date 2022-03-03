/*

    This file is part of Jellyfish.
    This work is dual-licensed under 3-Clause BSD License.

    Original code from https://github.com/gmarcais/Jellyfish.
    Modified on Friday 17th September, 2021.

*/

#include <stdlib.h>
#include <stdint.h>


static uint64_t word_reverse_complement(uint64_t w, int k) {
    w = ((w >> 2 )  & 0x3333333333333333) | ((w & 0x3333333333333333) << 2 );
    w = ((w >> 4 )  & 0x0F0F0F0F0F0F0F0F) | ((w & 0x0F0F0F0F0F0F0F0F) << 4 );
    w = ((w >> 8 )  & 0x00FF00FF00FF00FF) | ((w & 0x00FF00FF00FF00FF) << 8 );
    w = ((w >> 16)  & 0x0000FFFF0000FFFF) | ((w & 0x0000FFFF0000FFFF) << 16);
    w = ( w >> 32                       ) | ( w                       << 32);

    w = ~w;                 // Take reverse complement.
    w >>= 2 * (32 - k);     // Account for k-values less than 32.
    return w;
}

