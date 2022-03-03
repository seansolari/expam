// kmers.c

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

static void shift_kmer(uint64_t *kmer, const unsigned char base, int k) {
	uint64_t	next_bits, mask;
	int			chunks, next_chunk;
	int			i;

    chunks = k / 32 + !! (k % 32 != 0);
	mask = ~(3ULL);

	for (i=0; i < chunks; ++i) {
		kmer[i] <<= 2;

		next_chunk = (i == chunks - 2) ? k % 32 : 32;
		if (i < chunks - 1)
			// Hard-coded for 64-bit chunks.
			next_bits = kmer[i+1] >> (2*(next_chunk-1));

		else {
			switch(base) {
				case 'A':
					next_bits = 0ULL;
					break;

				case 'C':
					next_bits = 1ULL;
					break;

				case 'G':
					next_bits = 2ULL;
					break;

				case 'T':
					next_bits = 3ULL;
					break;

				default:
					next_bits = 0ULL;
					break;
			}

			/* Zero out non-2k bits */
			if (k % 32 != 0)
				mask &= ~(3ULL << (2*(k % 32)));
		}

		kmer[i]	= (kmer[i] & mask) | next_bits;
	}
}

static void update_kmer(uint64_t *kmer, const unsigned char *sequence, uint64_t index, int k) {
	uint64_t	kmer_value, next_bit;
	int			chunks, current_chunk;
	int			i, j;

    chunks = k / 32 + !! (k % 32 != 0);

	for (i=0; i < chunks; ++i){
		kmer_value = 0;

		if (i == chunks - 1)
			current_chunk = k % 32;
		else
			current_chunk = 32;

		for (j=0; j < current_chunk; ++j) {
			switch(sequence[index + i*32 + j]) {
				case 'A':
					next_bit = 0ULL;
					break;

				case 'C':
					next_bit = 1ULL << (2 * (current_chunk-(j+1)) );
					break;

				case 'G':
					next_bit = 2ULL << (2 * (current_chunk-(j+1)) );
					break;

				case 'T':
					next_bit = 3ULL << (2 * (current_chunk-(j+1)) );
					break;

				default:
					next_bit = 0ULL;
					break;
			}

			kmer_value |= next_bit;
		}

		kmer[i] = kmer_value;
	}
}

static inline uint64_t* init_kmer(const unsigned char *sequence, uint64_t index, int k) {
	int			chunks;
	uint64_t	*kmer;

	chunks = k / 32 + !!(k % 32 != 0);
	kmer = (uint64_t *)malloc(sizeof(uint64_t) * chunks);

	update_kmer(kmer, sequence, index, k);

	return kmer;
}

static uint64_t find_next_kmer(const unsigned char *sequence, uint64_t index, int k) {
	bool			can_start;
	uint64_t		seq_length, max_index;
    unsigned char	l;

    can_start = true;  // We are moving onto the next index.
    seq_length = max_index = strlen(sequence);
    max_index -= k - 1;

    while (index <= max_index)
    {
        for (int i=0; i<k; ++i)
        {
            l = sequence[index + i];
            if ((l!='A') && (l!='C') && (l!='G') && (l!='T'))
            {
                can_start = false;
                break;
            }
        }

        if (can_start)
        	break;
        else {
            can_start = true;
            ++index;
        }
    }

    if (index < max_index)
    	return index;
    return seq_length;
};

static uint64_t raw_kmers(const unsigned char *sequence, int k, uint64_t *arr) {
    uint64_t        max_index, sequence_length;
    uint64_t        i, j, x;
    uint64_t        *kmer;
    unsigned char   base;
    uint64_t        chunks;

    max_index = sequence_length = strlen(sequence);
    max_index -= k-1;

    i = x = 0;
	chunks = (uint64_t ) (k / 32) + !! (k % 32 != 0);

    kmer = (uint64_t *)malloc(chunks * sizeof(uint64_t));

    if (!kmer)
        return 0;

    /* Initialise kmer */
	i = find_next_kmer(sequence, i, k);

	if (i != sequence_length) {
		kmer = init_kmer(sequence, i, k);

		for (j=0; j<chunks; j++)
			*(arr + (j + 1)) = kmer[j];

		*(arr + 0) = i;

		++x, ++i;
	}

	/* Filling loop */
	while (i < max_index)
	{
		/* Find next kmer */
		base = sequence[i + k - 1];

		if ((base!='A') && (base!='C') && (base!='G') && (base!='T'))
		{
			i = find_next_kmer(sequence, i, k);

			if (i != sequence_length)
				update_kmer(kmer, sequence, i, k);
			else
				// This will quit filling loop.
				continue;
		}
		else
			shift_kmer(kmer, base, k);

		/* Write this kmer to array and its corresponding index */
		for (j=0; j<chunks; ++j)
			*(arr + x*(chunks + 1) + (j + 1)) = *(kmer + j);

		*(arr + x*(chunks + 1)) = i;

		++x, ++i;
	}

	free(kmer);

	return x;
}

static void fill_kmers(const unsigned char *sequence, int k, uint64_t n, uint64_t *arr) {
    uint64_t    *kmer;
    uint64_t	i, j, chunks;

    i = j = 0;
    chunks = k / 32 + !! (k % 32 != 0);

    kmer = (uint64_t *)malloc(chunks * sizeof(uint64_t));

    if (!kmer)
        return;

    for (i=0; i<n; ++i) {
        update_kmer(kmer, sequence, *(arr + i*(chunks + 1)), k);

        for (j=0; j<chunks; ++j)
            *(arr + i*(chunks + 1) + (j + 1)) = *(kmer + j);
    }

    free(kmer);
}
