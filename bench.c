
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "cvec.h"
#include "bench.h"

static inline
char *make_random_seq(
	int64_t len)
{
	char *buf = malloc(len + 1);

	char const *table = "ACGTURYSWKMBDHVN";
	for(int64_t i = 0; i < len; i++) {
		buf[i] = table[rand() % strlen(table)];
	}
	return(buf);
}

#define _b(x)	( (x) & 0x1f )
enum bases {
	A = 0x01, C = 0x02, G = 0x04, T = 0x08
};
char const table[] __attribute__(( aligned(16) )) = {
	[_b('A')] = A,
	[_b('C')] = C,
	[_b('G')] = G,
	[_b('T')] = T,
	[_b('U')] = T,
	[_b('R')] = A | G,
	[_b('Y')] = C | T,
	[_b('S')] = G | C,
	[_b('W')] = A | T,
	[_b('K')] = G | T,
	[_b('M')] = A | C,
	[_b('B')] = C | G | T,
	[_b('D')] = A | G | T,
	[_b('H')] = A | C | T,
	[_b('V')] = A | C | G,
	[_b('N')] = 0,		/* treat 'N' as a gap */
	[_b('_')] = 0		/* sentinel */
};

static inline
char conv_ascii_to_4bit(
	char c)
{
	return(table[_b(c)]);
}

static inline
cvec_t conv_ascii_to_4bit_vec(
	cvec_t cv)
{
	return(_conv_5t_cvec(_conv_a5_cvec(cv), table));
}

#undef _b

int main(void)
{
	int64_t len = 128 * 1024 * 1024;

	printf("nucleotide sequence conversion benchmark: len(%lld)\n", len);

	char *src = make_random_seq(len);
	char *sdst = malloc(len);
	char *vdst = malloc(len);

	/* scalar conversion (with conv table) */
	bench_t s;
	bench_init(s);
	bench_start(s);
	for(int64_t i = 0; i < len; i++) {
		sdst[i] = conv_ascii_to_4bit(src[i]);
	}
	bench_end(s);

	/* vector conversion */
	bench_t v;
	bench_init(v);
	bench_start(v);
	for(int64_t i = 0; i < len; i += sizeof(cvec_t)) {
		_store_cvec(&vdst[i], conv_ascii_to_4bit_vec(_load_cvec(&src[i])));
	}
	bench_end(v);

	/* result check */
	if(memcmp(sdst, vdst, len) != 0) {
		printf("results differ.\n");
	}

	printf("scalar:\t%lld\nvector:\t%lld\n", bench_get(s), bench_get(v));

	free(src);
	free(sdst);
	free(vdst);
	return(0);
}
