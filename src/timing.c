/** @license 2021 Neil Edelman, distributed under the terms of the
 [MIT License](https://opensource.org/licenses/MIT).

 Timings.

 @std C89/90 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include "orcish.h"


/** On-line numerically stable first-order statistics, <Welford, 1962, Note>. */
struct datum { size_t count; double mean, ssdm; };

static void datum_reset(struct datum *const d)
	{ assert(d); d->count = 0, d->mean = d->ssdm = 0; }

static void datum_add(struct datum *const d, const double replica) {
	const size_t n = ++d->count;
	const double delta = replica - d->mean;
	assert(d);
	d->mean += delta / n;
	d->ssdm += delta * (replica - d->mean);
}

static double datum_mean(const struct datum *const d)
	{ assert(d); return d->count ? d->mean : (double)NAN; }

static double datum_sample_variance(const struct datum *const d)
	{ assert(d); return d->count > 1 ? d->ssdm / (d->count - 1) : (double)NAN; }

static double datum_stddev(const struct datum *const d)
	{ return sqrt(datum_sample_variance(d)); }

/** Returns a time difference in microseconds from `then`. */
static double diff_us(clock_t then)
	{ return 1000000.0 / CLOCKS_PER_SEC * (clock() - then); }


/* The test data. */

static void word_to_str(const char *const *const w, char (*const a)[12])
	{ sprintf(*a, "%.11s", *w); }
static int word_cmp(const char *const *const a, const char *const *const b)
	{ return strcmp(*a, *b); }
static int word_cmp_void(const void *const a, const void *const b)
	{ return word_cmp(a, b); }


/* Load all the boxes. */

/* Trie trie-btree. */
#define TRIE_NAME word
#define TRIE_TO_STRING
#define TRIE_TEST
#include "../../trie-btree/src/trie.h"

/* Old versions require a different setup. */
#undef PT_
#undef T_
#undef TRIE_H

/* Fixed-size trie. */
static void fill_word(const char *const *a) { assert(a && (!a ^ 0)); }
#define TRIE_NAME fixed
#define TRIE_TO_STRING
#define TRIE_TEST &fill_word
#include "../../../trie/src/trie.h"

/* Array with `qsort` and `bsearch`. */
#define ARRAY_NAME word
#define ARRAY_TYPE const char *
#define ARRAY_EXPECT_TRAIT
#include "../../../array/src/array.h"
#define ARRAY_TO_STRING &word_to_str
#define ARRAY_EXPECT_TRAIT
#include "../../../array/src/array.h"
#define ARRAY_COMPARE &word_cmp
#include "../../../array/src/array.h"
/** Fills `strs` with `words` of size `words_size` from `words_start` to
 `words_chosen`, wrapping around. */
static int array_fill(struct word_array *const array,
	const char *const*const words, const size_t words_size,
	const size_t words_start, const size_t words_chosen) {
	assert(array && words && words_chosen && words_chosen <= words_size
		&& words_start < words_size);
	word_array_clear(array);
	if(!word_array_append(array, words_chosen)) return 0;
	if(words_start + words_chosen > words_size) {
		const size_t size_a = words_size - words_start,
		size_b = words_chosen - size_a;
		memcpy(array->data, words + words_start, sizeof *words * size_a);
		memcpy(array->data + size_a, words, sizeof *words * size_b);
	} else {
		memcpy(array->data, words + words_start, sizeof *words * words_chosen);
	}
	return 1;
}

/** Perform a 32 bit
 [Fowler/Noll/Vo FNV-1a](http://www.isthe.com/chongo/tech/comp/fnv/) hash on
 `str`. */
static unsigned fnv_32a_str(const char *const str) {
	const unsigned char *s = (const unsigned char *)str;
	/* 32 bit FNV-1 and FNV-1a non-zero initial basis, `FNV1_32A_INIT`. */
	unsigned hval = 0x811c9dc5;
	/* FNV magic prime `FNV_32_PRIME 0x01000193`. */
	while(*s) {
		hval ^= *s++;
		hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
	}
	return hval;
}
static int word_is_equal(const char *const a, const char *const b) {
	return !strcmp(a, b);
}
static void word_to_string(const char *const*const ps, char (*const a)[12]) {
	strncpy(*a, *ps, sizeof(*a) - 1);
	(*a)[sizeof(*a) - 1] = '\0';
}
#define SET_NAME word
#define SET_TYPE const char *
#define SET_HASH &fnv_32a_str
#define SET_IS_EQUAL &word_is_equal
#define SET_EXPECT_TRAIT
#include "../../../set/src/set.h"
#define SET_TO_STRING &word_to_string
#include "../../../set/src/set.h"
#define POOL_NAME word_node
#define POOL_TYPE struct word_set_node
#include "../../../pool/src/pool.h"


/* How many experiments is an X-macro. `gnuplot` doesn't like `_`. */
#define TESTS X(ARRAY), X(HASH), X(FIXEDTRIE), X(BTRIE)

static int timing_comparison(const char *const *const words,
	const size_t words_size) {
	size_t i, r, n, x, replicas = 5, wtf;
	clock_t t, t_total;
	int success = 0, is_full = 0;
	const char *const gnu_name = "experiment";
	const int esystem = -300;

	/* The containers. */
	struct word_array array = ARRAY_IDLE;
	struct word_set set = SET_IDLE;
	struct word_node_pool word_pool = POOL_IDLE;
	struct fixed_trie ftrie = TRIE_IDLE;
	struct word_trie btrie = TRIE_IDLE;

	/* How many files we open simultaneously. */
	enum { INIT, LOOK } act;
	const char *const act_str[] = { "init", "look" };
#define X(name) name
	enum { TESTS };
#undef X
#define X(name) { { #name, 0, { 0, 0, 0 } }, { #name, 0, { 0, 0, 0 } } }
	struct { const char *name; FILE *fp; struct datum d; }
		tests[][2] = { TESTS };
#undef X
	const size_t tests_size = sizeof tests / sizeof *tests;

	assert(words && words_size);

	/* Open all graphs for writing. */
	for(x = 0; x < tests_size; x++) {
		for(act = INIT; act <= LOOK; act++) {
			char fn[64];
			if(sprintf(fn, "graph/%s-%s.tsv",
				tests[x][act].name, act_str[act]) < 0
				|| !(tests[x][act].fp = fopen(fn, "w"))) goto catch;
			fprintf(tests[x][act].fp, "# %s\n"
				"# <items>\t<t (ms)>\t<sample error on t with %lu replicas>\n",
				tests[x][act].name, (unsigned long)replicas);
		}
	}
	for(n = 1; !is_full; n <<= 1) {
		if(n >= words_size) is_full = 1, n = words_size;
		for(x = 0; x < tests_size; x++)
			for(act = INIT; act <= LOOK; act++)
				datum_reset(&tests[x][act].d);
		for(r = 0; r < replicas; r++) {
			size_t start_i = (unsigned)rand() / (RAND_MAX / words_size + 1);
			printf("Replica %lu/%lu.\n",
				(unsigned long)(r + 1), (unsigned long)replicas);
			t_total = clock();

			/* Sorted array; pre-allocate for fair test. */
			array_fill(&array, words, words_size, start_i, n);
			t = clock();
			qsort(array.data, array.size, sizeof array.data,
				&word_cmp_void);
			word_array_unique(&array);
			datum_add(&tests[ARRAY][INIT].d, diff_us(t));
			printf("Added init array size %lu: %s.\n",
				(unsigned long)array.size, word_array_to_string(&array));
			t = clock();
			for(i = 0; i < n; i++) {
				const char *const word = words[(start_i + i) % words_size],
					**const key = bsearch(&word, array.data, array.size,
					sizeof array.data, &word_cmp_void);
				int cmp;
				assert(key);
				cmp = strcmp(word, *key);
				assert(!cmp);
				(void)cmp;
			}
			datum_add(&tests[ARRAY][LOOK].d, diff_us(t));
			printf("Added look array size %lu.\n", (unsigned long)array.size);

			/* Hash map. */
			word_set_clear(&set);
			word_node_pool_clear(&word_pool);
			t = clock();
			for(i = 0; i < n; i++) {
				struct word_set_node *elem = word_node_pool_new(&word_pool);
				elem->key = words[(start_i + i) % words_size];
				if(word_set_policy_put(&set, elem, 0))
					word_node_pool_remove(&word_pool, elem);
			}
			datum_add(&tests[HASH][INIT].d, diff_us(t));
			printf("Added init hash size %lu: %s.\n",
				(unsigned long)set.size, word_set_to_string(&set));
			t = clock();
			for(i = 0; i < n; i++) {
				const char *const word = words[(start_i + i) % words_size];
				const struct word_set_node *const elem
					= word_set_get(&set, word);
				int cmp;
				assert(elem);
				cmp = strcmp(word, elem->key);
				assert(!cmp);
			}
			datum_add(&tests[HASH][LOOK].d, diff_us(t));
			printf("Added look hash size %lu.\n",
				(unsigned long)set.size);

			/* Fixed trie. */
			t = clock();
			array_fill(&array, words, words_size, start_i, n);
			fixed_trie_clear(&ftrie);
			/*for(i = 0; i < n; i++)
				str_trie_add(&trie, array.data[i]); <- this is very slow! */
			fixed_trie_from_array(&ftrie, array.data, array.size);
			datum_add(&tests[FIXEDTRIE][INIT].d, diff_us(t));
			printf("Added init fixed trie size %lu: %s.\n",
				(unsigned long)fixed_trie_size(&ftrie),
				fixed_trie_to_string(&ftrie));
			t = clock();
			for(i = 0; i < n; i++) {
				const char *const word = array.data[i],
					*const key = fixed_trie_get(&ftrie, word);
				const int cmp = strcmp(word, key);
				(void)cmp, assert(key && !cmp);
			}
			datum_add(&tests[FIXEDTRIE][LOOK].d, diff_us(t));
			printf("Added look fixed trie size %lu.\n",
				(unsigned long)fixed_trie_size(&ftrie));

			/* B-Tree trie. */
			wtf = 0;
			array_fill(&array, words, words_size, start_i, n);
			t = clock();
			for(i = 0; i < n; i++) word_trie_add(&btrie, array.data[i]);
			datum_add(&tests[BTRIE][INIT].d, diff_us(t));
			{
				struct word_trie_iterator it;
				size_t size;
				word_trie_prefix(&btrie, "", &it);
				size = word_trie_size(&it);
				printf("Added init new B-trie, size %lu.\n", size);
				if(size < 10000)
					trie_word_no++, trie_word_graph(&btrie, "graph/trie-no.gv");
			}
			t = clock();
			for(i = 0; i < n; i++) {
				const char *const word = array.data[i],
					*const key = word_trie_get(&btrie, word);
				int cmp;
				if(!key) { wtf++; continue; }
				assert(key), cmp = strcmp(word, key);
				assert(!cmp);
			}
			datum_add(&tests[BTRIE][LOOK].d, diff_us(t));
			word_trie_(&btrie);
			printf("Added look B-trie.\n");
			if(wtf) printf("%lu UNEXPLAINED OCCURANCES!\n", wtf);

			/* Took took much time; decrease the replicas for next time. */
			if(replicas != 1
				&& 10.0 * (clock() - t_total) / CLOCKS_PER_SEC > 1.0 * replicas)
				replicas--;
		}
		for(x = 0; x < tests_size; x++) {
			for(act = INIT; act <= LOOK; act++) {
				double stddev = datum_stddev(&tests[x][act].d);
				if(stddev != stddev) stddev = 0; /* Is nan; happens. */
				fprintf(tests[x][act].fp, "%lu\t%f\t%f\n",
					(unsigned long)n, datum_mean(&tests[x][act].d), stddev);
			}
		}
		/* Crashes?
		 if(n == 512) trie_fixed_graph(&trie, "graph/example-thing.gv");*/
	}
	printf("Test passed.\n");
	{ success = 1; goto finally; }
catch:
	perror("graph");
	printf("Test failed.\n");
finally:
	word_array_(&array);
	word_set_(&set);
	word_node_pool_(&word_pool);
	fixed_trie_(&ftrie);
	word_trie_(&btrie);
	for(x = 0; x < tests_size; x++)
		for(act = INIT; act <= LOOK; act++)
		if(tests[x][act].fp && fclose(tests[x][act].fp))
		perror(tests[x][act].name);
	if(!success) return 0;

	/* Output a `gnuplot` script. */
	{
		char fn[64];
		for(act = INIT; act <= LOOK; act++) {
			FILE *gnu_fp = 0;
			if(sprintf(fn, "graph/%s-%s.gnu", gnu_name, act_str[act]) < 0
				|| !(gnu_fp = fopen(fn, "w"))) goto catch2;
			fprintf(gnu_fp,
				/*"set for [i=1:8] linetype i dashtype i\n"*/
				"set style line 1 lt 5 lw 3 lc rgb '#0072bd'\n"
				"set style line 2 lt 5 lw 3 lc rgb '#ff0000'\n" /* a2142f */
				"set style line 3 lt 5 lw 3 lc rgb '#00ac33'\n" /* 30ac77 */
				"set style line 4 lt 5 lw 3 lc rgb '#19d3f5'\n");
			fprintf(gnu_fp, "set term postscript eps enhanced color\n"
				/*"set encoding utf8\n" Doesn't work at all; {/Symbol m}. */
				"set output \"graph/%s-%s.eps\"\n"
				"set grid\n"
				"set xlabel \"elements\"\n"
				"set ylabel \"time per element, t (ns)\"\n"
				"set yrange [0:2000]\n"
				"set log x\n"
				"plot", gnu_name, act_str[act]);
			for(x = 0; x < tests_size; x++) fprintf(gnu_fp,
				"%s \\\n\"graph/%s-%s.tsv\" using 1:($2/$1*1000):($3/$1*1000) "
				"with errorlines title \"%s\" ls %d", x ? "," : "",
				tests[x][act].name, act_str[act], tests[x][act].name,
				(int)x + 1);
			fprintf(gnu_fp, "\n");
			if(fclose(gnu_fp)) goto catch2;
		}
	}
	{
		int result;
		char cmd[64];
		fprintf(stderr, "Running Gnuplot to get a graph "
			"(http://www.gnuplot.info/.)\n");
		if((result = system("/usr/local/bin/gnuplot --version")) == -1)
			goto catch2;
		else if(result != EXIT_SUCCESS) { errno = esystem; goto catch2; }
		for(act = INIT; act <= LOOK; act++) {
			if(sprintf(cmd, "/usr/local/bin/gnuplot graph/%s-%s.gnu",
				gnu_name, act_str[act]) < 0
				|| (result = system(cmd)) == -1) goto catch2;
			else if(result != EXIT_SUCCESS) { errno = esystem; goto catch2; }
		}
	}
	goto finally2;
catch2:
	errno == esystem
		? fprintf(stderr, "gnuplot: didn't do it.\n")
		: (perror(gnu_name), 0);
finally2:
	return 1;
}

int main(void) {
#if 0
	/* This is a dictionary defined in `parole_inglesi.c`. */
	extern const char *const parole[];
	extern const size_t parole_size;
	fprintf(stderr, "parole_size %lu\n", (unsigned long)parole_size);
	return timing_comparison(parole, parole_size) ? EXIT_SUCCESS : EXIT_FAILURE;
#else
	/* Generated data. */
	const char **words = 0;
	const size_t words_size = 1000000;
	struct orc { char name[64]; } *orcs = 0;
	size_t i;
	int success = EXIT_FAILURE;

	if(!(words = malloc(sizeof *words * words_size))
		|| !(orcs = malloc(sizeof *orcs * words_size))) goto catch;
	for(i = 0; i < words_size; i++)
		orcish(orcs[i].name, sizeof orcs[i].name),
		/*printf("%s\n", orcs[i].name), */words[i] = orcs[i].name;
	return timing_comparison(words, words_size) ? EXIT_SUCCESS : EXIT_FAILURE;
	{ success = EXIT_SUCCESS; goto finally; }
catch:
	perror("orcs");
finally:
	free(words);
	free(orcs);
	return success;
#endif
}
