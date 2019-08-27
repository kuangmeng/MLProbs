#include "deterministic_random.h"


int parkmiller(int seed) // 1<=seed<m
{
	int const a = 75; 
	int const m = RND_MAX + 1;
	seed = ((seed * a)) % m;

	return seed;
}