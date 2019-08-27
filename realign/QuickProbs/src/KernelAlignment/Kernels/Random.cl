#define RND_MAX 65536
#define RND_MAX_INV 0.000015298473212373405134167610072515f

int parkmiller(int seed) // 1<=seed<m
{
	int const a = 75; 
	int const m = RND_MAX + 1;
	seed = ((seed * a)) % m;

	return seed;
}