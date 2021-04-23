#include "../random.h"

void safe_rand(double* rand_num, int a, int b, int c)
{
	// get initial random value
    a += (int)9471948645;
    b -= (int)55187636245;
    c += (int)717952745;
    int d = (int)10598377245;

    int val = (a * rand()) ^ (b * c + d) ^ (a + (b << 5) + (c << 10) + d);

    if (val < 0) val *= -1;

    *rand_num = (double)val / INT_MAX;


}