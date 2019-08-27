#pragma once

#define _sqr(a) (a)*(a)

#define _max(a,b)		((a) >= (b) ? (a) : (b))
#define _min(a,b)		((a) <= (b) ? (a) : (b))
#define _min3(a,b,c)	_min(_min(a,b),c)
#define _max3(a,b,c)	_max(_max(a,b),c)
#define _min4(a,b,c,d)	_min(_min(a,b),_min(c,d))
#define _max4(a,b,c,d)	_max(_max(a,b),_max(c,d))

// divides a by b and gets ceil of result
#define _ceildiv(a,b) (((a) + (b) - 1) / (b))

// perform ceil rounding to the closest multiplicity of b
#define _ceilround(a,b) (((a) + (b) - 1) / (b) * (b))
