#ifndef _UNISTD_H
#define _UNISTD_H    1

#include <intrin.h>
#include <io.h>

#define F_OK    0       /* Test for existence.  */

#define access _access

inline int __builtin_ffsll(long long mask)
{
	unsigned long index;
	_BitScanForward64(&index, mask);
	return index+1;
}
#endif