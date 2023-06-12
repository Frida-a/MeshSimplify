#ifndef __TIMING_H__
#define __TIMING_H__

#include <string.h>
#include <sys/timeb.h>

#ifndef WIN32
#include <sys/time.h>
#endif // WIN32

inline double Time(void) {
#ifdef WIN32
	struct _timeb t;
	_ftime(&t);
	return double(t.time) + double(t.millitm) / 1000.0;
#else // WIN32
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec + double(t.tv_usec) / 1000000;
#endif // WIN32
}

#endif
