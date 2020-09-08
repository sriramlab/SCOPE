#ifndef PROPCA_HELPER_H_
#define PROPCA_HELPER_H_


#include "time.h"

#include <iostream>


extern struct timespec t0;

struct timespec elapsed() {
	struct timespec ts;
	clock_gettime(CLOCK_REALTIME, &ts);
	if (ts.tv_nsec < t0.tv_nsec) {
		ts.tv_nsec = 1000000000 + ts.tv_nsec - t0.tv_nsec;
		ts.tv_sec--;
	}
	ts.tv_sec -= t0.tv_sec;
	return (ts);
}

int timelog(const char* message) {
  struct timespec ts = elapsed();
  return (printf("[%06ld.%09ld] %s\n", ts.tv_sec, ts.tv_nsec, message));
}

void * malloc_double_align(size_t n, unsigned int a /*alignment*/, double * &output) {
    void *adres = NULL;
    void *adres2 = NULL;
    adres = malloc(n * sizeof(double) + a);
    size_t adr = (size_t) adres;
    size_t adr2 = adr + a - (adr & (a - 1u)); 	// a valid address for a alignment
    adres2 = reinterpret_cast<void *>(adr2);
    output = reinterpret_cast<double *>(adres2);
    return adres;                		// pointer to be used in free()
}

void print_timenl() {
	clock_t c = clock();
	double t = static_cast<double>(c) / CLOCKS_PER_SEC;
	std::cout << "Time = " << t << std::endl;
}

void print_time() {
	clock_t c = clock();
	double t = static_cast<double>(c) / CLOCKS_PER_SEC;
	std::cout << "Time = " << t  << " : ";
}

#endif  // PROPCA_HELPER_H_
