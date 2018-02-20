#ifndef ALGEBRA_H_INCLUDED
#define ALGEBRA_H_INCLUDED
#include <cassert>
#include <stdexcept>
#include <vector>
#include <array>
#include <ctgmath>

template <size_t N>
void dscal (std::array<double, N>& a, double b)
{		
	for(auto i = 0; i < N; i++) a[i] *= b;
}

template <size_t N>
double norm (const std::array<double, N>& a)
{		
	double norm = 0.;
	for(auto i = 0; i < N; i++)	norm += a[i]*a[i];
	return sqrt(norm);
}

template <size_t N>
double norm2 (const std::array<double, N>& a)
{		
	double norm2 = 0.;
	for(auto i = 0; i < N; i++)	norm2 += a[i]*a[i];
	return norm2;
}

template <size_t N>
double dist2 (const std::array<double, N>& a, const std::array<double, N>& b)
{		
	double dist2 = 0.;
	for(auto i = 0; i < N; i++)	dist2 += pow(a[i] - b[i], 2.);
	return dist2;
}

template <size_t N>
double dist (const std::array<double, N>& a, const std::array<double, N>& b)
{		
	double dist = 0.;
	for(auto i = 0; i < N; i++)	dist += pow(a[i] - b[i], 2.);
	return sqrt(dist);
}

template <size_t N>
std::array<double, N> mid_v (std::array<double, N> a, std::array<double, N> b)
{
	std::array<double, N> mid_v;
	for(auto i = 0; i < N; i++)	mid_v[i] = (a[i] + b[i])/2.;
	return mid_v;
}

#endif
