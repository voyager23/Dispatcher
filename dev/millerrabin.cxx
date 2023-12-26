/*
 * MillerRabin.hxx
 * 
 * Copyright 2023 mike <mike@fedora38-1.home>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef __MILLERRABIN__
#define __MILLERRABIN__
	#include <iostream>
	#include <cstdint>

	using u64 = uint64_t;
	using u128 = __uint128_t;
	using namespace std;

	u64 binpower(u64 base, u64 e, u64 mod) {
		u64 result = 1;
		base %= mod;
		while (e) {
			if (e & 1)
				result = (u128)result * base % mod;
			base = (u128)base * base % mod;
			e >>= 1;
		}
		return result;
	}

	bool check_composite(u64 n, u64 a, u64 d, int s) {
		u64 x = binpower(a, d, n);
		if (x == 1 || x == n - 1)
			return false;
		for (int r = 1; r < s; r++) {
			x = (u128)x * x % n;
			if (x == n - 1)
				return false;
		}
		return true;
	};

	bool MillerRabin(uint64_t n) { // returns true if n is prime, else returns false.
		if (n < 2)
			return false;

		int r = 0;
		uint64_t d = n - 1;
		while ((d & 1) == 0) {
			d >>= 1;
			r++;
		}

		for (uint64_t a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
			if (n == a)
				return true;
			if (check_composite(n, a, d, r))
				return false;
		}
		return true;
	}

#endif
