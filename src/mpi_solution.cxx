/*
 * asoln.cxx
 * 
 * Copyright 2023 mike <mike@Fedora37>
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
 * 492. 242586962923928
 */

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <array>
#include <unordered_map>
#include <climits>
#include <cstdint>

using namespace std;

// Declaration
bool MillerRabin(uint64_t n);
uint64_t finite_field(std::vector<uint64_t> &nodeprime);

// Constants common to all nodes
const int MAX_PROCS = 6 + 4 + 8 + (8 + 8 + 8 + 8);	// maximum number of cores available {50}
const uint64_t n = 1e15;		// maximum number of iteration to calc. Example 1 n=10^3
const uint64_t y = 1e7;
const uint64_t x = 1e9;
// Variable common to all nodes
int taskid;


// Definitions
vector<uint64_t> prime_modulus(uint64_t x, uint64_t y){
	// approx 25 seconds for x=10^9 and y=10^7 returns 482449 values
	std::vector<uint64_t> primes = {};
	if((x % 2)==0) x += 1; // test odd values
	for(uint64_t p = x; p <= x+y; p+=2) {
		if (MillerRabin(p)) {
			if(p <= UINT_MAX){
				primes.push_back(p);
			} else {
			cout << "Warning: prime_modulus overflow" << endl;
			}
		}
	}
	return primes;
}

uint64_t finite_field(std::vector<uint64_t> &nodeprime)
{
	// Version: 19.12.2023 14:09
	// Modify the reference point to a[1000]. This seems to allow all moduli
	// to return a cycle.
	
	uint64_t aRef, iRef, a, idx, local_b=0;
	const uint64_t block_size = 1000;	// record progress every 1000 values
	std::unordered_map<uint64_t,uint64_t> progress;
	std::unordered_map<uint64_t,uint64_t>::iterator iter;
	uint64_t processed = 0;
	for(uint64_t p : nodeprime) {
			
		//cout << "Modulus:" << p << endl;
		
		// Move fwd to a[block_size]
		a = 2359; idx = 3;	
		do {
			a = (6*a*a + 10*a + 3) % p;
			idx += 1;
		}while(idx != block_size);
		
		aRef = a;	//Reference value
		
		// Save first map entry here - map.emplace(1000, aRef)
		progress.clear();
		progress.emplace(3,2359);	// cover case where reset idx = 0
		progress.emplace(block_size, aRef);
		
		// Search forward for matching value to aRef
		
		do {
			a = (6*a*a + 10*a + 3) % p;
			if (((++idx) % block_size) == 0) { 
				progress.emplace(idx,a); //  progress has (4000, a) etc
			}
			
		} while ((a != aRef)&&(idx <= n));
		
		if(idx > n){
			std::cout<<"Error: idx > n."<< endl;			
			exit(1);
		} 
		else 
		{ // Match: a[idx] => aRef
			uint64_t order = idx - block_size;
			uint64_t residue = (n - block_size + 1) % order;	// corrected this line
			int64_t offset = residue - 1;
			if(offset < 0) offset += order;
			uint64_t req_idx = block_size + offset;
			
			// refer to progress map
			idx = (req_idx / block_size) * block_size;			
			if(idx==0) {	// special case
				idx = 3;	// progress has (3,2359)
			} else {
				iter = progress.find(idx);
				if(iter == progress.end()) {
					cout<<"entry not found in progress"<<endl;
					exit(1);
				}

			}
			// reset value of a from progress map			
			a = (*iter).second;			
			// ----------------
			while(idx != req_idx){
				a = (6*a*a + 10*a + 3) % p;
				++idx;
			}
			// a now has the value of a{n}
			local_b += a;
			// << "a[" << n << "] mod " << p << " = " << a << endl;
			processed += 1;
			if (processed % 100 == 0){
				cout << taskid << ": " << processed << "/" << nodeprime.size() << endl;
			}
		} // Match: a[idx] => a7
		
	} // for prime:nodename
	cout << taskid <<") local_b = " << local_b << endl;
	return local_b;			
}


//======================================================================

int main (int argc, char *argv[])
{
	
	int  numtasks, len, partner, message;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	
	int count[MAX_PROCS];
	int displace[MAX_PROCS];
	
	// setup vector of primes
	std::vector<uint64_t> primes;	// Referenced in Thread data block
	
	// Calculate the prime vector
	primes = prime_modulus(x,y);	// Expect 49 primes for example 2
	
	// setup scatterv
	vector<uint64_t> local_p;
	int c = primes.size() / numtasks;
	int r = primes.size() % numtasks;
	// count array
	for(auto i = 0; i != numtasks; ++i) count[i] = c;	// basic counts for each process
	for(auto i = 1, residue = r; residue != 0; ++i, --residue) count[i] += 1; // add residues
	// displace array;
	displace[0] = 0;
	for(auto i = 0, j = 1; i != numtasks; ++i,++j) displace[j] = count[i] + displace[j-1];
	
	if(taskid == 0) {
		
		// DEBUG
		cout << primes.size() << " primes." << endl;
		cout << "n = " << n << endl;
		cout << "Tasks/node: ";
		for(auto i=0; i != numtasks; ++i) cout << count[i] << " ";
		cout << endl;
		cout << "Displacements: ";
		for(auto i=0; i != numtasks; ++i) cout << displace[i] << " ";
		cout << endl;	
		// END DEBUG
		
		// Local result storage
		vector<array<uint64_t,3>> v_results;
		
		// Scatter
		local_p.resize(count[taskid]);
		int count_recv = local_p.size();
		MPI_Scatterv(primes.data(), count, displace, MPI_UINT64_T,
			local_p.data(), count_recv,
			MPI_UINT64_T, 0, MPI_COMM_WORLD);
			
		// master node process primes
		uint64_t root_B = 0;
	
		// ----------Cyclic Root search begins here----------
		root_B = finite_field(local_p);

		// Reduce_Sum the local sums
		uint64_t B = 0;

		MPI_Reduce( &root_B, &B, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
		
		// Show final sum B;
		cout << "Final result B = " << B << endl;
		
	} else { // Node
		// Nodes expect to receive a vector of primes which is stored in nodeprime
		vector<uint64_t> nodeprime;
		nodeprime.resize(count[taskid]);	// make space for the primes
		int count_recv = count[taskid];		// set to expected receive count
		
		// Receive primes
		MPI_Scatterv(primes.data(), count, displace, MPI_UINT64_T,
		nodeprime.data(), count_recv,
		MPI_UINT64_T, 0, MPI_COMM_WORLD);
		
		// Work Node Process primes
		uint64_t local_B = 0;
	
		// ----------Cyclic Node search begins here----------
		local_B = finite_field(nodeprime);
		
		// Reduce_Sum the result package for this node
		MPI_Reduce( &local_B, NULL, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
					
	} // else Node...
	
	// Required clean up.
	MPI_Finalize();
	// Required by C++ standard
	return 0;
}

