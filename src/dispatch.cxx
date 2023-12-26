/*
 * bdev.cxx
 * 
 * Copyright 2023 mike <mike@Fedora37>
 * 
 */

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <queue>
#include <unordered_map>
#include <set>

#include "../inc/MillerRabin.hxx"

using namespace std;

uint64_t finite_field(std::vector<uint64_t> &nodeprime, int taskid, uint64_t n = 1e15);
// Declaration
bool MillerRabin(uint64_t n);

// Constants common to all nodes
const int MAX_PROCS = 6 + 4 + 8 + (8 + 8 + 8 + 8);	// maximum number of cores available {50}
const uint64_t n = 1e15;		// maximum number of iteration to calc. Example 1 n=10^3
const uint64_t y = 1e4;			// 1e6 requires 27m 18s with 11 + 1 processes
const uint64_t x = 1e9;

uint64_t finite_field(std::vector<uint64_t> &nodeprime, int taskid, uint64_t n)
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
		
		// deal with part blocks
		if (p == 0) continue;
			
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
			// cout << "a[" << n << "] mod " << p << " = " << a << endl;
			processed += 1;
			//~ if (processed % 100 == 0){
				//~ cout << taskid << ": " << processed << "/" << nodeprime.size() << endl;
			//~ }
		} // Match: a[idx] => a[1000]		
	} // for prime:nodename
	// cout << taskid <<") local_b = " << local_b << endl;
	return local_b;			
}

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


//======================================================================

int main (int argc, char *argv[])
{	
	int  taskid, numtasks, len, partner, message ;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;	// used for either root or work nodes
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	unsigned stride;	// For this to work numtasks * stride <= primes.size()

	// Assume Root node(0) are dedicated to preparation, dispatch and receive
	// All other nodes are process nodes.
	
	if(taskid == 0) { // -----Root Node-----
		
		vector<uint64_t>primes = prime_modulus(x,y);
		stride = primes.size() / (numtasks*4);
		
		// divide the primes vector into blocks of 'stride' or less and save as discrete
		// vectors in a queue of blocks.		
		cout << "Stride: " << stride << endl;		
		queue<std::vector<uint64_t>> blocks;
		vector<uint64_t> new_blk;
		vector<uint64_t>::iterator i,j;
		// setup iterators and temp variable
		i = primes.begin();
		j = i + stride;
		new_blk.clear();
		
		// populate the blocks queue
		while(1) {
			do {
				new_blk.push_back(*i);
				if(++i == primes.end()) break;
			} while (i != j);
			// pad with zeros
			while(new_blk.size() < stride) new_blk.push_back(0);
			blocks.push(new_blk);
			if(i == primes.end()) break;
			// prepare next loop
			j = i + stride;
			new_blk.clear();
		}
		
		// Send/Recv Server
		uint64_t B = 0;
		uint64_t recvbuff; // expect 1 element
		int busy = 0;
		const vector<uint64_t> stop_block = {0,stride};
		
		// Send a data block to all nodes
		for(int task = 1; task != numtasks; ++task){
			MPI_Send((blocks.front()).data(), stride, MPI_UNSIGNED_LONG_LONG, task, 123, MPI_COMM_WORLD);
			busy += 1;
			blocks.pop();
		}
		// now cycle through remaining blocks
		cout << blocks.size() << " blocks remaining." << endl;
		while(1){
			MPI_Recv(&recvbuff, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, 321, MPI_COMM_WORLD, &status);
			B += recvbuff;
			if (blocks.empty() == true){
				// empty send zero block
				MPI_Send(stop_block.data(), stride, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, 123, MPI_COMM_WORLD);
				if(--busy == 0) break;
			} else {
				// next block
				MPI_Send((blocks.front()).data(), stride, MPI_UNSIGNED_LONG_LONG, status.MPI_SOURCE, 123, MPI_COMM_WORLD);
				blocks.pop();
				cout << blocks.size() << " blocks remaining." << endl;
			}
		}
		
		cout << "Final Sum:" << B << endl;
	}
	
	else // -----Work node-----
	
	{ 
		vector<uint64_t> workbuff;
		workbuff.resize(stride);	// expect 'stride' elements
		uint64_t local_b = 0;
		while(1){
			MPI_Recv(workbuff.data(), stride, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, 123, MPI_COMM_WORLD, &status);
			if(workbuff[0] == 0) break;	// stop block
				cout << taskid << ") working." << endl;
			local_b = finite_field(workbuff, taskid);
			MPI_Send(&local_b, 1, MPI_UNSIGNED_LONG_LONG, 0, 321, MPI_COMM_WORLD);
		}
	}	
	// Required clean up.
	MPI_Finalize();
	// Required by C++ standard
	return 0;
	
}

