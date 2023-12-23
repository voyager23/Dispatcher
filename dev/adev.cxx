/*
 * adev.cxx
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

using namespace std;

uint64_t finite_field(std::vector<uint64_t> &nodeprime, int taskid, uint64_t n = 1e15);

std::vector<uint64_t> primes = {	// 200 primes
1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,
1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,
1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,
1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,
1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,
1367,1373,1381,1399,1409,1423,1427,1429,1433,1439,
1459,1471,1481,1483,1487,1489,1493,1447,1451,1453, 
1499,1511,1523,1531,1543,1549,1553,1559,1567,1571,
1579,1583,1597,1601,1607,1609,1613,1619,1621,1627,
1637,1657,1663,1667,1669,1693,1697,1699,1709,1721,
1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,
1801,1811,1823,1831,1847,1861,1867,1871,1873,1877,
1879,1889,1901,1907,1913,1931,1933,1949,1951,1973,
1979,1987,1993,1997,1999,2003,2011,2017,2027,2029,
2039,2053,2063,2069,2081,2083,2087,2089,2099,2111,
2113,2129,2131,2137,2141,2143,2153,2161,2179,2203,
2207,2213,2221,2237,2239,2243,2251,2267,2269,2273,
2281,2287,2293,2297,2309,2311,2333,2339,2341,2347,
2351,2357,2371,2377,2381,2383,2389,2393,2399,2411,
2417,2423,2437,2441,2447,2459,2467,2473,2477,2503,1,2,3};

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
			if (processed % 100 == 0){
				cout << taskid << ": " << processed << "/" << nodeprime.size() << endl;
			}
		} // Match: a[idx] => a[1000]		
	} // for prime:nodename
	// cout << taskid <<") local_b = " << local_b << endl;
	return local_b;			
}



//======================================================================

int main (int argc, char *argv[])
{
	
	int  taskid, numtasks, len, partner, message ;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	// For this to work numtasks * stride <= primes.size()
	const unsigned stride = primes.size() / (numtasks*4);
	
	// Assume Root node(0) are dedicated to preparation, dispatch and receive
	// All other nodes are process nodes.
	
	if(taskid == 0) {
		// divide the primes vector into blocks of 'stride' or less and save as discrete
		// vectors in a queue of blocks.
		
		cout << "Stride: " << stride << endl;
		
		queue<std::vector<uint64_t>> blocks;
		vector<uint64_t> new_blk;
		vector<uint64_t>::iterator i,j;
		
		i = primes.begin();
		j = i + stride;
		new_blk.clear();
		
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
		 
		// Dispatch Mode
		// initial round of data disribution
		for (int n = 1; n != numtasks; ++n) {
			MPI_Send( (blocks.front()).data(), stride, MPI_UNSIGNED_LONG_LONG, n, 123, MPI_COMM_WORLD);
			blocks.pop();
		}
		// Recv Mode
		uint64_t B = 0;
		uint64_t recvbuff;
		// recv all the results and sum to B
		set<int>nodes;
		for(int n = 1; n != numtasks; ++n) nodes.insert(n);
		// listen for results
		do {
			MPI_Recv(&recvbuff, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE, 321, MPI_COMM_WORLD, &status);
			B += recvbuff;
			nodes.erase(status.MPI_SOURCE);
			cout << status.MPI_SOURCE << ") " << recvbuff << ":" << B << endl;
		} while(!nodes.empty());
		
		
		
		
	} 
	else 
	{ 
		// Work node
		vector<uint64_t> recvbuff;
		int count;
		recvbuff.resize(stride);
		
		MPI_Recv(recvbuff.data(), stride, MPI_UNSIGNED_LONG_LONG, 0, 123, MPI_COMM_WORLD, &status);
		MPI_Get_count( &status, MPI_UNSIGNED_LONG_LONG, &count );
		
		//~ for (auto p : recvbuff) cout << p << " ";
		//~ cout << endl;
		
		uint64_t local_b = finite_field(recvbuff, taskid);
		// send this back to root node
		MPI_Send(&local_b, 1, MPI_UNSIGNED_LONG_LONG, 0, 321, MPI_COMM_WORLD);
		
		
	}
	
	
	// Required clean up.
	MPI_Finalize();
	// Required by C++ standard
	return 0;
	
}

