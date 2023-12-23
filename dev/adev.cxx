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

using namespace std;

std::vector<uint64_t> primes = {	// 200 primes
1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,
1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,
1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,
1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,
1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,
1367,1373,1381,1399,1409,1423,1427,1429,1433,1439,
1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,
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


//======================================================================

int main (int argc, char *argv[])
{
	const unsigned stride = 7;
	
	int  taskid, numtasks, len, partner, message ;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	
	// Assume Root node(0) are dedicated to preparation, dispatch and receive
	// All other nodes are process nodes.
	
	if(taskid == 0) {
		// divide the primes vector into blocks of 'stride' or less and save as discrete
		// vectors in a queue.
		std::queue<std::vector<uint64_t>> blocks;
		std::vector<uint64_t> temp;
		std::vector<uint64_t>::iterator i,j;
		
		i = primes.begin();
		j = i + stride;
		temp.clear();
		
		while(1) {
			do {
				temp.push_back(*i);
				if(++i == primes.end()) break;
			} while (i != j);
			// pad with zeros
			while(temp.size() < stride) temp.push_back(0);
			blocks.push(temp);
			temp.clear();
			
			if(i == primes.end()) break;
			
			j = i + stride;
		}
		 
		// Dispatch Mode
		// initial round of data disribution
		for (int n = 1; n != numtasks; ++n) {
			MPI_Send( (blocks.front()).data(), stride, MPI_UNSIGNED_LONG_LONG, n, 123, MPI_COMM_WORLD);
			blocks.pop();
		}
		
	

	} else { 
		// Work node
		vector<uint64_t> recvbuff;
		int count;
		recvbuff.resize(stride);
		
		MPI_Recv(recvbuff.data(), stride, MPI_UNSIGNED_LONG_LONG, 0, 123, MPI_COMM_WORLD, &status);
		MPI_Get_count( &status, MPI_UNSIGNED_LONG_LONG, &count );
		
		cout << "Count: " << count << endl;
		
		for (auto p : recvbuff) cout << p << " ";
		cout << endl;
		
	}
	
	
	// Required clean up.
	MPI_Finalize();
	// Required by C++ standard
	return 0;
	
}

