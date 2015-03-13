#include "common.h"
#include "Point.h"
#include "MST.h"
#include "Minmatching/PerfectMatching.h"

/*
This project is a starter code and wrappers for CSE101W15 Implementation project.

point.h - uniform random pointset generator

MST.h - minimum spanning tree

PerfectMatching.h - interface to min cost perfect matching code

-------------------------------------
PerfectMatching is from the paper:

Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.

sourcecode : pub.ist.ac.at/~vnk/software/blossom5-v2.05.src.tar.gz

*/

// void LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, float** adjacentMatrix, int N) {
// 	int e = 0;
// 	node_num = N;
// 	edge_num = N*(N-1)/2 ; //complete graph

// 	edges = new int[2*edge_num];
// 	weights = new int[edge_num];

// 	for(int i = 0; i < N ; ++i) {
// 		for(int j = i+1 ; j< N ; ++j) {
// 			edges[2*e] = i;
// 			edges[2*e+1] = j;
// 			weights[e] = adjacentMatrix[i][j];
// 			e++;
// 		}
// 	}

// 	if (e != edge_num) {
// 		cout<<"the number of edge is wrong"<<endl;

// 		exit(1);
// 	}
// }

// void PrintMatching(int node_num, PerfectMatching* pm) {
// 	int i, j;

// 	for (i=0; i<node_num; i++) {
// 		j = pm->GetMatch(i);
// 		if (i < j) printf("%d %d\n", i, j);
// 	}
// }

int main() {
	set< pair<int,int> > generatedPointset;
	float** adjacentMatrix;
	int W, H, N;
	Point pointset;

	W = 5000;
	H = 7000;
	N = 100;


	cout << "Please Enter Number of Trials n = .... " << endl;
	int numTrials ;
	cin >> numTrials;

	//making a freaking for loop

	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl;

	int sumMST =0 , sumTSP2=0 , sumTSP1_5=0;
	int sqMST =0, sqTSP2 =0, sqTSP1_5=0;

	vector<int> MST_COSTS ;
	vector<int> TSP2_COSTS ;
	vector<int> TSP1_5_COSTS ;


	for(int i=0 ; i<numTrials ; i++) {

			pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation
			pointset.printPointset();

			generatedPointset = pointset.getPointset();
			adjacentMatrix = pointset.getAdjacentMatrix();

			//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal
			MST mst(adjacentMatrix, N);

			mst.makeTree();
			//mst.printMST();
			int mstCost = mst.mstCost();
			int tsp2Cost = mst.makeTSP2();
			int tsp15Cost = mst.makeTSP1_5();

			MST_COSTS.push_back(mstCost);
			// TSP2_COSTS.push_back(tsp2Cost);
			// TSP1_5_COSTS.push_back(tsp15Cost);

			sumMST += mstCost;
			sumTSP2 += tsp2Cost;
			sumTSP1_5 += tsp15Cost;

			//change this
			// sqMST += pow(mstCost , 2);
			// sqTSP2 += pow(tsp2Cost , 2);
			// sqTSP1_5 += pow(tsp15Cost , 2);
			cout << mstCost << " " << tsp2Cost << " " << tsp15Cost << endl;


	}


	// double sum = std::accumulate(MST_COSTS.begin(), MST_COSTS.end(), 0.0);
	// cout << "sum is " << sum << endl;
	// double mean2 = sum / MST_COSTS.size();
	// cout << "mean2 is " << mean2 << endl;
	// //
	// double sq_sum = std::inner_product(MST_COSTS.begin(), MST_COSTS.end(), MST_COSTS.begin(), 0.0);
	// cout << "sq_sum is " << sq_sum << endl;
	// double stdev = std::sqrt(sq_sum / MST_COSTS.size() - mean2 * mean2);
	// cout << "stdDev is " << stdev << endl;
	// cout << "THE STANDARD DEV FOR MST IS    " << stdev << endl;







	cout << "THE MEAN FOR MST IS .. " << sumMST/numTrials << endl;
	cout << "THE MEAN FOR TSP2 IS .. " << sumTSP2/numTrials << endl;
	cout << "THE MEAN FOR TSP1_5 IS .. " << sumTSP1_5/numTrials << endl;

	//cout << (sqMST) << endl;
	double mean = sumMST/numTrials;
	double sqSum = 0;
	// for(auto it = MST_COSTS.begin() ; it != MST_COSTS.end()-1 ; it++) {
	//
	// 	//sqSum += pow((*it) - mean , 2);
	// 	cout << pow((*it) - mean , 2) << endl;
	// }

	for(int i=0 ; i<numTrials ; i++) {
		cout << MST_COSTS[i] << endl;
	}

	double stdev = 0;

	//delete MST_COSTS;

// 	static long double sum = std::accumulate(std::begin(MST_COSTS), std::end(MST_COSTS), 0.0);
// static long double m =  sum / MST_COSTS.size();
//
//  static long double accum = 0.0;
// std::for_each (std::begin(MST_COSTS), std::end(MST_COSTS), [&](const long double d) {
//     accum += (d - m) * (d - m);
// });
//
// static long double stdev = sqrt(accum / (MST_COSTS.size()-1));

cout << "STD DEV IS---> " << stdev << endl;


	// int a = (sqMST/numTrials) ;
	// int b =  pow((((sumMST/numTrials))) , 2) ;
	// float std = sqrt(a-b);
	// cout << "STD Dev MST is .. " << std << endl;


	// mst.makeTSP2();
	// mst.makeTSP1_5();

	//Deliverable B: Find TSP2 path from the constructed MST
	//You won't need any wrappers for B.


	//Deliverable C: Find TSP1.5 path from the constructed MST

	//Find the perfect minimum-weight matching
	// struct PerfectMatching::Options options;
	// int i, e, node_num = N, edge_num = N*(N-1)/2;
	// int* edges;
	// int* weights;
	// PerfectMatching *pm = new PerfectMatching(node_num, edge_num);

	// LoadInput(node_num, edge_num, edges, weights, adjacentMatrix, N);

	// for (e=0; e<edge_num; e++) {
	// 	pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
	// }

	// pm->options = options;
	// pm->Solve();

	// double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
	// printf("Total cost of the perfect min-weight matching = %.1f\n", cost);

	// PrintMatching(node_num, pm);

	// delete pm;
	// delete [] edges;
	// delete [] weights;

	return 0;
}
