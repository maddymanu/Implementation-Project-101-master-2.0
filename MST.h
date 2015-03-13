#include "common.h"
#include <vector>
#include <list>

#include "Minmatching/PerfectMatching.h"


#pragma once

class MST {
public:
	float** adjacentMatrix;

	float** adjacentMatrixMST;
	int* parent; //Array to store constructed MST
	int* key; //Key values used to pick minimum weight edge in cut
	bool* mstSet; //To represent set of vertices not yet included in MST
	int N; //the size of pointset
	vector<int> vTSP2;
	vector<pair<int , int>> vTSP15;

	int heuristicTSP2Cost;
	list<int> *adjListMST;
	int numOfOdd = 0;
	int numOfOddVertices = 0;
	vector<int> oddArray;
	int numOfEdges=0;
	
	

	MST(float** adjacentMatrix, int size);
	~MST();

	//deliverable a
	void makeTree();
	void printMST();

	void LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, float** adjacentMatrix, int N);
	void PrintMatching(int node_num, PerfectMatching* pm);

	//deliverable b
	int mstCost();

	int makeTSP2();

	//deliverable c
	int makeTSP1_5();
	
	// float calMean(int option);
	// float calStd(int option);

private:
	// void minimumMatching();
	// void combine();
	int minKey(int key[], bool mstSet[]);

};