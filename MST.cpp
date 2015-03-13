#include "MST.h"
#include <list>
#include "Minmatching/PerfectMatching.h"




class Graph
{
    int V;   

    void Util(int v, bool visited[]);
public:
    Graph(int V);
    void addEd(int v, int w);
    void DFS();
    vector<int> vertciesForTSP2;
    list<int> *adj;
    int totalEulerCost = 0;
    vector<pair<int,int>> vertciesForTSP15;


    void removeEd(int u, int v);

    void printEulerTour();
    void printEulerUtil(int s);

    int DFSCt(int v, bool visited[]);


    bool isValid(int u, int v);
};


Graph::Graph(int V)
{
    this->V = V;
    adj = new list<int>[V];
}

void Graph::addEd(int v, int w)
{
    adj[v].push_back(w); // Add w to v’s list.
    adj[w].push_back(v);
}


void Graph::Util(int v, bool visited[])
{
    // Mark the current node as visited and print it
    visited[v] = true;
    //cout << v << " ";
    vertciesForTSP2.push_back(v);

    list<int>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i)
        if(!visited[*i])
            Util(*i, visited);
}


void Graph::DFS()
{
    bool *visited = new bool[V];
    for(int i = 0; i < V; i++)
        visited[i] = false;


    for(int i = 0; i < V; i++)
        if(visited[i] == false)
            Util(i, visited);
}






//For TSP1,5

void Graph::removeEd(int u, int v)
{
  // Find v in adjacency list of u and replace it with -1
  list<int>::iterator iv = find(adj[u].begin(), adj[u].end(), v);
  *iv = -1;

  // Find u in adjacency list of v and replace it with -1
  list<int>::iterator iu = find(adj[v].begin(), adj[v].end(), u);
  *iu = -1;
}


void Graph::printEulerTour()
{
  // Find a vertex with odd degree
  int u = 0;
  for (int i = 0; i < V; i++)
      if (adj[i].size() & 1)
        {   u = i; break;  }

  // Print tour starting from oddv
  printEulerUtil(u);
  cout << endl;
}



// Print Euler tour starting from vertex u
void Graph::printEulerUtil(int u)
{
  // Recur for all the vertices adjacent to this vertex
  list<int>::iterator i;
  for (i = adj[u].begin(); i != adj[u].end(); ++i)
  {
      int v = *i;

      // If edge u-v is not removed and it's a a valid next edge
      if (v != -1 && isValid(u, v))
      {
          //cout << u << "-" << v << "  ";
          vertciesForTSP15.push_back(make_pair (u,v));
          removeEd(u, v);
          printEulerUtil(v);
      }
  }
}


bool Graph::isValid(int u, int v)
{
  // The edge u-v is valid in one of the following two cases:

  // 1) If v is the only adjacent vertex of u
  int count = 0;  // To store count of adjacent vertices
  list<int>::iterator i;
  for (i = adj[u].begin(); i != adj[u].end(); ++i)
     if (*i != -1)
        count++;
  if (count == 1)
    return true;


  // 2) If there are multiple adjacents, then u-v is not a bridge
  // Do following steps to check if u-v is a bridge

  // 2.a) count of vertices reachable from u
  bool visited[V];
  memset(visited, false, V);
  int count1 = DFSCt(u, visited);

  // 2.b) Remove edge (u, v) and after removing the edge, count
  // vertices reachable from u
  removeEd(u, v);
  memset(visited, false, V);
  int count2 = DFSCt(u, visited);

  // 2.c) Add the edge back to the graph
  addEd(u, v);

  // 2.d) If count1 is greater, then edge (u, v) is a bridge
  return (count1 > count2)? false: true;
}



int Graph::DFSCt(int v, bool visited[])
{
  // Mark the current node as visited
  visited[v] = true;
  int count = 1;

  // Recur for all vertices adjacent to this vertex
  list<int>::iterator i;
  for (i = adj[v].begin(); i != adj[v].end(); ++i)
      if (*i != -1 && !visited[*i])
          count += DFSCt(*i, visited);

  return count;
}



// struct Edge
// {
//     int src, dest, weight;
// };


// struct Graph
// {
//     int V, E;

//     struct Edge* edge;
// };

// struct Graph* createGraph(int V, int E)
// {
//     struct Graph* graph = (struct Graph*) malloc( sizeof(struct Graph) );
//     graph->V = V;
//     graph->E = E;

//     graph->edge = (struct Edge*) malloc( graph->E * sizeof( struct Edge ) );

//     return graph;
// }


//use Prim's algorithm or Kruskal algorithm. Copied from 'http://www.geeksforgeeks.org/greedy-algorithms-set-5-prims-minimum-spanning-tree-mst-2/'
MST::MST(float** input, int size) {
    adjacentMatrix = input;
    key = new int[size];
    mstSet = new bool[size];
    parent = new int[size];

    N = size;
}

MST::~MST() {

}

//use Prim's algorithm or Kruskal algorithm. Copied from 'http://www.geeksforgeeks.org/greedy-algorithms-set-5-prims-minimum-spanning-tree-mst-2/'
void MST::makeTree() {
     // Initialize all keys as INFINITE
     for (int i = 0; i < N; i++)
        key[i] = INT_MAX, mstSet[i] = false;

     // Always include first 1st vertex in MST.
     key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
     parent[0] = -1; // First node is always root of MST

     // The MST will have V vertices
     for (int count = 0; count < N-1; count++)
     {
        // Pick thd minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet);

        // Add the picked vertex to the MST Set
        mstSet[u] = true;

        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < N; v++)
           // mstSet[v] is false for vertices not yet included in MST
           // Update the key only if adjacentMatrix[u][v] is smaller than key[v]
          if (adjacentMatrix[u][v] && mstSet[v] == false && adjacentMatrix[u][v] <  key[v])
             parent[v]  = u, key[v] = adjacentMatrix[u][v];
     }


}

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int MST::minKey(int key[], bool mstSet[])
{
   // Initialize min value
   int min = INT_MAX, min_index;

   for (int v = 0; v < N; v++)
     if (mstSet[v] == false && key[v] < min)
         min = key[v], min_index = v;

   return min_index;
}

// A utility function to print the constructed MST stored in parent[]
void MST::printMST() {

  int mstCost = 0;
	cout<<endl;
	cout<<"Minimum spanning tree from the adjacency matrix"<<endl;
	cout<<"Edge   Weight"<<endl;
	for (int i = 1; i < N; i++) {
		cout<<parent[i]<<" - "<<i<<"  "<<adjacentMatrix[i][parent[i]]<<endl;
        //adjacentMatrixMST[parent[i]][i] = adjacentMatrix[i][parent[i]];
    mstCost += adjacentMatrix[i][parent[i]];

	}

  cout << "MST COST IS ---------------------------------> " << mstCost << endl;

	// Graph g(N);
	// for (int i = 1; i < N; i++) {
	// 	g.addEdge(parent[i], i);
	// }

 //    adjListMST = g.adj;

	// g.DFS();



}


void MST::LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, float** adjacentMatrix, int N) {
    int e = 0;


    for (int i = 0; i < N; i++) {


        if(adjListMST[i].size() % 2 !=0) {


            numOfOddVertices++;
            //cout << "the size of " << i << " is odd!" << endl;
            oddArray.push_back(i);
            numOfOdd++;


        }

    }

    cout << " Odd Vertices are" << oddArray.size() << endl;

    node_num = numOfOdd;
    edge_num = numOfOdd*(numOfOdd-1)/2; //complete graph





    edges = new int[2*(numOfOdd*(numOfOdd-1)/2)];
    weights = new int[(numOfOdd*(numOfOdd-1)/2)];


    e=0;

    for(int i = 0; i < numOfOdd ; ++i) {
        for(int j = i+1 ; j< numOfOdd ; ++j) {

           // cout << "Adding at indeces " << i << " " << j;
          //  cout << " Vertices are " << oddArray[i] << " to " << oddArray[j] << endl;
            //cout << "Adding to EDGES array at " << 2*e << " and " << 2*e+1 << endl;
            edges[2*e] = i;
            edges[2*e+1] = j;
            weights[e] = adjacentMatrix[oddArray[i]][oddArray[j]];
            e++;
        }
    }




    if (e != edge_num) {
        cout<<"the number of edge is wrong"<<endl;

        exit(1);
    }
}

void MST::PrintMatching(int node_num, PerfectMatching* pm) {
    int i, j;

    for (i=0; i<node_num; i++) {
        j = pm->GetMatch(i);
        if (i < j) {
            //printf("%d %d\n", i, j);
            //cout << " Vertices are " << oddArray[i] << " to " << oddArray[j] << endl;

        }
    }
}

//calculate mean of all edges in the MST
// float MST::calMean(int option) {
// 	float mean = 0.0;

// 	if(option == MST_1) {
// 		//calculate
// 	}else if(option == TSP2) {

// 	} else if(option == TSP1_5) {

// 	}

// 	return mean;
// }

// //calculate standard deviation of all edges in the MST
// float MST::calStd(int option) {
// 	float std = 0.0;

// 	if(option == MST_1) {
// 		//calculate
// 	}else if(option == TSP2) {

// 	} else if(option == TSP1_5) {

// 	}

// 	return std;
// }

int MST::mstCost() {


  int mstCost = 0;

  for (int i = 1; i < N; i++) {
    //cout<<parent[i]<<" - "<<i<<"  "<<adjacentMatrix[i][parent[i]]<<endl;
        //adjacentMatrixMST[parent[i]][i] = adjacentMatrix[i][parent[i]];
    mstCost += adjacentMatrix[i][parent[i]];

  }



  return mstCost;
}

int MST::makeTSP2() {

    heuristicTSP2Cost = 0;

    Graph g(N);
    for (int i = 1; i < N; i++) {
        g.addEd(parent[i], i);
    }

    g.DFS();
    vTSP2 = g.vertciesForTSP2;
    adjListMST = g.adj;

    for(auto it = vTSP2.begin() ; it!=vTSP2.end() ; it++) {
        //cout << *it << endl;
        heuristicTSP2Cost+=adjacentMatrix[*it][*(it+1)];
        //cout << "The  cost from " << *it << " to " << *(it+1) << " is " <<  adjacentMatrix[*it][*(it+1)] << endl;
    }

    cout << "The total cost is " << heuristicTSP2Cost << endl;




    return heuristicTSP2Cost;

    // cout << "The  cost from 3 to 7 is " <<  adjacentMatrix[3][7] << endl;
	//make a Eulerian tour by DFS

	//add shortcuts if a vertex has no detours.

	//calculate heuristic TSP cost
}






int MST::makeTSP1_5() {

	//construct minimum-weight-matching for the given MST
	// minimumMatching();

	// //make all edges has even degree by combining mimimum-weight matching and MST
	// combine();





    struct PerfectMatching::Options options;
    int i, e, node_num = N-1, edge_num;
    int* edges;
    int* weights;


    LoadInput(node_num, edge_num, edges, weights, adjacentMatrix, N);
    node_num = numOfOdd;
    edge_num = node_num*(node_num-1)/2;


    cout << " The numOfOdd VERTICESis" << node_num << endl;
    cout << " The EDGE num is " << edge_num << endl;


    PerfectMatching *pm = new PerfectMatching(node_num, edge_num);


    for (e=0; e<edge_num; e++) {
        //cout << edges[2*e] << "-" << edges[2*e+1] << " weight = " << weights[e] << endl;
       pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
    }

    pm->options = options;
    pm->Solve();

    // double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
    // printf("Total cost of the perfect min-weight matching = %.1f\n", cost);

    PrintMatching(node_num, pm);

    //addung original MST
    Graph g(N);
    for (int i = 1; i < N; i++) {
        g.addEd(parent[i], i);
    }


    //ading PM edges
    int i2, j2;

    for (i2=0; i2<node_num; i2++) {
        j2 = pm->GetMatch(i2);
        if (i2 < j2) {
            //printf("%d %d\n", i2, j2);
            //cout << " Vertices are " << oddArray[i2] << " to " << oddArray[j2] << endl;
            g.addEd(oddArray[i2], oddArray[j2]);

        }
    }

    g.printEulerTour();
    vTSP15 = g.vertciesForTSP15;



    //now g should have all the edges upto PM+MST

    //complete euler tour on g.
    cout << "printing tour again" << endl;
    int totalTSP15 = 0;
    for (auto pos = vTSP15.begin() ; pos != vTSP15.end(); ++pos)
    {
      //cout << pos->first << " " << pos->second << endl;
      totalTSP15 += adjacentMatrix[pos->first][pos->second];
    }


    cout << "END COST TSP1.5 is ----- " << totalTSP15 << endl;


    //make new graph with MST and PM.

    delete pm;
    delete [] edges;
    delete [] weights;


    return totalTSP15;

	//calculate heuristic TSP cost
}

// void MST::minimumMatching() { //if you choose O(n^2)
// 	//find minimum-weight matching for the MST.

// 	//you should carefully choose a matching algorithm to optimize the TSP cost.
// }

// void MST::combine() {
// 	//combine minimum-weight matching with the MST to get a multigraph which has vertices with even degree
// }
