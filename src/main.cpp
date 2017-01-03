#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <numeric>
#include "POLRGraph.h"

using namespace std;

double getCostOfPath(vector<Vertex<int, double>* > in)
{
    double out = 0.0;
    for(int i = 0; i < in.size(); i++)
    {
        out += in[i]->adjList.find(in[i+1])->second;
    }
    return out;
}

int main()
{
    Graph<int, double> g;
    g.addNode(1, 1, 1);
    g.addNode(2, 4, 1);
    g.addNode(3, 8, 2);
    g.addNode(4, 5, 3);
    g.addNode(5, 2, 4);
    g.addNode(6, 4, 5);
    g.addNode(7, 3, 7);
    g.addNode(8, 6, 7);
    g.addNode(9, 7, 5);
    g.addNode(10, 8, 8);
    g.addNode(11, 5, 9);
    g.addNode(12, 8, 10);
    g.addNode(13, 11, 5);
    g.addNode(14, 11, 8);
    g.addNode(15, 17, 5);
    g.addNode(16, 14, 7);
    g.addNode(17, 20, 7);
    g.addNode(18, 13, 10);
    g.addNode(19, 11, 12);
    g.addNode(20, 14, 14);
    g.addNode(21, 17, 12);
    g.addNode(22, 17, 9);
    g.addNode(23, 20, 11);
    g.addNode(24, 20, 13);
    g.addEdge1(1, 2);
    g.addEdge1(1, 4);
    g.addEdge1(1, 5);
    g.addEdge1(2, 3);
    g.addEdge1(3, 9);
    g.addEdge1(4, 9);
    g.addEdge1(4, 6);
    g.addEdge1(5, 6);
    g.addEdge1(5, 7);
    g.addEdge1(6, 8);
    g.addEdge1(7, 8);
    g.addEdge1(9, 8);
    g.addEdge1(8, 10);
    g.addEdge1(9, 10);
    g.addEdge1(7, 11);
    g.addEdge1(10, 12);
    g.addEdge1(11, 12);
    g.addEdge1(11, 13);
    g.addEdge1(3, 13);
    g.addEdge1(9, 13);
    g.addEdge1(10, 13);
    g.addEdge1(10, 14);
    g.addEdge1(12, 19);
    g.addEdge1(12, 13); 
    g.addEdge1(10, 20);
    g.addEdge1(19, 20);
    g.addEdge1(13, 14);
    g.addEdge1(13, 15);
    g.addEdge1(13, 16);
    g.addEdge1(16, 15);
    g.addEdge1(16, 17);
    g.addEdge1(15, 17);
    g.addEdge1(14, 16);
    g.addEdge1(14, 18);
    g.addEdge1(14, 22);
    g.addEdge1(18, 22);
    g.addEdge1(22, 17);
    g.addEdge1(17, 23);
    g.addEdge1(22, 23);
    g.addEdge1(18, 21);
    g.addEdge1(20, 21);
    g.addEdge1(21, 23);
    g.addEdge1(23, 24);


    vector<double> aStarDist, POLRQuickDist, POLRDist, leastDeviationLDist, leastDeviationPosDist;
    vector<double> aStarTime, POLRQuickTime, POLRTime, leastDevLTime, leastDevPosTime;
    int startTime, stopTime;
    double tm;
    vector<Vertex<int, double>* > path ;
    vector<int> xS, yS;
    for(int i = 1; i <= 5; i++)
    {
        xS.emplace_back(i);
    }
    for(int i = 19; i <= 24; i++)
    {
        yS.emplace_back(i);
    }
    for(const auto& x: xS)
    {
        for(const auto& y: yS)
        {
            if(x != y)
            {
                startTime = clock();
                path = g.aStar(x,y);
                stopTime = clock();
                tm = ((double)(stopTime - startTime) / (double)CLOCKS_PER_SEC);
                aStarDist.push_back(getCostOfPath(path));
                aStarTime.push_back(tm);

                startTime = clock();
                path = g.createPath(x, y);
                stopTime = clock();
                tm = ((double)(stopTime - startTime) / (double)CLOCKS_PER_SEC);
                POLRDist.push_back(getCostOfPath(path));
                POLRTime.push_back(tm);


                startTime = clock();
                path = g.createPathStraightLine1(x, y);
                stopTime = clock();
                tm = ((double)(stopTime - startTime) / (double)CLOCKS_PER_SEC);
                leastDeviationLDist.push_back(getCostOfPath(path));
                leastDevLTime.push_back(tm);
            }
        }
    }
    cout << "paths found" << endl;
    double aD = accumulate(aStarDist.begin(), aStarDist.end(), 0.0)/aStarDist.size();
    double aT = accumulate(aStarTime.begin(), aStarTime.end(), 0.0)/aStarTime.size();
    double pD = accumulate(POLRDist.begin(), POLRDist.end(), 0.0)/POLRDist.size();
    double pT = accumulate(POLRTime.begin(), POLRTime.end(), 0.0)/POLRTime.size();
    double llD = accumulate(leastDeviationLDist.begin(), leastDeviationLDist.end(), 0.0)/leastDeviationLDist.size();
    double llT = accumulate(leastDevLTime.begin(), leastDevLTime.end(), 0.0)/leastDevLTime.size();
    vector<tuple<string, double, double> > stuffs;
    stuffs.push_back(make_tuple("AStar", aD, aT));
    stuffs.push_back(make_tuple("POLR", pD, pT));
    stuffs.push_back(make_tuple("LeastDev - length" , llD, llT));
    sort(stuffs.begin(), stuffs.end(), [](tuple<string, double, double> a, tuple<string, double, double> b) { return (get<1>(a) < get<1>(b));});
    for(const auto tp: stuffs)
    {
        cout << get<0>(tp) << " ran  in an average time of " << get<2>(tp) << " creating an average path length of " << get<1>(tp) << endl;
    }


    return 0;
}
