#ifndef _POLR_GRAPH_H
#define _POLR_GRAPH_H
#include <queue>
#include <functional>
#include <map>
#include <set>
#include <stack>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <future>


//try creating a PF alg that takes the straight line dist between the two, creates a line equation, and makes shortest path based on least deviation from that line 


using namespace std;


template <class Object, class CostType>
class Vertex
{
    typedef Vertex<Object, CostType>* VertPtr;
    typedef std::map<VertPtr, CostType> AdjacencyList;
    typedef Vertex<Object, CostType> Vert;

    public:
        static constexpr CostType INFINITY_VT = (CostType)UINT64_MAX;
        AdjacencyList adjList;
        AdjacencyList revAdjList;
        Object data;
        VertPtr forwardsParent, backwardsParent;
        CostType totalCost, estimatedCost, dist, x, y;
        double error;
        Vertex(const Object & obj = Object(), CostType x = CostType(), CostType y = CostType());
        void addToAdjList(VertPtr neighbor, CostType cost);
        void addToRevAdjList(VertPtr neighbor, CostType cost);
        bool operator<(const Vert& rhs) const;
        bool operator>(const Vert& rhs) const;
        const Vertex<Object, CostType>& operator=(const Vert& rhs);
        void showAdjList();

        static int nSortKey;
        static stack<int> keyStack;
        static enum { SORT_BY_DATA, SORT_BY_DIST, SORT_BY_ERR } eSortType;
        static bool setNSortType( int whichType );
        static void pushSortKey() { keyStack.push(nSortKey); }
        static void popSortKey() { nSortKey = keyStack.top(); keyStack.pop(); }
};

template <class Object, typename CostType>
int Vertex<Object, CostType>::nSortKey = Vertex<Object, CostType>::SORT_BY_DATA;

template <class Object, typename CostType>
stack<int> Vertex<Object, CostType>::keyStack;

template <class Object, class CostType>
Vertex<Object, CostType>::Vertex(const Object& obj, CostType x, CostType y): data(obj), x(x), y(y), forwardsParent(NULL), backwardsParent(NULL), dist(INFINITY_VT), estimatedCost(INFINITY_VT), error((double)INFINITY_VT)
{

}



template <class Object, typename CostType>
bool Vertex<Object, CostType>::setNSortType( int whichType )
{
   switch (whichType)
   {
   case SORT_BY_DATA:
   case SORT_BY_DIST:
   case SORT_BY_ERR:
      nSortKey = whichType;
      return true;
   default:
      return false;
   }
}





template <class Object, class CostType>
void Vertex<Object, CostType>::addToAdjList(VertPtr neighbor, CostType cost)
{
    adjList[neighbor] = cost;
}


template <class Object, class CostType>
void Vertex<Object, CostType>::addToRevAdjList(VertPtr neighbor, CostType cost)
{
    revAdjList[neighbor] = cost;
}



template <class Object, class CostType>
bool Vertex<Object, CostType>::operator<(const Vert& rhs) const
{
    switch(nSortKey)
    {
        case SORT_BY_DATA:
            return (data < rhs.data);
        case SORT_BY_DIST:
            return (totalCost < rhs.totalCost);//TODO: maybe edit, depends on how we want make_heap to work
        case SORT_BY_ERR:
            return (error < rhs.error);
        default:
            return (data < rhs.data);
    }
}



template <class Object, class CostType>
bool Vertex<Object, CostType>::operator>(const Vert& rhs) const
{
    switch(nSortKey)
    {
        case SORT_BY_DATA:
            return (data > rhs.data);
        case SORT_BY_DIST:
            return (totalCost > rhs.totalCost);//TODO: maybe edit, depends on how we want make_heap to work
        case SORT_BY_ERR:
            return (error > rhs.error);
        default:
            return (data > rhs.data);
    }

}



template <class Object, class CostType>
const Vertex<Object, CostType>& Vertex<Object, CostType>::operator=(const Vert& rhs)
{   
    adjList = rhs.adjList;
    revAdjList = rhs.revAdjList;
    data = rhs.data;
    forwardsParent = rhs.forwardsParent;
    backwardsParent = rhs.backwardsParent;
    totalCost = rhs.totalCost;
    estimatedCost = rhs.estimatedCost;
    dist = rhs.dist;
    x = rhs.x;
    y = rhs.y;
    return *this;
}




template <class Object, class CostType>
class Graph
{
    typedef Vertex<Object, CostType> Vert;
    typedef Vertex<Object, CostType>* VertPtr;
    typedef map<VertPtr, CostType> AdjacencyList;
    typedef set<Vert> VertexSet;
    typedef set<VertPtr> VertPtrSet;

    private:
        VertPtrSet vertPtrSet;
        VertexSet vertexSet;
        VertPtr mid;
    public:
        Graph(){}
        Vertex<Object, CostType>* addNode(const Object& object, CostType x, CostType y);
        void addEdge(const Object& source, const Object& dest, CostType cost);
        void addEdge1(const Object& source, const Object& dest);//automatically calculates distance formula
        VertPtr addToVertexSet(const Object& object);
        void showAdjTable();
        VertPtrSet getVertPtrSet() const {return vertPtrSet;}
        void clear();
        vector<VertPtr> createPath(const Object& source, const Object& dest);//TODO: find a way to elegantly fail if path doesn't exist
        vector<VertPtr> createPathStraightLine1(const Object& source, const Object& dest);//pretty damn good shit right here 
        void resetVertices();
        vector<VertPtr> aStar(const Object& source, const Object& dest);

    private:
        static CostType sqDist(VertPtr first, VertPtr sec){ CostType x = first->x - sec->x; CostType y = first->y - sec->y; return ((x*x) + (y*y));}
        VertPtr getVertexWithThisData(const Object& x);
        CostType findMinCost(VertPtr v, vector<VertPtr>& oppFront);
        vector<VertPtr> getPath();
        vector<VertPtr> getAStarPath(VertPtr e);
        double getDeviation(pair<double, double>& ln, VertPtr Vert, VertPtr s);
        double theoreticalDist(pair<double, double>& ln, CostType x, VertPtr s);
        pair<double, double> createLine(VertPtr f, VertPtr s);
        double getDeviation1(pair<double, double>& ln, VertPtr Vert, VertPtr s);
        double theoreticalDist1(pair<double, double>& ln, CostType x, VertPtr s);
};

template <class Object, class CostType>
Vertex<Object, CostType>* Graph<Object, CostType>::addNode(const Object& object, CostType x, CostType y)
{
    pair<typename VertexSet::iterator, bool> retVal;
    VertPtr vPtr;
    Vert::pushSortKey();
    Vert::setNSortType(Vert::SORT_BY_DATA);
    retVal = vertexSet.insert( Vert(object, x, y));;
    vPtr = (VertPtr)&*retVal.first;
    vertPtrSet.insert(vPtr);
    Vert::popSortKey();
    return vPtr;
}

template <class Object, class CostType>
CostType Graph<Object, CostType>::findMinCost(VertPtr v, vector<VertPtr>& oppFront)
{
    CostType minDist = Vert::INFINITY_VT;
    for(int i = 0; i < oppFront.size(); i++)
    {   
        minDist = min(minDist, sqDist(v, oppFront[i]));
    }
    return minDist;
}



template <class Object, class CostType>
void Graph<Object, CostType>::addEdge(const Object& source, const Object& dest, CostType cost)
{
    VertPtr src, dst;

    // put both source and dest into vertex list(s) if not already there
    src = addToVertexSet(source);
    dst = addToVertexSet(dest);

    // add dest to source's adjacency list
    src->addToAdjList(dst, cost); 
    dst->addToRevAdjList(src, cost); 
}
template <class Object, class CostType>
void Graph<Object, CostType>::addEdge1(const Object& source, const Object& dest)
{
    VertPtr src, dst;
    

    // put both source and dest into vertex list(s) if not already there
    src = addToVertexSet(source);
    dst = addToVertexSet(dest);
    auto cost = sqrt(sqDist(src, dst));

    // add dest to source's adjacency list
    src->addToAdjList(dst, cost); 
    dst->addToRevAdjList(src, cost); 
}

template <class Object, class CostType>
Vertex<Object, CostType>* Graph<Object, CostType>::addToVertexSet(const Object& object)
{
    pair<typename VertexSet::iterator, bool> retVal;
    VertPtr vPtr; 
    Vert::pushSortKey();
    Vert::setNSortType(Vert::SORT_BY_DATA);
    retVal = vertexSet.insert( Vert(object) );

    // get pointer to this vertex and put into vert pointer list
    vPtr = (VertPtr)&*retVal.first;
    vertPtrSet.insert(vPtr);
    Vert::popSortKey();
    return vPtr;
}


template <class Object, class CostType>
void Graph<Object, CostType>::clear()
{
    vertexSet.clear();
    vertPtrSet.clear();
}


template <class Object, class CostType>
Vertex<Object, CostType>* Graph<Object, CostType>::getVertexWithThisData(const Object& x)
{
    //TODO: add back in key mechanics so we can find the object based on data not totalDistance
    typename VertexSet::iterator findResult;
    Vert vert(x);

    Vert::pushSortKey();
    Vert::setNSortType(Vert::SORT_BY_DATA);
    findResult = vertexSet.find(vert);
    Vert::popSortKey();

    if ( findResult == vertexSet.end() )
        return NULL;
    return  (VertPtr)&*findResult;
}




template <class Object, class CostType>
vector<Vertex<Object, CostType>* > Graph<Object, CostType>::createPath(const Object& source, const Object& dest)
{
    //TODO: maybe make another priority_queue or revise algorithm to also take into account deviation from straight line distance between the two points
    //auto 
    auto rst = std::async(launch::async, [this]()mutable {resetVertices(); return true;});
    auto cmp = [](VertPtr a, VertPtr b){ return (*a > *b);};
    vector<VertPtr> qC, rC;
    vector<VertPtr> fF, bF;
    set<VertPtr> fE, bE;
    vector<VertPtr> intersect(1);//set to size of 1
    auto it = intersect.begin();
    rst.get();
    VertPtr s = getVertexWithThisData(source);
    s->dist = 0;
    VertPtr e = getVertexWithThisData(dest);
    e->dist = 0;
    VertPtr q = s;
    VertPtr r = e;



    fF.push_back(s);
    bF.push_back(e);
    Vert::pushSortKey();
    Vert::setNSortType(Vert::SORT_BY_DIST);
    while(distance(intersect.begin(), it) == 0)
    {
        /*
        Vert::pushSortKey();
        Vert::setNSortType(Vert::SORT_BY_DATA);
        if((fE.find(e) != fE.end()) && (bE.find(s) != bE.end()))
        {
            vector<VertPtr> v{s};
            return v;
        }
        Vert::popSortKey();
        */
        q = fF.front();
        pop_heap(fF.begin(), fF.end(), cmp);
        fF.pop_back();
        r = bF.front();
        pop_heap(bF.begin(), bF.end(), cmp);
        bF.pop_back();

        for(auto& pr: q->adjList)
        {
            
            auto vt = pr.first;
            auto val = std::async([this, &vt, &bF]() mutable {return findMinCost(vt, bF);});
            vt->dist = q->dist + pr.second;
            vt->forwardsParent = q;
            vt->estimatedCost = val.get();//findMinCost(vt, fF);
            vt->totalCost = vt->dist + vt->estimatedCost;
            fF.push_back(vt);
            push_heap(fF.begin(), fF.end(), cmp); 
            
        }
        for(auto& pr: r->revAdjList)
        {
            
            auto vt = pr.first;
            auto val = std::async([this, &vt, &fF]() mutable {return findMinCost(vt, fF);});
            vt->dist = r->dist + pr.second;
            vt->backwardsParent = r;
            vt->estimatedCost = val.get();//findMinCost(vt, fF);
            vt->totalCost = vt->dist + vt->estimatedCost;
            rC.push_back(vt);
            bF.push_back(vt);
            push_heap(bF.begin(), bF.end(), cmp); 
            
        }

        for(auto& vt: fF)
        {
            vt->estimatedCost = min(vt->estimatedCost, findMinCost(vt, rC));
            vt->totalCost = vt->estimatedCost + vt->dist;
        }
        qC.clear();
        
        make_heap(fF.begin(), fF.end(), cmp);
        Vert::pushSortKey();
        Vert::setNSortType(Vert::SORT_BY_DATA);
        fE.insert(q);
        bE.insert(r);
        it = set_intersection(fE.begin(), fE.end(), bE.begin(), bE.end(), intersect.begin());
        Vert::popSortKey();
    }

    mid = intersect.front();
    return getPath();
}
template <class Object, class CostType>
vector<Vertex<Object, CostType>* > Graph<Object, CostType>::getPath()
{
    stack<VertPtr> bk;
    VertPtr it = mid;
    while(it->forwardsParent != NULL)
    {
        bk.push((it->forwardsParent));
        it = it->forwardsParent;
    }
    vector<VertPtr> ret;
    while(!bk.empty())
    {
        ret.push_back(bk.top());
        bk.pop();
    }
    ret.push_back(mid);
    it = mid;
    while(it->backwardsParent != NULL)
    {
        ret.push_back((it->backwardsParent));
        it = it->backwardsParent;
    }
    return ret;
}

template <class Object, class CostType>
void Graph<Object, CostType>::resetVertices()
{
    for(auto& v: vertPtrSet)
    {
        v->forwardsParent = NULL;
        v->backwardsParent = NULL;
        v->dist = Vert::INFINITY_VT;
        v->estimatedCost = Vert::INFINITY_VT;
        v->totalCost = Vert::INFINITY_VT;
    }
}

template <class Object, class CostType>

vector<Vertex<Object, CostType>* > Graph<Object, CostType>::aStar(const Object& source, const Object& dest)
{
    auto cmp = [](VertPtr a, VertPtr b){ return (*a > *b);};
    resetVertices();
    VertPtr s = getVertexWithThisData(source);
    s->dist = 0;
    VertPtr e = getVertexWithThisData(dest);
    VertPtr q = s;
    q->totalCost = 0;
    set<VertPtr> closed;
    vector<VertPtr> open;
    open.push_back(s);
    Vert::pushSortKey();
    Vert::setNSortType(Vert::SORT_BY_DIST);
    while(!open.empty())
    {
        q = open.front();
        pop_heap(open.begin(), open.end(), cmp);
        open.pop_back();
        for(auto& pr: q->adjList)
        {
            auto& nd = pr.first;
            if(nd == e)
            {
                nd->forwardsParent = q;
                return getAStarPath(e);
            }
            auto potentialDist = q->dist + pr.second;
            auto potentialCost = sqDist(nd, e);
            auto potentialTotal = potentialDist + potentialCost;
            Vert::pushSortKey();
            Vert::setNSortType(Vert::SORT_BY_DATA);
            auto it = find(closed.begin(), closed.end(), nd);
            Vert::popSortKey();
            if(it != closed.end())
            {
                if(potentialTotal < (*it)->totalCost)
                {
                    nd->dist = potentialDist;
                    nd->estimatedCost = potentialCost;
                    nd->totalCost = potentialTotal;
                    nd->forwardsParent = q;
                    open.push_back(nd);
                    push_heap(open.begin(), open.end(), cmp);
                }
            }
            else
            {
                nd->dist = potentialDist;
                nd->estimatedCost = potentialCost;
                nd->totalCost = potentialTotal;
                nd->forwardsParent = q;
                open.push_back(nd);
                push_heap(open.begin(), open.end(), cmp);
            }
            Vert::pushSortKey();
            Vert::setNSortType(Vert::SORT_BY_DATA);
            auto it1 = find(open.begin(), open.end(), nd);
            Vert::popSortKey();
            if(it1 != open.end())
            {
                if(potentialTotal < (*it1)->totalCost)
                {
                    nd->dist = potentialDist;
                    nd->estimatedCost = potentialCost;
                    nd->totalCost = potentialTotal;
                    nd->forwardsParent = q;
                    open.push_back(nd);
                    push_heap(open.begin(), open.end(), cmp);
                }
            }
            else
            {
                nd->dist = potentialDist;
                nd->estimatedCost = potentialCost;
                nd->totalCost = potentialTotal;
                nd->forwardsParent = q;
                open.push_back(nd);
                push_heap(open.begin(), open.end(), cmp);
            }
            Vert::pushSortKey();
            Vert::setNSortType(Vert::SORT_BY_DATA);
            closed.insert(q);
            Vert::popSortKey();
        }
        
    }
    return getAStarPath(e);
}


template <class Object, class CostType>
vector<Vertex<Object, CostType>* > Graph<Object, CostType>::getAStarPath(VertPtr e)
{
    stack<VertPtr> s;
    vector<VertPtr> v;
    VertPtr it = e;
    s.push(e);
    while(it->forwardsParent != NULL)
    {
        s.push((it->forwardsParent));
        it = it->forwardsParent;
    }
    while(!s.empty())
    {
        v.push_back(s.top());
        s.pop();
    }
    return v;
}

template <class Object, class CostType>
pair<double, double> Graph<Object, CostType>::createLine(VertPtr f, VertPtr s)
{
    double m = (f->y - s->y)/(f->x - s->x);
    double b = (f->y + static_cast<double>(m*(f->x)));
    return pair<double, double>(m, b);
}

template <class Object, class CostType>
double Graph<Object, CostType>::theoreticalDist(pair<double, double>& ln, CostType x, VertPtr s)
{
    //This is probably better in general
    /*
    double X = static_cast<double>(x);
    double xDif = (X - static_cast<double>(s->x));
    double funcY = static_cast<double>( (ln.first*(double)x + ln.second));
    double yDif = (funcY- static_cast<double>(s->y));
    return sqrt( (xDif*xDif) + (yDif*yDif));
    */
    //this is better when working with straight line distanecs between points
    double y = static_cast<double>(ln.first*(double)x + ln.second);
    return y;
}

template <class Object, class CostType>
double Graph<Object, CostType>::getDeviation(pair<double, double>& ln, VertPtr vert, VertPtr s)
{
    double theoretical = theoreticalDist(ln, vert->x, s);
    //dist is better in general, but y is better in straight line graphs
    double top = static_cast<double>(vert->y/*dist*/ - theoretical);
    return abs(top/theoretical);
}

template <class Object, class CostType>
double Graph<Object, CostType>::theoreticalDist1(pair<double, double>& ln, CostType x, VertPtr s)
{
    //This is probably better in general
    double X = static_cast<double>(x);
    double xDif = (X - static_cast<double>(s->x));
    double funcY = static_cast<double>( (ln.first*(double)x + ln.second));
    double yDif = (funcY- static_cast<double>(s->y));
    return sqrt( (xDif*xDif) + (yDif*yDif));
}

template <class Object, class CostType>
double Graph<Object, CostType>::getDeviation1(pair<double, double>& ln, VertPtr vert, VertPtr s)
{
    double theoretical = theoreticalDist1(ln, vert->x, s);
    //dist is better in general, but y is better in straight line graphs
    double top = static_cast<double>(vert->dist - theoretical);
    return abs(top/theoretical);
}

template <class Object, class CostType>
vector<Vertex<Object, CostType>* > Graph<Object, CostType>::createPathStraightLine1(const Object& source, const Object& dest)//workhorse right here
{
    auto cmp = [](VertPtr a, VertPtr b){return (*a > *b);}; 
    auto s = getVertexWithThisData(source);
    auto e = getVertexWithThisData(dest);
    resetVertices();
    auto ln = createLine(s, e);
    VertPtr q = s;
    //either make queue out of vertices with least deviation from main line at point, or least deviation of distance from start to vert as from start to point on line 
    vector<VertPtr> frontier(0);
    set<VertPtr> explored;
    frontier.push_back(s);
    Vert::pushSortKey();
    Vert::setNSortType(Vert::SORT_BY_ERR);
    while(!frontier.empty())
    {
    
        q = frontier.front();
        pop_heap(frontier.begin(), frontier.end(), cmp);
        frontier.pop_back();
        if(q == e)
        {
            break;
        }
        for(auto& pr: q->adjList)
        {
            auto& nd = pr.first;
            Vert::pushSortKey();
            Vert::setNSortType(Vert::SORT_BY_DATA);
            auto it1 = find(explored.begin(), explored.end(), nd);
            Vert::popSortKey();          
            auto newDist = q->dist + pr.second;
            auto dev = getDeviation1(ln, nd, s);
            if(it1 != explored.end())
            {
                if(newDist < nd->dist)
                {
                    nd->forwardsParent = q;
                    nd->dist = newDist;
                    nd->error = dev;
                    Vert::pushSortKey();
                    Vert::setNSortType(Vert::SORT_BY_DATA);
                    explored.erase(nd);
                    Vert::popSortKey();
                    frontier.push_back(nd);
                    push_heap(frontier.begin(), frontier.end(), cmp);
                }
            }
            else
            { 
                Vert::pushSortKey();
                Vert::setNSortType(Vert::SORT_BY_DATA);
                auto it2 = find(frontier.begin(), frontier.end(), nd);
                Vert::popSortKey();
                if(it2 != frontier.end())
                {
                    if(newDist < nd->dist)
                    {
                        nd->forwardsParent = q;
                        nd->dist = newDist;
                        nd->error = dev;
                        make_heap(frontier.begin(), frontier.end(), cmp);
                    }
                }
                else
                {
                    nd->forwardsParent = q;
                    nd->dist = newDist;
                    nd->error = dev;
                    frontier.push_back(nd);
                    push_heap(frontier.begin(), frontier.end(), cmp);
                }
            }
        }
        Vert::pushSortKey();
        Vert::setNSortType(Vert::SORT_BY_DATA);
        explored.insert(q);
        Vert::popSortKey();
    }
    Vert::popSortKey();
    return getAStarPath(e);
}
#endif
