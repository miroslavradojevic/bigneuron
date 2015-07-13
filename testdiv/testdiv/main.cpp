#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <queue>

using namespace std;

class Body {
public:
    float height;
    Body(float h);
    Body();
    ~Body();
};

static bool compare(Body b1, Body b2){
    return (b1.height>b2.height);
}

Body::Body(float h) {height=h;}
Body::Body(){height=0;}
Body::~Body(){}

// possible with c++11 lambdas using 1 method only - this way it's class + operator
class CompareIndicesByNodeVectorCorrValues {
    std::vector<Body>* _bdies;
public:
    CompareIndicesByNodeVectorCorrValues(std::vector<Body>* values) : _bdies(values) {}
public:
    bool operator() (const int& a, const int& b) const { return (*_bdies)[a].height > (*_bdies)[b].height; }
};

template<typename T>
class BfsQueue // queue used for the breadth first search
{
public:
    std::queue<T> kew;

    BfsQueue() {}

    void enqueue(T item) {kew.push(item);}

    T dequeue() {T output = kew.front(); kew.pop(); return output;}

    int size() {return kew.size();}

    bool hasItems(){return !kew.empty();}
};

vector<int> exampleF(int in){

    vector<int> out;
    out.clear();

    if (in==0) {out.push_back(-1);}
    else if (in==1) {out.push_back(-2);}
    else {out.push_back(1); out.push_back(2);}

    return out;

}

int main()
{

    cout << "test return funciton"<<endl;

    cout << (exampleF(1)[0]==-2) << endl;

    if (1) return 0;

    cout << "count number of falses in std::vector<bool>" << endl;



    cout << " find index of value in std::vector" << endl;
    vector<int> vec1;
    vec1.push_back(5);
    vec1.push_back(9);
    vec1.push_back(7);
    vec1.push_back(1);

    cout << vec1.size() << endl;

    int pos = find(vec1.begin(), vec1.end(), 9) - vec1.begin();

    cout << pos << endl;

    if (1) return 0;

    cout << "test array of std::vector" << endl;
    vector<int> testarrayvec[10];
    testarrayvec[0].push_back(3);
    testarrayvec[0].push_back(4);
    testarrayvec[0].push_back(8);

    testarrayvec[5].push_back(9);
    testarrayvec[5].push_back(99);

    //
    for (int i = 0; i < 10; ++i) {
        cout << i << " : \t\t";
        for (int j = 0; j < testarrayvec[i].size(); ++j) {
            cout << testarrayvec[i][j] << "\t";
        }
        cout << endl;
    }

    if (1) return 0;

    vector<int> tstvec;
    tstvec.push_back(4);
    cout << tstvec.size() << " elements" << endl;
    tstvec.push_back(8);
    cout << tstvec.size() << " elements" << endl;
    tstvec.push_back(2);
    cout << tstvec.size() << " elements" << endl;
    if (1) return 0;

    srand(time(NULL)); // randomization initialize
    cout << "testing queue (fifo) " << endl;

    BfsQueue< vector<int> > boob; // initialize the queue

    for (int var = 0; var < 10; ++var) { // push 10 vector<int> elements (integer arrays)
        vector<int> lnk(2);
        lnk[0] = round(rand()/(RAND_MAX/50));
        lnk[1] = round(rand()/(RAND_MAX/50));
        cout << " ->(" << lnk[0] << ", " << lnk[1] << ")<- \t" << flush;
        boob.enqueue(lnk);
    }

    cout << "\nqueue size " << boob.size() << "\n-----\n" << endl;

    for (int var = 0; var < 5; ++var) { // dequeue 5 elements
        vector<int> takeit = boob.dequeue();
        cout << "[" << takeit[0] << ", " << takeit[1] << "]\t" << flush;
    }

    cout << endl;
    for (int var = 0; var < 10; ++var) { // push 10 more vector<int> elements (integer arrays)
        vector<int> lnk(2);
        lnk[0] = round(rand()/(RAND_MAX/50));
        lnk[1] = round(rand()/(RAND_MAX/50));
        cout << " ->(" << lnk[0] << ", " << lnk[1] << ")<- \t" << flush;
        boob.enqueue(lnk);
    }

    cout << endl;
    while(boob.hasItems()) {
        vector<int> takeit = boob.dequeue();
        cout << "[" << takeit[0] << ", " << takeit[1] << "]\t" << flush;
    }

    cout << endl;
    if (1) return 0;

    queue< vector<int> > my_queue;// = new queue<vector<int>>;

    for (int var = 0; var < 10; ++var) {
        vector<int> lnk(2);
        lnk[0] = round(rand()/(RAND_MAX/50));// * 50;
        lnk[1] = round(rand()/(RAND_MAX/50));// * 50;
        cout << var << " -push- " << lnk[0] << ", " << lnk[1] << endl;
        my_queue.push(lnk);
    }

    cout << "---" << endl;

//    cout << my_queue.front()[0] << ", " << my_queue.front()[1] << endl;
//    vector<int> aa = my_queue.pop();
//    cout << my_queue.front()[0] << ", " << my_queue.front()[1] << endl;
//    my_queue.pop();
//    cout << my_queue.front()[0] << ", " << my_queue.front()[1] << endl;
//    cout << my_queue.back()[0]  << ", " << my_queue.back()[1] << endl;
//    cout << my_queue.size() << " elements" << endl;
//    for (int var = 0; var < my_queue.size(); ++var) {
//        my_queue.
//
//        cout << my_queue.back()[0]  << ", " << my_queue.back()[1] << endl;
////
//        cout << " --- " << endl;
//    }

    cout << "finished." << endl;

    if (1) return 0;


    std::vector<Body> nodl = std::vector<Body>();

    for (int ii = 0; ii < 20; ++ii) {
        Body bb(rand()/(RAND_MAX/50.0));
        nodl.push_back(bb);
    }

    std::vector<int> indices(nodl.size());
//    std::iota(indices.begin(), indices.end(), 0); // c++11

    cout << "before:"<<endl;
    for (int i = 0; i < nodl.size(); ++i) {
        indices.at(i) = i;
        cout << indices.at(i) << " : " << nodl.at(indices.at(i)).height << endl;
    }

    std::sort(indices.begin(), indices.end(), CompareIndicesByNodeVectorCorrValues(&nodl));

    cout << "after:"<<endl;
    for (int i = 0; i < indices.size(); ++i) {
        cout << indices.at(i) << " : " << nodl.at(indices.at(i)).height << endl;
    }

    std::vector<Body> tt(10); // sort vector of classes example

    cout << "before\t";
    for (int i = 0; i < 10; ++i) {
        float hh = rand()/(RAND_MAX/50.0);
        Body bdy(hh);
        tt.at(i).height = hh; // tt.push_back(bdy);
        cout << tt.at(i).height << "\t";
    }
    cout << endl;

    //std::sort(tt.begin(), tt.end(), compare);
//    for (auto i: sort_indexes(v)) {
//      cout << v[i] << endl;
//    }

    cout << "after\t";
    for (int i = 0; i < 10; ++i) {
        cout << tt.at(i).height << "\t";
    }
    cout << endl;


    cout << "some test" << endl;

    float wtst = -FLT_MAX;
    wtst = 3.4;

    cout << "is it equal? " << (wtst==-FLT_MAX) << endl;

    int a = 5;
    int b = 9;
    int c = a * b * b * b;

    cout << "test..." << endl;
    for (int i = 0; i < c; ++i) {
        if (c/b/b != c/(b*b)) cout << "happened!" << endl;
    }
    cout << "done" << endl;


    cout << "one more test" << endl;

    int r0 = (int)round(pow(a,0));
    int r1 = (int)round(pow(a,2));

    b = 6;

    cout << b/r0 << "\t" << b/r1 << endl;

    return 0;
}

