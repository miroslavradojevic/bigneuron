//#include <QCoreApplication>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{

    cout << " testing..."<<endl;

    long N = 512;
    long M = 512;
    long P = 8;
    long size = N*M*P;

//    size = 55;

    cout << size << endl;

//    int t1[size];

    int * t2 = new int[size];

    for (long var = 0; var < size; ++var) {
        t2[var] = 0;
    }

    cout << "done." << endl;

}
