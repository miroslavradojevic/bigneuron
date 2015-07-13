#include <iostream>
#include <cstdlib>
#include <stdio.h>

using namespace std;

//////////////////////////////

//template<typename T>
//void    PrintArray(T* array, int n);

template<typename T>
void    quicksort(T* array, int* idxs, int startIndex, int endIndex);

template<typename T>
int     quicksort_partition(T* array, int* indxs, T pivot, int startIndex, int endIndex);

template<typename T>
void swp(T &a, T &b) {T temp; temp = a; a = b; b = temp;}

template<typename T>
void quickselect(T* G, int* Gidxs, int first, int last, int k, T &Gout, int &Gidx);

template<typename T>
int partition(T* G, int* Gidxs, int first, int last);

template<typename T>
void getKhighest(T* G, int* Gidxs, int Glen, int K, T* Kvals, int* Kidxs);

//////////////////////////////

//template<typename T>
//void PrintArray(T* arr, int n) {
//    for(int i = 0; i < n; i++) cout<<arr[i]<<'\t';
//    cout<<""<<endl;
//}

template<typename T>
void quickselect(T* G, int* Gidxs, int first, int last, int k, T &Gout, int &Gidx) {

    if (first <= last) {

        int pivot = partition(G, Gidxs, first, last);

        if (pivot == k) {
            Gout=G[k];
            Gidx=Gidxs[k];
        }
        else if (pivot > k) {
            quickselect(G, Gidxs, first,    pivot-1,    k, Gout, Gidx);
        }
        else {
            quickselect(G, Gidxs, pivot+1,  last,       k, Gout, Gidx);
        }

    }
    else { // illegal case
        Gout = -1;
        Gidx = -1;
    }

}

template<typename T>
int partition(T* G, int* Gidxs, int first, int last) {

    int pivot = first + rand()%(last - first + 1); // random number between 0 and (last - first + 1) excluded

    swp(G[last], G[pivot]);
    swp(Gidxs[last], Gidxs[pivot]);

    for (int i = first; i < last; i++) {
        if (G[i] > G[last]) {

            swp(G[i], G[first]);
            swp(Gidxs[i],Gidxs[first]);

            first++;
        }
    }

    swp(G[first], G[last]);
    swp(Gidxs[first], Gidxs[last]);

    return first;
}

template<typename T>
void quicksort(T* array, int* idxs, int startIndex, int endIndex) {

    // pivot element - value and index
    T pivot = array[startIndex];
    int pivot_idx = idxs[startIndex];

    int splitPoint;
    if(endIndex > startIndex) {

        splitPoint = quicksort_partition(array, idxs, pivot, startIndex, endIndex);    // returns the pivot position

        array[splitPoint] = pivot;
        idxs[splitPoint] = pivot_idx;

        quicksort(array, idxs, startIndex, splitPoint-1);                     // Quick sort first half
        quicksort(array, idxs, splitPoint+1, endIndex);                       // Quick sort second half
    }

}

template<typename T>
int quicksort_partition(T* array, int* indxs,    T pivot,    int startIndex, int endIndex) {

    int lBdry = startIndex;
    int rBdry = endIndex;

    while(lBdry < rBdry) {

        while( pivot > array[rBdry] && rBdry > lBdry) rBdry--;

        swp(array[lBdry], array[rBdry]);
        swp(indxs[lBdry], indxs[rBdry]);

        while( pivot <= array[lBdry] && lBdry < rBdry) lBdry++;

        swp(array[rBdry], array[lBdry]);
        swp(indxs[rBdry], indxs[lBdry]);

    }

    return lBdry;

}

template<typename T>
void getKhighest(T* G, int* Gidxs, int Glen, int K, T* Kvals, int* Kidxs) {

    // quick select Kth highest -> pivot, pivot_idx
    T pivot;
    int pivot_idx;
    quickselect(G, Gidxs, 0, Glen-1, K, pivot, pivot_idx); // quick select K-th highest

    cout<<"--------->\t"<<K<<"th highest: pivot "<<"["<<pivot_idx<<"] "<<pivot<<endl;

    // use pivot, pivot_idx to partition the array around the pivot value
    // (this stage might not add anything - thet are partitioned already, just pick top K from G)
    // it will additionally swap values in G and Gidxs... probably redundant
    int splitP = quicksort_partition(G, Gidxs, pivot, 0, Glen-1);
    swp(G[splitP],  G[0]);  // 0 was the left index
    swp(Gidxs[splitP], Gidxs[0]);

    for (int kk = 0; kk < K; ++kk) {
        Kvals[kk] = G[kk];
        Kidxs[kk] = Gidxs[kk];
    }

}

#define SZ 20
#define K 10
#define VAL_LIMIT 100
#define PRINT_MARGIN 10

int main(void) {

    srand ( time(NULL) );

    float arr1[SZ];
    int arr1i[SZ];
    for (int i = 0; i < SZ; ++i) {
        arr1[i] = rand()%VAL_LIMIT;
        arr1i[i] = i;
    }

    cout<<"\noriginal"<<endl;
    for (int i = 0; i < SZ; ++i) {
        printf("[%3d] %4.0f\t", arr1i[i], arr1[i]);
        if ((i+1)%PRINT_MARGIN==0) cout<<endl;
    }
    cout<<endl;

    // test quicksort
    float arr2[SZ];
    int arr2i[SZ];
    for (int i = 0; i < SZ; ++i) {arr2[i] = arr1[i]; arr2i[i]=i;}
    quicksort(arr2, arr2i, 0, SZ-1); // quicksort
    cout<<"\n---\nquicksort:"<<endl;
    for (int i = 0; i < SZ; ++i) {
        printf("[%3d] %4.0f\t", arr2i[i], arr2[i]);
        if ((i+1)%PRINT_MARGIN==0) cout<<endl;
    }
    cout<<"---\n"<<endl;

    // test find top K values with indexes
    float arr3[SZ]; int arr3i[SZ]; // copy for test
    for (int i = 0; i < SZ; ++i) {arr3[i] = arr1[i]; arr3i[i]=i;}
    float kvals[K]; int kidxs[K];

    getKhighest(arr3, arr3i, SZ, K, kvals, kidxs);

    cout<<"\n***\nK highest:"<<endl;
    for (int i = 0; i < K; ++i) {
        printf("[%3d] %4.0f\t", kidxs[i], kvals[i]);
        if ((i+1)%PRINT_MARGIN==0) cout<<endl;
    }
    cout<<"\n=================================================================="<<endl;

    // quick select K-th highest
    float pivot;
    int pivot_idx;
    cout<<"\nquickselect"<<endl;
    quickselect(arr1, arr1i, 0, SZ-1, K, pivot, pivot_idx); // find pivot and it's index as k-th largest

    for (int i = 0; i < SZ; ++i) {
        printf("[%3d] %4.0f\t", arr1i[i], arr1[i]);
        if ((i+1)%PRINT_MARGIN==0) cout<<endl;
    }
    cout<<endl;

    cout<<"--------->\t"<<K<<"th highest: pivot "<<"["<<pivot_idx<<"] "<<pivot<<endl;

    // use pivot to partition the array
    int splitP = quicksort_partition(arr1, arr1i, pivot, 0, SZ-1);
    swp(arr1[splitP],  arr1[0]);  // 0 was the left index
    swp(arr1i[splitP], arr1i[0]);

    cout<<"\npartition:" <<endl;
    for (int i = 0; i < SZ; ++i) {
        printf("[%3d] %4.0f\t", arr1i[i], arr1[i]);
        if ((i+1)%PRINT_MARGIN==0) cout<<endl;
    }
    cout<<endl;

    // top K
    int topKidxs[K];
    float topKvals[K];
    cout << "\n***\nK highest:" << endl;
    for (int kk = 0; kk < K; ++kk) {
        topKvals[kk] = arr1[kk];
        topKidxs[kk] = arr1i[kk];
        printf("[%3d] %4.0f\t", topKidxs[kk], topKvals[kk]); if ((kk+1)%PRINT_MARGIN==0) cout<<endl;
    }

}
