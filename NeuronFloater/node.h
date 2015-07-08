#ifndef NODE_H
#define NODE_H

#include <vector>

class Node
{
public:

    float x;
    float y;
    float z;
    float r;
    int type;

    std::vector<int> nbr = std::vector<int>(); // list of neighbouring nodes

    // types as in neuromorpho.org description
    static int NOTHING;
    static int SOMA;
    static int AXON;
    static int BASAL_DENDRITE;
    static int APICAL_DENDRITE;
    static int FORK;
    static int END;
    static int UNDEFINED;

    Node(float xn, float yn, float zn, float rn);
    Node(float xn, float yn, float zn, float rn, int typ);
    Node(float xn, float yn, float zn, float rn, int typ, int pred);

    ~Node();

};

#endif // NODE_H
