#include "node.h"

int Node::NOTHING = 0;
int Node::SOMA = 1;
int Node::AXON = 2;
int Node::BASAL_DENDRITE = 3;
int Node::APICAL_DENDRITE = 4;
int Node::FORK = 5;
int Node::END = 6;
int Node::UNDEFINED = 7;

Node::Node(float xn, float yn, float zn, float rn)
{
    x = xn;
    y = yn;
    z = zn;
    r = rn;
    type = NOTHING;
//    nbr = std::vector<int>();
}

Node::Node(float xn, float yn, float zn, float rn, int typ)
{
    x = xn;
    y = yn;
    z = zn;
    r = rn;
    type = typ;
//    nbr = std::vector<int>();
}

Node::Node(float xn, float yn, float zn, float rn, int typ, int pred)
{
    x = xn;
    y = yn;
    z = zn;
    r = rn;
    type = typ;
//    nbr = std::vector<int>();
    nbr.push_back(pred);
}

Node::~Node()
{
}
