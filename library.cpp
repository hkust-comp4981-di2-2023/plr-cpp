#include "library.h"

#include <iostream>

int main() {
    GreedyPLR<uint64_t, double> s {0.5f};
    auto a = s.process(Point<double>(0,0));
    return 0;
}