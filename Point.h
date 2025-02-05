#pragma once

#include <iostream>

using namespace std;

class Point {
    public:
        double x, y;

        long locmm;

        Point();
        Point(double _x, double _y, long locmm );
};