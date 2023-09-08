#include <stdbool.h>

#include "epsilon-equal.h"

#include "point.h"



/* Operations on points with real coordinates.  It is not orthogonal,
   but more convenient, to have the subtraction operator return a
   vector, and the addition operator return a point.
*/



Point
point_make(float const x,
           float const y,
           float const z) {

    Point retval;

    retval.x = x;
    retval.y = y;
    retval.z = z;

    return retval;
}



bool
point_equal(Point const comparand,
            Point const comparator) {

    return
        epsilon_equal(comparand.x, comparator.x)
        &&
        epsilon_equal(comparand.y, comparator.y)
        &&
        epsilon_equal(comparand.z, comparator.z)
        ;
}



Point
point_sum(Point const coord1,
          Point const coord2) {

    Point retval;

    retval.x = coord1.x + coord2.x;
    retval.y = coord1.y + coord2.y;
    retval.z = coord1.z + coord2.z;

    return retval;
}



Point
point_scaled(Point const coord,
             float const r) {

    Point retval;

    retval.x = coord.x * r;
    retval.y = coord.y * r;
    retval.z = coord.z * r;

    return retval;
}


