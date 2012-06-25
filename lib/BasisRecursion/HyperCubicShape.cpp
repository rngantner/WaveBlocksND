#include <Eigen/Core>
//#include "boost/python.hpp"
#include "HyperCubicShape.h"
#include "convenienceFunctions.h"
#include <string>
#include <vector>

#include <boost/python.hpp>
using namespace boost::python;

/*
#ifdef PYTHONMODULE
    #include <boost/python.hpp>
    using namespace boost::python;
    //#include <Python.h>
    //#include <numpy/arrayobject.h>
#endif
*/

// remove after debugging:
#include <iostream>

//
// HyperCubicShape
//

/**
 * Test if a tuple is contained in the HyperCubicShape.
 * TODO: optimize this! (is in inner loop! don't use python dict lookup)
 */
bool HyperCubicShape::contains(Eigen::VectorXi o){
    tuple op = toTuple(o);
    return _lima.has_key(op);
}

/**
 * Return list of neighbors of the index (vector) k.
 * \param k Eigen vector containing the index
 * \param selection Specifies whether to look forward, backward or give all neighbours
 * \param direction Gives the direction in which to look. If it is -1, return all possibilities.
 * \return std::vector of neighbours (Eigen::VectorXi instances)
 */
std::vector<Eigen::VectorXi>
HyperCubicShape::get_neighbours(Eigen::VectorXi k, std::string selection, int direction) {
    std::vector<Eigen::VectorXi> neighbours;
    // first look at all possibilities
    if (direction != -1){
        // only look in the given direction
        Eigen::VectorXi e(k);
        e.setZero();
        e[direction] = 1;
        if (selection == "backward" || selection == "all")
            neighbours.push_back(k-e);
        if (selection == "forward" || selection == "all")
            neighbours.push_back(k+e);
    } else {
        // look in all directions
        Eigen::VectorXi e(k);
        e.setZero();
        for (size_t d=0; d < this->D; d++) {
            e[d] = 1;
            if (selection == "backward" || selection == "all") {
                Eigen::VectorXi knew = k-e;
                if (contains(knew));
                    neighbours.push_back(knew);
            }
            if (selection == "forward" || selection == "all") {
                Eigen::VectorXi knew = k+e;
                if (contains(knew));
                    neighbours.push_back(knew);
            }
            e[d] = 0;
        }
    }
    return neighbours;
}


/**
 * Returns a new Chain Iterator
 */
HyperCubicShape::iterator* HyperCubicShape::get_index_iterator_chain(size_t direction)const {
    return new HyperCubicShape::iterator(_limits,D,direction);
}


//
// boost::python stuff
//

#ifdef PYTHONMODULE
BOOST_PYTHON_MODULE(HyperCubicShape) {
// need init call here, or bp will assume a default constructor exists!
class_<HyperCubicShape>("HyperCubicShape",init<size_t,tuple,dict,dict>())
    //.def(init<>())  // default constructor
    //.def(init<HyperCubicShape>()) // copy constructor
    //.def(init<size_t,tuple,dict,dict>()) // see above!
    .def("contains", &HyperCubicShape::contains)
    .def("get_neighbours", &HyperCubicShape::get_neighbours)
    //.def("get_index_iterator_chain", &HyperCubicShape::get_index_iterator_chain)
    .def("__iter__",iterator<HyperCubicShape>())
    ;
}

#else // no python module (PYTHONMODULE not defined)
int main(int argc, const char *argv[])
{
    Py_Initialize();
    std::cout << "creating tuple" << std::endl;
    // create variables
    tuple limtuple = make_tuple(2,2); // 0,1 in each dim
    std::cout << "tuple created" << std::endl;
    //std::cout << limtuple[0] << std::endl;
    list indices;
    //Eigen::VectorXi a,b,c,d;
    tuple a,b,c,d;
    std::cout << "creating indices" << std::endl;
    //a << 0,0; indices.append(a);
    a = make_tuple(0,0); indices.append(a);
    b = make_tuple(1,0); indices.append(b);
    c = make_tuple(0,1); indices.append(c);
    d = make_tuple(1,1); indices.append(d);
    std::cout << "created indices" << std::endl;
    // example dicts
    std::cout << "creating lima, lima_inv dicts" << std::endl;
    dict lima,lima_inv;
    try {
    lima[a] = 0; lima_inv[0] = a;
    lima[b] = 1; lima_inv[1] = b;
    lima[c] = 2; lima_inv[2] = c;
    lima[d] = 3; lima_inv[3] = d;
    } catch (...) {
        PyObject *ptype, *pvalue, *ptraceback;
        PyErr_Fetch(&ptype, &pvalue, &ptraceback);
        std::string error = extract<std::string>(pvalue);
        std::cout << "ERROR: " << error << std::endl;
    }
    std::cout << "created lima, lima_inv dicts" << std::endl;
    std::cout << "instantiating HCS" << std::endl;
    HyperCubicShape hc(2,limtuple,lima,lima_inv);
    std::cout << "getting HCS iterator" << std::endl;
    HyperCubicShape::iterator* it = hc.get_index_iterator_chain();
    std::cout << it->index << std::endl;
    return 0;
}

#endif

