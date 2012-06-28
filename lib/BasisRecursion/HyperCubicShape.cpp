#include <Eigen/Core>
//#include "boost/python.hpp"
#include "HyperCubicShape.h"
#include "convenienceFunctions.h"
#include <boost/python.hpp>


//
// boost::python stuff
//

#ifdef PYTHONMODULE
using namespace boost::python;
BOOST_PYTHON_MODULE(HyperCubicShape) {
// need init call here, or bp will assume a default constructor exists!
class_<HyperCubicShape<PyIndexIterator> >("HyperCubicShape",init<size_t,boost::python::tuple,boost::python::dict,boost::python::dict>())
    //.def(init<>())  // default constructor
    //.def(init<HyperCubicShape>()) // copy constructor
    .def(init<size_t,tuple,boost::python::dict,boost::python::dict>()) // see above!
    .def("contains", &HyperCubicShape<PyIndexIterator>::contains_py)
    .def("get_neighbours", &HyperCubicShape<PyIndexIterator>::get_neighbours)
    //.def("get_index_iterator_chain", &HyperCubicShape::get_index_iterator_chain)
    .def("__iter__",iterator<HyperCubicShape<PyIndexIterator> >())
    ;
}

#else // no python module (PYTHONMODULE not defined)
#include <iostream>
#include <vector>
#include <string>
int main(int argc, const char *argv[])
{
    Py_Initialize();
    std::cout << "creating tuple" << std::endl;
    // create variables
    boost::python::tuple limtuple = boost::python::make_tuple(2,3); // 0,1 in each dim
    std::cout << "tuple created" << std::endl;
    //std::cout << limtuple[0] << std::endl;
    boost::python::list indices;
    //Eigen::VectorXi a,b,c,d;
    boost::python::tuple a,b,c,d,e,f;
    std::cout << "creating indices" << std::endl;
    //a << 0,0; indices.append(a);
    a = boost::python::make_tuple(0,0); indices.append(a);
    b = boost::python::make_tuple(1,0); indices.append(b);
    c = boost::python::make_tuple(0,1); indices.append(c);
    d = boost::python::make_tuple(1,1); indices.append(d);
    e = boost::python::make_tuple(0,2); indices.append(e);
    f = boost::python::make_tuple(1,2); indices.append(f);
    std::cout << "created indices" << std::endl;
    // example dicts
    std::cout << "creating lima, lima_inv dicts" << std::endl;
    boost::python::dict lima,lima_inv;
    lima[a] = 0; lima_inv[0] = a;
    lima[b] = 1; lima_inv[1] = b;
    lima[c] = 2; lima_inv[2] = c;
    lima[d] = 3; lima_inv[3] = d;
    lima[e] = 3; lima_inv[3] = e;
    lima[f] = 3; lima_inv[3] = f;
    std::cout << "created lima, lima_inv dicts" << std::endl;
    std::cout << "instantiating HCS" << std::endl;
    HyperCubicShape<EigIndexIterator> hc(2,limtuple,lima,lima_inv);
    std::cout << "getting neighbours of 1,1" << std::endl;
    Eigen::VectorXi vec(2);
    vec << 1,1;
    std::vector<Eigen::VectorXi> n;
    std::cout << "vec size: " << vec.size() << std::endl;
    try {
    n = hc.get_neighbours(vec,"all",0);
    } catch (...) {
        PyObject *ptype, *pvalue, *ptraceback;
        PyErr_Fetch(&ptype, &pvalue, &ptraceback);
        std::string error = boost::python::extract<std::string>(pvalue);
        std::cout << "ERROR: " << error << std::endl;
    }
    for (size_t i=0; i < n.size(); i++)
        std::cout << "i: " << i << " n.T= " << n[i].transpose() << std::endl;
    std::cout << "getting HCS iterator" << std::endl;
    //HyperCubicShape<EigIndexIterator>::iterator it; // doesn't work- no default constructor (TODO?)
    HyperCubicShape<EigIndexIterator>::iterator it = hc.get_index_iterator_chain(1);
    std::cout << "got HCS iterator" << std::endl;
    for (it = hc.begin(); it != hc.end(); it++)
        std::cout << "index: " << (*it).transpose() << std::endl;

    return 0;
}

#endif

