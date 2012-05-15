#include <boost/python.hpp>
#include <Eigen/Core>
using namespace boost::python;

class HyperCubicShape;

//
// Iterator
//
template< class ValueType, class NodeType=Eigen::VectorXi >
struct IndexIterator : std::iterator< std::forward_iterator_tag, ValueType > {

    ValueType& operator*() { return index; }

    template< class VT2, class NT2 >
    friend bool operator==( IndexIterator const &lhs, IndexIterator< VT2, NT2 > const &rhs ){ return lhs.index == rhs.index; }
    friend bool operator!=( IndexIterator const &lhs, IndexIterator< VT2, NT2 > const &rhs ){ return lhs.index != rhs.index; }
    NodeType operator*(){ return index; }
    IndexIterator& operator++();
    IndexIterator operator++(int);

private:
    NodeType index;
    bool done;
    NodeType limits;

    friend class HyperCubicShape;
    IndexIterator(size_t dim, bool begin ); // private constructor for begin, end
};



//class myiterator : public iterator<input_iterator_tag, int> {
//  int* p;
//public:
//  myiterator(int* x) :p(x) {}
//  myiterator(const myiterator& mit) : p(mit.p) {}
//  myiterator& operator++() {++p;return *this;}
//  myiterator operator++(int) {myiterator tmp(*this); operator++(); return tmp;}
//  bool operator==(const myiterator& rhs) {return p==rhs.p;}
//  bool operator!=(const myiterator& rhs) {return p!=rhs.p;}
//  int& operator*() {return *p;}
//};


/**
 * Provides a few access functions given the data of a Python instance of HyperCubicShape.
 * Would be slow if C++ needed to construct an own Python object and call its functions.
 */
class HyperCubicShape {
public:
    HyperCubicShape (size_t dimension, tuple limits, dict lima, dict lima_inv) :
        D(dimension), _limits(limits), _lima(lima), _lima_inv(lima_inv){}
    //virtual ~HyperCubicShape ();

    //typedef IndexIterator< T, my_node< T > > iterator;
    //typedef IndexIterator< T const, my_node< T > const > const_iterator;
    typedef IndexIterator< T, Eigen::VectorXi > iterator;
    typedef IndexIterator< T const, Eigen::VectorXi const > const_iterator;
    IndexIterator begin(){ return IndexIterator(_limits,D,true); }
    IndexIterator end(){ return IndexIterator(_limits,D,false); }
    bool contains(tuple o);
    list get_neighbours(tuple k, bool forward=true); // forward=false is backward
    IndexIterator* get_index_iterator_chain(size_t direction=0)const;
    size_t getD()const {return D;}

private:
    size_t D;
    tuple _limits;
    dict _lima, _lima_inv;
};

namespace boost
{
    namespace python
    {
        template <>
        struct iterators< HyperCubicShape >
        {
            typedef HyperCubicShape::IndexIterator iterator;
            static iterator begin( Polygon& x) { return x.vertices_begin(); }
            static iterator end( Polygon& x) { return x.vertices_end(); }
        };
    }
}


//class IndexIterator {
//    const HyperCubicShape* _hcs;
//    size_t index, dim;
//    Eigen::VectorXi limits,z;
//    bool done;
//public:
//    IndexIterator(const HyperCubicShape* hcs, tuple _limits, size_t dimension);
//    void next();
//    bool operator()(){ return done; } //isDone?
//    Eigen::VectorXi& operator *(){ return z; }
//};
