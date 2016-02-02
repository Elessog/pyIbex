#include "ibex_IntervalVector.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include "pyIbex_to_python_converter.h"

using namespace boost;
using namespace boost::python;
using namespace ibex;
namespace py = boost::python;


enum IBOOL { IN = 0 , OUT = 1 ,MAYBE =2 , UNK = 3 , EMPTY = 4 , UNK2 = 5 };



int testDist(const IntervalVector X, const IntervalVector m, const double rang2, const bool efficient ) {
    Interval reach = Interval(0,rang2);
    Interval Xm = max(Interval(0), sign( (X[0]-m[0].ub())*(X[0]-m[0].lb()))) * 
              min(sqr(X[0]-m[0].lb()),sqr(X[0]-m[0].ub()) ) 
            + max(Interval(0), sign( (X[1]-m[1].ub())*(X[1]-m[1].lb()))) 
              * min(sqr(X[1]-m[1].lb()),sqr(X[1]-m[1].ub()) ); // f-(||[X]-m||²)

    Interval Xp = max(sqr(X[0]-m[0].lb()),sqr(X[0]-m[0].ub())) +
                  max(sqr(X[1]-m[1].lb()),sqr(X[1]-m[1].ub()));//f+(||[X]-m||²)
      
    Interval Xub = Xm | Xp;
    if (reach.is_disjoint(Xub))
      return IBOOL::OUT;
    else if (Xub.is_subset(reach))
      return IBOOL::IN;
      
    bool b1 = ( Xm - reach.ub()).is_subset(Interval(-1000, 0));
    Interval B2 = ( reach.ub() - Xp);
    if  (B2.ub() < 0)
    {
        if (b1 )
          return IBOOL::MAYBE;
        if (efficient)
           return IBOOL::OUT;
	return IBOOL::UNK2;
    }
    return IBOOL::UNK;
}

bool find(std::vector<int> res,int val){
       for ( auto it = res.cbegin(); it != res.cend(); ++it ){
            if (*it==val)
              return true;
       }
       return false;
    }

int testR(const IntervalVector& X,const py::list& lst,double rang2,bool efficient){
    boost::python::ssize_t n = boost::python::len(lst);
    std::vector<IntervalVector> m;
    for (int i=0;i<n;i++){
        IntervalVector vect = IntervalVector(2);
        vect[0] = to_std_vector<Interval>(lst[i])[0];
        vect[1] = to_std_vector<Interval>(lst[i])[1];
        /*std::cout << "m : [" << vect[0].lb() << " , " << vect[0].ub() << " ]"  << " , [" << vect[1].lb() << " , " << vect[1].ub() << " ]"<< std::endl;*/
        m.push_back(vect);
    }
    //=to_std_vector<IntervalVector>(lst);
    std::vector<int> res;
    for ( auto it = m.cbegin(); it != m.cend(); ++it )
      res.push_back(testDist(X,*it,rang2,efficient));
    //std::cout << res[0] << std::endl;
    if (find(res,IBOOL::IN))
      return IBOOL::IN;
    if (find(res,IBOOL::UNK))
      return IBOOL::UNK;
    if (find(res,IBOOL::MAYBE))
      return IBOOL::MAYBE;
    bool test = true;
    for ( auto it = res.cbegin(); it != res.cend(); ++it )
      test &= (*it == IBOOL::OUT);
    if (test)
      return IBOOL::OUT;
    return IBOOL::UNK;
}



void export_testFunction(){
        def("testR",&testR);
}

