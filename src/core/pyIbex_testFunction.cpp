#include "ibex_IntervalVector.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include "pyIbex_to_python_converter.h"
#include<stack>
#if defined VIBES_COMP_ 
#include "vibes/vibes.h"
#endif

using namespace boost;
using namespace boost::python;
using namespace ibex;
namespace py = boost::python;


enum IBOOL { IN = 0 , OUT = 1 ,MAYBE =2 , UNK = 3 , EMPTY = 4 , UNK2 = 5 };

void printIntVect(const IntervalVector vect){
   std::cout << "m : [" << vect[0].lb() << " , " << vect[0].ub() << " ]"  << " , [" 
             << vect[1].lb() << " , " << vect[1].ub() << " ]"<< std::endl;
}

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
    int i;
    for (i=0;i<n;i++){
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


int testRC(const IntervalVector& X,const std::vector<IntervalVector> m,double rang2,bool efficient){
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


std::vector<IntervalVector> rSIVIA(const IntervalVector& X,const py::list& lst,double rang2,double eps=0.1,bool efficient=true){
    boost::python::ssize_t n = boost::python::len(lst);
    std::vector<IntervalVector> m;
    #if defined VIBES_COMP_
    vibes::beginDrawing ();
    vibes::newFigure("lab2");
    #endif
    int i;
    for (i=0;i<n;i++){
        IntervalVector vect = IntervalVector(2);
        vect[0] = to_std_vector<Interval>(lst[i])[0];
        vect[1] = to_std_vector<Interval>(lst[i])[1];
        m.push_back(vect);
    }
    std::vector<IntervalVector> AllIn;
    std::stack<IntervalVector> s;
    //printIntVect(X);
    s.push(X);
    while (!s.empty()) 
    {
		IntervalVector box=s.top();
		s.pop();
                //printIntVect(box);
		int res=testRC(box,m,rang2,efficient);

		if (res==IBOOL::IN)
                {
                        #if defined VIBES_COMP_
			vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[r]");
                        #endif
                        AllIn.push_back(box);
                }
		else if (res==IBOOL::OUT)
                {
                        #if defined VIBES_COMP_
			vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[o]");
                        #endif
                        continue;
                }
                else if (res==IBOOL::MAYBE)
                {
                        #if defined VIBES_COMP_
			vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[b]");
                        #endif
                        continue;
                }
		else if (box.max_diam()>eps) {
			// get the index of the dimension of maximal size (false stands for "max")
			i = box.extr_diam_index(false);
			std::pair<IntervalVector,IntervalVector> p=box.bisect(i);
			s.push(p.first);
			s.push(p.second);
		}
                #if defined VIBES_COMP_
                else
                        vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[y]");
                #endif
    }

    #if defined VIBES_COMP_
    vibes::endDrawing();
    #endif
    return AllIn;
}

/*std::vector<IntervalVector> rSIVIA(const IntervalVector& X0,std::vector<IntervalVector> m,double rang2,bool efficient){
    std::vector<IntervalVector> X;
    X	

}*/

void export_testFunction(){
        python::to_python_converter<std::vector<IntervalVector>, std_vector_to_list<IntervalVector> >();
        def("testR",&testR);
        def("SIVIAtest",&rSIVIA);
}

