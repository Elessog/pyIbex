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
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include "opencv2/videoio/videoio.hpp"
#include <pthread.h>

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
    s.push(X);
    while (!s.empty()) 
    {
		IntervalVector box=s.top();
		s.pop();
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

py::list aSIVIA(const IntervalVector& X,const py::list& lst,double rang2,double eps=0.1,bool efficient=true){
   std::vector<IntervalVector> in;
   std::vector<IntervalVector> out;
   std::vector<IntervalVector> maybe;
   std::vector<IntervalVector> border;
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
   std::stack<IntervalVector> s;
   s.push(X);
   while (!s.empty()) 
   {
       IntervalVector box=s.top();
       s.pop();
       int res=testRC(box,m,rang2,efficient);

       if (res==IBOOL::IN)
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[r]");
          #endif
          in.push_back(box);
       }
       else if (res==IBOOL::OUT)
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[o]");
          #endif
          out.push_back(box);
       }
       else if (res==IBOOL::MAYBE)
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[b]");
          #endif
          maybe.push_back(box);
       }
       else if (box.max_diam()>eps)
       {
       // get the index of the dimension of maximal size (false stands for "max")
           i = box.extr_diam_index(false);
           std::pair<IntervalVector,IntervalVector> p=box.bisect(i);
           s.push(p.first);
           s.push(p.second);
       }
       else
       {
           #if defined VIBES_COMP_
           vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[y]");
           #endif
           border.push_back(box);
       }
   }

   #if defined VIBES_COMP_
   vibes::endDrawing();
   #endif
   py::list list;
   list.append(in);
   list.append(out);
   list.append(maybe);
   list.append(border);
   py::incref(list.ptr());
   return list;
}


   //Method for the transformation of the SIVIA scale to the image (pixel) scale:
void toPixels(const float echellePixel,int resp[2],const int x,const int y,const int i0,const int j0){
   //get the value in the matrix imgIntegral at the position x, y in the interval world
   resp[0] = round(x/echellePixel) + i0;
   resp[1] = -round(y/echellePixel) + j0;
}
    //Method for calculate the sum and the size of the box in the integral image:
void computeBox(const float &echellePixel,int resp[2],const IntervalVector& X,const std::vector<std::vector<int>> &l_imgIntegral,const int &i0,const int &j0){
     //Calcul of the points A,B,C and D:
     int A[2];
     toPixels(echellePixel,A,X[0].lb(),X[1].ub(),i0,j0);
     int C[2];
     toPixels(echellePixel,C,X[0].ub(),X[1].lb(),i0,j0);
     //Calcul of the sum:
     resp[0] = l_imgIntegral[A[1]][A[0]] + l_imgIntegral[C[1]][C[0]] - l_imgIntegral[A[1]][C[0]] - l_imgIntegral[C[1]][A[0]]; //Sum=A + C - B - D  
     //Calcul of the size of the box:
     resp[1] = (C[0] - A[0]) * (C[1] - A[1]);
}

//Method which create the pavage with SIVIA:
int testIm(const float echellePixel,const IntervalVector& X,const std::vector<std::vector<int>>& image,const int i0,const int j0){
        int resp[2];
        computeBox(echellePixel,resp,X,image,i0,j0);

        if (resp[0] == 0)
            return IBOOL::OUT;
        if (resp[0] == resp[1])
            return IBOOL::IN;
        else
            return IBOOL::UNK;
}

py::list imSIVIA(const IntervalVector& X,const py::list& lst,double rang2,const py::list& image,const float echellePixel,
                 const int i0,const int j0,double eps=0.1,bool efficient=true){
   std::vector<IntervalVector> in;
   std::vector<IntervalVector> out;
   std::vector<IntervalVector> maybe;
   std::vector<IntervalVector> border;
   boost::python::ssize_t n = boost::python::len(lst);
   boost::python::ssize_t width = boost::python::len(image[0]);
   boost::python::ssize_t height = boost::python::len(image);
   std::vector<IntervalVector> m;
   std::vector<std::vector<int>> imageC;
   int i;
   int j;
   for (i=0;i<n;i++){
       IntervalVector vect = IntervalVector(2);
       vect[0] = to_std_vector<Interval>(lst[i])[0];
       vect[1] = to_std_vector<Interval>(lst[i])[1];
       m.push_back(vect);
   }
   for (j=0;j<height;j++){
       std::vector<int> line;
       for (i=0;i<width;i++){
        line.push_back(py::extract<int>(image[j][i]));
       }
       imageC.push_back(line);       
   }
   #if defined VIBES_COMP_
   vibes::beginDrawing ();
   vibes::newFigure("lab2");
   #endif
   std::stack<IntervalVector> s;
   s.push(X);
   while (!s.empty()) 
   {
       IntervalVector box=s.top();
       s.pop();
       int res=testIm(echellePixel,box,imageC,i0,j0);
       int res2=testRC(box,m,rang2,efficient);
       if (res==IBOOL::IN || res2 ==IBOOL::IN)
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[r]");
          #endif
          in.push_back(box);
       }
       else if ((res == IBOOL::OUT && res2 == IBOOL::OUT ) || (res2 == IBOOL::OUT &&  res != IBOOL::UNK))
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[o]");
          #endif
          out.push_back(box);
       }
       else if (res2==IBOOL::MAYBE)
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[b]");
          #endif
          maybe.push_back(box);
       }
       else if (box.max_diam()>eps)
       {
       // get the index of the dimension of maximal size (false stands for "max")
           i = box.extr_diam_index(false);
           std::pair<IntervalVector,IntervalVector> p=box.bisect(i);
           s.push(p.first);
           s.push(p.second);
       }
       else
       {
           #if defined VIBES_COMP_
           vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[y]");
           #endif
           border.push_back(box);
       }
   }

   #if defined VIBES_COMP_
   vibes::endDrawing();
   #endif
   py::list list;
   list.append(in);
   list.append(out);
   list.append(maybe);
   list.append(border);
   py::incref(list.ptr());
   return list;
}


class clSIVIA
{
  public:
      clSIVIA(const py::list& image,int kernelSize = 3 ,int erosion_elem = 0,int iter = 1) {
       boost::python::ssize_t width = boost::python::len(image[0]);
       boost::python::ssize_t height = boost::python::len(image);
       int i;
       int j;
       int erosion_type;
       iteration = iter;
       if( erosion_elem == 0 ){ erosion_type = cv::MORPH_RECT; }
       else if( erosion_elem == 1 ){ erosion_type = cv::MORPH_CROSS; }
       else if( erosion_elem == 2) { erosion_type = cv::MORPH_ELLIPSE; }

       kernelTrail = cv::getStructuringElement( erosion_type,
                                       cv::Size( 2*kernelSize + 1, 2*kernelSize+1 ),
                                       cv::Point(kernelSize, kernelSize) );
       is_recording = false;
       is_rect = false;
       do_erode = true;
       show = false;
       width_ = width;
       height_ = height;
       imageCBin = cv::Mat(height_,width_,CV_8UC1);
       imageC = cv::Mat(height_+1,width_+1,CV_32SC1);
       maskColor = cv::Mat(height_,width_,CV_8UC3);
       maskColor.setTo(cv::Scalar(51,153,255));
       for (j=0;j<height_;j++){
          for (i=0;i<width_;i++)
		imageCBin.at<uint8_t>(j,i)=py::extract<int>(image[j][i]);
       }
       cv::integral(imageCBin,imageC,CV_32SC1);
       imageC = imageC(cv::Rect(1,1,width_,height_));
       std::cout << "First Image saved.." << std::endl;
       past = cv::Mat(height_,width_,CV_8UC1);
       imageColor = cv::Mat(height_,width_,CV_8UC3);
       imageColor1 = cv::Mat(height_,width_,CV_8UC3);
       imageColor2 = cv::Mat(height_,width_,CV_8UC3);
       pastInt = cv::Mat(height_,width_,CV_32SC1);
       past.setTo(cv::Scalar(0,0,0));
       pastInt.setTo(cv::Scalar(0,0,0));
       imageColor1.setTo(cv::Scalar(0,0,0));
       imageColor2.setTo(cv::Scalar(0,0,0));
       std::cout << "clSIVIA initialised !" << std::endl;
      }
      py::list imSIVIA(const IntervalVector& X,const py::list& lst,double rang2,const float echellePixel,
                 const int i0,const int j0,double eps,bool efficient);
      void setRecord(std::string const & msg,float fps,bool rectangle);
      void setErode(bool do_erod);
     void stopScreening();
      ~clSIVIA(){
        show = false;
        if (is_recording)
           pthread_join(thread,NULL);
      }
  private:
     pthread_mutex_t mutex;
     pthread_t thread;
     bool change;
     cv::Mat imageCBin;
     cv::Mat imageC;
     cv::Mat past;
     cv::Mat pastInt;
     cv::Mat imageColor;
     cv::Mat imageColor1;
     cv::Mat imageColor2;
     cv::VideoWriter record;
     cv::Mat maskColor;
     cv::Mat kernelTrail;
     bool is_recording;
     bool is_rect;
     bool do_erode;
     int iteration;
     int width_;
     int height_;
     bool show;
     int testIm(const float echellePixel,const IntervalVector& X,const int i0,const int j0);
     void computeBox(const float &echellePixel,int resp[2],const IntervalVector& X,const int &i0,const int &j0);
     int testIm3(const float echellePixel,const IntervalVector& X,const int i0,const int j0);
     void computeBox2(const float &echellePixel,int resp[2],const IntervalVector& X,const int &i0,const int &j0);
     void drawBox(const IntervalVector& X,const float echellePixel,const int i0,const int j0);
     void drawReal(const IntervalVector& X,const float echellePixel,const int i0,const int j0,cv::Scalar color);
     void screening();
     static void *screening_helper(void *context)
     {
        ((clSIVIA *)context)->screening();
        return NULL;
     }
};

void clSIVIA::setRecord(std::string const & msg,float fps,bool rectangle = false){
   record = cv::VideoWriter(msg,CV_FOURCC('D','I','V','X'), fps,cv::Size(width_,height_), true);
   is_rect = rectangle;
   is_recording = true;
   mutex = PTHREAD_MUTEX_INITIALIZER;
   if (pthread_create(&thread, NULL,&clSIVIA::screening_helper,this) != 0)
             exit(EXIT_FAILURE);
}


void clSIVIA::screening(){
     show = true;
     while(show){
        pthread_mutex_lock(&mutex);
        if (change)
          imageColor1.copyTo(imageColor2);
        change = false;
        pthread_mutex_unlock(&mutex);
        cv::imshow("Boxes",imageColor2);
        cv::waitKey(20);
     }
}

void clSIVIA::stopScreening(){
     show = false;
}

void clSIVIA::setErode(bool do_erod){
    do_erode = do_erod;
}

void clSIVIA::computeBox(const float &echellePixel,int resp[2],const IntervalVector& X,const int &i0,const int &j0){
     //Calcul of the points A,B,C and D:
     int A[2];
     toPixels(echellePixel,A,X[0].lb(),X[1].ub(),i0,j0);
     int C[2];
     toPixels(echellePixel,C,X[0].ub(),X[1].lb(),i0,j0);
     //Calcul of the sum:
     resp[0] = imageC.at<int>(A[1],A[0]) + imageC.at<int>(C[1],C[0]) - imageC.at<int>(A[1],C[0]) - imageC.at<int>(C[1],A[0]); //Sum=A + C - B - D  
     //Calcul of the size of the box:
     resp[1] = (C[0] - A[0]) * (C[1] - A[1]);
}

//Method which create the pavage with SIVIA:
int clSIVIA::testIm(const float echellePixel,const IntervalVector& X,const int i0,const int j0){
        int resp[2];
        clSIVIA::computeBox(echellePixel,resp,X,i0,j0);

        if (resp[0] == 0)
            return IBOOL::OUT;
        if (resp[0] == resp[1])
            return IBOOL::IN;
        else
            return IBOOL::UNK;
}

void clSIVIA::computeBox2(const float &echellePixel,int resp[2],const IntervalVector& X,const int &i0,const int &j0){
     //Calcul of the points A,B,C and D:
     int A[2];
     toPixels(echellePixel,A,X[0].lb(),X[1].ub(),i0,j0);
     int C[2];
     toPixels(echellePixel,C,X[0].ub(),X[1].lb(),i0,j0);
     //Calcul of the sum:
     int sign[4] = {1,1,-1,-1};
     resp[0] = sign[0]*pastInt.at<int>(A[1],A[0]) + sign[1]*pastInt.at<int>(C[1],C[0]) + sign[2]*pastInt.at<int>(A[1],C[0]) + sign[3]*pastInt.at<int>(C[1],A[0]); //Sum=A + C - B - D  

     //Calcul of the size of the box:
     resp[1] = (C[0] - A[0]) * (C[1] - A[1]);
     if (resp[0]<0)
        std::cout<<"error A:"<< A[0] <<" "<<A[1] <<"  C :" << C[0] <<" " <<C[1] <<" "<<resp[1] <<std::endl;
}

//Method which create the pavage with SIVIA:
int clSIVIA::testIm3(const float echellePixel,const IntervalVector& X,const int i0,const int j0){
        int resp[2];
        clSIVIA::computeBox2(echellePixel,resp,X,i0,j0);

        if (resp[0] == 0)
            return IBOOL::OUT;
        if (resp[0] == resp[1]){
               return IBOOL::IN;
        }
        else
            return IBOOL::UNK;
}

void clSIVIA::drawBox(const IntervalVector& X,const float echellePixel,const int i0,const int j0){
     int A[2];
     toPixels(echellePixel,A,X[0].lb(),X[1].ub(),i0,j0);
     int C[2];
     toPixels(echellePixel,C,X[0].ub(),X[1].lb(),i0,j0);
     std::vector<cv::Point> fillContSingle;
     fillContSingle.push_back(cv::Point(A[0],A[1]));
     fillContSingle.push_back(cv::Point(C[0],A[1]));
     fillContSingle.push_back(cv::Point(C[0],C[1]));
     fillContSingle.push_back(cv::Point(A[0],C[1]));
     std::vector<std::vector<cv::Point> > fillContAll;
     fillContAll.push_back(fillContSingle);
     cv::fillPoly( past, fillContAll, cv::Scalar(1,1,1));
}

void clSIVIA::drawReal(const IntervalVector& X,const float echellePixel,const int i0,const int j0,cv::Scalar color){
     int A[2];
     toPixels(echellePixel,A,X[0].lb(),X[1].ub(),i0,j0);
     int C[2];
     toPixels(echellePixel,C,X[0].ub(),X[1].lb(),i0,j0);
     std::vector<cv::Point> fillContSingle;
     fillContSingle.push_back(cv::Point(A[0],A[1]));
     fillContSingle.push_back(cv::Point(C[0],A[1]));
     fillContSingle.push_back(cv::Point(C[0],C[1]));
     fillContSingle.push_back(cv::Point(A[0],C[1]));
     std::vector<std::vector<cv::Point> > fillContAll;
     fillContAll.push_back(fillContSingle);
     cv::fillPoly( imageColor, fillContAll, color);
     if (is_rect)
       cv::rectangle(imageColor,cv::Point(A[0],A[1]),cv::Point(C[0],C[1]), cv::Scalar(0,0,0));
}

py::list clSIVIA::imSIVIA(const IntervalVector& X,const py::list& lst,double rang2,const float echellePixel,
                 const int i0,const int j0,double eps=0.1,bool efficient=true){
   std::vector<IntervalVector> in;
   std::vector<IntervalVector> out;
   std::vector<IntervalVector> maybe;
   std::vector<IntervalVector> border;
   boost::python::ssize_t n = boost::python::len(lst);
   std::vector<IntervalVector> m;
   int i;
   for (i=0;i<n;i++){
       IntervalVector vect = IntervalVector(2);
       vect[0] = to_std_vector<Interval>(lst[i])[0];
       vect[1] = to_std_vector<Interval>(lst[i])[1];
       m.push_back(vect);
   }
   #if defined VIBES_COMP_
   vibes::beginDrawing ();
   vibes::newFigure("lab2");
   #endif
   std::stack<IntervalVector> s;
   s.push(X);
   std::cout << "Initialisation Finished" << std::endl; 
   while (!s.empty()) 
   {
       IntervalVector box=s.top();
       s.pop();

       int res=clSIVIA::testIm(echellePixel,box,i0,j0);
       int res3=clSIVIA::testIm3(echellePixel,box,i0,j0);
       int res2=testRC(box,m,rang2,efficient);

       if (res==IBOOL::IN || res2 ==IBOOL::IN || res3 == IBOOL::IN )
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[r]");
          #endif
          in.push_back(box);
          clSIVIA::drawBox(box,echellePixel,i0,j0);
          clSIVIA::drawReal(box,echellePixel,i0,j0,cv::Scalar(0,204,0));
       }
       else if (res == IBOOL::OUT && res2 == IBOOL::OUT && res3 == IBOOL::OUT)
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[o]");
          #endif
          out.push_back(box);
          if (is_recording)
             clSIVIA::drawReal(box,echellePixel,i0,j0,cv::Scalar(255,128,0));
       }
       else if (res2==IBOOL::MAYBE)
       {
          #if defined VIBES_COMP_
          vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[b]");
          #endif
          maybe.push_back(box);
       }
       else if (box.max_diam()>eps)
       {
       // get the index of the dimension of maximal size (false stands for "max")
           i = box.extr_diam_index(false);
           std::pair<IntervalVector,IntervalVector> p=box.bisect(i);
           s.push(p.first);
           s.push(p.second);
       }
       else
       {
           #if defined VIBES_COMP_
           vibes::drawBox(box[0].lb(), box[0].ub(), box[1].lb(), box[1].ub(), "k[y]");
           #endif
           border.push_back(box);
       }
   }
   #if defined VIBES_COMP_
   vibes::endDrawing();
   #endif
   pastInt.setTo(cv::Scalar(0,0,0));
   if (do_erode)
     cv::erode(past,past,kernelTrail,cv::Point(-1,-1),iteration);
   cv::integral(past,pastInt,CV_32SC1);
   pastInt = pastInt(cv::Rect(1,1,width_,height_));
   //reset imgBin
   past.setTo(cv::Scalar(0,0,0));
   
   if (is_recording){
     maskColor.copyTo(imageColor,imageCBin);
     for(auto it = m.cbegin(); it!=m.cend();it++){
        int C[2];
        toPixels(echellePixel,C,(*it)[0].mid(),(*it)[1].mid(),i0,j0);
        circle(imageColor,cv::Point(C[0],C[1]), 2, cv::Scalar(0,0,255), -1);
     }
     record << imageColor;
     pthread_mutex_lock(&mutex);
     imageColor.copyTo(imageColor1);
     change = true;
     pthread_mutex_unlock(&mutex);
     imageColor.setTo(cv::Scalar(0,0,0));
   }
   
   py::list list;
   list.append(in);
   list.append(out);
   list.append(maybe);
   list.append(border);
   py::incref(list.ptr());
   return list;
}

void export_testFunction(){
        python::to_python_converter<std::vector<IntervalVector>, std_vector_to_list<IntervalVector> >();
        def("testR",&testR);
        def("SIVIAtest",&rSIVIA);
        def("fSIVIAtest",&aSIVIA);
        def("imSIVIAtest",&imSIVIA);
        class_<clSIVIA>("clSIVIA", init<const py::list&,int,int,int>())
           .def("imSIVIA", &clSIVIA::imSIVIA)  // Add a regular member function.
           .def("setRecord", &clSIVIA::setRecord)
           .def("setErode", &clSIVIA::setErode)
           .def("stopScreening",&clSIVIA::stopScreening)
        ;
}

