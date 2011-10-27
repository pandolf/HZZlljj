#include <iostream>
#include <math.h>

#include "RooRelBW.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

using namespace RooFit;

 RooRelBW::RooRelBW(){}

 RooRelBW::RooRelBW(const char *name, const char *title, 
		    RooAbsReal& _x,
		    RooAbsReal& _mean,
		    RooAbsReal& _width,
		    RooAbsReal& _n
		    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   n("n","n",this,_n)
 { 
 } 


 RooRelBW::RooRelBW(const RooRelBW& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   n("n",this,other.n)

 { 
 } 



 double RooRelBW::evaluate() const 
 { 
   return pow(x*x,n)/((x*x-mean*mean)*(x*x-mean*mean)+pow(x*x/(mean*mean),2*n)*mean*mean*width*width);
 } 

