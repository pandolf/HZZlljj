#include <iostream>
#include <math.h>

#include "RooDoubleGauss.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

 ClassImp(RooDoubleGauss) 

 RooDoubleGauss::RooDoubleGauss(const char *name, const char *title, 
				RooAbsReal& _x,
				RooAbsReal& _mean,
				RooAbsReal& _shift,
				RooAbsReal& _width,
				RooAbsReal& _alpha,
				RooAbsReal& _f
		    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   shift("shift","shift",this,_shift),
   width("width","width",this,_width),
   alpha("alpha","alpha",this,_alpha),
   f("f","f",this,_f)
 { 
 } 


 RooDoubleGauss::RooDoubleGauss(const RooDoubleGauss& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   shift("shift",this,other.shift),
   width("width",this,other.width),
   alpha("alpha",this,other.alpha),
   f("f",this,other.f)

 { 
 } 



 double RooDoubleGauss::evaluate() const 
 { 
   double pi=3.1415;
   double gauss1 = (1/sqrt(2*pi*width*width))*exp(-(x-mean)*(x-mean)/(2*width*width));
   double gauss2 = (1/sqrt(2*pi*width*width*alpha*alpha))*exp(-(x-mean-shift)*(x-mean-shift)/(2*width*width*alpha*alpha));
   return gauss1+f*gauss2;
 } 

