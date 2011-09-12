#include <iostream>
#include <math.h>
#include "TMath.h"
#include "RooCB.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

using namespace RooFit;

 ClassImp(RooCB) 

 RooCB::RooCB(){}

 RooCB::RooCB(const char *name, const char *title, 
		    RooAbsReal& _x,
		    RooAbsReal& _mean,
		    RooAbsReal& _width,
		    RooAbsReal& _alpha,
	      RooAbsReal& _n,
	      RooAbsReal& _theta	      
		    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   alpha("alpha","alpha",this,_alpha),
   n("n","n",this,_n),
   theta("theta","theta",this,_theta)
 { 
 } 


 RooCB::RooCB(const RooCB& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   alpha("alpha",this,other.alpha),
   n("n",this,other.n),
   theta("theta",this,other.theta)
 { 
 } 

 double RooCB::evaluate() const 
 { 
   double a = cos(theta)*alpha - sin(theta)*width;
   double w = sin(theta)*alpha + cos(theta)*width;

   double t = (x-mean)/w;
   if(a<0) t = -t;
   
   double absa = fabs((double)a);
   
   double A = TMath::Power(n/absa,n)*exp(-0.5*absa*absa);
   double B = n/absa-absa;
   
   if(t >= -absa){
     return exp(-0.5*t*t);
   }else{
     return A/TMath::Power(B-t,n);
   }   
 } 

