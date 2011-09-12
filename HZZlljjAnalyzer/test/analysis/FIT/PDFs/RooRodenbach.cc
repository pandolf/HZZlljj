#include <iostream>
#include <math.h>

#include "RooRodenbach.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

using namespace RooFit;

ClassImp(RooRodenbach) 
  
  RooRodenbach::RooRodenbach(){}

RooRodenbach::RooRodenbach(const char *name, const char *title, 
			   RooAbsReal& _x,
			   RooAbsReal& _mean,
			   RooAbsReal& _width,
			   RooAbsReal& _alpha,
			    RooAbsReal& _theta
			   ) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  mean("mean","mean",this,_mean),
  width("width","width",this,_width),
  alpha("alpha","alpha",this,_alpha),
  theta("theta","theta",this,_theta)
{ 
} 


RooRodenbach::RooRodenbach(const RooRodenbach& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  mean("mean",this,other.mean),
  width("width",this,other.width),
  alpha("alpha",this,other.alpha),
  theta("theta",this,other.theta)
  
{ 
} 

double RooRodenbach::evaluate() const 
{ 
  double a = cos(theta)*alpha - sin(theta)*width;
  double w = sin(theta)*alpha + cos(theta)*width;
  
  double beta=(a-mean)/(w*w);
  double A=exp(-(a-mean)*(a-mean)/(2*w*w))/exp(-beta*a);
  
  if(x<a)
    return exp(-(x-mean)*(x-mean)/(2*w*w));
  else
    return A*exp(-beta*x);
  
} 

