#include <iostream>
#include <math.h>

#include "RooFermi.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

using namespace RooFit;


 RooFermi::RooFermi(){}

 RooFermi::RooFermi(const char *name, const char *title, 
		      RooAbsReal& _x,
		      RooAbsReal& _cutOff,
		    RooAbsReal& _beta
		    ) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   cutOff("cutOff","cutOff",this,_cutOff),
   beta("beta","beta",this,_beta)
 { 
 } 


 RooFermi::RooFermi(const RooFermi& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   cutOff("cutOff",this,other.cutOff),
   beta("beta",this,other.beta)

 { 
 } 



 double RooFermi::evaluate() const 
 { 
   return 1.0/(exp((cutOff-x)/beta)+1);
 } 

