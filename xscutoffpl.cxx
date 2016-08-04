#include <XSFunctions/functionMap.h>
#include <XSUtil/Numerics/IncGamma.h>

void cutoffPowerLaw (const RealArray& energyArray, const RealArray& params, 
                int spectrumNumber, RealArray& fluxArray, 
                RealArray& fluxErrArray, const string& initString);

extern "C" double incgamma(double a, double x)
{
  Numerics::IncGamma incGamma;
  Real y(x);
  Real b(a);
  return double(incGamma(b,y)); 
}

extern "C" void cutoffpl(double *ear, const int ne, double *param, double *photar)

{ int i;
  RealArray energyArray(ear, ne+1);
  RealArray params(param, 2);
  RealArray fluxArray(ne);
  RealArray fluxErrArray;
  
  cutoffPowerLaw(energyArray, params, 0, fluxArray, fluxErrArray, "");
  for(i=1;i<=ne;i++)photar[i-1]=fluxArray[i-1];
  return;
}