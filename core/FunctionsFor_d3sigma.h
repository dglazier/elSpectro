#pragma once
//Copied from jpacPhoto::inclusive_kinematics

namespace d3sigma{
  
  inline double Kallen(double x, double y, double z)
  {
    return x*x + y*y + z*z - 2. * (x*y + x*z + y*z);
  };
  
  // Also useful is M2 as a function of X and T
  inline double M2fromTX(double t, double x)
  {
    double lami = Kallen(_s, 0., _mT2);
    double lamf = Kallen(_s, _mX2, _minM2);
    double num = _mT2 * _mX2 + _mT2 * _s + _mX2 * _s - _s*_s - 2.*_s*t + sqrt(lami * lamf) * x;
    return num / (_mT2 - _s);
  };

}
