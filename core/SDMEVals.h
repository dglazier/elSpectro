#include <vector>

namespace elSpectro{


//so rho^alpha_x,y  -> sdme_alphaxy[alpha][x][y]

using  sdme_y = std::vector<double>;
using  sdme_xy = std::vector<sdme_y>;
using  sdme_alpha = std::vector<sdme_xy>;

class SDME {

 SDME(uint J=1,uint alphaMax=1) : _J{J}, _alphaMax{alphaMax} {
    _elements.resize(alphaMax, sdme_xy(2J+1, sdme_y(2J+1,0) ) );
  }

  double Val(uint alpha, int x, int y)const noexcept {return _elements[alpha][_J+x][_J+y]; }

  void SetAlpha(uint alpha,const sdme_xy& xy){
    _elements[alpha]=xy;
  }
  void SetElement(uint alpha,int x,int y,double val){
    _elements[alpha][_J+x][_J+y] = val;
  }
  
 private:
  
  uint _J={0};
  uint _alphaMax={0};
  
  sdme_alpha _elements;
  

};

}
