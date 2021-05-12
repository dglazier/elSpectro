#pragma once

#include <vector>
#include <complex>

namespace elSpectro{


  //so rho^alpha_x,y  -> sdme_alphaxy[alpha][x][y]

  using  sdme_y = std::vector<std::complex<double>>;
  using  sdme_xy = std::vector<sdme_y>;
  using  sdme_alpha = std::vector<sdme_xy>;
 
  class SDME {

  public:

    
  SDME(uint J=1,uint alphaMax=1) : _J{J}, _alphaMax{alphaMax}{
      
      _elements.resize((int)(alphaMax), sdme_xy((int)(J+1), sdme_y((int)(2*J+1)) ) );
      
    }

    std::complex<double> Val(uint alpha, int x, int y)const noexcept {
      return _elements[alpha][x][_J+y];
    }
    double Re(uint alpha, int x, int y)const noexcept {
      return std::real(_elements[alpha][x][_J+y]);
    }
    double Im(uint alpha, int x, int y)const noexcept {
      return std::imag(_elements[alpha][x][_J+y]);
    }
    double Abs(uint alpha, int x, int y)const noexcept {
      return std::abs(_elements[alpha][x][_J+y]);
    }

    void SetAlpha(uint alpha,const sdme_xy& xy){
      _elements[alpha]=xy;
    }
    void SetElement(uint alpha,int x,int y,std::complex<double>  val){
      _elements[alpha][x][_J+y] = val;
    }
    
    uint Spin()const {return _J;}
  private:
  
    uint _J={0};
    uint _alphaMax={0};
  
    sdme_alpha _elements;
  

  };

}
