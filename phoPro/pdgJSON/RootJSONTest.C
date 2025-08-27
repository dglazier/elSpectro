std::string string_first_numbers(std::string const & str)
    {
      char const* digits = "0123456789.";
      std::size_t const n = str.find_first_of(digits);
      if (n != std::string::npos)
      {
        std::size_t const m = str.find_first_not_of(digits, n);
        return str.substr(n, m != std::string::npos ? m-n : m);
      }
      return std::string();
    }

bool string_contains(const std::string&  s1,const std::string& s2){

  if (s1.find(s2) != std::string::npos) return true;
  return false;
}

string string_between(const std::string& str,const char  s1,const char s2){
  int p1 = str.find(s1);
  int p2 = str.find(s2);
  auto result = str.substr(p1 + 1, p2 - p1 - 1);
  //remove whitespace
  result.erase (std::remove (result.begin(), result.end(), ' '), result.end()); 
  return result; 
}

double convert_unit(const std::string& str){
  if ( string_contains(str,"\"eV\"") ) return 1E-6;
  if ( string_contains(str,"\"MeV\"") ) return 1E-3;
  if ( string_contains(str,"\"GeV\"") ) return 1.;
  if ( string_contains(str,"\"TeV\"") ) return 1E3;
  return 1.;
}

class PDGJsonParticle {

public :
  
  string CheckForString( const std::string& type, const std::string& sline){
    if ( string_contains(sline,type) ){
     return string_between(sline,':',',');
    }
    return "";
  }
  
  void CheckForName(const std::string& sline){
    if(_takeNextValueForName==false) return;
    
    _name = CheckForString("\"description\"",sline);
    if( _name!=string() )_takeNextValueForName=false;
  }
  
  void CheckForMass(const std::string& sline){
    if(_takeNextValueForMass==-2) return;
    
    if( string_contains(sline,"MASS") &&  _takeNextValueForMass==false){
      _takeNextValueForMass=true;
    }
    if( string_contains(sline,"\"value\"") && _takeNextValueForMass==true ){
      _mass = std::strtod( string_first_numbers(sline).data(), nullptr );
      _takeNextValueForMass=-1;//take first value only
    }
    if( string_contains(sline,"\"unit\"") && _takeNextValueForMass==-1 ){
      _mass*= convert_unit(sline); //convert mass to GeV
      _takeNextValueForMass=-2;//take first value only
    }

  }
  void CheckForWidth(const std::string& sline){
    if(_takeNextValueForWidth==-2) return;
    
    if( string_contains(sline,"WIDTH") &&  _takeNextValueForWidth==false){
      _takeNextValueForWidth=true;
    }
    if( string_contains(sline,"\"value\"") && _takeNextValueForWidth==true ){
      _width = std::strtod( string_first_numbers(sline).data(), nullptr );
      //cout<<std::setprecision(15)<<Mass<<endl;
      _takeNextValueForWidth=-1;//take first value only
    }
    if( string_contains(sline,"\"unit\"") && _takeNextValueForWidth==-1 ){
      _width*= convert_unit(sline); //convert mass to GeV
      _takeNextValueForWidth=-2;//take first value only
    }


  }

  void Print(){
    std::cout<<"Name = "<<_name<<std::endl;
    std::cout<<"\t Mass = "<<_mass<<std::endl;
    std::cout<<"\t Width = "<<_width<<std::endl;
    
  }
  
private:
  string _name;
  bool _takeNextValueForName=true;
    
  double _mass=0;
  short _takeNextValueForMass=false;
    
  double _width=0;
  short _takeNextValueForWidth=false;

};
  
void RootJSONTest(){
  TMacro m("M001.json");
  //m.Print();
  auto lines = m.GetListOfLines();
  TString json;

  PDGJsonParticle particle;
  
  for(auto* line:*lines){
    std::string sline = line->GetName();
    particle.CheckForName(sline);
    particle.CheckForMass(sline);
    particle.CheckForWidth(sline);
    
  }
  particle.Print();
  // cout<<json<<endl;
}
