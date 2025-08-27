#include <nlohmann/json.hpp> // nlohmann::json::parse
#include <fstream>           // std::ifstream
#include <iostream>           // std::ifstream

bool string_contains(const std::string&  s1,const std::string& s2){

  if (s1.find(s2) != std::string::npos) return true;
  return false;
}
double convert_unit(const std::string& str){
  if ( str=="eV" ) return 1E-6;
  if ( str=="MeV" ) return 1E-3;
  if ( str=="GeV" ) return 1.;
  if ( str=="TeV" ) return 1E3;
  return 1.;
}
double string_first_numbers(std::string const & str)
    {
      char const* digits = "0123456789.";
      std::size_t const n = str.find_first_of(digits);
      if (n != std::string::npos)
      {
        std::size_t const m = str.find_first_not_of(digits, n);
        return std::strtod( str.substr(n, m != std::string::npos ? m-n : m).data(), nullptr);
      }
      return 0.;
    }

void RootJSONTest2(){
  const nlohmann::ordered_json fullData = nlohmann::ordered_json::parse(std::ifstream("M001.json"));
  auto name = fullData.find("description").value();
  std::cout<<name<<" "<<std::endl;
  auto summaries = fullData.find("summaries");
  auto properties = summaries->find("properties");
  double mass = 0;
  double width = 0;
  
  for(auto& prop:*properties){
    auto getVal = [&prop](double& val){
      auto pdg = prop.find("pdg_values");
       for(auto& p:*pdg){
	 //mass value
	 val= p.find("value").value().get<double>();
	 auto unit = convert_unit(p.find("unit").value().get<string>() );
	 val*=unit;
	 break;
       }
    };
    
    if( string_contains(prop.find("description").value(),"MASS") ){
      getVal(mass);
    }
    if( string_contains(prop.find("description").value(),"WIDTH") ){
      getVal(width);
    }
    
  }
  auto branching = summaries->find("branching_ratio");


  cout<<name <<" mass = "<<mass<<" width "<<width<<endl;
}
