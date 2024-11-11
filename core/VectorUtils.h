#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <TMath.h>

namespace vector_utils{

  ////////////////print contents of vector
  template <typename T>
    void print(const std::vector<T>& vec){
    for(const auto& item:vec){
      std::cout<<" "<<item;
    }
    std::cout<<" "<<std::endl;
  }
  
  ///////////////check if vector contains item  
  template <typename T>
    bool contains(const std::vector<T>& vec,T item){
    return( std::find(vec.begin(),vec.end(),item)!=vec.end() );
  }

  ////////////// set nans to 0
  template <typename T>
    void zero_nans(std::vector<T>& vec){
    for(auto& val:vec){
      if(TMath::IsNaN(val))val=0;
    }
  }

    
  //for simulataneous sorting
  //https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
  
/*   vector<MyObject> vectorA; */
/* vector<int> vectorB; */

/* auto p = sort_permutation(vectorA, */
/*     [](T const& a, T const& b){ /\*some comparison*\/ }); */

/* apply_permutation_in_place(vectorA, p); */
/* apply_permutation_in_place(vectorB, p); */

  
  template <typename T, typename Compare>
    std::vector<std::size_t> sort_permutation(
					      const std::vector<T>& vec,
					      const Compare& compare)
  {
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
	      [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
    return p;
  }

  template <typename T>
    void apply_permutation_in_place(
				    std::vector<T>& vec,
				    const std::vector<std::size_t>& p)
    {
      std::vector<bool> done(vec.size());
      for (std::size_t i = 0; i < vec.size(); ++i)
	{
	  if (done[i])
	    {
	      continue;
	    }
	  done[i] = true;
	  std::size_t prev_j = i;
	  std::size_t j = p[i];
	  while (i != j)
	    {
	      std::swap(vec[prev_j], vec[j]);
	      done[j] = true;
	      prev_j = j;
	      j = p[j];
	    }
	}
    }

  //sort v1,v2,v3 in assending order of v1
  template <typename T1,typename T2,typename T3>
    void synched_sort(std::vector<T1>& v1,std::vector<T2>& v2,std::vector<T3>& v3){
    
    auto p = vector_utils::sort_permutation(v1,
					   [](double const& a, double const& b){ return a<b; });
    
    vector_utils::apply_permutation_in_place(v1, p);
    vector_utils::apply_permutation_in_place(v2, p);
    vector_utils::apply_permutation_in_place(v3, p);
    
  }
  //sort v1,v2 in assending order of v1
  template <typename T1,typename T2>
    void synched_sort(std::vector<T1>& v1,std::vector<T2>& v2){
    
    auto p = vector_utils::sort_permutation(v1,
					   [](double const& a, double const& b){ return a<b; });
    
    vector_utils::apply_permutation_in_place(v1, p);
    vector_utils::apply_permutation_in_place(v2, p);
    
  }
  
  template <typename T>
  double max_consecutive_diff(const std::vector<T>& vec){
    size_t n = vec.size();
    double ans=0.;
    for(size_t i=1;i<n;i++){
      auto diff = std::abs(vec[i] - vec[i-1]);
      ans = (diff>ans) ? diff : ans;
    }
    return ans;
  }

   template <typename T>
   void take_max_of_neighbours(std::vector<T>& vec){
     if(vec.size()<2) return;
     auto vec2 = vec; //take a copy
     size_t n = vec.size();
     //do first and last which only have 1 neighbour
     vec[0] = vec2[1] > vec2[0] ? vec2[1] : vec2[0];
    if(vec.size()>2)
      vec[n-1] = vec2[n-2] > vec2[n-1] ? vec2[n-2] : vec2[n-1];

     for(int i=1;i<n-1;++i){
       vec[i] = vec2[i-1] > vec2[i] ? vec2[i-1] : vec2[i];
       vec[i] = vec2[i+1] > vec2[i] ? vec2[i+1] : vec2[i];
     }
   }

}
