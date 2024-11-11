#pragma once


namespace cpp{
  /// Take an object of a derived class and return
  /// a shared_ptr of the base class
  /// usage : auto sh = MakeBaseShared<TObject>(derived_obj);
  
  template <class Base, class Derive >
    std::shared_ptr<Base> MakeBaseShared(const Derive& d){
    std::shared_ptr<Base>  v{new Derive(d)};
    return v;
  }

}
