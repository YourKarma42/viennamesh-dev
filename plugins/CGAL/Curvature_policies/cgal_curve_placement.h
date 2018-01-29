#include <CGAL/Surface_mesh_simplification/Detail/Common.h>


namespace CGAL {

namespace Surface_mesh_simplification  
{

  template<class ECM_>
class Curvature_placement
{
public:
    
  typedef ECM_ ECM ;
  
public:

  Curvature_placement() {}
     
  template <typename Profile> 
  //optional<typename Profile::Point> //probably need a profile for later
  optional<typename Profile::Point> operator()( Profile const& aProfile ) const
  {
    //probably delete old vertices from the hastable
    return optional<typename Profile::Point>(midpoint(aProfile.p0(),aProfile.p1())) ;
  }
  
/*private:

  LindstromTurk_params mParams ;    */

};

} // namespace Surface_mesh_simplification

} //namespace CGAL

//TODO:  #defines machen wie in cgal