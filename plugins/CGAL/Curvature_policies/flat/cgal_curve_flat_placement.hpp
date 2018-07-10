#include <CGAL/Surface_mesh_simplification/Detail/Common.h>





namespace CGAL {

namespace Surface_mesh_simplification  
{

  template<class ECM_>
class Curvature_flat_placement
{
public:
    
  typedef ECM_ ECM ;

  typedef typename ECM::Point_3 Point_3;



public:

  Curvature_flat_placement()
  {}
     
  template <typename Profile> 
  //optional<typename Profile::Point> //probably need a profile for later
  optional<typename Profile::Point> operator()( Profile const& aProfile ) const
  {


    return optional<typename Profile::Point>(midpoint(aProfile.p0(),aProfile.p1())) ;
  }


};

} // namespace Surface_mesh_simplification

} //namespace CGAL

//TODO:  #defines machen wie in cgal