#include <CGAL/Surface_mesh_simplification/Detail/Common.h>

//hash table to store curvatures
#include <unordered_map>


#include <iostream>
#include <fstream>

//analytics
#include "../../analizing/cgal_mesh_analytics.hpp"


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


  viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;

 
  Curvature_flat_placement(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a) 
  : analytics(a) {}
     
  template <typename Profile> 
  //optional<typename Profile::Point> //probably need a profile for later
  optional<typename Profile::Point> operator()( Profile const& aProfile ) const
  {


    return optional<typename Profile::Point>(midpoint(aProfile.p0(),aProfile.p1())) ;
  }
  


};

} // namespace Surface_mesh_simplification

} //namespace CGAL















