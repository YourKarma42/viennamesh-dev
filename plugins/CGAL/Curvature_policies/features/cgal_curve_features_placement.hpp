#include <CGAL/Surface_mesh_simplification/Detail/Common.h>

//hash table to store curvatures
#include <unordered_map>


#include <iostream>
#include <fstream>

//analytics
#include "../../analizing/cgal_mesh_analytics.hpp"


//#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core.h>


namespace CGAL {

namespace Surface_mesh_simplification  
{

  template<class ECM_>
class Curvature_features_placement
{
public:
    
  typedef ECM_ ECM ;

  typedef typename ECM::Point_3 Point_3;





public:


  viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;

 
  Curvature_features_placement(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a, LindstromTurk_params const& aParams = LindstromTurk_params()) 
  : analytics(a), mParams(aParams){}
     
  template <typename Profile> 
  //optional<typename Profile::Point> //probably need a profile for later
  optional<typename Profile::Point> operator()( Profile const& aProfile ) const
  {

    int area_max = 20;

    /*if(analytics.get_transition_area(aProfile.v0()) != -1 && analytics.get_transition_area(aProfile.v1()) != -1){

      if (analytics.get_transition_area(aProfile.v0()) == area_max){

        analytics.set_transition_area(aProfile.v1(), 20);

      }else if(analytics.get_transition_area(aProfile.v1()) == area_max){
        analytics.set_transition_area(aProfile.v0(), 20);

      }
    }*/

    return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_placement() ;
  }
  
  private:

  LindstromTurk_params mParams ;    


};

} // namespace Surface_mesh_simplification

} //namespace CGAL















