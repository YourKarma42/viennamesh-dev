#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core.h>


#include <iostream>
#include <fstream>

#include<chrono>

//hash table to store curvatures (delete)
#include <unordered_map>


#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

//analytics
#include "../../analizing/cgal_mesh_analytics.hpp"



namespace CGAL {

namespace Surface_mesh_simplification
{


  template<class ECM_>
class Curvature_features_cost
{
  
public:

  //LindstromTurk_cost( LindstromTurk_params const& aParams = LindstromTurk_params() ) : mParams(aParams) {}

  typedef ECM_ ECM;

  typedef typename ECM::Vertex::Halfedge_around_vertex_circulator Vertex_circulator;

  typedef typename ECM::Point_3 Point_3;

  typedef typename ECM::Facet_handle Facet_handle;

  typedef typename viennamesh::cgal::Vector_3 Vector_3;



 

public:

 int b;

   double max_size;

  Curvature_features_cost(double ms,viennamesh::cgal::cgal_mesh_analytics<ECM>  & a,  LindstromTurk_params const& aParams = LindstromTurk_params()) 
  : max_size(ms), analytics(a),  mParams(aParams) {}

  //Curvature_cost(int a) : b(a) {}

  //Curvature_cost() {}

  template <typename Profile> 
  optional<typename Profile::FT> 
  operator()( Profile const& aProfile, optional<typename Profile::Point> const& aPlacement ) const
  {
    typedef optional<typename Profile::FT> result_type;

    //for the moment a dummy value
    double max_value = 150;



    int area_v0 = analytics.get_transition_area(aProfile.v0());

    int area_v1 = analytics.get_transition_area(aProfile.v1());


//return CGAL::squared_distance(aProfile.p0(), aProfile.p1());



   /* if((area_v0 != -1 && area_v1 != -1))

    {
        //if(area_v0 < area_max && area_v1 < area_max){

          //if(analytics.get_transition_distance(aProfile.v0()) >= min_size && analytics.get_transition_distance(aProfile.v1()) >= min_size) 
            return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement);
       // }

    }*/

    //if(analytics.get_feature(aProfile.v0()) == -1 && analytics.get_feature(aProfile.v1()) == -1) 
    //    return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement);


    if(analytics.get_transition_distance(aProfile.v0()) >= max_size && analytics.get_transition_distance(aProfile.v1()) >= max_size) {
      return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement);
      //return CGAL::squared_distance(aProfile.p0(), aProfile.p1());
    }

 
    return result_type(max_value);


  }

private:
  LindstromTurk_params mParams ;



  viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;
  
};


} // namespace Surface_mesh_simplification


} //namespace CGAL

 

