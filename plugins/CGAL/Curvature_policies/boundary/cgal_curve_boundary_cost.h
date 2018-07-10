#include <CGAL/Surface_mesh_simplification/Detail/Common.h>


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
class Curvature_boundary_cost
{
  
public:

  typedef ECM_ ECM;

  typedef typename ECM::Vertex::Halfedge_around_vertex_circulator Vertex_circulator;

  typedef typename ECM::Point_3 Point_3;

  viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;

 

public:

 int b;

  Curvature_boundary_cost(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a) 
  : analytics(a) {}

  //Curvature_cost(int a) : b(a) {}

  //Curvature_cost() {}

  template <typename Profile, typename T> 
  optional<typename Profile::FT> operator()( Profile const& aProfile, T const& /*aPlacement*/ ) const
  {
    typedef optional<typename Profile::FT> result_type;

    //for the moment a dummy value
    double max_value = 100;


    //if < 0 then the vertex is not part of the boundary
    if(analytics.get_color(aProfile.v0()) != -1 && analytics.get_color(aProfile.v1()) != -1){

      if(analytics.get_color(aProfile.v0()) == analytics.get_color(aProfile.v1())){

        if(analytics.get_color(aProfile.v0()) % 2 == 1)
          return result_type(50);

      }
        

    }
 
    return result_type(max_value);


  }
  
};


} // namespace Surface_mesh_simplification


} //namespace CGAL

 


