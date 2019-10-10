
#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
//#include <CGAL/cgal_statistic_funktions.cpp>
//#include "cgal_mesh.hpp"




#include <iostream>
#include <fstream>

#include<chrono>


#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

//analytics
#include "../../analizing/cgal_mesh_analytics.hpp"



namespace CGAL {

namespace Surface_mesh_simplification
{


template<class ECM_>
class Curvature_flat_cost
{
  
public:

  typedef ECM_ ECM;

  typedef typename ECM::Vertex::Halfedge_around_vertex_circulator Vertex_circulator;

  typedef typename ECM::Point_3 Point_3;


//TODO: nicht gut überlegen ob mans anders machen kann
  typedef CGAL::Simple_cartesian<double> Kernel;
  typedef Kernel::Vector_3 Vector_3;




public:

 int b;

  Curvature_flat_cost(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a, LindstromTurk_params const& aParams = LindstromTurk_params()) 
  :  analytics(a) , mParams(aParams) {}


  template <typename Profile, typename T> 
  optional<typename Profile::FT> operator()( Profile const& aProfile, T const& aPlacement ) const
  {
    typedef optional<typename Profile::FT> result_type;



    double flat_boundary = 0.001;

    double c1;

    double c2;

  
    if(analytics.get_transition_distance(aProfile.v0()) >= max_size && analytics.get_transition_distance(aProfile.v1()) >= max_size) {
        return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement);
    }
    

    return result_type(100);


  }

private:
    viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;

    LindstromTurk_params mParams ; 

    double max_size;

  
};


} // namespace Surface_mesh_simplification


} //namespace CGAL

 


