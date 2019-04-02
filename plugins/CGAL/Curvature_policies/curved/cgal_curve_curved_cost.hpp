
#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Detail/Lindstrom_Turk_core.h>
//#include <CGAL/cgal_statistic_funktions.cpp>
//#include "cgal_mesh.hpp"




#include <iostream>
#include <fstream>

#include<chrono>


#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

//analytics
#include "../../analizing/cgal_mesh_analytics.hpp"
//for testing
#include "../../curvature_calculator.hpp"



namespace CGAL {

namespace Surface_mesh_simplification
{


template<class ECM_>
class Curvature_curved_cost
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

  Curvature_curved_cost(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a, LindstromTurk_params const& aParams = LindstromTurk_params()) 
  : analytics(a) , mParams(aParams) {}


  template <typename Profile, typename T> 
  optional<typename Profile::FT> operator()( Profile const& aProfile, T const& aPlacement ) const
  {
    typedef optional<typename Profile::FT> result_type;

    bool coars_features = true;



    // values only for testing

    //std::cout<< analytics.get_mean(aProfile.v0()) << std::endl;

    //viennamesh::cgal::Curves c1 = viennamesh::cgal::calc_curvatures(*aProfile.v0());
    
    //viennamesh::cgal::Curves c2 = viennamesh::cgal::calc_curvatures(*aProfile.v1());

    //if(analytics.get_curved_area(aProfile.v0())!= -1 && analytics.get_curved_area(aProfile.v1())!= -1){
    // for testing use feature
   /* if((analytics.get_feature(aProfile.v0())!= -1 && analytics.get_feature(aProfile.v1())!= -1)){
       //std::cout << LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement) << std::endl;
        return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement) ; 
    }

    if(coars_features){
  
      if((analytics.get_feature(aProfile.v0())!= -1 && analytics.get_transition_area(aProfile.v1())!= -1) ||
      (analytics.get_feature(aProfile.v1())!= -1 && analytics.get_transition_area(aProfile.v0())!= -1)){
        return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement) ; 
      }
    }*/

    return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement) ; 

    return result_type(1000);


  }

  private:

  LindstromTurk_params mParams ; 

  viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;
  
};


} // namespace Surface_mesh_simplification


} //namespace CGAL

 


