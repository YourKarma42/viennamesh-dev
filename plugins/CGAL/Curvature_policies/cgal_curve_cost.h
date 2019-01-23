
#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
//#include <CGAL/cgal_statistic_funktions.cpp>
//#include "cgal_mesh.hpp"

//#include "../cgal_curvature_functions.hpp"


#include <iostream>
#include <fstream>

#include<chrono>

//hash table to store curvatures (delete)
#include <unordered_map>


#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

//analytics
#include "../analizing/cgal_mesh_analytics.hpp"



namespace CGAL {

namespace Surface_mesh_simplification
{


  template<class ECM_>
class Curvature_cost
{
  
public:

  typedef ECM_ ECM;

  typedef typename ECM::Vertex::Halfedge_around_vertex_circulator Vertex_circulator;

  typedef typename ECM::Point_3 Point_3;


//TODO: nicht gut überlegen ob mans anders machen kann
  typedef CGAL::Simple_cartesian<double> Kernel;
  typedef Kernel::Vector_3 Vector_3;



  viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;

    LindstromTurk_params mParams ;  



public:

 int b;

  Curvature_cost(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a, 
  LindstromTurk_params const& aParams = LindstromTurk_params()) 
  : analytics(a),  mParams(aParams) {}

  //Curvature_cost(int a) : b(a) {}

  //Curvature_cost() {}

  template <typename Profile, typename T> 
  optional<typename Profile::FT> operator()( Profile const& aProfile, T const& aPlacement ) const
  {
    typedef optional<typename Profile::FT> result_type;


    //costs File is deleted in cgal_simplify_curve.cpp

    if((analytics.get_feature(aProfile.v0())!= -1 || analytics.get_feature(aProfile.v1())!= -1)){
       //std::cout << "blub" << std::endl;
        return LindstromTurkCore<ECM,Profile>(mParams,aProfile).compute_cost(aPlacement) ; 
    }


  return result_type(1000);

  }
  
};


} // namespace Surface_mesh_simplification


} //namespace CGAL

 


