
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



  viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;


public:

 int b;

  Curvature_flat_cost(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a) 
  : analytics(a) {}


  template <typename Profile, typename T> 
  optional<typename Profile::FT> operator()( Profile const& aProfile, T const& /*aPlacement*/ ) const
  {
    typedef optional<typename Profile::FT> result_type;



    double flat_boundary = 0.001;

    double c1;

    double c2;

    // values only for testing

    //std::cout<< analytics.get_mean(aProfile.v0()) << std::endl;
    

    if(analytics.has_curvature(aProfile.v0())){
        c1 = analytics.get_mean(aProfile.v0());
        //std::cout<< "blubb" << std::endl;
    }else{
        c1=0;
        std::cout<< "blubb" << std::endl;
    }

    if(analytics.has_curvature(aProfile.v1())){
        c2 = analytics.get_mean(aProfile.v1());
        //std::cout<< "blubb" << std::endl;
    }else{
        c2=0;
        std::cout<< "blubb" << std::endl;
    }

    if(c1 < flat_boundary && c2 < flat_boundary){
       //std::cout << CGAL::squared_distance(aProfile.p0(), aProfile.p1()) << std::endl;
        return CGAL::squared_distance(aProfile.p0(), aProfile.p1());

    }
      


    return result_type(100);


  }
  
};


} // namespace Surface_mesh_simplification


} //namespace CGAL

 


