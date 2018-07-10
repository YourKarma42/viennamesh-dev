
#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
//#include <CGAL/cgal_statistic_funktions.cpp>
//#include "cgal_mesh.hpp"

#include "../cgal_curvature_functions.hpp"


#include <iostream>
#include <fstream>

#include<chrono>

//hash table to store curvatures
#include <unordered_map>

#include <CGAL/Curvature_policies/cgal_hash_points.hpp>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>



namespace CGAL {

namespace Surface_mesh_simplification
{


  template<class ECM_>
class my_Curvature_cost
{
  
public:

  typedef ECM_ ECM;

  typedef typename ECM::Vertex::Halfedge_around_vertex_circulator Vertex_circulator;

  typedef typename ECM::Point_3 Point_3;

  typedef typename ECM::Vertex Vertex;


//TODO: nicht gut überlegen ob mans anders machen kann
  typedef CGAL::Simple_cartesian<double> Kernel;
  typedef Kernel::Vector_3 Vector_3;



  std::unordered_map<Vertex, double>   & curvatures;

  std::vector<double> & times;


public:

 int b;

  my_Curvature_cost(std::unordered_map<Point_3, std::unordered_map<Vertex, double>> & c, 
  std::vector<double> &t,
  int & ee) 
  : curvatures(c), times(t){}

  //Curvature_cost(int a) : b(a) {}

  //Curvature_cost() {}

  template <typename Profile, typename T> 
  optional<typename Profile::FT> operator()( Profile const& aProfile, T const& /*aPlacement*/ ) const
  {
    typedef optional<typename Profile::FT> result_type;


    //costs File is deleted in cgal_simplify_curve.cpp


    double flat_boundary = 0.00001;

    double curve_p0 = 0.0;
    double curve_p1 = 0.0;

    double mod = 0.0;


    double curves_p0[2];
    double curves_p1[2];

/*
   
    //kann man eventuell wegtun
    ECM& mesh = aProfile.surface_mesh();


    // nach dem collapse ändern sich curvatures  um den collapsten punkt
    // cgal berücksichtigt das
  

    //auto start = std::chrono::high_resolution_clock::now();


    auto current_point = curvatures.find(aProfile.p0());
    if(current_point != curvatures.end()){
        //curvature was already calculated

        

        curves_p0[0] = (current_point->second).curvature_1;

        curves_p0[1] = (current_point->second).curvature_2; 

        curve_p0 =  curves_p0[0]*curves_p0[1]; 

    }
    else{
        //calculate curvature and save it in the hashtable

        times.push_back(viennamesh::cgal::principal_curvatures(*(aProfile.v0()), mesh, curves_p0));

        //viennamesh::cgal::principal_curvatures(*(aProfile.v0()), mesh, curves_p0);

        curve_p0=curves_p0[0]*curves_p0[1];

        curvatures[aProfile.p0()].curvature_1 = curves_p0[0];
        curvatures[aProfile.p0()].curvature_2 = curves_p0[1];

    }  
    
    current_point = curvatures.find(aProfile.p1());
    if(current_point != curvatures.end()){
        //curvature was already calculated

        curves_p1[0] = (current_point->second).curvature_1;

        curves_p1[1] = (current_point->second).curvature_2;

        curve_p1=curves_p1[0]*curves_p1[1];

    }
    else{
        //calculate curvature and save it in the hashtable

        times.push_back(viennamesh::cgal::principal_curvatures(*(aProfile.v1()), mesh, curves_p1));

        //viennamesh::cgal::principal_curvatures(*(aProfile.v1()), mesh, curves_p1);

        curve_p1=curves_p1[0]*curves_p1[1];

        curvatures[aProfile.p1()].curvature_1 = curves_p1[0];
        curvatures[aProfile.p1()].curvature_2 = curves_p1[1];

      }  
      

      */

    //auto finish = std::chrono::high_resolution_clock::now();

    //times.push_back(std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count());



    //the smaller the cost the earlier its coursed

    double c1 = viennamesh::cgal::mean_curvature_my(*(aProfile.v0()));

    double c2 = viennamesh::cgal::mean_curvature_my(*(aProfile.v1()));

    if ((c1 <= flat_boundary) && (c2 <= flat_boundary)){

    }

    return result_type(c1+c2);


  }
  
};


} // namespace Surface_mesh_simplification


} //namespace CGAL

 


