#include <CGAL/Surface_mesh_simplification/Detail/Common.h>

//hash table to store curvatures
#include <unordered_map>


#include <iostream>
#include <fstream>

#include "../cgal_curvature_functions.hpp"

#include <CGAL/Curvature_policies/cgal_hash_points.hpp>


namespace CGAL {

namespace Surface_mesh_simplification  
{

  template<class ECM_>
class Curvature_placement
{
public:
    
  typedef ECM_ ECM ;

  typedef typename ECM::Point_3 Point_3;



//TODO: nicht gut überlegen ob mans anders machen kann
  typedef CGAL::Simple_cartesian<double> Kernel;
  typedef Kernel::Vector_3 Vector_3;




  std::unordered_map<Point_3, viennamesh::cgal::vertex_statistics, viennamesh::cgal::Point_3_Hash, viennamesh::cgal::Point_3_Equal>   & curvatures;
  
  int & ecke;

public:



  Curvature_placement(std::unordered_map<Point_3, viennamesh::cgal::vertex_statistics, viennamesh::cgal::Point_3_Hash, viennamesh::cgal::Point_3_Equal> & c, 
    int & e)
  : curvatures(c), ecke(e) {}
     
  template <typename Profile> 
  //optional<typename Profile::Point> //probably need a profile for later
  optional<typename Profile::Point> operator()( Profile const& aProfile ) const
  {
    //probably delete old vertices from the hastable


    double flat_boundary = 0.001;

  
 
    auto v0 = curvatures.find(aProfile.p0());

    auto v1 = curvatures.find(aProfile.p1());


//___________CGAL calculates placements in the collection phase ???____________________
//___________CGAL dosnt calcuclate costs of new Vectors???_____________________________

//hack not nice!

    double curves_p0[2];
    double curves_p1[2];

    if(v0 == curvatures.end()){

        viennamesh::cgal::principal_curvatures(*(aProfile.v0()),  aProfile.surface_mesh(), curves_p0);




        curvatures[aProfile.p0()].curvature_1 = curves_p0[0];
        curvatures[aProfile.p0()].curvature_2 = curves_p0[1];

        v0 = curvatures.find(aProfile.p0());

               //return optional<typename Profile::Point>(Point_3(0,0,0)) ;

    }

    if(v1 == curvatures.end()){
        viennamesh::cgal::principal_curvatures(*(aProfile.v1()),  aProfile.surface_mesh(), curves_p1);



        curvatures[aProfile.p1()].curvature_1 = curves_p1[0];
        curvatures[aProfile.p1()].curvature_2 = curves_p1[1];

        v1 = curvatures.find(aProfile.p1());

                 //return optional<typename Profile::Point>(Point_3(0,0,0)) ;

    }
//______________________________________________________________________________________

    /*if(fabs((v0->second).curvature_1) > flat_boundary)
      ecke++; 

    if(fabs((v0->second).curvature_2> flat_boundary))
      ecke++;

    if(fabs((v1->second).curvature_1> flat_boundary))
      ecke++;

    if(fabs((v1->second).curvature_2> flat_boundary))
      ecke++;*/

      //ecke++;

            //__DEBUG___

        
          /*std::ofstream myfile;
          myfile.open ("costs.txt", std::ios_base::app);
          if (myfile.is_open()){

            myfile <<  

            "Punkt p0: " <<  aProfile.p0() << std::endl <<

            "Punkt p1: " <<  aProfile.p1() << std::endl << std::endl;
  
            myfile.close();
          }*/

    if (fabs((v0->second).curvature_1) < flat_boundary && fabs((v0->second).curvature_2) < flat_boundary){

      if (fabs((v1->second).curvature_1) < flat_boundary && fabs((v1->second).curvature_2) < flat_boundary){
        // "flat" edge so midpoint

        //ecke++;
        return optional<typename Profile::Point>(midpoint(aProfile.p0(),aProfile.p1())) ;

      }else{

        //ecke++;

        // v1 is not flat
        return optional<typename Profile::Point>(aProfile.p1()) ;

      }
    }else{

      if (fabs((v1->second).curvature_1) < flat_boundary && fabs((v1->second).curvature_2) < flat_boundary){

        //ecke++;

        // v0 is not flat
        return optional<typename Profile::Point>(aProfile.p0()) ;


      }else{
        
        //ecke++;
        // edge is part of the boundary
        
        //__DEBUG___


        //Face test = face(aProfile.v0_v1(), aProfile.surface_mesh());

        double epsilon_normal = 0.9;

        //(aProfile.v0_v1())->normal()


        Vector_3 normal_left_face = CGAL::normal((aProfile.v0_v1())->vertex()->point(), 
                                                (aProfile.v0_v1())->next()->vertex()->point(), 
                                                (aProfile.v0_v1())->next()->next()->vertex()->point());

        Vector_3 normal_right_face = CGAL::normal((aProfile.v1_v0())->vertex()->point(), 
                                                 (aProfile.v1_v0())->next()->vertex()->point(), 
                                                 (aProfile.v1_v0())->next()->next()->vertex()->point()); 

                    Vector_3 a = normal_right_face;

            Vector_3 b = normal_left_face;

        //double ang = CGAL::angle(normal_right_face, normal_left_face);                                        


        // Abfrage genau überlegen

        if(((a.x()*b.x()+a.y()*b.y()+a.z()*b.z()) /
           (std::sqrt(a.x()*a.x()+a.y()*a.y()+a.z()*a.z())*std::sqrt(b.x()*b.x()+b.y()*b.y()+b.z()*b.z())))
           > epsilon_normal){

          auto e0 = curvatures.find(aProfile.v0_v1()->next()->vertex()->point());
          auto e1 = curvatures.find(aProfile.v1_v0()->next()->vertex()->point());
          //___________CGAL calculates placements in the collection phase ???____________________
          //___________CGAL dosnt calcuclate costs of new Vectors???_____________________________

          //hack not nice!

              double curves_e0[2];
              double curves_e1[2];

              if(e0 == curvatures.end()){

                  viennamesh::cgal::principal_curvatures(*(aProfile.v0_v1()->next()->vertex()),  aProfile.surface_mesh(), curves_e0);

                  curvatures[aProfile.p0()].curvature_1 = curves_e0[0];
                  curvatures[aProfile.p0()].curvature_2 = curves_e0[1];

                  e0 = curvatures.find(aProfile.p0());

              }

              if(e1 == curvatures.end()){
                  viennamesh::cgal::principal_curvatures(*(aProfile.v1_v0()->next()->vertex()),  aProfile.surface_mesh(), curves_e1);

                  curvatures[aProfile.p1()].curvature_1 = curves_e1[0];
                  curvatures[aProfile.p1()].curvature_2 = curves_e1[1];

                  e1 = curvatures.find(aProfile.p1());


              }
          //______________________________________________________________________________________


         
          /*std::ofstream myfile;
          myfile.open ("costs.txt", std::ios_base::app);
          if (myfile.is_open()){

            myfile <<  

            "Punkt p0: " <<  aProfile.p0() << std::endl <<

            "Punkt p1: " <<  aProfile.p1() << std::endl <<

            "ECKE v0_v1:" <<  aProfile.v0_v1()->next()->vertex()->point() << std::endl <<

            "ECKE v1_v0:" <<  aProfile.v1_v0()->next()->vertex()->point() << std::endl <<
            
            "normalvektore: " <<  normal_left_face << std::endl <<
            
            "Krümmungen p0: " << fabs((e0->second).curvature_1) << "  " << fabs((e0->second).curvature_2) << std::endl <<

            "Krümmungen p1: " << fabs((e1->second).curvature_1) << "  " << fabs((e1->second).curvature_2) << std::endl;



            myfile.close();
          }*/


          if(fabs((e0->second).curvature_1) < flat_boundary && fabs((e0->second).curvature_2) < flat_boundary){
            return optional<typename Profile::Point>(aProfile.v1_v0()->next()->vertex()->point()) ;
          }else{
            return optional<typename Profile::Point>(aProfile.v0_v1()->next()->vertex()->point()) ;
          }

         
          
        }
        //__DEBUG___

        
          std::ofstream myfile;
          myfile.open ("costs.txt", std::ios_base::app);
          if (myfile.is_open()){

            Vector_3 a = normal_right_face;

            Vector_3 b = normal_left_face;


            myfile <<  

            "Punkt p0: " <<  aProfile.p0() << std::endl <<

            "Punkt p1: " <<  aProfile.p1() << std::endl <<

            "ECKE v0_v1:" << aProfile.v0_v1()->next()->vertex()->point() << " " 
                          << aProfile.v0_v1()->next()->next()->vertex()->point() << " "
                          << aProfile.v0_v1()->next()->next()->next()->vertex()->point()  << std::endl <<

            "ECKE v1_v0:" << aProfile.v0_v1()->opposite()->next()->vertex()->point() << " " 
                          << aProfile.v0_v1()->opposite()->next()->next()->vertex()->point() << " "
                          << aProfile.v0_v1()->opposite()->next()->next()->next()->vertex()->point()<< std::endl <<


            "neuer punkt: " << midpoint(aProfile.p0(),aProfile.p1()) << std::endl <<

            "normalvektore right: " << normal_right_face << std::endl <<
            
            "normalvektore left: " <<  normal_left_face << std::endl;

            
            
            
            myfile.close();
          }

          

        
        
        return optional<typename Profile::Point>(midpoint(aProfile.p0(),aProfile.p1())) ;


      }

    }
    



    //return optional<typename Profile::Point>(midpoint(aProfile.p0(),aProfile.p1())) ;
  }
  
      // delete just tmp
  double tmp_skp(Vector_3 a, Vector_3 b){

    return (a.x()*b.x()+a.y()*b.y()+a.z()*b.z())/
    (sqrt(a.x()+a.y()+a.z())*sqrt(b.x()+b.y()+b.z()));
  }

/*private:

  LindstromTurk_params mParams ;    */

};

} // namespace Surface_mesh_simplification

} //namespace CGAL

//TODO:  #defines machen wie in cgal