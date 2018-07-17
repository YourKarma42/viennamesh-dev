

/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                ViennaMesh - The Vienna Meshing Framework
                            -----------------

                    http://viennamesh.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

//#include "viennameshpp/plugin.hpp"

#ifndef CURVE_FUNC_T
#define CURVE_FUNC_T


#define CGAL_EIGEN3_ENABLED
#include <CGAL/Monge_via_jet_fitting.h>

#include<chrono>

#include <sstream>

#include <unordered_map>

#include <CGAL/Curvature_policies/cgal_hash_points.hpp>

//curvature testfunction
#include <CGAL/curve_calc_dev/curve_test_1.hpp>

//TODO: Rerwrite with template types!


namespace viennamesh
{
  namespace cgal
  {
    typedef CGAL::Simple_cartesian<double> Kernel;
    typedef CGAL::Polyhedron_3<Kernel> mesh_t;
    typedef mesh_t::Point_3 Point_3;
    typedef Kernel::Vector_3 Vector_3;
    typedef mesh_t::Face_handle Face;
    typedef mesh_t::Facet Facet_t;
    typedef CGAL::Monge_via_jet_fitting<Kernel>::Monge_form Monge_form;

    Vector_3 exproduct(Vector_3 a,Vector_3 b)
    {
        return Vector_3(a[1]*b[2]-a[2]*b[1],-a[0]*b[2]+a[2]*b[0],a[0]*b[1]-a[1]*b[0]);
    }

    double norm(Vector_3 e)
    {
        return std::sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
    }

    Vector_3 get_surface_normal(mesh_t::Vertex& vertex)
    {
        // überdenken wegen hash funktion
        // irgendwie normalen am rand definieren? annehmen das es eben weitergeht?
        Vector_3 out(0,0,0),normal;
        mesh_t::Vertex::Halfedge_around_vertex_circulator at=vertex.vertex_begin(),end=at;
        
        do
        {
            if(! at->is_border_edge())
            {
                Vector_3 a=at->opposite()->vertex()->point() - at->vertex()->point();
                Vector_3 b=at->next()->opposite()->vertex()->point() - at->next()->vertex()->point();
                out=out+exproduct(a,b);
            }
            ++at;
        }while(at!=end);
    
        if(norm(out)==0)
            return out;     //if it is not defined just return (0,0,0)
    
        return out/norm(out);
    }   

//_____________________My Stuff____________________________________________________________________________________________________

    //TODO: ANSCHAUN UND GENAU verstehen
    //gathers the closest points by collecting the neighbors and there neighbors and so on if nessesary
    std::vector<Point_3> closest_points(long amount, mesh_t::Vertex& vertex,mesh_t& mesh)
    {
      std::vector<Point_3> out,tmp;
      std::vector<mesh_t::Vertex> new_vertexes;
      //new_vertexes.reserve(40);
      long id=0,id_end;
    
      //std::cout << "test\n";
      out.push_back(vertex.point());
      new_vertexes.push_back(vertex);

      while(out.size()<amount)
      {

          //muss ein long "iterator"
          id_end=new_vertexes.size();
          while(id!=id_end)
          {
              mesh_t::Vertex::Halfedge_around_vertex_circulator at=new_vertexes.at(id).vertex_begin(),end=at;
              do
              {
                
                  if(std::find(out.begin(),out.end(),at->opposite()->vertex()->point())==out.end())
                  {
                      out.push_back(at->opposite()->vertex()->point());
                      //std::cout << "now adding the vertex\n";
                      new_vertexes.push_back(* (at->opposite()->vertex()));
                  }
              
                  ++at;
              }while(at!=end);
              ++id;
          }
      }

      return out;
    }    



    Monge_form get_monge_form (mesh_t::Vertex& vertex ,mesh_t& mesh)
    {
      unsigned int d_fitting = 2;
      unsigned int d_monge = 2;
      unsigned int needed = 6;
      
      Vector_3  normal=get_surface_normal(vertex);

      std::vector<Point_3> points=closest_points(needed,vertex,mesh);              //<---- These Points
      if(points.size()<needed)
      {
          std::cerr << "not enough points: " << points.size() <<"\n";
          // return 0;
      }

      CGAL::Monge_via_jet_fitting<Kernel> monge_fit;
      //monge fit anschaun wirklich points?
      CGAL::Monge_via_jet_fitting<Kernel>::Monge_form monge_form 
              = monge_fit(points.begin(),points.end(),d_fitting, d_monge);
              
      monge_form.comply_wrt_given_normal(normal);

      return monge_form;
    }



    //my

    //needs array with 2 fields of floating type
    double principal_curvatures_my(mesh_t::Vertex& vertex, double  curvatures[]){

        auto start = std::chrono::high_resolution_clock::now();

        
        calc_principle_curveatures(vertex, curvatures);

        //curvatures[0] = 0;  // max
        //curvatures[1] = 0; // min

        auto finish = std::chrono::high_resolution_clock::now();



        return std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();

       // return 0;

    }

    double mean_curvature_my(mesh_t::Vertex& vertex){


        return mean_test(vertex);

        

       // return 0;

    }

    //for testint delte
    double gauss_curvatures_my(mesh_t::Vertex& vertex){

        return gauss_test(vertex);

       // return 0;

    }


    //generic types verwenden
    //needs array with 2 fields of floating type
    double principal_curvatures_cgal(mesh_t::Vertex& vertex, mesh_t& mesh, double  curvatures[]){

        auto start = std::chrono::high_resolution_clock::now();

        Monge_form monge_form=get_monge_form (vertex,mesh);


        curvatures[0] = monge_form.principal_curvatures(0); // max
        curvatures[1] = monge_form.principal_curvatures(1); // min

        auto finish = std::chrono::high_resolution_clock::now();



        return std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();

       // return 0;

    }


    double principal_curvatures(mesh_t::Vertex& vertex, mesh_t& mesh, double  curvatures[]){

         //return principal_curvatures_cgal(vertex, mesh, curvatures);

         return principal_curvatures_my(vertex,  curvatures);

    }

    std::string curvature_output(std::unordered_map<Point_3, viennamesh::cgal::vertex_statistics, viennamesh::cgal::Point_3_Hash, viennamesh::cgal::Point_3_Equal> & m){
        std::stringstream  out;

        double max_1 = 0.0;
        double max_2 = 0.0;

        double min_1 = 0.0;
        double min_2 = 0.0;

        for ( auto &it : m ){
            if((it.second).curvature_1 > max_1)
                max_1 = (it.second).curvature_1;
            
            if((it.second).curvature_2 > max_2)
                max_2 = (it.second).curvature_2;

            if((it.second).curvature_1 < min_1)
                min_1 = (it.second).curvature_1;
            
            if((it.second).curvature_2 < min_2)
                min_2 = (it.second).curvature_2;
        }

        out << std::endl << "----------------++++++++++++++++++ curvatures ++++++++++++++++++----------------" << std::endl;

        out << "Max abs curve 1: " << max_1 << std::endl;

        out << "Max abs curve 2: " << max_2 << std::endl;

        out << "Min abs curve 1: " << min_1 << std::endl;

        out << "Min abs curve 2: " << min_2 << std::endl;


        out << "--------------------------------------------------------------------------------" << std::endl;


        return out.str();

    }

//_____________________________________ END MY STUFF _________________________________________________________________

   




    mesh_t::Vertex& convert_Point_to_Vertex(Point_3 point,mesh_t & mesh)
    {
        //TODO: anschaun eventuell ändern iteriert ganzen mesh?
        for(mesh_t::Vertex_iterator at=mesh.vertices_begin(),end=mesh.vertices_end();at!=end;++at)
            if(at->point()==point)
                return *at;
        
        std::cerr << "nothing found maby try a diffrent mesh\n";
        static mesh_t::Vertex v(Point_3(0,0,0));
        return v;
    }  

    


    //gathers the closest points by collecting the neighbors and there neighbors and so on if nessesary
    std::vector<Point_3> closest_points_old(long amount,Point_3 point,mesh_t& mesh)
    {
      mesh_t::Vertex vertex=convert_Point_to_Vertex(point,mesh);
      std::vector<Point_3> out,tmp;
      std::vector<mesh_t::Vertex> new_vertexes;
      //new_vertexes.reserve(40);
      long id=0,id_end;
    
      //std::cout << "test\n";
      out.push_back(point);
      new_vertexes.push_back(vertex);

      while(out.size()<amount)
      {

          //muss ein long "iterator"
          id_end=new_vertexes.size();
          while(id!=id_end)
          {
              mesh_t::Vertex::Halfedge_around_vertex_circulator at=new_vertexes.at(id).vertex_begin(),end=at;
              do
              {
                
                  if(std::find(out.begin(),out.end(),at->opposite()->vertex()->point())==out.end())
                  {
                      out.push_back(at->opposite()->vertex()->point());
                      //std::cout << "now adding the vertex\n";
                      new_vertexes.push_back(* (at->opposite()->vertex()));
                  }
              
                  ++at;
              }while(at!=end);
              ++id;
          }
      }

      return out;
    }
    

    Monge_form get_monge_form_old (Point_3 point,mesh_t& mesh)
    {
      unsigned int d_fitting = 2;
      unsigned int d_monge = 2;
      unsigned int needed = 6;
      
      Vector_3  normal=get_surface_normal(convert_Point_to_Vertex(point,mesh));

      std::vector<Point_3> points=closest_points_old(needed,point,mesh);              //<---- These Points
      if(points.size()<needed)
      {
          std::cerr << "not enough points: " << points.size() <<"\n";
          // return 0;
      }

      CGAL::Monge_via_jet_fitting<Kernel> monge_fit;
      CGAL::Monge_via_jet_fitting<Kernel>::Monge_form monge_form 
              = monge_fit(points.begin(),points.end(),d_fitting, d_monge);
              
      monge_form.comply_wrt_given_normal(normal);

      return monge_form;
    }


    double mean_curvature(Point_3 point,mesh_t& mesh)
    {
      Monge_form monge_form=get_monge_form_old (point,mesh);
      //monge_form.principal_curvatures(0 is maximum 1 is minimum)
      return (monge_form.principal_curvatures(0)+monge_form.principal_curvatures(1))/2;
    }

    double gaus_curvature(Point_3 point,mesh_t& mesh)        
    {
      Monge_form monge_form=get_monge_form_old (point,mesh);
      return monge_form.principal_curvatures(0)*monge_form.principal_curvatures(1);
    }

    
    double max_curvature(Point_3 point,mesh_t& mesh)
    {
      Monge_form monge_form=get_monge_form_old (point,mesh);
      return monge_form.principal_curvatures(0);
    }
   
    double min_curvature(Point_3 point,mesh_t& mesh)
    {
      Monge_form monge_form=get_monge_form_old (point,mesh);
      return monge_form.principal_curvatures(1);
    }

    double t()
    {
      return 1.0;
    }
   
  }

}

#endif /* CURVE_FUNC_T */


