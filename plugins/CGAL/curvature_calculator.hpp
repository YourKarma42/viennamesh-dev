#ifndef CURVATURE_CALCULATOR
#define CURVATURE_CALCULATOR

#include "cgal_mesh.hpp"

namespace viennamesh
{
    namespace cgal
    {

        typedef typename viennamesh::cgal::polyhedron_surface_mesh::Vertex Vertex;
        typedef typename viennamesh::cgal::polyhedron_surface_mesh::Vertex::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

        struct Triangle{
            double alpha;
            double beta;
            double theta;

            double len_at_v;
            double len_at_prev;
            double len_prev_v;

            Vector_3 at_v;
            Vector_3 at_prev;
            Vector_3 prev_v;

            bool theta_border;
        };

        struct Curves{
            double gauss;
            double mean;
            double p1;
            double p2;
        };


        double const edge_curvature = 100.0;



        // for testing
        std::string my_curvature_test_print(Vertex& v);

        //test with sincos 
        double cot(double angle);

        std::list<Triangle> calc_triangle_angles_sides(Vertex& v);

        double triangle_area(Triangle &triangle);

        double mixed_area(std::list<Triangle> &triangles);


        double mean_curvature(std::list<Triangle> &triangles,  double area  );

        double gauss_curvature(std::list<Triangle> &triangles, double area );

        double k1(double gauss_curve, double mean_curve);

        double k2(double gauss_curve, double mean_curve);


        Curves calc_curvatures(Vertex& v);

        void calc_principle_curveatures(Vertex& v, double curvatures[]);


        // for testing
        std::string my_curvature_test_print(Vertex& v);

    }
}




#endif /*CURVATURE_CALCULATOR*/