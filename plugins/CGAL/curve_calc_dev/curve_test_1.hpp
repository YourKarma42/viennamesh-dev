




#ifndef CURVE_TEST
#define CURVE_TEST

#include <math.h>

#include <sstream>

namespace viennamesh
{
  namespace cgal
  {
    typedef CGAL::Simple_cartesian<double> Kernel;
    typedef CGAL::Polyhedron_3<Kernel> mesh_t;
    typedef mesh_t::Point_3 Point_3;
    typedef Kernel::Vector_3 Vector_3;
    typedef mesh_t::Vertex& Vertex;
    typedef mesh_t::Face_handle Face;
    typedef mesh_t::Facet Facet_t;


    struct Triangle{
        double alpha;
        double beta;
        double theta;

        double len_at_v;
        double len_at_prev;
        double len_prev_v;

        /*Vector_3 at_v;
        Vector_3 at_prev;
        Vector_3 prev_v;*/

    };


    double mixed_area(std::list<Triangle> &triangles);

    double mean_curvature(std::list<Triangle> &triangles, Vertex& v, double area  );

    double gauss_curvature(std::list<Triangle> &triangles, double area );

    double triangle_area(Triangle &triangle);

    double k1(double gauss_curve, double mean_curve);

    double k2(double gauss_curve, double mean_curve);


    double triangle_area(Point_3 p1,Point_3 p2,Point_3 p3){

        return 0.0;
    }

    //Triangle calc_angles(Vertex v1, );



/*
		    v		
		    /\
	       /  \
	      /    \
	     /theta \
	    /        \
	   /          \
	  /            \
	 / alpha    beta\  
	/________________\ 
   at               prev
*/


    std::string calc_curvature(Vertex& v){

        std::stringstream  out;

        out << std::endl << "----------------++++++++++++++++++ curvatures ++++++++++++++++++----------------" << std::endl;

        //eventuell Vektor
        std::list<Triangle> one_ring_neightbors;

        Vector_3 vec1, vec2, vec3;

        
        mesh_t::Vertex::Halfedge_around_vertex_circulator at=v.vertex_begin(), end = at;
        mesh_t::Vertex::Halfedge_around_vertex_circulator prev = at;

        vec1 = Vector_3(at->opposite()->vertex()->point(), v.point());

        ++at;
        
        do
        {
            vec2 = Vector_3(prev->opposite()->vertex()->point(), at->opposite()->vertex()->point());

            double norm1 = std::sqrt(vec1.squared_length());

            double norm2 = std::sqrt(vec2.squared_length());

            double tmp_beta = (vec1*(1/norm1))*(vec2*(1/norm2));
          
            vec2 = vec2*(-1);

            vec3 = Vector_3(at->opposite()->vertex()->point(), v.point());

            double norm3 = std::sqrt(vec3.squared_length());

            double tmp_alpha = (vec2*(1/norm2))*(vec3*(1/norm3));

            vec3 = vec3*(-1);

            vec1 = vec1*(-1);

            double tmp_theta = (vec1*(1/norm1))*(vec3*(1/norm3));

            Triangle newTriangle;
            newTriangle.alpha = tmp_alpha;
            newTriangle.beta = tmp_beta;
            newTriangle.theta = tmp_theta;
            newTriangle.len_prev_v = norm1;
            newTriangle.len_at_prev = norm2;
            newTriangle.len_at_v = norm3;
            

            one_ring_neightbors.push_back(newTriangle);

            prev = at;

            vec1 = vec3*(-1);
              
            ++at;

        }while(at != end);
     
        vec2 = Vector_3(prev->opposite()->vertex()->point(), at->opposite()->vertex()->point());

        double norm1 = std::sqrt(vec1.squared_length());

        double norm2 = std::sqrt(vec2.squared_length());

        double tmp_beta = (vec1*(1/norm1))*(vec2*(1/norm2));


        vec2 = vec2*(-1);

        vec3 = Vector_3(at->opposite()->vertex()->point(), v.point());

        double norm3 = std::sqrt(vec3.squared_length());

        double tmp_alpha = (vec2*(1/norm2))*(vec3*(1/norm3));

        
        vec3 = vec3*(-1);

        vec1 = vec1*(-1);


        double tmp_theta = (vec1*(1/norm1))*(vec3*(1/norm3));

        Triangle newTriangle;
        newTriangle.alpha = tmp_alpha;
        newTriangle.beta = tmp_beta;
        newTriangle.theta = tmp_theta;
        newTriangle.len_prev_v = norm1;
        newTriangle.len_at_prev = norm2;
        newTriangle.len_at_v = norm3;
        
        

        one_ring_neightbors.push_back(newTriangle);

       

        //out << v.point() << std::endl;

        out << one_ring_neightbors.size() << std::endl;

        for(auto t : one_ring_neightbors){
            out << t.alpha << std::endl;
            out << t.beta << std::endl;
            out << t.theta << std::endl;
            out << t.len_prev_v<< std::endl;
            out << t.len_at_prev<< std::endl;
            out << t.len_at_v<< std::endl;

            out <<  std::endl;

        }

        double area = mixed_area(one_ring_neightbors);

        double mean = mean_curvature(one_ring_neightbors, v, area);

        double gauss = gauss_curvature(one_ring_neightbors, area);

        out << "area: " << area << std::endl;

        out << "gaus: " << gauss << std::endl;

        out << "mean: " << mean << std::endl;


//******************************Testing*************************


        //paper formular (10)
        double k_1 = mean + std::sqrt(mean*mean - gauss);
        //paper formular (11)
        double k_2 = mean - std::sqrt(mean*mean - gauss);

        out << std::sqrt(mean*mean - gauss) << std::endl;

        out << "k1:  " << k_1 << std::endl;

        out << "k2:  " << k_2 << std::endl;




       /* Vector_3 mean_curve_normal_op;

        mesh_t::Vertex::Halfedge_around_vertex_circulator at_a=v.vertex_begin(), end_a = at_a;

        std::list<Triangle>::const_iterator tri_it = one_ring_neightbors.begin();

        Triangle pre = one_ring_neightbors.back();

        do {

            out <<
            tan(M_PI/2 -(*tri_it).alpha) << std::endl <<
            tan(M_PI/2 -pre.beta) << std::endl
            << std::endl;

            //paper formular (8)
            mean_curve_normal_op = mean_curve_normal_op +  (1/(2*area)) *
            (tan(M_PI/2 -(*tri_it).alpha) + tan(M_PI/2 -pre.beta)) * Vector_3(at_a->opposite()->vertex()->point(), v.point());

            pre = *tri_it;

            at_a ++;
            tri_it++;

        }while(at_a != end_a);

        out << mean_curve_normal_op << std::endl;

        out << std::sqrt(mean_curve_normal_op.squared_length())/2 << std::endl;*/

//***************************************************************


        return out.str();
    }

    double k1(double gauss_curve, double mean_curve){
        return 0;

    }

    double k2(double gauss_curve, double mean_curve){
        return 0;
    }


    double mean_curvature(std::list<Triangle> &triangles, Vertex& v, double area ){

        Vector_3 mean_curve_normal_op;

        mesh_t::Vertex::Halfedge_around_vertex_circulator at_a=v.vertex_begin(), end_a = at_a;

        std::list<Triangle>::const_iterator tri_it = triangles.begin();

        Triangle pre = triangles.back();

        do {

            //paper formular (8)
            mean_curve_normal_op = mean_curve_normal_op +  (1/(2*area)) *
            (tan(M_PI/2 -(*tri_it).alpha) + tan(M_PI/2 -pre.beta)) * Vector_3(at_a->opposite()->vertex()->point(), v.point());

            pre = *tri_it;

            at_a ++;
            tri_it++;

        }while(at_a != end_a);

        return std::sqrt(mean_curve_normal_op.squared_length())/2;
    }


    double gauss_curvature(std::list<Triangle> &triangles, double area ){

        double sum_theta=0;

        for(auto t : triangles)
            sum_theta +=  t.theta;

        //paper Formula (9)

        return (2*M_PI - sum_theta)/area;


    }

   
    double mixed_area(std::list<Triangle> &triangles){

        double area=0;

        for(auto t: triangles){

            if(t.theta <= 0){ //angle at v obtuse

                area += triangle_area(t)/2;

            }else if(t.alpha <= 0 || t.beta <= 0){ //angle at alpha or beta obtuse

                area += triangle_area(t)/4;

            }else{ // nonobtuse trinangle

                //Voronoi formula
                area += 0.125 * (tan(M_PI/2 - t.alpha)*t.len_prev_v*t.len_prev_v + tan(M_PI/2 - t.beta)*t.len_at_v*t.len_at_v);


            }           

        }

        return area;

    }

    //Herons Formular numericaly stable
    //https://en.wikipedia.org/wiki/Heron%27s_formula
    double triangle_area(Triangle &triangle){

        double a, b, c, tmp;

        a=triangle.len_at_v;
        b=triangle.len_at_prev;
        c=triangle.len_prev_v;

        if(a<b){
            tmp=b;
            b=a;
            a=tmp;
        }
        if(b<c){
            tmp=c;
            c=b;
            b=tmp;            
        }
        if(a<b){
            tmp=b;
            b=a;
            a=tmp;
        }

        return 0.25*std::sqrt((a+(b+c))*(c-(a-b))*(c+(a-b))*(a+(b-c)));
    }



  }
}


#endif /* CURVE_TEST */