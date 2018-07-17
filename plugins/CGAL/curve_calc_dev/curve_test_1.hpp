




#ifndef CURVE_TEST
#define CURVE_TEST

#include <math.h>

#include <sstream>

namespace viennamesh
{
  namespace cgal
  {
      //mit dem schon definierten Mesh machen
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

        Vector_3 at_v;
        Vector_3 at_prev;
        Vector_3 prev_v;

        bool theta_border;
    };

    double const edge_curvature = 100.0;

    //test with sincos 
    double cot(double angle){


        return cos(angle)/sin(angle);
    }

    std::list<Triangle> calc_triangle_angles_sides(Vertex& v);

    double triangle_area(Triangle &triangle);

    double mixed_area(std::list<Triangle> &triangles);

    double mean_curvature(std::list<Triangle> &triangles,  double area  );

    double gauss_curvature(std::list<Triangle> &triangles, double area );

    double k1(double gauss_curve, double mean_curve);

    double k2(double gauss_curve, double mean_curve);





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


    struct Curves{
        double gauss;
        double mean;
        double p1;
        double p2;
    };

    Curves calc_curvatures(Vertex& v){

        std::list<Triangle> one_ring_neightbors;
        one_ring_neightbors = calc_triangle_angles_sides(v);


        double area = mixed_area(one_ring_neightbors);

        double mean = mean_curvature(one_ring_neightbors, area);

        double gauss = gauss_curvature(one_ring_neightbors, area);

        double k_1 = k1(gauss, mean);

        double k_2 = k2(gauss, mean);

        return Curves{gauss, mean, k_1, k_2};




    }



    void calc_principle_curveatures(Vertex& v, double curvatures[]){

        std::list<Triangle> one_ring_neightbors;
        one_ring_neightbors = calc_triangle_angles_sides(v);


        double area = mixed_area(one_ring_neightbors);

        double mean = mean_curvature(one_ring_neightbors, area);

        double gauss = gauss_curvature(one_ring_neightbors, area);

        double k_1 = k1(gauss, mean);

        double k_2 = k2(gauss, mean);

        curvatures[0] = k_1;

        curvatures[1] = k_2;


    }

    //for testing delte
    double gauss_test(Vertex& v){

        std::list<Triangle> one_ring_neightbors;
        one_ring_neightbors = calc_triangle_angles_sides(v);

        double area = mixed_area(one_ring_neightbors);

        double gauss = gauss_curvature(one_ring_neightbors, area);

        return gauss;


    }

        //for testint delte
    double mean_test(Vertex& v){

        std::list<Triangle> one_ring_neightbors;
        one_ring_neightbors = calc_triangle_angles_sides(v);

        double area = mixed_area(one_ring_neightbors);

        double mean = mean_curvature(one_ring_neightbors, area);

        return mean;
    }


    double mean_curvature(std::list<Triangle> &triangles, double area ){

        Vector_3 mean_curve_normal_op=Vector_3(0.0,0.0,0.0);

        
        Triangle pre = triangles.back();


        for(auto tri: triangles){

            //besprechen was tun 
            //TODO: Border erkennen mithilfe der Facets der halfedge

            /*if(tri.theta >= 3.1){
                std::cout << tri.theta_border << std::endl;
                return edge_curvature;
            }*/

            if(tri.theta_border){
                return edge_curvature;
            }


            mean_curve_normal_op = mean_curve_normal_op +  
            (cot((tri).alpha) + cot(pre.beta)) *tri.prev_v;

            pre = tri;
        }

        mean_curve_normal_op = (1/(2*area)) * mean_curve_normal_op;

        
        return std::sqrt(mean_curve_normal_op.squared_length())/2;

    }


    // for testing
    std::string my_curvature_test_print(Vertex& v){

        std::stringstream  out;

        out << std::endl << "----------------++++++++++++++++++ my_curvatures ++++++++++++++++++----------------" << std::endl;

        //eventuell Vektor
        std::list<Triangle> one_ring_neightbors;


        one_ring_neightbors = calc_triangle_angles_sides(v);

        mesh_t::Vertex::Halfedge_around_vertex_circulator at=v.vertex_begin(), end = at;


        do {
            out << "points around Triangle: " << at->opposite()->vertex()->point() << std::endl;
            at++;
        }while(at != end);

     


        //out << v.point() << std::endl;

        //out << one_ring_neightbors.size() << std::endl;


        for(auto t : one_ring_neightbors){
            //out << t.alpha << std::endl;
            //out << t.beta << std::endl;
            //out << t.theta << std::endl;
            /*out << t.len_prev_v<< std::endl;
            out << t.len_at_prev<< std::endl;
            out << t.len_at_v<< std::endl;

            out << "at_v: " <<t.at_v << std::endl;
            out << "prev_v: " <<t.prev_v << std::endl; */

  
            out << std::endl;


        } 

        double area = mixed_area(one_ring_neightbors);

        double mean = mean_curvature(one_ring_neightbors, area);

        double gauss = gauss_curvature(one_ring_neightbors, area);

        double k_1 = k1(gauss, mean);

        double k_2 = k2(gauss, mean);


         Vector_3 mean_curve_normal_op=Vector_3(0.0,0.0,0.0);

        
        Triangle pre = one_ring_neightbors.back();


        for(auto tri: one_ring_neightbors){

            //besprechen was tun 

           /* if(tri.theta >= 3.1){
                // eventuell Randkonstante einführen

                out <<"error"<< std::endl;
                //return 100.0;
            }*/

            mean_curve_normal_op = mean_curve_normal_op +  
            (cot((tri).alpha) + cot(pre.beta)) *tri.prev_v;

           /* out << "angles" << std::endl;
            out << (tri).alpha << "    "<< pre.beta << "    "<< std::endl;

            out << "at vectors " << std::endl;
            out << tri.len_at_prev<< std::endl;
            out << tri.len_at_v<< std::endl;

            out << "prev vectors " << std::endl;
            out << pre.len_at_prev<< std::endl;
            out << pre.len_at_v<< std::endl;*/

            out << "cot alpha " << cot((tri).alpha) << std::endl;

            out << "cot beta " << cot(pre.beta) << std::endl;

            out << (cot((tri).alpha) + cot(pre.beta)) *tri.prev_v << std::endl;

            pre = tri;
        }

        mean_curve_normal_op = (1/(2*area)) * mean_curve_normal_op;

        
        //return std::sqrt(mean_curve_normal_op.squared_length())/2;

       // out << "area: " << area << std::endl;

        out << "gaus: " << gauss << std::endl;

        out << "mean: " << mean << std::endl;

      /*  out << "k1:  " << k_1 << std::endl;

        out << "k2:  " << k_2 << std::endl;*/

        out << std::endl << "----------------++++++++++++++++++ my_curvatures end ++++++++++++++++++----------------" << std::endl;

        return out.str();
    }


    std::list<Triangle> calc_triangle_angles_sides(Vertex& v){


        //initalisation later
        Vector_3 vec1, vec2, vec3;

        std::list<Triangle> one_ring_neightbors;

        
        mesh_t::Vertex::Halfedge_around_vertex_circulator at=v.vertex_begin(), end = at;
        mesh_t::Vertex::Halfedge_around_vertex_circulator prev = at;

        vec1 = Vector_3(at->opposite()->vertex()->point(), v.point());

        ++at;
        
        do
        {
            Vector_3 vec2 = Vector_3(prev->opposite()->vertex()->point(), at->opposite()->vertex()->point());

            double norm1 = std::sqrt(vec1.squared_length());

            double norm2 = std::sqrt(vec2.squared_length());

            double tmp_beta = acos((vec1*(1/norm1))*(vec2*(1/norm2)));
          
            vec2 = vec2*(-1);

            vec3 = Vector_3(at->opposite()->vertex()->point(), v.point());

            double norm3 = std::sqrt(vec3.squared_length());

            double tmp_alpha = acos((vec2*(1/norm2))*(vec3*(1/norm3)));

            vec3 = vec3*(-1);

            vec1 = vec1*(-1);

            double tmp_theta = acos((vec1*(1/norm1))*(vec3*(1/norm3)));

            //check vec1 and vec3 are part of the border

            bool is_border = false;

            if(at->is_border_edge() && prev->is_border_edge())
                is_border = true;
            



            Triangle newTriangle;
            newTriangle.alpha = tmp_alpha;
            newTriangle.beta = tmp_beta;
            newTriangle.theta = tmp_theta;

            newTriangle.len_prev_v = norm1;
            newTriangle.len_at_prev = norm2;
            newTriangle.len_at_v = norm3;

            newTriangle.at_v = vec3;
            newTriangle.at_prev = vec2;
            newTriangle.prev_v = vec1;

            newTriangle.theta_border = is_border;
            

            one_ring_neightbors.push_back(newTriangle);

            prev = at;

            vec1 = vec3*(-1);
              
            ++at;

        }while(at != end);
     
        vec2 = Vector_3(prev->opposite()->vertex()->point(), at->opposite()->vertex()->point());

        double norm1 = std::sqrt(vec1.squared_length());

        double norm2 = std::sqrt(vec2.squared_length());

        double tmp_beta = acos((vec1*(1/norm1))*(vec2*(1/norm2)));


        vec2 = vec2*(-1);

        vec3 = Vector_3(at->opposite()->vertex()->point(), v.point());

        double norm3 = std::sqrt(vec3.squared_length());

        double tmp_alpha = acos((vec2*(1/norm2))*(vec3*(1/norm3)));

        
        vec3 = vec3*(-1);

        vec1 = vec1*(-1);


        double tmp_theta = acos((vec1*(1/norm1))*(vec3*(1/norm3)));

        bool is_border = false;

        if(at->is_border_edge() && prev->is_border_edge())
            is_border = true;

        Triangle newTriangle;
        newTriangle.alpha = tmp_alpha;
        newTriangle.beta = tmp_beta;
        newTriangle.theta = tmp_theta;
        newTriangle.len_prev_v = norm1;
        newTriangle.len_at_prev = norm2;
        newTriangle.len_at_v = norm3;

        newTriangle.at_v = vec3;
        newTriangle.at_prev = vec2;
        newTriangle.prev_v = vec1;

        newTriangle.theta_border = is_border;
        
        

        one_ring_neightbors.push_back(newTriangle);

        return one_ring_neightbors;


    }

    double k1(double gauss_curve, double mean_curve){

        //paper formular (10)
        if( (mean_curve*mean_curve - gauss_curve)>0 )
            return mean_curve + std::sqrt(mean_curve*mean_curve - gauss_curve);

        return mean_curve;

    }

    double k2(double gauss_curve, double mean_curve){

       //paper formular (11)
        if( (mean_curve*mean_curve - gauss_curve)>0 )
            return mean_curve - std::sqrt(mean_curve*mean_curve - gauss_curve);

        return mean_curve;
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

            if(t.theta >= M_PI_2){ //angle at v obtuse

                area += triangle_area(t)/2;

            }else if(t.alpha >= M_PI_2 || t.beta >= M_PI_2){ //angle at alpha or beta obtuse

                area += triangle_area(t)/4;

            }else{ // nonobtuse trinangle

                //Voronoi formula
                area += 0.125 * (cot(t.alpha)*t.len_prev_v*t.len_prev_v + cot(t.beta)*t.len_at_v*t.len_at_v);


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