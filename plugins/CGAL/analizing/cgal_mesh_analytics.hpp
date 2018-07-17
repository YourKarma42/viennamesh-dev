
#ifndef MESH_ANALYTICS
#define MESH_ANALYTICS

//hash table to store curvatures
#include <unordered_map>

#include "../cgal_curvature_functions.hpp"

#include "../curve_calc_dev/curve_test_1.hpp"


namespace viennamesh
{
    namespace cgal
    {
        template<class ECM_>    
        class cgal_mesh_analytics
        {
        public:

            typedef ECM_ ECM;

            typedef typename ECM::Vertex Vertex;

            typedef typename ECM::Vertex_handle Vertex_handle;

            typedef typename ECM::Halfedge_handle Halfedge_handle;

            typedef typename ECM::Facet_handle Facet_handle;

            typedef typename ECM::Facet::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;

            typedef typename ECM::Facet_iterator Facet_iterator;

            typedef typename ECM::Vertex::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

            typedef typename ECM::Point_3 Point_3;

            //cgal_mesh.hpp auch pushen
            typedef typename viennamesh::cgal::Vector_3 Vector_3;


            //std::vector<Vertex_handle> get_trace(Vertex_handle v, Vertex_handle prev);


          //  cgal_mesh_analytics(double f, double b, double c):flat_edges(f), between_edges(b), curved_edges(c)
            //{}

            cgal_mesh_analytics(ECM & mesh_h) : mesh(mesh_h){


                double flat_boundary = 0.00001;

                
                for (Facet_iterator facet = mesh.facets_begin(); facet != mesh.facets_end(); ++facet)
                {

                    //dosnt work
                    Vector_3 face_normal = CGAL::normal(facet->facet_begin()->vertex()->point(),
                                 facet->facet_begin()->next()->vertex()->point(),
                                 facet->facet_begin()->opposite()->vertex()->point());


                    facet_normals[facet] = std::sqrt(face_normal.squared_length()) * face_normal;
                

                }



                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=mesh.vertices_begin(),end=mesh.vertices_end();at!=end;++at)
                {
                    Vertex_handle t = at;
                    curvatures[at] = calc_curvatures(*at);
                }

               
                std::unordered_map<Vertex_handle, Vertex_handle> trace_startingpoints;

                
 
                int i=0;
                for(cgal::polyhedron_surface_mesh::Halfedge_iterator at = mesh.halfedges_begin(), end = mesh.halfedges_end(); at !=end; ++at){


                    // boundary edge
                    if( curvatures[at->vertex()].mean > flat_boundary && curvatures[at->opposite()->vertex()].mean > flat_boundary){

                        //check if we have an important feature

                        if(curvatures[at->vertex()].p2 < 0.01 && curvatures[at->opposite()->vertex()].p2 < 0.01){



                                                    
                            Facet_handle facet = at->face();
                            Facet_handle facet2 = at->opposite()->face();

                            double angle=0;

                            //check if edge is part of the border

                            if(facet != 0 && facet2 != 0){

                            
                                Vector_3 face_normal = CGAL::normal(facet->facet_begin()->vertex()->point(),
                                    facet->facet_begin()->next()->vertex()->point(),
                                    facet->facet_begin()->opposite()->vertex()->point());


                                Vector_3 face_normal2 = CGAL::normal(facet2->facet_begin()->vertex()->point(),
                                    facet2->facet_begin()->next()->vertex()->point(),
                                    facet2->facet_begin()->opposite()->vertex()->point());

                                double skp = (face_normal*(1/std::sqrt(face_normal.squared_length())))*
                                                        (face_normal2*(1/std::sqrt(face_normal2.squared_length())));
                                
                                //had rounding errors in some Meshes
                                if (skp > 1.0) skp = 1.0;

                                angle = acos(skp);

                                

                                /*if(std::isnan(angle)){
                                    std::cout << "sqrt len: " << std::sqrt(face_normal.squared_length()) << std::endl;
                                    std::cout << "1/sqrt len: " << 1/std::sqrt(face_normal.squared_length()) << std::endl;

                                    std::cout << "sqrt len: " << std::sqrt(face_normal2.squared_length()) << std::endl;
                                    std::cout << "1/sqrt len: " << 1/std::sqrt(face_normal2.squared_length()) << std::endl;

                                    double blub = (face_normal*(1/std::sqrt(face_normal.squared_length())))*
                                                        (face_normal2*(1/std::sqrt(face_normal2.squared_length())));

                                    std::cout << "SKP: " <<  (blub > 1.0) << std::endl;

                                    std::cout << "SKP: " << acos(blub)<< std::endl;

                                    std::cout <<  std::endl;
                                }*/

                                
                            }else{

                                angle = 1.5707;
                                
                            }



                           
                         
                            //std::cout <<  std::endl;

                            //check if edge intersects the mesh
                            if(angle < 0.001){
                                //std::cout << "angle: " << angle << std::endl;

                                if(facet != 0 && facet2 != 0){

                                    Vertex_handle v1, v2;

                                    

                                    Halfedge_around_facet_circulator begin =facet->facet_begin(), f_end = begin; 
                                    
                                    do{
                                        if(begin->vertex() != at->vertex() && begin->vertex() != at->opposite()->vertex()){
                                            v1 = begin->vertex();
                                            break;
                                        }
                                        begin++;
                                    }while(begin != f_end);

                                    begin =facet2->facet_begin(), f_end = begin; 

                                    do{

                                        if(begin->vertex() != at->vertex() && begin->vertex() != at->opposite()->vertex()){
                                            v2 = begin->vertex();
                                            break;
                                        }
                                        begin++;
                                    }while(begin != f_end);

                                    //check if the edge is in the middle of the mesh (not part of a border)
                                    if(curvatures[v1].mean > flat_boundary && curvatures[v2].mean > flat_boundary ){

                                        //nicht color 3 Color 4 (fürs testen 3)

                                        curvatures[at->vertex()].gauss = 35;

                                        curvatures[at->opposite()->vertex()].gauss = 35;

                                        colors[at->vertex()]=3;

                                        colors[at->opposite()->vertex()]=3;

                                    }else{
                                        
                                        //only color if it isnt colored jet
                                        if(colors.find(at->vertex()) == colors.end() && colors.find(at->opposite()->vertex()) == colors.end()){

                                            curvatures[at->vertex()].gauss = 10;

                                            curvatures[at->opposite()->vertex()].gauss = 10;

                                            colors[at->vertex()]=3;

                                            colors[at->opposite()->vertex()]=3;
                                        }

                                    }




                                }else{

                                    //only color if it isnt colored jet
                                    if(colors.find(at->vertex()) == colors.end() && colors.find(at->opposite()->vertex()) == colors.end()){

                                        curvatures[at->vertex()].gauss = 10;

                                        curvatures[at->opposite()->vertex()].gauss = 10;

                                        colors[at->vertex()]=3;

                                        colors[at->opposite()->vertex()]=3;
                                    }
                                }

                                
                            }



                            //check if one Vertex is already colored
                            //if(!((colors.find(at->vertex()) != colors.end()) && (colors.find(at->opposite()->vertex()) != colors.end() ))){
                            //if(colors.find(at->vertex()) == colors.end() && colors.find(at->opposite()->vertex()) == colors.end()){ 
                                
                                //std::cout << << std::endl;        

                            if(colors.find(at->vertex()) == colors.end()){
                                colors[at->vertex()]=1;
                                curvatures[at->vertex()].gauss = 25;
                               // std::cout << at->vertex()->point() << "   " << at->opposite()->vertex()->point() << std::endl; 
                                //std::cout << angle << std::endl;  
                            }

                            if(colors.find(at->opposite()->vertex()) == colors.end()){
                                curvatures[at->opposite()->vertex()].gauss = 25;
                                colors[at->opposite()->vertex()]=1;
                                //std::cout << at->vertex()->point() << "   " << at->opposite()->vertex()->point() << std::endl; 
                                //std::cout << angle << std::endl;
                            }

                            //std::cout << std::endl;
                                
                           // }

                        } else{

                            //check wich vertex of the edge is not flat

                            //eventuell umdenken eventuell zuviele reinholen und dann später raushauen
                            if(curvatures[at->vertex()].p2 > 0.001){

                                colors[at->vertex()]=0;

                                curvatures[at->vertex()].gauss = 1;

                                if (trace_startingpoints.find(at->vertex()) == trace_startingpoints.end())
                                    trace_startingpoints[at->vertex()]=at->vertex();

                                
                            }else{
                                colors[at->opposite()->vertex()]=0;

                                curvatures[at->opposite()->vertex()].gauss = 1;

                                if (trace_startingpoints.find(at->opposite()->vertex()) == trace_startingpoints.end())
                                    trace_startingpoints[at->opposite()->vertex()]=at->opposite()->vertex();
                                }


                                
                            }

                            bool flat_1=true, flat_2 = true;
                            
                            //check the curvatures of each vertex around the the edge
                            Halfedge_around_vertex_circulator v1=at->vertex()->vertex_begin(), v1_end = v1;
                            do{      
                                
                                if(curvatures[v1->opposite()->vertex()].mean < flat_boundary){
                                    flat_1=false;
                                    break;                              
                                }
                                v1++;
                            }while(v1 != v1_end);

                            Halfedge_around_vertex_circulator v2=at->opposite()->vertex()->vertex_begin(), v2_end = v2;
                            do{      
                                
                                if(curvatures[v2->opposite()->vertex()].mean < flat_boundary){
                                    flat_2=false;
                                    break;                              
                                }
                                v2++;
                            }while(v2 != v2_end);

                            if(flat_1){
                                curvatures[at->vertex()].gauss = 50;                           
                            }

                            if(flat_2){
                                curvatures[at->opposite()->vertex()].gauss = 50;
                            }


                       i++;
                   }
                }

                

                /*std::cout << "size: " << trace_startingpoints.size() << std::endl;
*/

                /*for(auto v: trace_startingpoints){

                    std::cout << "Point: " << v.second->point() << std::endl;

                }*/
                

                std::vector<std::vector<Vertex_handle>> traces;


                for(auto start: trace_startingpoints){
                    
                    Vertex_handle v = start.second;
                    
                    Vertex_handle prev = start.second;

                    std::vector<Vertex_handle> vertices;

                    std::vector<Halfedge_handle> he;


                    Halfedge_around_vertex_circulator at=v->vertex_begin(), end = at;
                    do{
                        if(colors.find(at->opposite()->vertex()) != colors.end())                                                
                            he.push_back(at);
                        
                        at++;
                    }while(at != end);

                    std::cout << "start: " << v->point() << std::endl;

                    std::cout << "size: " << he.size() << std::endl;

                    for(auto he_n: he){

                        if(colors[he_n->opposite()->vertex()] == 4){

                            if(get_color(he_n->next()->vertex()) == -1 || get_color(he_n->opposite()->next()->vertex()) == -1 ){

                                //std::cout << "thats the point: " << he_n->opposite()->vertex()->point() << std::endl;
                                //curvatures[he_n->opposite()->vertex()].gauss = 50;

                                //traces.push_back(get_trace_test(he_n->opposite()->vertex(), prev));



                            }



                        }

                        //std::cout << "  point: " << v1->point() << std::endl;

                        traces.push_back(get_trace(he_n->opposite()->vertex(), prev));

                        //traces.push_back(get_trace_test(v1, prev));

                        //std::cout <<  std::endl;

                    }

                                               
                        //std::cout << "colors " << colors[at->opposite()->vertex()] << std::endl;

                }
                    
                colors.clear();

                int next = 1;

                int col_output = 50;

                for(auto t: traces){

                    //std::cout << "size: " << t.size() << std::endl;

                    int i = 1;
                   
                    int color = (next*2)-1;

                    //col_output = 50;

                    

                    for(auto v: t){
                        //curvatures[v].gauss = col_output;
                        //col_output += 10;
                       // std::cout << "point: " << v->point() << std::endl;

                        //std::cout << " color " << color << std::endl;

                        //std::cout << " i " << i << "   color " << color << std::endl;

                        colors[v] = color;

                        //curvatures[v].gauss = col_output;

                        if(i == 10){
                            if(color == ((next*2)-1)){
                                col_output = 100;
                                color++;
                                i=0;

                            }else{
                                col_output = 50;
                                color--;
                                i=0;
                            }
                        }

                        i++;


                    }
                    //std::cout <<  std::endl;

                    next++;

                }



                std::cout << "edges: " << i << std::endl;

            
            }

            std::vector<Vertex_handle> get_trace(Vertex_handle v, Vertex_handle prev){

                std::vector<Vertex_handle> vertices;

                std::vector<Vertex_handle> current_trace;

                std::vector<Halfedge_handle> he_around_v;

                current_trace.push_back(v);

                int i = 0;

                    //find all colored neighbors

                    do{
                        vertices.clear();
                        he_around_v.clear();

                        //std::cout << "  v    " << v->point() << std::endl;
                        //std::cout << "  prev " << prev->point() << std::endl;

                        Halfedge_around_vertex_circulator at=v->vertex_begin(), end = at;
                        do{                            
                            if(colors.find(at->opposite()->vertex()) != colors.end()){
                                //found previous vertex
                                if(prev != at->opposite()->vertex()){

                                    //found vertex that corsses the mesh
                                    if(colors[v] != 3 || colors[at->opposite()->vertex()] != 3){
                                        vertices.push_back(at->opposite()->vertex());
                                        he_around_v.push_back(at);
                                    }                                   
                                }
                            }
                            at++;
                        }while(at != end);

                        //std::cout << vertices.size() << std::endl;


                        if(vertices.size() > 0){

                            if(vertices.size() == 1){
                                if(colors[vertices[0]] != 0){

                                    current_trace.push_back(vertices[0]);

                                    prev = v;

                                    colors.erase(v);

                                    v = vertices[0];

                                }else {
                                    //?
                                    //current_trace.push_back(vertices[0]);
                                    
                                    colors.erase(v);
                                    return current_trace;

                                }
                            }else{

                                /*if(current_trace.empty()){
                                    std::cout << "start: " << v->point() << std::endl;                                   
                                }

                                //std::cout << "at: " << v->point() << std::endl;

                                //funkt irgendwie seltsamerweise...

                                for(auto he: he_around_v){
                                    // eventuell direkt auf colors zugreifen
                                    if(get_color(he->next()->vertex()) == -1 || get_color(he->opposite()->next()->vertex()) == -1 ){

                                        if(colors[he->opposite()->vertex()] != 0){

                                            curvatures[he->opposite()->vertex()].gauss = 50;

                                            std::cout <<he->opposite()->vertex()->point() <<std::endl;

                                            current_trace.push_back(he->opposite()->vertex());

                                            prev = v;

                                            colors.erase(v);

                                            v = he->opposite()->vertex();


                                        }else {
                                                //?
                                                //current_trace.push_back(vertices[0]);
                                                
                                                colors.erase(v);
                                                return current_trace;
                                        }
                                    }

                                }

                                if(i>10){
                                    
                                    std::cout << std::endl;
                                    
                                }

                                i++;*/

                                


                                //case basteln sodass neue 0 eingeführt wird

                                return current_trace;
                                
                            }

                        }
                    }while(vertices.size() > 0);

                    


                return current_trace;


            }


            std::vector<Vertex_handle> get_trace_test(Vertex_handle v, Vertex_handle prev){

                std::vector<Vertex_handle> vertices;

                std::vector<Vertex_handle> current_trace;

                std::vector<Halfedge_handle> he_around_v;

                current_trace.push_back(v);

                int i = 0;

                    //find all colored neighbors

                    do{
                        vertices.clear();
                        he_around_v.clear();

                        //std::cout << "  v    " << v->point() << std::endl;
                        //std::cout << "  prev " << prev->point() << std::endl;

                        Halfedge_around_vertex_circulator at=v->vertex_begin(), end = at;
                        do{                            
                            if(colors.find(at->opposite()->vertex()) != colors.end()){
                                //found previous vertex
                                if(prev != at->opposite()->vertex()){

                                    //found vertex that corsses the mesh
                                    if(colors[v] != 3 || colors[at->opposite()->vertex()] != 3){
                                        vertices.push_back(at->opposite()->vertex());
                                        he_around_v.push_back(at);
                                    }                                   
                                }
                            }
                            at++;
                        }while(at != end);

                        //std::cout << vertices.size() << std::endl;


                        if(vertices.size() > 0){

                            //std::cout << "start: ";

                            if(vertices.size() > 1){

                                std::cout << "start: " << v->point() << std::endl;                                   
                                
                                //std::cout << "at: " << v->point() << std::endl;

                                for(auto he: he_around_v){
                                    std::cout << he->opposite()->vertex()->point() << std::endl;

                                    std::cout << colors[he->opposite()->vertex()] << std::endl;

                                    std::cout << colors[he->opposite()->vertex()] << std::endl;
                                    // eventuell direkt auf colors zugreifen
                                    if(get_color(he->next()->vertex()) == -1 || get_color(he->opposite()->next()->vertex()) == -1 ){

                                        if(colors[he->opposite()->vertex()] != 0){

                                            curvatures[he->opposite()->vertex()].gauss = 50;

                                           /* std::cout << he->opposite()->vertex()->point() <<std::endl;
                                           
                                             std::cout << colors[he->next()->vertex()] << std::endl;

                                             std::cout << colors[he->opposite()->next()->vertex()] << std::endl;*/

                                            current_trace.push_back(he->opposite()->vertex());

                                            prev = v;

                                            colors.erase(v);

                                            v = he->opposite()->vertex();

                                            i++;
                                                                             

                                            break;


                                        }else {
                                                //?
                                                //current_trace.push_back(vertices[0]);
                                                
                                                colors.erase(v);
                                                return current_trace;
                                        }
                                    }

                                }

                                if(i>10){
                                    
                                    std::cout << std::endl;
                                    return current_trace;
                                }

                                

                            }

                        }else{
                            return current_trace;
                        }
                        
                    }while(vertices.size() > 0);

                    


                return current_trace;


            }



            

            void set_new_point(Point_3 p){
                new_point = p;
            }

            bool is_new_point(Point_3 p){

                if(new_point == p)
                    return true;


                return false;

            }

            Vertex_handle search_vertex(){



                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=mesh.vertices_begin(),end=mesh.vertices_end();at!=end;++at)
                {
                    if(at->point() == new_point)
                        return at;
                }

            }




            bool has_curvature(Vertex_handle v){
                return (curvatures.find(v) != curvatures.end());
            }

            double get_mean(Vertex_handle v){
                return curvatures[v].mean;
            }

            double get_gauss(Vertex_handle v){
                return curvatures[v].gauss;
            }

            double get_p1(Vertex_handle v){
                return curvatures[v].p1;
            }

            double get_p2(Vertex_handle v){
                return curvatures[v].p2;
            }

            int get_color(Vertex_handle v){

                if(colors.find(v) != colors.end()){
                    return colors[v];
                }

                return -1;
            }

            int set_color(Vertex_handle v, int color){

                return colors[v] = color;
            }

            //test vars

            int num=0;

            bool collection_over = false;


        private:

            ECM & mesh;

            Point_3 new_point;
            
            double flat_edges;

            double between_edges;

            double curved_edges;

            std::unordered_map<Vertex_handle, Curves> curvatures;

            std::unordered_map<Facet_handle, Vector_3> facet_normals;
            
            std::unordered_map<Vertex_handle, int> colors; 
        };
    }
}

#endif /*MESH_ANALYTICS*/