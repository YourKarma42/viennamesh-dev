
#ifndef MESH_ANALYTICS
#define MESH_ANALYTICS

//hash table to store curvatures
#include <unordered_map>

//#include "../cgal_curvature_functions.hpp"

//#include "../curve_calc_dev/curve_test_1.hpp"

#include "../curvature_calculator.hpp"


namespace viennamesh
{
    namespace cgal
    {
        template<class ECM_>    
        class cgal_mesh_analytics
        {

            typedef ECM_ ECM;

            typedef typename viennamesh::cgal::polyhedron_surface_mesh::Vertex Vertex;

            //typedef typename ECM::Vertex Vertex;
            typedef typename ECM::Vertex_handle Vertex_handle;
            typedef typename ECM::Vertex::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;

            typedef typename ECM::Halfedge_handle Halfedge_handle;

            typedef typename ECM::Facet_handle Facet_handle;
            typedef typename ECM::Facet::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
            typedef typename ECM::Facet_iterator Facet_iterator;

            typedef typename ECM::Point_3 Point_3;

            //cgal_mesh.hpp auch pushen
            typedef typename viennamesh::cgal::Vector_3 Vector_3;

        private:

            ECM & _mesh;

            Point_3 new_point;
            
            int _curved_edges=0;
            int _flat_edges=0;
            int _feature_edges=0;

            //0.2 macht gute ergebnisse
            double _flat_boundary = 0.13;



            double _feature_mean_boundary = 8.0;

            double const _edge_curvature = 100;

            std::unordered_map<Vertex_handle, Curves> _curvatures;

            std::unordered_map<Facet_handle, Vector_3> _facet_normals;
            
            std::unordered_map<Vertex_handle, int> _transition_area; 

            std::unordered_map<Vertex_handle, Vertex_handle> _trace_startingpoints;

            std::unordered_map<Vertex_handle, int> _border;

            std::unordered_map<Vertex_handle, int> _features;

            std::unordered_map<Vertex_handle, int> _curved_area;

            std::unordered_map<Vertex_handle, int> _visited_area;

            std::unordered_map<Vertex_handle, double> _transition_distance;

        
        public:



           

            cgal_mesh_analytics(ECM & mesh_h) : _mesh(mesh_h){

                calculate_curvatures();
                
                new_extrect_features(); 

                test_fill();

                //TODO weg für performance
                //calculate_metrics();

                //calculate_transition_areas();

                //extract_features(); 

                std::cout << "SIZE OF FEATURE: " << _features.size() << std::endl;

           
            }

            bool calculate_transition_areas(){

                _transition_area.clear();

                std::list<Vertex_handle> transition_area;

                std::list<Vertex_handle> transition_area_fill;            

                int transition_zone_size = 100;

                double transition_distance_max = 5;

                for(cgal::polyhedron_surface_mesh::Halfedge_iterator at = _mesh.halfedges_begin(), end = _mesh.halfedges_end(); at !=end; ++at){

                    //find borders between zones

                    if((get_feature(at->vertex()) == -1) ^ (get_feature(at->opposite()->vertex()) == -1)){
                        if(get_feature(at->vertex()) == -1){
                            _transition_area[at->vertex()] = 10;
                            transition_area.push_back(at->vertex());
                            transition_area_fill.push_back(at->vertex());

                            _transition_distance[at->vertex()]=0.0;
                        }else{
                            _transition_area[at->opposite()->vertex()] = 10;
                            transition_area.push_back(at->opposite()->vertex());
                            transition_area_fill.push_back(at->vertex());

                            _transition_distance[at->vertex()]=0.0;
                        }
                    }
                }


                //create the transition area
                while(!transition_area.empty()){


                    Vertex_handle v = transition_area.front();


                    if(_transition_area[v] == (10 + transition_zone_size) || _transition_distance[v] >= transition_distance_max){
                        //std::cout << "blub " << std::endl;
                        transition_area.pop_front();
                        continue;
                    }

                    Halfedge_around_vertex_circulator begin_h=v->vertex_begin(),end_h=begin_h;                    
                    
                    //add new vertices
                    do{
                        //std::cout << "blub1 " << std::endl;
                        if(get_transition_area(begin_h->opposite()->vertex()) == -1 && get_feature(begin_h->opposite()->vertex()) == -1){
                            //std::cout << "blub " << std::endl;
                            transition_area.push_back(begin_h->opposite()->vertex());
                            _transition_area[begin_h->opposite()->vertex()] = _transition_area[begin_h->vertex()]+1;


                            double new_distance = std::sqrt(CGAL::squared_distance(begin_h->opposite()->vertex()->point(), begin_h->vertex()->point()));

                            if(_transition_distance.find(begin_h->opposite()->vertex()) != _transition_distance.end()){

                                if(_transition_distance[begin_h->opposite()->vertex()] < _transition_distance[begin_h->vertex()] + new_distance){
                                    _transition_distance[begin_h->opposite()->vertex()] = 
                                        _transition_distance[begin_h->vertex()] + new_distance;
                                }
                            } else {
                                _transition_distance[begin_h->opposite()->vertex()] = 
                                    _transition_distance[begin_h->vertex()] + new_distance;
                            }
                        }
            
                        begin_h++;
                    }while(begin_h != end_h);

                    transition_area.pop_front();
                }

                //delete small paths in transition area
                while(!transition_area_fill.empty()){

                    Vertex_handle v = transition_area_fill.front();

                    Halfedge_around_vertex_circulator begin_h=v->vertex_begin(),end_h=begin_h;

                    bool recolor = true;                    
                    
                    //add new vertices
                    do{

                        if(get_transition_area(begin_h->opposite()->vertex()) == -1 && get_feature(begin_h->opposite()->vertex()) == -1){
                            recolor = false;
                            
                            break;
                        }

                        if(get_transition_area(begin_h->opposite()->vertex()) != -1 && get_transition_area(begin_h->opposite()->vertex()) > 10) {
                            recolor = false;
                             
                            break; 
                        }
            
                        begin_h++;
                    }while(begin_h != end_h);

                    if(recolor == true){
                        _transition_area.erase(v);
                        _transition_distance.erase(v);
                        _features[v] = 10;
                    }

                    transition_area_fill.pop_front();


                }

                /*for(cgal::polyhedron_surface_mesh::Vertex_iterator at=_mesh.vertices_begin(),end=_mesh.vertices_end();at!=end;++at)
                {
                    
                    if(get_feature(at) == -1 && get_transition_area(at) == -1){

                        //fill_area_transition(80, at);
                        
                    }
                }*/
                

            }

            bool calculate_metrics(){

                int total_number=0;

                double min_length=10000000.0;
                double max_length=0.0;

                _curved_edges=0;
                _flat_edges=0;
                _feature_edges=0;


                for(cgal::polyhedron_surface_mesh::Halfedge_iterator at = _mesh.halfedges_begin(), end = _mesh.halfedges_end(); at !=end; ++at){
                    
                    double length = CGAL::squared_distance(at->vertex()->point(), at->opposite()->vertex()->point());
                    if(length < min_length){
                        min_length = length;
                    }

                    if(length > max_length){
                        max_length = length;
                    }

                    //calculate number of edges in each region
                    if(get_feature(at->vertex()) ==-1 && get_feature(at->opposite()->vertex()) == -1){
                        if(get_curved_area(at->vertex()) == -1 && get_curved_area(at->opposite()->vertex()) == -1){
                            _flat_edges++;
                        }else{
                            _curved_edges++;
                        }
                    }else{
                        _feature_edges++;

                    }
                    total_number++;

                }

                _curved_edges = _curved_edges/2;
                _flat_edges = _flat_edges/2;
                _feature_edges = _feature_edges/2;

                std::cout << "flat: "<< _flat_edges << std::endl;

                std::cout << "curved: "<< _curved_edges << std::endl;

                std::cout << "feature: "<< _feature_edges << std::endl;

                std::cout << "total number of edges: "<< total_number/2 << std::endl;

                std::cout << "MAX length of edge: "<< std::sqrt(max_length) << std::endl;

                std::cout << "MIN length of edge: "<< std::sqrt(min_length) << std::endl;

            }

            double average_length_curved_edges(){

                int curved_edges_count=0;

                double sum_length = 0.0;


                for(cgal::polyhedron_surface_mesh::Halfedge_iterator at = _mesh.halfedges_begin(), end = _mesh.halfedges_end(); at !=end; ++at){

                    //calculate number of edges in each region
                    if(get_feature(at->vertex()) ==-1 && get_feature(at->opposite()->vertex()) == -1){
  
                    }else{

                        sum_length = sum_length + (CGAL::squared_distance(at->vertex()->point(), at->opposite()->vertex()->point())/2);
                        curved_edges_count++;

                    }


                }

                std::cout << "anz "<< ((double)curved_edges_count/2) << std::endl;
                return (sum_length/((double)curved_edges_count/2));

                
                //return (sum_length);

            }

            void extend_transition_area(){

                int extension_val = 30;


                for(cgal::polyhedron_surface_mesh::Halfedge_iterator at = _mesh.halfedges_begin(), end = _mesh.halfedges_end(); at !=end; ++at){

                    if(get_transition_area(at->vertex()) !=-1 && get_transition_area(at->opposite()->vertex()) != -1){
                        continue;
                    }
                    
                    
                    if((get_transition_area(at->vertex()) !=-1 || get_transition_area(at->opposite()->vertex()) != -1) &&
                     (get_feature(at->vertex()) ==-1 || get_feature(at->opposite()->vertex()) == -1)){

                        //das if eventuell weg noch überlegen was ich damit mach!
                        if(get_feature(at->vertex()) ==-1 && get_feature(at->opposite()->vertex()) == -1){

                            if((get_transition_area(at->vertex()) != extension_val && get_transition_area(at->opposite()->vertex()) != extension_val)){

                                if(get_transition_area(at->vertex()) ==-1){
                                    set_transition_area(at->vertex(), extension_val);
                                }else{
                                    set_transition_area(at->opposite()->vertex(), extension_val);
                                }

                            }
                        }


                    }


                }


            }

            bool fill_area_transition(int max_size, Vertex_handle v){
                
                //std::unordered_map<Vertex_handle, int> tmp_features = _features;

                std::unordered_map<Vertex_handle, int> tmp_transition;

                int node_counter=1;

                tmp_transition[v] = 10;

                std::list<Vertex_handle> que;

                std::list<Vertex_handle> to_color;

                to_color.push_back(v);

                while(node_counter <= max_size){

                    tmp_transition[v] = 10;

                    Halfedge_around_vertex_circulator begin_h=v->vertex_begin(),end_h=begin_h;                    

                    do{
                        if(tmp_transition.find(begin_h->opposite()->vertex()) == tmp_transition.end()
                           && _transition_area.find(begin_h->opposite()->vertex()) == _transition_area.end()){


                            tmp_transition[begin_h->opposite()->vertex()] = 5;

                            que.push_back(begin_h->opposite()->vertex());

                            to_color.push_back(begin_h->opposite()->vertex());


                            node_counter++;

                        }
                                        
                        begin_h++;
                    }while(begin_h != end_h);


                    //if no vertices are in the que we have found an area that is small enough to fill
                    if(que.empty() ){
                        for(auto vert: to_color){      


                            _transition_area[vert] = 10;
                            //machen
                            _transition_distance[vert]=1;
                        }
                        return true;
                    }
                    
                    //pop next from que
                    v = que.front();
                    que.pop_front();

                    

                }        

                return false;

            }

            bool fill_area(int max_size, Vertex_handle v){
                
                //std::unordered_map<Vertex_handle, int> tmp_features = _features;

                std::unordered_map<Vertex_handle, int> tmp_features;

                int node_counter=1;

                tmp_features[v] = 10;

                std::list<Vertex_handle> que;

                std::list<Vertex_handle> to_color;

                to_color.push_back(v);

                while(node_counter <= max_size){

                    tmp_features[v] = 10;

                    Halfedge_around_vertex_circulator begin_h=v->vertex_begin(),end_h=begin_h;                    

                    do{
                        if(tmp_features.find(begin_h->opposite()->vertex()) == tmp_features.end()
                           && _features.find(begin_h->opposite()->vertex()) == _features.end()){


                            tmp_features[begin_h->opposite()->vertex()] = 5;

                            que.push_back(begin_h->opposite()->vertex());

                            to_color.push_back(begin_h->opposite()->vertex());


                            node_counter++;

                        }
                                        
                        begin_h++;
                    }while(begin_h != end_h);


                    //if no vertices are in the que we have found an area that is small enough to fill
                    if(que.empty() ){
                        for(auto vert: to_color)
                            _features[vert] = 10;
                        return true;
                    }
                    
                    //pop next from que
                    v = que.front();
                    que.pop_front();

                    

                }        

                return false;

            }

            void test_fill(){               

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=_mesh.vertices_begin(),end=_mesh.vertices_end();at!=end;++at)
                {
                    
                    if(_features.find(at) == _features.end()){

                        fill_area(5, at);
                        
                    }
                }
            }

            void new_extrect_features(){

                _features.clear();


                std::list<Vertex_handle> border_vertices;


                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=_mesh.vertices_begin(),end=_mesh.vertices_end();at!=end;++at)
                {

                    //check if vertex is part of the border
                    //maby extract the border of the mesh
                    if(_curvatures[at].gauss == _edge_curvature && _curvatures[at].mean == _edge_curvature){

                        border_vertices.push_back(at);
                        continue;
                    }
                    


                    //vertex inside the mesh
                    if(fabs(_curvatures[at].p1) > _flat_boundary || fabs(_curvatures[at].p2) > _flat_boundary){

                        //vertex is part of a curved area

                        double diff_in_curv =  1;

                        if(fabs(_curvatures[at].gauss) >= 10.0 || fabs(_curvatures[at].mean) >=10.0){

                            //vertex is part of a feature

                            _curved_area[at] = 10;
                            _features[at] = 10;



                        }else{

                            //vertex is part of a curved area

                            _features[at] = 10; 

                        }

                    }

                }
                
                //add feature designation to boarder
                for(auto vert: border_vertices){


                    
                    Halfedge_around_vertex_circulator begin_h=vert->vertex_begin(),end_h=begin_h;

                    int num_features=0;
                    //ignore first 
                    begin_h++;

                    do{
                        if(get_feature(begin_h->opposite()->vertex())==10){
                            if(!(_curvatures[begin_h->opposite()->vertex()].gauss == _edge_curvature && 
                               _curvatures[begin_h->opposite()->vertex()].mean == _edge_curvature)){                          
                                                    
                                _features[vert] = 10;
                                goto differing_features_break;
                                
                            }
                        }

                        begin_h++;

                    }while(begin_h != end_h);

                differing_features_break:;
                }

            }

            void clean_curved_area(){

                //dont iterate over all vertices save veature vertices in a map

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=_mesh.vertices_begin(),end=_mesh.vertices_end();at!=end;++at)
                {
                    if(get_feature != -1){
                        Halfedge_around_vertex_circulator begin_h=at->vertex_begin(),end_h=begin_h;
                        do{


                            begin_h++;
                        }while(begin_h != end_h);
                    }


                }

                //umgebung aller features betrachten abhängig vom Grad des vertex entscheiden ob es sich wirklich um ein feature handelt

            }

            void clean_features(){

                //überlegen wie man features verbinden kann

            }


            void extract_features(){
                int i=0;
                for(cgal::polyhedron_surface_mesh::Halfedge_iterator at = _mesh.halfedges_begin(), end = _mesh.halfedges_end(); at !=end; ++at){


                    // boundary edge
                    if( _curvatures[at->vertex()].mean > _flat_boundary && _curvatures[at->opposite()->vertex()].mean > _flat_boundary){

                        //check if we have an important feature

                        if(_curvatures[at->vertex()].p2 < 0.01 && _curvatures[at->opposite()->vertex()].p2 < 0.01){



                                                    
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
                                
                            }else{

                                angle = 1.5707;
                                
                            }


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
                                    if(_curvatures[v1].mean > _flat_boundary && _curvatures[v2].mean > _flat_boundary ){

                                        //nicht color 3 Color 4 (fürs testen 3)

                                        _features[at->vertex()] = 35;

                                        _features[at->opposite()->vertex()] = 35;

                                        _transition_area[at->vertex()]=3;

                                        _transition_area[at->opposite()->vertex()]=3;

                                        //_transition_area[at->vertex()]=0;

                                        //_transition_area[at->opposite()->vertex()]=0;

                                    }else{
                                        
                                        //only color if it isnt colored jet
                                        if(_transition_area.find(at->vertex()) == _transition_area.end() && _transition_area.find(at->opposite()->vertex()) == _transition_area.end()){

                                            _features[at->vertex()] = 10;

                                            _features[at->opposite()->vertex()] = 10;

                                            _transition_area[at->vertex()]=3;

                                            _transition_area[at->opposite()->vertex()]=3;

                                            //_transition_area[at->vertex()]=0;

                                            //_transition_area[at->opposite()->vertex()]=0;
                                        }

                                    }




                                }else{

                                    //only color if it isnt colored jet
                                    if(_transition_area.find(at->vertex()) == _transition_area.end() && _transition_area.find(at->opposite()->vertex()) == _transition_area.end()){

                                        _features[at->vertex()] = 10;

                                        _features[at->opposite()->vertex()] = 10;

                                        _transition_area[at->vertex()]=3;

                                        _transition_area[at->opposite()->vertex()]=3;
                                    }
                                }

                                
                            }



                            //check if one Vertex is already colored
                            //if(!((colors.find(at->vertex()) != colors.end()) && (colors.find(at->opposite()->vertex()) != colors.end() ))){
                            //if(colors.find(at->vertex()) == colors.end() && colors.find(at->opposite()->vertex()) == colors.end()){ 
                                
                                //std::cout << << std::endl;        

                            if(_transition_area.find(at->vertex()) == _transition_area.end()){
                                _transition_area[at->vertex()]=1;
                                _features[at->vertex()] = 25;
                               // std::cout << at->vertex()->point() << "   " << at->opposite()->vertex()->point() << std::endl; 
                                //std::cout << angle << std::endl;  
                            }

                            if(_transition_area.find(at->opposite()->vertex()) == _transition_area.end()){
                                _features[at->opposite()->vertex()] = 25;
                                _transition_area[at->opposite()->vertex()]=1;
                                //std::cout << at->vertex()->point() << "   " << at->opposite()->vertex()->point() << std::endl; 
                                //std::cout << angle << std::endl;
                            }

                            //std::cout << std::endl;
                                
                           // }

                        } else{

                            //check wich vertex of the edge is not flat

                            //eventuell umdenken eventuell zuviele reinholen und dann später raushauen
                            if(_curvatures[at->vertex()].p2 > 0.001){

                                _transition_area[at->vertex()]=0;

                                _features[at->vertex()] = 1;

                                if (_trace_startingpoints.find(at->vertex()) == _trace_startingpoints.end())
                                    _trace_startingpoints[at->vertex()]=at->vertex();

                                
                            }else{
                                _transition_area[at->opposite()->vertex()]=0;

                                _features[at->opposite()->vertex()] = 1;

                                if (_trace_startingpoints.find(at->opposite()->vertex()) == _trace_startingpoints.end())
                                    _trace_startingpoints[at->opposite()->vertex()]=at->opposite()->vertex();
                                }


                                
                            }

                            bool flat_1=true, flat_2 = true;
                            
                            //check the curvatures of each vertex around the the edge
                            Halfedge_around_vertex_circulator v1=at->vertex()->vertex_begin(), v1_end = v1;
                            do{      
                                
                                if(_curvatures[v1->opposite()->vertex()].mean < _flat_boundary){
                                    flat_1=false;
                                    break;                              
                                }
                                v1++;
                            }while(v1 != v1_end);

                            Halfedge_around_vertex_circulator v2=at->opposite()->vertex()->vertex_begin(), v2_end = v2;
                            do{      
                                
                                if(_curvatures[v2->opposite()->vertex()].mean < _flat_boundary){
                                    flat_2=false;
                                    break;                              
                                }
                                v2++;
                            }while(v2 != v2_end);

                            if(flat_1){
                                _curved_area[at->vertex()] = 50;                           
                            }

                            if(flat_2){
                                _curved_area[at->opposite()->vertex()] = 50;
                            }


                       i++;
                   }
                }
            }
            
            void create_traces(){
                std::vector<std::vector<Vertex_handle>> traces;


                for(auto start: _trace_startingpoints){
                    
                    Vertex_handle v = start.second;
                    
                    Vertex_handle prev = start.second;

                    std::vector<Vertex_handle> vertices;

                    std::vector<Halfedge_handle> he;


                    Halfedge_around_vertex_circulator at=v->vertex_begin(), end = at;
                    do{
                        if(_transition_area.find(at->opposite()->vertex()) != _transition_area.end())                                                
                            he.push_back(at);
                        
                        at++;
                    }while(at != end);

                    //std::cout << "start: " << v->point() << std::endl;

                    //std::cout << "size: " << he.size() << std::endl;

                    for(auto he_n: he){


                        //std::cout << "  point: " << v1->point() << std::endl;

                        traces.push_back(get_trace(he_n->opposite()->vertex(), prev));

                        //traces.push_back(get_trace_test(v1, prev));

                        //std::cout <<  std::endl;

                    }

                                               
                        //std::cout << "colors " << colors[at->opposite()->vertex()] << std::endl;

                }
                    
                _transition_area.clear();

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

                        _transition_area[v] = color;

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
                            if(_transition_area.find(at->opposite()->vertex()) != _transition_area.end()){
                                //found previous vertex
                                if(prev != at->opposite()->vertex()){

                                    //found vertex that corsses the mesh
                                    if(_transition_area[v] != 3 || _transition_area[at->opposite()->vertex()] != 3){
                                        vertices.push_back(at->opposite()->vertex());
                                        he_around_v.push_back(at);
                                    }                                   
                                }
                            }
                            at++;
                        }while(at != end);

                        if(fabs((*v).point().x() + 1.0) <= 0.00001 &&
                            fabs((*v).point().y() - 0.5125 ) <= 0.00001 && /*1.275*/
                            fabs((*v).point().z() + 1.3) <= 0.00001) {   

                            std::cout << std::endl << std::endl;

                                
                            std::cout << "BLUB" << std::endl;

                            //std::cout << he_around_v[0]->opposite()->vertex()->point() << std::endl;
                            //std::cout << he_around_v[1]->opposite()->vertex()->point() << std::endl;

                           // std::cout << _transition_area[he_around_v[0]->opposite()->vertex()] << std::endl;
                            //std::cout << _transition_area[he_around_v[1]->opposite()->vertex()] << std::endl;

                            /*for(auto vert: current_trace){
                                std::cout << vert->point() << std::endl;
                            }*/
                            

                            std::cout<< std::endl << std::endl;
                  
                        }   

                        //std::cout << vertices.size() << std::endl;


                        if(vertices.size() > 0){

                            if(vertices.size() == 1){
                                if(_transition_area[vertices[0]] != 0){

                                    current_trace.push_back(vertices[0]);

                                    prev = v;

                                    _transition_area.erase(v);

                                    v = vertices[0];

                                }else {
                                    //?
                                    //current_trace.push_back(vertices[0]);
                                    
                                    _transition_area.erase(v);
                                    return current_trace;

                                }
                            }else{

                                //temporary fix for exp_conf

                                if(vertices.size() == 2){
                                    if(_transition_area[prev]==0){
                                        if(_transition_area[he_around_v[0]->opposite()->vertex()] != 0){

                                            current_trace.push_back(vertices[0]);

                                            prev = v;

                                            _transition_area.erase(v);

                                            v = vertices[0];

                                            continue;

                                        }
                                        if(_transition_area[he_around_v[1]->opposite()->vertex()] != 0){
                                            current_trace.push_back(vertices[1]);

                                            prev = v;

                                            _transition_area.erase(v);

                                            v = vertices[1];

                                            continue;
                                            
                                        }
                                    }

                                }

                               

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
                            if(_transition_area.find(at->opposite()->vertex()) != _transition_area.end()){
                                //found previous vertex
                                if(prev != at->opposite()->vertex()){

                                    //found vertex that corsses the mesh
                                    if(_transition_area[v] != 3 || _transition_area[at->opposite()->vertex()] != 3){
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

                                    std::cout << _transition_area[he->opposite()->vertex()] << std::endl;

                                    std::cout << _transition_area[he->opposite()->vertex()] << std::endl;
                                    // eventuell direkt auf colors zugreifen
                                    if(get_transition_area(he->next()->vertex()) == -1 || get_transition_area(he->opposite()->next()->vertex()) == -1 ){

                                        if(_transition_area[he->opposite()->vertex()] != 0){

                                            

                                           /* std::cout << he->opposite()->vertex()->point() <<std::endl;
                                           
                                             std::cout << colors[he->next()->vertex()] << std::endl;

                                             std::cout << colors[he->opposite()->next()->vertex()] << std::endl;*/

                                            current_trace.push_back(he->opposite()->vertex());

                                            prev = v;

                                            _transition_area.erase(v);

                                            v = he->opposite()->vertex();

                                            i++;
                                                                             

                                            break;


                                        }else {
                                                //?
                                                //current_trace.push_back(vertices[0]);
                                                
                                                _transition_area.erase(v);
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


            //normalized facet normals
            void calculate_facet_normals(){

                for (Facet_iterator facet = _mesh.facets_begin(); facet != _mesh.facets_end(); ++facet)
                {

                    //dosnt work
                    Vector_3 face_normal = CGAL::normal(facet->facet_begin()->vertex()->point(),
                                 facet->facet_begin()->next()->vertex()->point(),
                                 facet->facet_begin()->opposite()->vertex()->point());


                    _facet_normals[facet] = std::sqrt(face_normal.squared_length()) * face_normal;
                

                }

            }

            void calculate_curvatures(){

                _curvatures.clear();

                std::list<int> sv;

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=_mesh.vertices_begin(),end=_mesh.vertices_end();at!=end;++at)
                {
                    /*double blub[2];
                    principal_curvatures_cgal(*at, _mesh, blub);*/

                    Vertex_handle t = at;
                    //in curvature calculator
                    _curvatures[at] = calc_curvatures(*at);



                    Halfedge_around_vertex_circulator at1=at->vertex_begin(), end1 = at1;
                    int anz=0;
                        do{     

                            anz++;                       
    
                            at1++;
                        }while(at1 != end1);

                    sv.push_back(anz);



                }

            char cwd[PATH_MAX];
            std::string output_filename_my = getcwd(cwd, sizeof(cwd));
            output_filename_my += "/oneringneig.csv";

            std::ofstream csv_my;
            csv_my.open(output_filename_my.c_str(),  std::ios::app);

            for(auto a: sv){
                csv_my << a << ", "; 
            }
            


            //close csv file
            csv_my.close();






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

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=_mesh.vertices_begin(),end=_mesh.vertices_end();at!=end;++at)
                {
                    if(at->point() == new_point)
                        return at;
                }

            }




            int get_num_curved_edges(){
                return _curved_edges;
            }

            int get_num_feature_edges(){
                return _feature_edges;
            }

            int get_num_flat_edges(){
                return _flat_edges;
            }

            bool reduce_flat_edges(int rem){
                _flat_edges = _flat_edges - rem;
                return true;
            }



            bool has_curvature(Vertex_handle v){
                return (_curvatures.find(v) != _curvatures.end());
            }

            double get_mean(Vertex_handle v){
                return _curvatures[v].mean;
            }

            double get_gauss(Vertex_handle v){
                return _curvatures[v].gauss;
            }

            double get_p1(Vertex_handle v){
                return _curvatures[v].p1;
            }

            double get_p2(Vertex_handle v){
                return _curvatures[v].p2;
            }

            int get_transition_area(Vertex_handle v){

                if(_transition_area.find(v) != _transition_area.end()){
                    return _transition_area[v];
                }

                return -1;
            }


            double get_transition_distance(Vertex_handle v){

                if(_transition_distance.find(v) != _transition_distance.end()){
                    return _transition_distance[v];
                }

                return -1;
            }



            int get_feature(Vertex_handle v){
                if(_features.find(v) != _features.end()){
                    return _features[v];
                }

                return -1;                
            }

            int get_curved_area(Vertex_handle v){
                if(_curved_area.find(v) != _curved_area.end()){
                    return _curved_area[v];
                }

                return -1;
            }



            int set_transition_area(Vertex_handle v, int color){

                return _transition_area[v] = color;
            }

            //test vars

            int num=0;

            bool collection_over = false;

        };
    }
}

#endif /*MESH_ANALYTICS*/