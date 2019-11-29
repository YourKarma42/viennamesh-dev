
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

            viennagrid_numeric const _flat_boundary = 2.0;

            double _feature_mean_boundary = 8.0;

            double const _edge_curvature = 100;

            std::unordered_map<Vertex_handle, Curves> _curvatures;

            std::unordered_map<Facet_handle, Vector_3> _facet_normals;
            
            std::unordered_map<Vertex_handle, int> _transition_area; 
            
            //remove in time
            //std::unordered_map<Vertex_handle, Vertex_handle> _trace_startingpoints;
            
            //remove in time
            //std::unordered_map<Vertex_handle, int> _border;

            std::unordered_map<Vertex_handle, int> _features;

            std::unordered_map<Vertex_handle, int> _curved_area;

            //remove in time
            //std::unordered_map<Vertex_handle, int> _visited_area;

            std::unordered_map<Vertex_handle, double> _transition_distance;



        
        public:


            //TODO: make private!
            double sum_features_edges= 0.0;
            double anz_features_edges=0.0;
            double max = 0.0;

            cgal_mesh_analytics(ECM & mesh_h, double flat_boundary) : _mesh(mesh_h),_flat_boundary(flat_boundary){

                //resereve space for containers

                long num_verts = _mesh.size_of_vertices(); 

                _curvatures.reserve(num_verts);
                _transition_area.reserve(num_verts);
                _features.reserve(num_verts);
                _curved_area.reserve(num_verts);
                _transition_distance.reserve(num_verts);
             

                calculate_curvatures();
                
                new_extrect_features(); 

                calculate_transition_areas();

                


                double new_distance = 0.0;
                for(cgal::polyhedron_surface_mesh::Edge_iterator at = _mesh.edges_begin(), end = _mesh.edges_end(); at !=end; ++at){

                    if((get_feature(at->vertex()) != -1) && (get_feature(at->opposite()->vertex()) != -1)){
                        new_distance = std::sqrt(CGAL::squared_distance(at->opposite()->vertex()->point(), at->vertex()->point()));
                        anz_features_edges++;
                        sum_features_edges += new_distance;
                        if(new_distance > max)
                            max = new_distance;
                        
                    }
                }

                std::cout << "max length: " << max << std::endl;
                std::cout << "avg length: " << (sum_features_edges/anz_features_edges) << std::endl;
           
            }

            void recalculate(){

                calculate_curvatures();
                
                new_extrect_features(); 

                calculate_transition_areas();

            }

            std::list<Vertex_handle> _global_transition_area;


            bool recalculate_transition_areas(){

                std::list<Vertex_handle> transition_area;// = _global_transition_area;

                for(cgal::polyhedron_surface_mesh::Halfedge_iterator at = _mesh.halfedges_begin(), end = _mesh.halfedges_end(); at !=end; ++at){

                    if((get_feature(at->vertex()) == -1) ^ (get_feature(at->opposite()->vertex()) == -1)){
                        if(get_feature(at->vertex()) == -1){
                            _transition_area[at->vertex()] = 10;
                            transition_area.push_back(at->vertex());
                            

                            _transition_distance[at->vertex()]=0.0;
                        }else{
                            _transition_area[at->opposite()->vertex()] = 10;
                            transition_area.push_back(at->opposite()->vertex());
                            

                            _transition_distance[at->vertex()]=0.0;
                        }
                    }
                }


 
                while(!transition_area.empty()){

                    Vertex_handle v = transition_area.front();

                    Halfedge_around_vertex_circulator begin_h=v->vertex_begin(),end_h=begin_h;                    
                    
                    //add new vertices
                    do{
                        if(get_transition_area(begin_h->opposite()->vertex()) == -1 && get_feature(begin_h->opposite()->vertex()) == -1){

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
               
                return true;
            }


            


            bool calculate_transition_areas(){

                //std::cout << _transition_area.capacity() <<std::endl;

                _transition_area.clear();

                //std::cout << _transition_area.capacity() <<std::endl;

                //check list v vector speed

                std::list<Vertex_handle> transition_area;

                //transition_area.

                std::list<Vertex_handle> transition_area_fill;            

                int transition_zone_size = 100;

                double transition_distance_max = 5;

                int i=0;

                




                //for(cgal::polyhedron_surface_mesh::Halfedge_iterator at = _mesh.halfedges_begin(), end = _mesh.halfedges_end(); at !=end; ++at){
                for(cgal::polyhedron_surface_mesh::Edge_iterator at = _mesh.edges_begin(), end = _mesh.edges_end(); at !=end; ++at){


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

                _global_transition_area = transition_area;


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

                //_curvatures.clear();

                //TODO: reserve space for the vector

                std::list<int> sv;

  
                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=_mesh.vertices_begin(),end=_mesh.vertices_end();at!=end;++at)
                {

                    Vertex_handle t = at;
                    //in curvature calculator
                    _curvatures[at] = calc_curvatures(*at);
                    //Curves test = calc_curvatures(*at);



                    /*Halfedge_around_vertex_circulator at1=at->vertex_begin(), end1 = at1;
                    int anz=0;
                        do{     

                            anz++;                       
    
                            at1++;
                        }while(at1 != end1);

                    sv.push_back(anz);*/



                }

                /*char cwd[PATH_MAX];
                std::string output_filename_my = getcwd(cwd, sizeof(cwd));
                output_filename_my += "/oneringneig.csv";

                std::ofstream csv_my;
                csv_my.open(output_filename_my.c_str(),  std::ios::app);

                for(auto a: sv){
                    csv_my << a << ", "; 
                }
            

                //close csv file
                csv_my.close();*/

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

            Vector_3 get_normal_vec(Vertex_handle v){
                return _curvatures[v].normal;
            }

            double get_area(Vertex_handle v){
                return _curvatures[v].area;
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