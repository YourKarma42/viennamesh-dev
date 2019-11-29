/* ============================================================================
   Copyright (c) 2011-2016, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                ViennaMesh - The Vienna Meshing Framework
                            -----------------

                    http://viennamesh.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */


#include "cgal_mesh.hpp"
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include<chrono>

#include <limits>


namespace viennamesh
{
    //ViennaGrid typedefs
    typedef viennagrid::mesh                                                        MeshType;

    typedef viennagrid::result_of::element<MeshType>::type                          VertexType;
    typedef viennagrid::result_of::element<MeshType>::type                          FacetType;
    typedef viennagrid::result_of::element<MeshType>::type                          TriangleType;
    typedef viennagrid::result_of::element<MeshType>::type ElementType;

    typedef viennagrid::result_of::const_cell_range<MeshType>::type                 ConstTriangleRange;
    typedef viennagrid::result_of::iterator<ConstTriangleRange>::type               ConstTriangleIterator;
    typedef viennagrid::result_of::const_neighbor_range<MeshType, 1, 2>::type       ConstNeighborRange;
    typedef viennagrid::result_of::iterator<ConstNeighborRange>::type               ConstNeighborIterator;

    typedef viennagrid::result_of::const_vertex_range<MeshType>::type               ConstVertexRangeType;
    typedef viennagrid::result_of::iterator<ConstVertexRangeType>::type             ConstVertexIteratorType;

    typedef viennagrid::result_of::const_element_range<MeshType, 2>::type           ConstCellRangeType;

    typedef viennagrid::result_of::const_facet_range<TriangleType>::type            ConstEdgeOnTriangleRange;
    typedef viennagrid::result_of::iterator<ConstEdgeOnTriangleRange>::type         ConstEdgeOnTriangleIterator;

    typedef viennagrid::result_of::const_vertex_range<TriangleType>::type           ConstVertexOnTriangleRange;
    typedef viennagrid::result_of::iterator<ConstVertexOnTriangleRange>::type       ConstVertexOnTriangleIterator;

       

    /*
     * Class to provide a functor for building a CGAL triangular mesh. GAL uses the half edge data structure.
    */
    template <class HDS>
    class Build_triangle_mesh : public CGAL::Modifier_base<HDS>
    {
    public:

        /*
         * points ... Coordinates of vertices (as they are stored in ViennaGrid data structure) are given in the format:
         *      point1_x, point1_y, point1_z, point2_x, point2_y, point2_z, point3_x, ...
         *
         * facets ... Contains the vertex id of the vertices which form a triangle, format:
         *      triangle1_vertex1_id, triangle1_vertex2_id, triangle1_vertex3_id, triangle2_vertex1_id, ...
         * 
         * scale ... Optional Parameter that scales the mesh by the supplied parameter
        */
        Build_triangle_mesh(const std::vector<double> points, const std::vector<int> facets) : points(points), facets(facets) {}
        Build_triangle_mesh(const std::vector<double> points, const std::vector<int> facets, double scale) : points(points), facets(facets), scale(scale) {}
        void operator()( HDS& hds)
        {
            typedef cgal::Point_3  Point;

            // Postcondition: hds is a valid polyhedral surface.
            CGAL::Polyhedron_incremental_builder_3<HDS> Builder( hds, true);

            //The triangle mesh formes one surface/segment, which is initialized here
            Builder.begin_surface( points.size() / 3, facets.size() / 3);

            //________________________________________Development area_______________________________________________________________________________________________________________________


            for(std::size_t i = 0; i < points.size()/3; ++i)
            {
                Builder.add_vertex( Point(points[3*i]*scale, points[3*i+1]*scale, points[3*i+2]*scale));

            }


            //________________________________________Development area_______________________________________________________________________________________________________________________


            //adding vertex coordinates (from ViennaGrid) to CGAL vertex list
            /*for(std::size_t i = 0; i < points.size()/3; ++i)
            {
                Builder.add_vertex( Point(points[3*i], points[3*i+1], points[3*i+2]));

            }*/

            
            std::vector<std::tuple<int, int, int>> failed_triangles;


            /*number of triangle, which failed to be added because none of the both possible orientation ensures a valid half
             * edge data structure (This is mostly due to topological or numeric problems.*/
            int failed_triangle_cnt = 0;

            /* Triangles are added to the surface. If a triangle's orientation does not fit to the already present half edge
             * data structre it is flipped by choosing an anticyclic permutation of the vertices.*/
            for(std::vector<int>::const_iterator it = facets.begin(); it < facets.end(); it = it+3)
            {

                if(Builder.test_facet(it, it+3) == true) //CGAL function that tests the current triangle orientation
                {
                    Builder.add_facet (it,it+3);
                }
                else //choose anticyclic vertex permutation
                {
                    int tmp_facet[3];

                    tmp_facet[0] = *it;
                    tmp_facet[1] = *(it+2);
                    tmp_facet[2] = *(it+1);

                    if(Builder.test_facet(tmp_facet, tmp_facet + 3) == false) //test again
                    {
                        error(1) <<  "Triangle insertion failed!" << std::endl;

                        failed_triangles.push_back(std::make_tuple(*it, *(it+1), *(it+2)));


                        error(1) <<  points[3*tmp_facet[0]] << " " << points[3*tmp_facet[0]+1] << " " << points[3*tmp_facet[0]+2] << std::endl;
                        error(1) <<  points[3*tmp_facet[1]] << " " << points[3*tmp_facet[1]+1] << " " << points[3*tmp_facet[1]+2] << std::endl;
                        error(1) <<  points[3*tmp_facet[2]] << " " << points[3*tmp_facet[2]+1] << " " << points[3*tmp_facet[2]+2] << std::endl;


                        ++failed_triangle_cnt;
                    }
                    else
                    {
                        Builder.add_facet (tmp_facet, tmp_facet + 3);
                    }
                }
            }

            Builder.end_surface();

            if(failed_triangle_cnt > 0)
            {
                info(5) << "CGAL half edge data structure built, but "<< failed_triangle_cnt << " triangles were not included (maybe due to topological reasons)."  << std::endl;
            }
            else
            {
                info(5) << "CGAL half edge data structure successfully built" << std::endl;
            }


        }
    private:
        const std::vector<double> points;
        const std::vector<int> facets;
        const double scale = 1.0;
        

    };


    /* Creates a vector of triangles with the property that a triangle has at least one edge in common with any of its predecessors.
     * Uses recursion, potential STACK OVERFLOW! Use only for small meshes.
     *
     * Currently not used for cgal mesh creation.
     * */
    void build_triangles_by_neighbor(const MeshType& mesh, TriangleType& tri, std::vector<TriangleType>& triangles, std::set<int>& triangle_visited)
    {

        triangles.push_back(tri);
        triangle_visited.insert(tri.id().index());


        ConstNeighborRange cr(mesh, tri);

        for(ConstNeighborIterator ci = cr.begin(); ci != cr.end(); ++ci)
        {
            if(triangle_visited.find((*ci).id().index()) == triangle_visited.end()) //not found
            {
                TriangleType triangle = *ci;
                build_triangles_by_neighbor(mesh, triangle, triangles, triangle_visited);
            }
        }
    }

    /*bool order_triangles(TriangleType t1, TriangleType t2){

        //atm only x coordinate
       

        ConstVertexOnTriangleRange votr(t1);

        double min=viennagrid::get_point(mesh, votr[0]);

        for(std::size_t i = 1; i < 3; ++i)
        {
            if(min > viennagrid::get_point(mesh, votr[i])){

                min = viennagrid::get_point(mesh, votr[i]);

            }
        }

        ConstVertexOnTriangleRange votr(t2);

        for(std::size_t i = 0; i < 3; ++i)
        {
            if(min > viennagrid::get_point(mesh, votr[i])){
               return false;
            }
        }

        return true;



    }*/

    //TODO: DELTE
    void build_triangles_by_iteration_new(const MeshType& mesh,  std::vector<TriangleType>& triangles, std::vector<TriangleType> & non_manifold_triangles)
    {
        //make has map ??
        std::set<int> visited_edges; // edges that form parts of the boundary of already visited triangles
        std::vector<TriangleType> remaining_triangles; //triangles that have to be inserted into the triangle list

        remaining_triangles.reserve(triangles.size());
        int t_i;
        //first iteration through all triangles in mesh as provided by ViennaGrid data structure and delete degenerated triangles
        ConstTriangleRange tr(mesh);
        for(ConstTriangleIterator it = tr.begin(); it != tr.end(); ++it)
        {
            //check if triangle is degenerate
            typedef viennagrid::result_of::point<MeshType>::type PointType;

            //   for a specific triangle get all vertices
            PointType points[3];
            ConstVertexOnTriangleRange votr(*it);

            for(std::size_t i = 0; i < 3; ++i)
            {
                points[i] = viennagrid::get_point(mesh, votr[i]);
            }

            if((points[0] == points[1]) || (points[0] == points[2]) || (points[1] == points[2]))
            {
                //triangle is degenerate, ignore it by not adding it to triangles vector, print it to info(5)?
                //info(5) << "Degerate Triangle encountered: " << points[0] << ", " << points[1] << ", " << points[2] << "--> ignored for triangle list" << std::endl;
                continue;
            }
            //degeneracy check end

            //________________________________________Development area_______________________________________________________________________________________________________________________
            bool found = false;
            for(auto tri: non_manifold_triangles){
                if(tri.id().index() == (*it).id().index()){
                    t_i++;
                    found = true;
                    break;
                }

            }

            if(found)
                continue;



            //________________________________________Development area_______________________________________________________________________________________________________________________

      

            remaining_triangles.push_back(*it);
            
        }



        info(5) << "Finished building triangle list." << std::endl;

       /* info(5) << "Sorting triangle list..." << std::endl;

        std::sort(remaining_triangles.begin(), remaining_triangles.end(), order_triangles);

        info(5) << "Finished sorting triangles list." << std::endl;*/

        info(5) << "Ordering Triangle list ..." << std::endl;

        //add the first triangle to the mesh 
        triangles.push_back(remaining_triangles.front());

        ConstEdgeOnTriangleRange er(remaining_triangles.front());
        for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
        {
            visited_edges.insert((*eit).id().index());
        } 

        remaining_triangles.erase(remaining_triangles.begin());



        int triangles_prev_iteration = 0;

        std::reverse(remaining_triangles.begin(), remaining_triangles.end());



        //the remaining triangles are added by growing the mesh out from the already added triangles
        while(!remaining_triangles.empty()) //not empty
        {
            //think how to initialize
            std::vector<ConstEdgeOnTriangleRange> inserted;

            std::cout << remaining_triangles.size() << std::endl;

            //if in one iteration no new triangles are connected to the triangle list
            if(remaining_triangles.size() == triangles_prev_iteration){

                info(5) << "Warning: The provided mesh has " << triangles_prev_iteration << " triangles that have no common edge with the other triangles." << std::endl;
                break;
            }
             
            triangles_prev_iteration = remaining_triangles.size();
            



            //iterator increment is provided at the end because vector elements are possibly erased
            for(std::vector<TriangleType>::const_iterator vit = remaining_triangles.begin(); vit != remaining_triangles.end();)
            {
                bool is_bounding_edge = false; //is edge part of the boundary of any triangle already transversed?

                ConstEdgeOnTriangleRange er(*vit);
                for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
                {
                    if((visited_edges.find((*eit).id().index()) != visited_edges.end()) || visited_edges.empty())
                    {
                        is_bounding_edge = true;

                        triangles.push_back(*vit);
                        break;
                    }
                }


                if(is_bounding_edge == true)
                {

                    inserted.push_back(er);
                    vit = remaining_triangles.erase(vit); //erase returns iterator to the element that follows the one to erase
                }
                else
                {
                    ++vit;
                }
            }

            //std::cout << "inserted into visited : " << inserted.size()<< std::endl;
            //std::cout << "size remaining tri    : " << remaining_triangles.size() << std::endl;
            //std::cout << "size remaining tri var: " << triangles_prev_iteration << std::endl;
            //tmp insert all edges into visited edges
            for(auto er: inserted){

                for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
                {
                        visited_edges.insert((*eit).id().index());
                }
            }

            //if(inserted.size() > 800)
             //   break;

            inserted.clear();
            
        }
        





    }

    



    /* Creates a vector of triangles with the property that a triangle has at least one edge in common with any of its predecessors.
     * Carried out using iterations only.
    */
    void build_triangles_by_iteration(const MeshType& mesh,  std::vector<TriangleType>& triangles, std::vector<TriangleType> & non_manifold_triangles)
    {
        std::set<int> visited_edges; // edges that form parts of the boundary of already visited triangles
        std::vector<TriangleType> remaining_triangles; //triangles which do not fulfill the above given requirement are collected here

        remaining_triangles.reserve(triangles.size()/10);

        //std::sort (remaining_triangles.begin(), remaining_triangles.end());


        //first iteration through all triangles in mesh as provided by ViennaGrid data structure
        ConstTriangleRange tr(mesh);
        for(ConstTriangleIterator it = tr.begin(); it != tr.end(); ++it)
        {

            //check if triangle is degenerate
            typedef viennagrid::result_of::point<MeshType>::type PointType;

            //   for a specific triangle get all vertices
            PointType points[3];
            ConstVertexOnTriangleRange votr(*it);

            for(std::size_t i = 0; i < 3; ++i)
            {
                points[i] = viennagrid::get_point(mesh, votr[i]);
            }

            if((points[0] == points[1]) || (points[0] == points[2]) || (points[1] == points[2]))
            {
                //triangle is degenerate, ignore it by not adding it to triangles vector, print it to info(5)?
                //info(5) << "Degerate Triangle encountered: " << points[0] << ", " << points[1] << ", " << points[2] << "--> ignored for triangle list" << std::endl;
                continue;
            }
            //degeneracy check end

            //________________________________________Development area_______________________________________________________________________________________________________________________

            bool found = false;
            for(auto tri: non_manifold_triangles){
                if(tri.id().index() == (*it).id().index()){
                    found = true;
                    break;
                }

            }

            if(found)
                continue;

            //________________________________________Development area_______________________________________________________________________________________________________________________

            


            bool is_bounding_edge = false; //is edge part of the boundary of any triangle already transversed?

            ConstEdgeOnTriangleRange er(*it);
            for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
            {
                if((visited_edges.find((*eit).id().index()) != visited_edges.end()) || visited_edges.empty() ) //suitable new triangle found
                {
                    is_bounding_edge = true;

                    triangles.push_back(*it);
                    break;
                }
            }

            if(is_bounding_edge == true)
            {

                //add all boundary edges of the suitable triangle to the visited edges
                for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
                {
                    
                    visited_edges.insert((*eit).id().index());
                }
            }
            else
            {
                remaining_triangles.push_back(*it);
            }
        }

        info(5) << "Finished building triangle list." << std::endl;

        info(5) << "Insert triangles that did not have a connecting edge..." << std::endl;


        int num_triangles_prev = 0;


        //the remaining triangles are added by repeatedly iterating through remaining_triangles
        while(!remaining_triangles.empty()) //not empty
        {

//________________________________________Development area_______________________________________________________________________________________________________________________



            std::cout << remaining_triangles.size() << std::endl;

            if(remaining_triangles.size() == num_triangles_prev){

                //if this happens there are Triangles that are not connected to the processed mesh

                for(std::vector<TriangleType>::const_iterator vit = remaining_triangles.begin(); vit != remaining_triangles.end();vit++){

                    std::cout << "Triangle Points:";

                    ConstVertexOnTriangleRange votr(*vit);

                    std::cout.precision(std::numeric_limits< double >::max_digits10);

                    for(std::size_t i = 0; i < 3; ++i)
                    {
                        std::cout << "  " << viennagrid::get_point(mesh, votr[i]);
                    }

                    std::cout << std::endl;
                }              
                //info(5) << "Error: " << num_triangles_prev << " could not be inserted into the half edge datastructure" << std::endl;
                break;

            }
//________________________________________Development area_______________________________________________________________________________________________________________________
               
            num_triangles_prev = remaining_triangles.size();
            
            //iterator increment is provided at the end because vector elements are possibly erased
            for(std::vector<TriangleType>::const_iterator vit = remaining_triangles.begin(); vit != remaining_triangles.end();)
            {
                bool is_bounding_edge = false; //is edge part of the boundary of any triangle already transversed?

                ConstEdgeOnTriangleRange er(*vit);
                for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
                {
                    if((visited_edges.find((*eit).id().index()) != visited_edges.end()) || visited_edges.empty())
                    {
                        is_bounding_edge = true;

                        triangles.push_back(*vit);
                        break;
                    }
                }


                if(is_bounding_edge == true)
                {
                    for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
                    {
                        visited_edges.insert((*eit).id().index());
                    }

                    vit = remaining_triangles.erase(vit); //erase returns iterator to the element that follows the one to erase
                }
                else
                {
                    ++vit;
                }
            }
        }
        
    }

    //________________________________________Development area_______________________________________________________________________________________________________________________


    
    void analyze_viennagrid_mesh(viennagrid::mesh const & input, std::vector<TriangleType> & non_manifold_triangles){


        std::vector<int> edges;

        ConstCellRangeType ccr(input);

        //TODO RESERVES FAR TOO MUCH SPACE
        edges.resize(3*ccr.size());

        ConstTriangleRange tr(input);

        for(ConstTriangleIterator it = tr.begin(); it != tr.end(); ++it)
        {
            ConstEdgeOnTriangleRange er(*it);
            for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
            {
               
               edges[(*eit).id().index()]++;

            }

        }

        

        std::vector<TriangleType> test;

        //get non manifold triangles
        for(ConstTriangleIterator it = tr.begin(); it != tr.end(); ++it)
        {
            ConstEdgeOnTriangleRange er(*it);
            for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
            {

                //is edge non manifold edge
                if(edges[(*eit).id().index()] > 2){

                    //non_manifold_triangles.push_back((*it));

                    std::cout << "number of edges: " << edges[(*eit).id().index()] << std::endl;

                    test.push_back((*it));

                    break;
                }            
            }           
        }


        std::cout << "\nnonmanifold triangles"<< test.size() << std::endl;
            

         for(auto t: test){

            ConstEdgeOnTriangleRange er(t);
            bool found = false;
            for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
            {
                if(edges[(*eit).id().index()] == 1){
                    non_manifold_triangles.push_back(t);
                    //std::cout << t.id() << " ";
                    found = true;
                }
            }

            ConstVertexOnTriangleRange verts(t);

            
                std::cout << "Triangle Points:";

                std::cout.precision(std::numeric_limits< double >::max_digits10);

                for(std::size_t i = 0; i < 3; ++i)
                {
                    std::cout << "  " << viennagrid::get_point(input, verts[i]);
                }

                    std::cout << std::endl << std::endl;
            

            //if(found)
                //std::cout << std::endl;
            
        } 

     /*   for(auto t: test){
            ConstVertexOnTriangleRange votr(t);
    
            std::cout << "x=np.append (x, [" << viennagrid::get_point(input, votr[0])[0] << "])" << std::endl;
            std::cout << "y=np.append (y, [" << viennagrid::get_point(input, votr[0])[1] << "])" << std::endl;
            std::cout << "z=np.append (z, [" << viennagrid::get_point(input, votr[0])[2] << "])" << std::endl;

            std::cout << "x=np.append (x, [" << viennagrid::get_point(input, votr[1])[0] << "])" << std::endl;
            std::cout << "y=np.append (y, [" << viennagrid::get_point(input, votr[1])[1] << "])" << std::endl;
            std::cout << "z=np.append (z, [" << viennagrid::get_point(input, votr[1])[2] << "])" << std::endl;   

            std::cout << "x=np.append (x, [" << viennagrid::get_point(input, votr[2])[0] << "])" << std::endl;
            std::cout << "y=np.append (y, [" << viennagrid::get_point(input, votr[2])[1] << "])" << std::endl;
            std::cout << "z=np.append (z, [" << viennagrid::get_point(input, votr[2])[2] << "])" << std::endl;      
        }*/

        std::vector<int> border_count;

        //see if we can exclude some of the non manifold triangles
       /* for(ConstTriangleIterator it = tr.begin(); it != tr.end(); ++it)
        {
      

            ConstEdgeOnTriangleRange er(*it);
            for(ConstEdgeOnTriangleIterator eit = er.begin(); eit != er.end(); ++eit)
            {
                for(auto e: border_non_manifold_edges){

                    if((*eit).id().index() == e]){


                    }
                }
            }
        }*/
    }



    //________________________________________Development area_______________________________________________________________________________________________________________________

    /* Converts a triangular surface mesh given in ViennaGrid data structure to CGAL half edge data structure.
     * In order to achieve a valid half edge data structure degenerated triangles are ignored and for each triangle orientation
     * (which is defined by the order of vertices) is set appropriately.
     *
     * Conversion involves 3 steps:
     *  1. Collect all points in mesh --> vector points
     *  2. Transverse all triangles in mesh in such a way that the following property is fulfilled:
     *       the next triangle added has at least one edge in common with any of its predecessors. --> build_triangles_by_iteration()
     *  3. Build CGAL half edge data structure by iterating through the triangles as given in 2 and change orientation if necessary
     *          --> functor of type Build_triangle_mesh
     * */
    viennamesh_error convert(viennagrid::mesh const & input, cgal::polyhedron_surface_mesh & output)
    {

        std::vector<TriangleType> non_manifold_triangles;

      

        analyze_viennagrid_mesh(input, non_manifold_triangles);




        //return VIENNAMESH_ERROR_CONVERSION_FAILED;


        typedef cgal::polyhedron_surface_mesh::HalfedgeDS HalfedgeDS;

        //points in mesh, sequence is given by ViennaGrid data structure; format: point1_x, point1_y, point1_z, point2_x, ...
        std::vector<double> points;

        ConstVertexRangeType vertices(input);
        points.reserve(vertices.size()*3);


        int index=0;
        for (ConstVertexIteratorType vit = vertices.begin(); vit != vertices.end(); ++vit, ++index)
        {

                   
            points.push_back(viennagrid::get_point(input, *vit)[0]);
            points.push_back(viennagrid::get_point(input, *vit)[1]);
            points.push_back(viennagrid::get_point(input, *vit)[2]);
        
        }


        // vector of triangles with the property that a triangle has at least one edge in common with any of its predecessors.
        std::vector<TriangleType> triangles;

        // vector of vertex IDs; format: triangle1_vertex1, triangle1_vertex2, triangle1_vertex3, triangle2_vertex1, ...
        std::vector<int> faces;


        ConstCellRangeType ccr(input);
        triangles.reserve(ccr.size());



        info(5) << "\nBuilding triangle list for CGAL half edge data structure creation" << std::endl;

        auto start = std::chrono::high_resolution_clock::now();

        auto finish = std::chrono::high_resolution_clock::now();
            
        start = std::chrono::high_resolution_clock::now();
           
        build_triangles_by_iteration(input, triangles, non_manifold_triangles);

        finish = std::chrono::high_resolution_clock::now();

        info(1) << "runtime mesh conversion: " << std::chrono::duration_cast<std::chrono::milliseconds>(finish-start).count() << " ms" << std::endl;

        faces.reserve(triangles.size()*3);

        //In order to build CGAL half edge data structure for all triangles the IDs of its vertices are collected
        for(std::size_t i= 0; i < triangles.size(); ++i)
        {
            ConstVertexOnTriangleRange boundary_vertices( triangles[i]);

            for (ConstVertexOnTriangleIterator vit = boundary_vertices.begin(); vit != boundary_vertices.end(); ++vit)
            {
                faces.push_back((*vit).id().index());
            }
        }

        info(5) << "Generate CGAL half edge data structure" << std::endl;

        Build_triangle_mesh<HalfedgeDS> triangle_mesh(points, faces);

        output.delegate(triangle_mesh); //triangle mesh validity is checked (postcondition)

        if(!output.is_valid()){
            info(5) << "Witty message that mesh is not valid" << std::endl;
        }

        return VIENNAMESH_SUCCESS;
    }

    //Lagacy function (Christoph)
    /*int get_id_of_vertex(cgal::Point_3 p, cgal::polyhedron_surface_mesh const & input)
    {
        typedef cgal::polyhedron_surface_mesh::Vertex_const_iterator Vertex_iterator;
        int pos =0;
        {
            Vertex_iterator begin = input.vertices_begin();
            for ( ; begin != input.vertices_end(); ++begin, ++pos)
                if(begin->point()==p)
                    return pos;
        }
        return -1;
    }*/

    //Unordered map to store the vertex id
    typedef cgal::polyhedron_surface_mesh::Vertex_const_handle Vertex_const_handle_cgal;
    std::unordered_map<Vertex_const_handle_cgal, int> vertex_ids; 

    viennamesh_error convert(cgal::polyhedron_surface_mesh const & input, viennagrid::mesh & output)
    {
        
        //info(5) << "\nConverting CGAL half edge data structure to viennagrid" << std::endl;

        typedef cgal::polyhedron_surface_mesh::Vertex_const_iterator	Vertex_iterator;
        typedef cgal::polyhedron_surface_mesh::Facet_const_iterator		Facet_iterator;
        int numberofpoints =0;

        numberofpoints = input.size_of_vertices();


        std::vector<VertexType> vertex_handles(numberofpoints);

        vertex_ids.reserve(numberofpoints);


        // iterate through all CGAL vertices and store the coordinates in the viennagrid vertices for each
        {
            int pos = 0;
            Vertex_iterator v = input.vertices_begin();
            for (int i = 0; i < numberofpoints; ++i, ++v)
            {
                vertex_handles[i] = viennagrid::make_vertex( output,
                                    viennagrid::make_point(v->point().x(),v->point().y(),v->point().z())             
                                    );
            
                vertex_ids[v] = pos;
                pos++;


            }
        }

        // interate through all CGAL Facets and create a Viennagrid triangle each time
        Facet_iterator f = input.facets_begin();
        for ( ; f != input.facets_end(); ++f)
        {


            viennagrid::make_triangle(
                output,
                //lagacy code remove? (Christoph)
                //vertex_handles[CGAL_ID(begin->facet_begin()->vertex()] // ToDo (Florian)
                //vertex_handles[get_id_of_vertex(f->facet_begin()->vertex()->point(),input)],
                //vertex_handles[get_id_of_vertex(f->facet_begin()->next()->vertex()->point(),input)],
                //vertex_handles[get_id_of_vertex(f->facet_begin()->opposite()->vertex()->point(),input)]

                vertex_handles[vertex_ids[f->facet_begin()->vertex()]],
                vertex_handles[vertex_ids[f->facet_begin()->next()->vertex()]],
                vertex_handles[vertex_ids[f->facet_begin()->opposite()->vertex()]]
            );
        }

        //------------------------------------------------------------
        //---------------------------- END ---------------------------
        //------------------------------------------------------------

        return VIENNAMESH_SUCCESS;
    }

    template<>
    viennamesh_error internal_convert<viennagrid_mesh, cgal::polyhedron_surface_mesh>(viennagrid_mesh const & input, cgal::polyhedron_surface_mesh & output)
    {
        return convert( input, output );
    }

    template<>
    viennamesh_error internal_convert<cgal::polyhedron_surface_mesh, viennagrid_mesh>(cgal::polyhedron_surface_mesh const & input, viennagrid_mesh & output)
    {
        viennagrid::mesh output_pp(output);
        return convert( input, output_pp );
    }
}
