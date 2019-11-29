#ifndef VIENNAGRID_SILVACO_STR_READER_HPP
#define VIENNAGRID_SILVACO_STR_READER_HPP

/* =======================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                     ViennaGrid - The Vienna Grid Library
                            -----------------

   License:      MIT (X11), see file LICENSE in the base directory
======================================================================= */


#include <fstream>
#include <iostream>
#include <assert.h>

#include "viennagrid/mesh/mesh.hpp"




/** @file viennagrid/io/tetgen_poly_reader.hpp
    @brief Provides a reader for Tetgen .poly files. See http://wias-berlin.de/software/tetgen/fformats.poly.html
*/



namespace viennagrid
{
  namespace io
  {
          template<typename stream_type>
          bool _my_get_valid_line( stream_type  & stream, std::string & line, char comment_line)
          {
            std::string tmp;

            while (true)
            {
              if (!std::getline(stream, tmp))
                return false;

              std::size_t pos = tmp.find(comment_line);
              if (pos != std::string::npos)
                tmp = tmp.substr(0, pos);

              for (std::size_t i = 0; i != tmp.size(); ++i)
              {
                if ( tmp[i] != ' ' )
                {
                  line = tmp.substr(i, std::string::npos);
                  return true;
                }
              }
            }
          }



    /** @brief Reader for Tetgen .poly files.
      *
      * See http://wias-berlin.de/software/tetgen/fformats.poly.html for a description of the file format
      */
    struct silvaco_str_reader
    {
      //public:
      /** @brief The functor interface triggering the read operation. Segmentations are not supported in this version.
       *
       * @param mesh_obj      The mesh where the file content is written to
       * @param filename      Name of the file
       */

      template <typename MeshT>
      void operator()(MeshT const & mesh, std::string const & filename) const
      {
        typedef typename viennagrid::result_of::point<MeshT>::type           PointType;


        static const std::size_t point_dim = 3;//viennagrid::result_of::static_size<PointType>::value;

        //typedef typename result_of::vertex<MeshT>::type         VertexType;

          //typedef typename viennagrid::result_of::point<MeshT>::type PointType;
        typedef typename viennagrid::result_of::element<MeshT>::type VertexType;
        typedef typename VertexType::id_type VertexIDType;

        std::ifstream reader(filename.c_str());

        #if defined VIENNAGRID_DEBUG_STATUS || defined VIENNAGRID_DEBUG_IO
        std::cout << "* silvaco_str_reader::operator(): Reading file " << filename << std::endl;
        #endif

        if (!reader)
        {
          throw cannot_open_file_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": Cannot open file!");
        }

        if (!reader.good())
        {
          throw bad_file_format_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": File is empty.");
        }







        std::string tmp;
        std::string local_tmp;
        std::string type;
        std::istringstream current_line;

        long vertex_count;
        long triangle_count;

/*        while (true)
        {
          if (!_my_get_valid_line(reader, tmp, '#'))
            throw bad_file_format_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": EOF encountered when reading information");

          if (tmp.find("silvaco.file binary") == 0)
          {
            std::cout << "binary file found!" << std::endl;
            break;
          }

          if (tmp.find("silvaco.file ASCII") == 0)
          {
            std::cout << "ASCII file found!" << std::endl;
            break;
          }
        }

        throw bad_file_format_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": EOF encountered when reading information");
*/
        
        while (true)
        {
          if (!_my_get_valid_line(reader, tmp, '#'))
            throw bad_file_format_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": EOF encountered when reading information");

          if (tmp.find("k 3") == 0)
          {
            current_line.str(tmp); current_line.clear();
            current_line >> local_tmp >> local_tmp >> vertex_count >> local_tmp >> triangle_count;

            break;
          }
        }
        
        //typedef typename viennagrid::result_of::vertex_handle<MeshT>::type VertexHandleType;
        //std::map<long, VertexHandleType> vertex_handles;

        std::vector<VertexType> vertices = std::vector<VertexType>();
        vertices.resize(vertex_count+1);

        //read vertices from .str        
        for (long i = 0; i < vertex_count; ++i)
        {
          if (!_my_get_valid_line(reader, tmp, '#'))
            throw bad_file_format_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": EOF encountered when reading information");

          PointType point(point_dim);

          long id;

          current_line.str(tmp); current_line.clear();
          current_line >> type >> id;

          if (type != "c")
            throw bad_file_format_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": expected point, got \"" + tmp + "\"");



          for (std::size_t j = 0; j < point_dim; ++j)
            current_line >> point[j];
          
          vertices[id] = viennagrid::make_vertex(mesh, point);

        }

        

        std::vector<std::tuple<long,long,long>> triangle_list;
        triangle_list.resize(triangle_count);


        //read triangles and material from .str
        for (long i = 0; i < triangle_count; ++i)
        {
          if (!_my_get_valid_line(reader, tmp, '#'))
            throw bad_file_format_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": EOF encountered when reading information");

          long id;
          long material_id;
          long vertex_id0;
          long vertex_id1;
          long vertex_id2;

          current_line.str(tmp); current_line.clear();
          current_line >> type >> id >> material_id >> vertex_id0 >> vertex_id1 >> vertex_id2;

          //store triangle list for later maby delete later
          triangle_list[i] = std::make_tuple(vertex_id0, vertex_id1, vertex_id2);

          if (type != "t")
            throw bad_file_format_exception("* ViennaGrid: silvaco_str_reader::operator(): File " + filename + ": expected triangle, got \"" + tmp + "\"");

          //TODO: read all materials / only surface / save surface information  

          if(material_id == 1){

          viennagrid::make_triangle(mesh, vertices[vertex_id0],
                                          vertices[vertex_id1],
                                          vertices[vertex_id2]);
          }
        }

        

        std::vector<long> triangles;
        std::string blub;

       /* do{

          _my_get_valid_line(reader, tmp, '$');
          current_line.str(tmp); current_line.clear();
          current_line >> type >> local_tmp >> blub;

          if(type == "w"){

            long triangle_id;

            std::cout << local_tmp << blub << std::endl;

            _my_get_valid_line(reader, tmp, '$');
            current_line.str(tmp); current_line.clear();
            current_line >> type >> triangle_id;

            while(type == "b"){

              triangles.push_back(triangle_id);

              viennagrid::make_triangle(mesh, vertices[std::get<0>(triangle_list)],
                                          vertices[std::get<1>(triangle_list)],
                                          vertices[std::get<2>(triangle_list)]);


              _my_get_valid_line(reader, tmp, '$');
              current_line.str(tmp); current_line.clear();
              current_line >> type >> triangle_id;

            }
            break;

          }
            
        }while(type != "EOF");*/

        


        //__________TESTCODE_________________________________________________

        /*for(auto t: triangles){
          long o =0;
          for(auto t_i: triangles){
            if(t == t_i){
              
              if(o != 0)
                std::cout << "it is double!!!!!" << std::endl;

              o++;
            }  
          }
        }*/



      } //operator()

    }; //class tetgen_poly_reader


  } //namespace io
} //namespace viennagrid

#endif
