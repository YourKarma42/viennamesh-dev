



#include <unordered_map>

#include "../cgal_mesh.hpp"

#include <sstream>

#ifndef HASH_FUNC_H
#define HASH_FUNC_H

//Alles anders machen nur tests

namespace viennamesh
{
  namespace cgal
  {

    // vertex handels sollten standardmäßig hasbar sein siehe
    //https://doc.cgal.org/latest/Polyhedron/classCGAL_1_1Polyhedron__3.html
    //All handles are model of LessThanComparable and Hashable, that is they can be used as keys in containers such as std::map and boost::unordered_map. 

    typedef polyhedron_surface_mesh::Point_3 Point_3;

    struct vertex_statistics {

        double curvature_1;

        double curvature_2;

        double gauss_curvature;

    };


    typedef double FT;


    struct Point_3_Hash {
    
        std::size_t operator()(const Point_3& p) const {
            //testen ob hash besser wenn die letzen kommastellen als haswert genommen werden
            // hash(x.yz) = z ?

            double table_size = 10.0;
            return std::hash<double>()(p.x() * table_size + p.y()*table_size + p.z()*table_size);
        }
        
    };

    struct Point_3_Equal {
        
        bool operator()(const Point_3& p1, const Point_3& p2) const {
            return (p1 == p2);
        }
    
    };

    std::string Hashtable_Statistics(std::unordered_map<Point_3, viennamesh::cgal::vertex_statistics, viennamesh::cgal::Point_3_Hash, viennamesh::cgal::Point_3_Equal> & m){
        std::stringstream  out;

        int buckets=m.bucket_count();

        int max_size_bucket=0;

        int not_empty = 0;

        int elements = 0;

        out << std::endl << "----------------++++++++++++++++++ HASH Statisticts ++++++++++++++++++----------------" << std::endl;

        out << "BUCKETS: " << buckets << std::endl;



        for(int i = 0; i < buckets; i++){
            int size=0;
            for ( auto local_it = m.begin(i); local_it!= m.end(i); ++local_it )
                size++;

            elements += size;

            if(size > max_size_bucket)
                max_size_bucket = size;

            if(size > 0)
                not_empty++;
        }


        out << "Elemente: " << elements << std::endl;

        out << "Größtes Bucket: " << max_size_bucket << std::endl;

        out << "Nicht leere Buckets: " << not_empty << std::endl;


        out << "--------------------------------------------------------------------------------" << std::endl;


        return out.str();

    }

  }
}

#endif /* HASH_FUNC_H */