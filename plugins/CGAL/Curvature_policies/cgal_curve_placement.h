#include <CGAL/Surface_mesh_simplification/Detail/Common.h>

//hash table to store curvatures
#include <unordered_map>

#include <CGAL/Curvature_policies/cgal_hash_points.hpp>


namespace CGAL {

namespace Surface_mesh_simplification  
{

  template<class ECM_>
class Curvature_placement
{
public:
    
  typedef ECM_ ECM ;

  typedef typename ECM::Point_3 Point_3;

  std::unordered_map<Point_3, viennamesh::cgal::vertex_statistics, viennamesh::cgal::Point_3_Hash, viennamesh::cgal::Point_3_Equal>   & curvatures;
  
public:

  Curvature_placement(std::unordered_map<Point_3, viennamesh::cgal::vertex_statistics, viennamesh::cgal::Point_3_Hash, viennamesh::cgal::Point_3_Equal> & c)
  : curvatures(c) {}
     
  template <typename Profile> 
  //optional<typename Profile::Point> //probably need a profile for later
  optional<typename Profile::Point> operator()( Profile const& aProfile ) const
  {
    //probably delete old vertices from the hastable


    double flat_boundary = 0.001;

    /*if(curvatures(aProfile.p0()) < flat_boundary && curvatures(aProfile.p1()) < flat_boundary){

      return optional<typename Profile::Point>(midpoint(aProfile.p0(),aProfile.p1()));
    } 
    else {
      //if p0 is flat then return the other point as ne point
      if(curvatures(aProfile.p0()) < flat_boundary){
        return optional<typename Profile::Point>(aProfile.p1()) ;
      }
      else{
        return optional<typename Profile::Point>(aProfile.p0()) ;
      }

    }*/

    



    return optional<typename Profile::Point>(midpoint(aProfile.p0(),aProfile.p1())) ;
  }
  
/*private:

  LindstromTurk_params mParams ;    */

};

} // namespace Surface_mesh_simplification

} //namespace CGAL

//TODO:  #defines machen wie in cgal