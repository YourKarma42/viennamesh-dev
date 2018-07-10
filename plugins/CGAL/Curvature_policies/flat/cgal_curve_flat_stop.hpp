#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>



namespace CGAL {

namespace Surface_mesh_simplification 
{

//*******************************************************************************************************************
//                                -= stopping condition predicate =-
//
// Determines whether the simplification has finished.
// The arguments are (current_cost,vertex,vertex,is_edge,initial_pair_count,current_pair_count,surface) and the result is bool
//
//*******************************************************************************************************************

// 
// Stops when the ratio of initial to current vertex pairs is below some value.
//
template<class ECM_>    
class Curvature_flat_stop
{
public:

  typedef ECM_ ECM ;
  
  typedef Edge_profile<ECM> Profile ;
  
  typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
  typedef typename boost::graph_traits<ECM>::edges_size_type size_type ;
    
  Curvature_flat_stop( double edge_length_threshold ) : m_edge_sq_length_threshold(edge_length_threshold*edge_length_threshold) {}
  //already implemented in higher cgal versions this version is 4.8.1
  template <typename F> 
  bool operator()( F const&       // aCurrentCost
                 , Profile const&    aProfile
                 , size_type         aInitialCount
                 , size_type         aCurrentCount
                 ) const 
  {
    // probably another hard break if edge length = penalty
  /*   std::cout << CGAL::squared_distance(aProfile.p0(), aProfile.p1())  << std::endl;
      std::cout <<  m_edge_sq_length_threshold << std::endl;
    std::cout << (CGAL::squared_distance(aProfile.p0(), aProfile.p1()) > m_edge_sq_length_threshold) << std::endl;
    */

    //as long as edges arnt too long
    return  CGAL::squared_distance(aProfile.p0(), aProfile.p1()) > m_edge_sq_length_threshold;
    
  }
  
private:
  
  double m_edge_sq_length_threshold ;
};    

} // namespace Surface_mesh_simplification

} //namespace CGAL


