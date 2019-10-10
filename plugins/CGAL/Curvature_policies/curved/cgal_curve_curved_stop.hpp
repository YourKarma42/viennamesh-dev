#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>


//analytics
#include "../../analizing/cgal_mesh_analytics.hpp"

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
class Curvature_curved_stop
{
public:

  typedef ECM_ ECM ;
  
  typedef Edge_profile<ECM> Profile ;

  
  typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
  typedef typename boost::graph_traits<ECM>::edges_size_type size_type ;
    
  Curvature_curved_stop( double curved_stop_ratio, viennamesh::cgal::cgal_mesh_analytics<ECM>  & a ) 
  : m_curved_stop_ratio(curved_stop_ratio), analytics(a)
  {}

  template <typename F> 
  bool operator()( F const&       // aCurrentCost
                 , Profile const&    aProfile
                 , size_type         aInitialCount
                 , size_type         aCurrentCount
                 ) const 
  {

    //int init_num_curved_edges = aInitialCount - analytics.get_num_flat_edges() - analytics.get_num_curved_edges();



    return std::sqrt(CGAL::squared_distance(aProfile.p0(), aProfile.p1())) > m_curved_stop_ratio;


  }
  
private:

  viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;
  
  double m_curved_stop_ratio ;
};    

} // namespace Surface_mesh_simplification

} //namespace CGAL


