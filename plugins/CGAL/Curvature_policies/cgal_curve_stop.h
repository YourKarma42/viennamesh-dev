#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>


//#include "../cgal_curvature_functions.hpp"

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
class Curvature_stop_predicate
{
public:

  typedef ECM_ ECM ;
  
  typedef Edge_profile<ECM> Profile ;
  
  typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
  typedef typename boost::graph_traits<ECM>::edges_size_type size_type ;
    
  Curvature_stop_predicate( double aRatio ) : mRatio(aRatio) {}
  
  template <typename F> 
  bool operator()( F const&       // aCurrentCost
                 , Profile const&    aProfile
                 , size_type         aInitialCount
                 , size_type         aCurrentCount
                 ) const 
  {
   /* double flat_boundary = 0.00001;

    double c1 = viennamesh::cgal::mean_curvature_my(*(aProfile.v0()));

    double c2 = viennamesh::cgal::mean_curvature_my(*(aProfile.v1()));

    if(c1 <= flat_boundary && c2 <= flat_boundary){
      if((static_cast<double>(aCurrentCount) / static_cast<double>(aInitialCount))  < mRatio){
        return true;
      }else{
        return false;
      }
        
    }else{
      return true;
    }*/
    
    //return ( static_cast<double>(aCurrentCount) / static_cast<double>(aInitialCount) ) < mRatio ;

    return CGAL::squared_distance(aProfile.p0(), aProfile.p1()) > mRatio;


  }
  
private:
  
  double mRatio ;
};    

} // namespace Surface_mesh_simplification

} //namespace CGAL


