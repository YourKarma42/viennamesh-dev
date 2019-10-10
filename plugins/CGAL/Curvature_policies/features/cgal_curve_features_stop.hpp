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

template<class ECM_>    
class Curvature_features_stop
{
public:

    typedef ECM_ ECM ;
  
    typedef Edge_profile<ECM> Profile ;
  
    typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
    typedef typename boost::graph_traits<ECM>::edges_size_type size_type ;
    
    viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;

    Curvature_features_stop(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a, double el) 
    : analytics(a), mEdgeLength(el) {}

  
    template <typename F> 
    bool operator()( F const&        aCurrentCost
                    , Profile const&    aProfile
                    , size_type         aInitialCount
                    , size_type         aCurrentCount
                    ) const 
    {

        if((analytics.get_feature(aProfile.v0()) != 30 && analytics.get_feature(aProfile.v1()) != 30)){
            return std::sqrt(CGAL::squared_distance(aProfile.p0(), aProfile.p1())) > mEdgeLength;
        }else{
            return true;
        }

        
            

    }
  
    private:
    
    double mEdgeLength ;


};    

} // namespace Surface_mesh_simplification

} //namespace CGAL





