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
class Curvature_boundary_stop
{
public:

    typedef ECM_ ECM ;
  
    typedef Edge_profile<ECM> Profile ;
  
    typedef typename boost::graph_traits<ECM>::edge_descriptor edge_descriptor ;
    typedef typename boost::graph_traits<ECM>::edges_size_type size_type ;
    
    viennamesh::cgal::cgal_mesh_analytics<ECM>  & analytics;

    Curvature_boundary_stop(viennamesh::cgal::cgal_mesh_analytics<ECM>  & a, double r) 
    : analytics(a), mRatio(r) {}

  
    template <typename F> 
    bool operator()( F const&        aCurrentCost
                    , Profile const&    aProfile
                    , size_type         aInitialCount
                    , size_type         aCurrentCount
                    ) const 
    {


            if(aCurrentCost == 100 )
                return true;

            return false;
           /* std::cout << analytics.get_color(aProfile.v0()) << " " << analytics.get_color(aProfile.v1()) << std::endl;
            std::cout << aProfile.p0() << " " << aProfile.p1() << std::endl;
            analytics.num++;*/

        

        //return ( static_cast<double>(aCurrentCount) / static_cast<double>(aInitialCount) ) < mRatio ;
    }
  
    private:
    
    double mRatio ;


};    

} // namespace Surface_mesh_simplification

} //namespace CGAL


