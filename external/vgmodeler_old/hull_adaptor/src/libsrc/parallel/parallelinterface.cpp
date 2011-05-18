#ifdef PARALLEL
#include <mystdlib.h>
#include <parallel.hpp>

#include <meshing.hpp>

#include "parallelinterface.hpp"

namespace netgen
{
  extern AutoPtr<Mesh> mesh;



using namespace netgen;

  int NgPar_GetLoc2Glob_VolEl ( const int & locnum )
  { return mesh -> GetParallelTopology().GetLoc2Glob_VolEl ( locnum+1) -1; }

int NgPar_Glob2Loc_SurfEl ( const int & locnum ) 
{ return mesh->GetParallelTopology().Glob2Loc_SurfEl(locnum+1) - 1; }

int NgPar_Glob2Loc_VolEl ( const int & locnum ) 
{ return mesh->GetParallelTopology().Glob2Loc_VolEl(locnum+1) - 1; }

 int NgPar_Glob2Loc_Segm (   const int & locnum )   
{ return mesh->GetParallelTopology().Glob2Loc_Segm(locnum+1) - 1; }

int NgPar_Glob2Loc_Vert ( const int & locnum )  
{ return mesh->GetParallelTopology().Glob2Loc_Vert(locnum+1) -1; }


int NgPar_GetDistantPNum ( const int & proc, const int & locpnum )  
{ 
  return mesh->GetParallelTopology().GetDistantPNum( proc, locpnum+1) - 1; 
}


int NgPar_GetDistantEdgeNum ( const int & proc, const int & locpnum )  
{ 
  return mesh->GetParallelTopology().GetDistantEdgeNum( proc, locpnum+1) - 1; 
}

int NgPar_GetDistantFaceNum ( const int & proc, const int & locpnum )  
{ 
  return mesh->GetParallelTopology().GetDistantFaceNum (proc, locpnum+1 ) - 1; 
}

int NgPar_GetDistantElNum ( const int & proc, const int & locelnum )
{ 
  return mesh->GetParallelTopology().GetDistantElNum (proc, locelnum+1 ) - 1; 
}

bool NgPar_IsExchangeFace ( const int & fnr ) 
{ return mesh->GetParallelTopology().IsExchangeFace ( fnr+1 ); }

  bool NgPar_IsExchangeVert ( const int & vnum )
{ return mesh->GetParallelTopology().IsExchangeVert ( vnum+1 ); }

bool NgPar_IsExchangeEdge ( const int & ednum )  
{ return mesh->GetParallelTopology().IsExchangeEdge ( ednum+1 ); }


  void NgPar_PrintParallelMeshTopology ()
  { mesh -> GetParallelTopology().Print (); }
 
bool NgPar_IsElementInPartition ( const int elnum, const int dest )
{
  return mesh -> GetParallelTopology().IsElementInPartition ( elnum+1, dest ); 
}



bool NgPar_IsGhostFace ( const int facenum )
{ return mesh -> GetParallelTopology().IsGhostFace ( facenum+1); }

bool NgPar_IsGhostEdge ( const int edgenum )
{ return mesh -> GetParallelTopology().IsGhostEdge ( edgenum+1); }

}

#endif