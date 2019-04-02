SET(TARGET_H
	curvature_calculator.hpp
)


SET(TARGET_SRC
	curvature_calculator.cpp
)

#add_library(library STATIC curvature_calculator.hpp)
#add_library(library curvature_calculator.cpp curvature_calculator.hpp)

#Eigen include path
message(STATUS "Including Eigen from System path")
INCLUDE_DIRECTORIES ( "/usr/include/eigen3" )


set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )


if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "" FORCE)
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()
endif()

include( ${CGAL_USE_FILE} )

find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()
endif()




VIENNAMESH_ADD_PLUGIN(cgal_module plugin.cpp
 		      cgal_mesh.cpp
          cgal_mesh_simplification.cpp
		      cgal_automatic_mesh_simplification.cpp
					cgal_curve.cpp
					cgal_remove_short_edges.cpp
					cgal_compare_curvature_calculations.cpp 
					cgal_simplify_curve.cpp
					#needs to be added for .hpp
					curvature_calculator.cpp
					cgal_data_creator.cpp
					)
