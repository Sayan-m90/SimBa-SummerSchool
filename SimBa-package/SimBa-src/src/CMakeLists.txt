#project( SimBa )

#cmake_minimum_required( VERSION 2.8 FATAL_ERROR )
#cmake_policy( SET CMP0017 NEW )
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


include_directories( ${CMAKE_SOURCE_DIR}/src "${CMAKE_SOURCE_DIR}/src/include" )
link_directories(${CMAKE_SOURCE_DIR}/src/lib)
# Boost
find_package( Boost 1.4 COMPONENTS program_options system REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

include_directories( "." ".." ) 

add_executable( SimBa
	ANNSearchSampling.h
	LongVector.h 
	PointSet.h
	SimpleGraph.h
	SimplexNode.h
	SimplicialComplex.h
	AnnotationMatrix.cpp
	ANNSearchSampling.cpp
	BCmain.cpp
	PointSet.cpp
	SimpleGraph.cpp
	UnionFindDeletion.cpp
	) 
target_link_libraries( SimBa ann.lib )
#target_link_libraries( SimBa ${Boost_LIBRARIES} libANN.a)