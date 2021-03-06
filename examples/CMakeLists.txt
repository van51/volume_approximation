project( GeomRandWalks_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

find_package(CGAL QUIET COMPONENTS )
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

#FIND_PACKAGE( OpenMP REQUIRED)
#if(OPENMP_FOUND)
#message("OPENMP FOUND")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  #include_directories (BEFORE ../external/extreme_points_d/include)
  include_directories (BEFORE ../../nanoflann/include)
  include_directories (BEFORE /home/aneva/Libraries/json/src)
  include_directories (BEFORE ../external/eigen_3_3)
  include_directories (BEFORE ../include)
  include_directories (BEFORE ../external/FALCONN/src/include)
  include_directories (BEFORE ../external/LSHBOX/include)
  
  #include_directories (BEFORE ../patches/include)
  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")
  #add_definitions(${CMAKE_CXX_FLAGS} "-std=c++0x")
  add_definitions(${CMAKE_CXX_FLAGS} "-pthread")
  add_definitions(${CMAKE_CXX_FLAGS} "-DRKD_PAR")
  add_definitions(${CMAKE_CXX_FLAGS} "-O3")
  #add_definitions(${CMAKE_CXX_FLAGS} "-g")

  #link_libraries("proc")
  #create_single_source_cgal_program("opt.cpp")
  create_single_source_cgal_program("vol.cpp")
  create_single_source_cgal_program("boundary_test.cpp")
  target_link_libraries(boundary_test boost_program_options)
  #  create_single_source_cgal_program("har_acc.cpp")
  #  target_link_libraries(har_acc boost_program_options)
  #create_single_source_cgal_program("min_ellipsoid.cpp")
  #create_single_source_cgal_program("polytope_generators/cross.cpp")
  #create_single_source_cgal_program("polytope_generators/skinny_cube.cpp")
  create_single_source_cgal_program("polytope_generators/rh_d_m.cpp")
  #create_single_source_cgal_program("polytope_generators/rrh_d_m.cpp")
  #create_single_source_cgal_program("polytope_generators/simplex.cpp")
  #create_single_source_cgal_program("polytope_generators/birk.cpp")
  #create_single_source_cgal_program("NN_cgal.cpp")
  #create_single_source_cgal_program("mink_sum_vol.cpp")

else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
