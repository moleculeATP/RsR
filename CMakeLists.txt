cmake_minimum_required(VERSION 3.10)
project(RsR)
set(CMAKE_CXX_STANDARD 20)

# --- Build type (pour éviter le warning CGAL) ---
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- Dépendances ---
find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(${Boost_INCLUDE_DIR})
include_directories($ENV{CONDA_PREFIX}/include)
link_directories($ENV{CONDA_PREFIX}/lib)

# --- Dossiers d'include ---
set(SOURCE_DIR "src")
include_directories(${SOURCE_DIR})
include_directories(3rd/SING_3D/src)
include(3rd/cmake/polyscope.cmake)

# --- Fichiers du projet principal ---
file(GLOB SOURCES "${SOURCE_DIR}/*.cpp" "${SOURCE_DIR}/*.h")

# --- Fichiers du sous-module sing_3D ---
set(SING_3D_SOURCES
    3rd/SING_3D/src/distances.h
    3rd/SING_3D/src/distances.cpp
    3rd/SING_3D/src/io_sing.h
    3rd/SING_3D/src/io_sing.cpp
)

# --- Point d’entrée ---
set(ENTRY_POINT "${SOURCE_DIR}/Recon.cpp")


# --- Exécutable ---
add_executable(${PROJECT_NAME}
    ${ENTRY_POINT}
    ${SOURCES}
    ${SING_3D_SOURCES}
)

# --- Liens ---
target_link_libraries(${PROJECT_NAME} PRIVATE
    CGAL::CGAL
    Eigen3::Eigen
    boost_filesystem
    boost_system
    polyscope::polyscope
)








# cmake_minimum_required(VERSION 3.10)
# set(CMAKE_CXX_STANDARD 20)

# # set the project name
# project(RsR)


# include_directories(${Boost_INCLUDE_DIR})
# find_package(CGAL CONFIG REQUIRED)
# find_package(Eigen3 CONFIG REQUIRED)
# include_directories(${EIGEN3_INCLUDE_DIR})
# # find_package(Open3D REQUIRED)
# # find_package(TBB CONFIG REQUIRED)

# include_directories($ENV{CONDA_PREFIX}/include)
# link_directories($ENV{CONDA_PREFIX}/lib)

# # add_subdirectory(src)
# add_subdirectory(3rd/SING_3D)

# # Set the directory containing the .h and .cpp files
# set(SOURCE_DIR "src")

# # Include all .cpp and .h files from the source directory
# file(GLOB SOURCES "${SOURCE_DIR}/*.cpp" "${SOURCE_DIR}/*.h")

# # Specify the entry point file
# set(ENTRY_POINT "${SOURCE_DIR}/Recon.cpp")

# # add the executable
# add_executable(RsR ${ENTRY_POINT} ${SOURCES})


# target_link_libraries(${PROJECT_NAME} PRIVATE
#     CGAL::CGAL
#     Eigen3::Eigen
#     sing_3D_lib
#     boost_filesystem
#     boost_system
# )
