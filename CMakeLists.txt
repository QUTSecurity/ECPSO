cmake_minimum_required(VERSION 3.10)
project(ECPSO_Project CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find MPI. This is crucial for cross-platform MPI support
if(WIN32)
    set(MPI_CXX_INCLUDE_PATH "C:/Program Files (x86)/Microsoft SDKs/MPI/Include")
    set(MPI_CXX_LIBRARIES "C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64/msmpi.lib")
    include_directories(${MPI_CXX_INCLUDE_PATH})
else()
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(MPI_CXX_LIBRARIES MPI::MPI_CXX)
endif()

# ECPSO Executable
add_executable(ecpso 
    ECPSO/algorithms.cpp
    ECPSO/fitness.cpp
    ECPSO/latin.cpp
    ECPSO/main.cpp
    ECPSO/variable.cpp
)
target_include_directories(ecpso PRIVATE ECPSO)
target_link_libraries(ecpso PRIVATE ${MPI_CXX_LIBRARIES})

