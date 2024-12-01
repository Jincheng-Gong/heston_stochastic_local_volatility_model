# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.20)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Set the project name
project(HestonSLVMC)

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/Bins.cpp
    src/Blackscholesformulas.cpp
    src/DupireLocalVolatilitySurface.cpp
    src/Gaussienne.cpp 
    src/HestonSolver.cpp
    src/ImpliedVolatilitySurface.cpp
    src/Modeles.cpp 
    src/PathSimulator.cpp
    src/Payoff.cpp 
    src/Pricer.cpp 
    src/main.cpp 
    src/ThomasSolver.cpp
)

# Add an executable with the above sources
add_executable(HestonSLVMC ${SOURCES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(HestonSLVMC
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)
