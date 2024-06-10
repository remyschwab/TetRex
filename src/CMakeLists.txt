cmake_minimum_required (VERSION 3.8)

# A interface to reuse common properties.
# You can add more external include paths of other projects that are needed for your project.
add_library ("${PROJECT_NAME}_interface" INTERFACE)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries("${PROJECT_NAME}_interface" INTERFACE OpenMP::OpenMP_CXX)
endif()
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE seqan3::seqan3 sharg::sharg)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE ../include ../lib/simde ../lib/lemon/include)
target_include_directories ("${PROJECT_NAME}_interface" SYSTEM INTERFACE ../lib/simde ../lib/lemon/include)
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra" "-fopenmp")


# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib" STATIC utils.cpp index_base.cpp query.cpp inspect_idx.cpp construct_nfa.cpp)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC "${PROJECT_NAME}_interface" re2::re2 seqan::hibf)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")
