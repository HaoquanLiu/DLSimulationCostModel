project(Geom)

get_filename_component(IncludeParentDir ${PROJECT_SOURCE_DIR} DIRECTORY)
include_directories(${IncludeParentDir}/Necessity)

set(LinkLibs libNecessity.so)
link_directories(${IncludeParentDir}/build)

set(SOURCES ${PROJECT_SOURCE_DIR}/GeomDigging.cpp ${PROJECT_SOURCE_DIR}/GeomDumping.cpp)
 
add_library(Geom SHARED ${SOURCES})
target_link_libraries(Geom ${LinkLibs})