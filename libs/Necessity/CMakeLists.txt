project(Necessity)

include_directories(${PROJECT_SOURCE_DIR})
set(SOURCES ${PROJECT_SOURCE_DIR}/IO.cpp ${PROJECT_SOURCE_DIR}/DataOperation.cpp ${PROJECT_SOURCE_DIR}/DataStruct.cpp ${PROJECT_SOURCE_DIR}/OptDataComputation.cpp ${PROJECT_SOURCE_DIR}/Cost.cpp)
 
add_library(Necessity SHARED ${SOURCES})
