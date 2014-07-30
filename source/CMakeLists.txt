set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Arnold
find_package(Arnold 4.0 REQUIRED)
include_directories(${ARNOLD_INCLUDE_DIR})

add_library(jf_nested_dielectric SHARED jf_nested_dielectric.cpp)

set_target_properties(jf_nested_dielectric PROPERTIES PREFIX "")
add_definitions(-w -O2 -fPIC -DLINUX -D_BOOL -ffast-math)

target_link_libraries(jf_nested_dielectric ${ARNOLD_LIBRARY})

install(TARGETS jf_nested_dielectric
    LIBRARY DESTINATION plugins)

install(FILES mtoa/jf_nested_dielectric.mtd
	DESTINATION plugins)

install(FILES mtoa/jf_nested_dielectricTemplate.py
	DESTINATION mtoa)

install(FILES sitoa/jf_nested_dielectric.spdl
	DESTINATION sitoa)

