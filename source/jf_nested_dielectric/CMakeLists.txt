set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Arnold
find_package(Arnold 4.0 REQUIRED)
include_directories(${ARNOLD_INCLUDE_DIR})

set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_SKIP_RPATH TRUE)

# Fix to stop crash on OSX < 10.9
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.8)

# Compiler flags
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    # Disable MSVC warnings
    add_definitions( "-W3 -D_CRT_SECURE_NO_WARNINGS -wd4005 -wd4996 -wd4305 -wd4244 -nologo" )
else()
    # -fPIC should be added automatically for shared libraries
	add_definitions(-w -D_BOOL -ffast-math)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fvisibility=hidden")
endif()

add_library(jf_nested_dielectric SHARED jf_nested_dielectric.cpp)

set_target_properties(jf_nested_dielectric PROPERTIES PREFIX "")

target_link_libraries(jf_nested_dielectric ${ARNOLD_LIBRARY})

# install(TARGETS jf_nested_dielectric
#     LIBRARY DESTINATION plugins)

install(FILES mtoa/jf_nested_dielectric.mtd
	DESTINATION plugins)

install(FILES mtoa/jf_nested_dielectricTemplate.py
	DESTINATION mtoa)

install(FILES sitoa/jf_nested_dielectric.spdl
	DESTINATION sitoa)