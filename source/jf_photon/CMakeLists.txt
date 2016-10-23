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

add_library(jf_photon SHARED jf_photon.cpp)

set_target_properties(jf_photon PROPERTIES PREFIX "")

target_link_libraries(jf_photon ${ARNOLD_LIBRARY})

# install(TARGETS jf_photon
#     LIBRARY DESTINATION plugins)

install(FILES mtoa/jf_photon.mtd
	DESTINATION plugins)

install(FILES mtoa/jf_photonTemplate.py
	DESTINATION mtoa)

install(FILES sitoa/jf_photon.spdl
	DESTINATION sitoa)