# this script is to install hwloc, to identify host CPU paramters for gemini3d.run
#
# cmake -P install_hwloc.cmake
# will install hwloc under the user's home directory.

cmake_minimum_required(VERSION 3.18...${CMAKE_VERSION})

set(CMAKE_TLS_VERIFY true)

if(NOT prefix)
  get_filename_component(prefix ~ ABSOLUTE)
endif()

set(version 2.4.1)

set(host https://download.open-mpi.org/release/hwloc/v2.4)

if(APPLE)
  find_program(brew
    NAMES brew
    PATHS /usr/local /opt/homeebrew
    PATH_SUFFIXES bin)

  if(brew)
    execute_process(COMMAND ${brew} install hwloc)
  endif()

  return()
endif()

if(WIN32)
  set(stem hwloc-win64-build-${version})
  set(name ${stem}.zip)
else()
  set(stem hwloc-${version})
  set(name ${stem}.tar.bz2)
endif()

set(url ${host}/${name})

get_filename_component(prefix ${prefix} ABSOLUTE)
set(path ${prefix}/${stem})

message(STATUS "installing hwloc ${version} to ${path}")

set(archive ${path}/${name})

if(EXISTS ${archive})
  file(SIZE ${archive} fsize)
  if(fsize LESS 100000)
    file(REMOVE ${archive})
  endif()
endif()

if(NOT EXISTS ${archive})
  message(STATUS "download ${url}")
  file(DOWNLOAD ${url} ${archive})
endif()

message(STATUS "extracting to ${path}")
file(ARCHIVE_EXTRACT INPUT ${archive} DESTINATION ${path})

if(WIN32)
  find_program(lstopo NAMES lstopo PATHS ${path}/${stem} PATH_SUFFIXES bin NO_DEFAULT_PATH)
  if(lstopo)
    get_filename_component(pathbin ${lstopo} DIRECTORY)
    message(STATUS "add to environment variable PATH ${pathbin}")
    message(STATUS "add environment variable HWLOC_ROOT ${path}/${stem}")
  else()
    message(FATAL_ERROR "failed to install from ${archive}")
  endif()

  return()
endif()

message(STATUS "Building HWLOC")
if(NOT EXISTS ${path}/${stem}/Makefile)
  execute_process(COMMAND ./configure --prefix=${path} WORKING_DIRECTORY ${path}/${stem})
endif()
execute_process(COMMAND make -j -C ${path}/${stem})
execute_process(COMMAND make install -C ${path}/${stem} RESULT_VARIABLE err)

if(NOT err EQUAL 0)
  message(FATAL_ERROR "Hwloc failed to build via autotools")
endif()

find_program(lstopo NAMES lstopo PATHS ${path} PATH_SUFFIXES bin NO_DEFAULT_PATH)
  if(lstopo)
    get_filename_component(pathbin ${lstopo} DIRECTORY)
    message(STATUS "add to environment variable PATH ${pathbin}")
    message(STATUS "add environment variable HWLOC_ROOT ${path}")
  else()
    message(FATAL_ERROR "failed to install from ${archive}")
  endif()
  return()
