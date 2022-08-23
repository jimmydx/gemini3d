include(cmake/cpu_count.cmake)

cmake_host_system_information(RESULT host_ramMB QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information(RESULT host_cpu QUERY PROCESSOR_DESCRIPTION)
math(EXPR host_ramGB "${host_ramMB} / 1000")
message(STATUS "Gemini3D: ${host_ramGB} GB RAM detected on ${CMAKE_HOST_SYSTEM_NAME} with ${host_cpu}.  Detected ${Ncpu} CPU cores.")
if(host_ramGB LESS 2)
  message(WARNING "Minimum RAM is about 2 GB--some tests or simulations may fail due to small memory (RAM)")
endif()


if(realbits EQUAL 32)
  message(VERBOSE " 32-bit real precision")
  set(arith s)
else()
  message(VERBOSE " 64-bit real precision")
  set(realbits 64)
  set(arith d)
endif()

option(dev "Gemini developer mode")

option(cpp "also build Gemini3D C++ frontend prototype" on)

option(glow "use NCAR GLOW airglow / aurora model" on)

option(hwm14 "use HWM14 neutral winds model")
option(msis2 "use MSIS 2.x neutral atmosphere model")

option(python "PyGemini checks")
# Matlab checks take much longer than Python, and Python covers much more
option(matlab "Matlab checks")

set(CMAKE_TLS_VERIFY true)  # for Git and Downloads

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

# append .debug to debug libraries, because the computation speed penalty is so great
set(CMAKE_DEBUG_POSTFIX .debug)

# to make Gemini3D more usable by external programs, put all Fortran .mod generated module files in a single directory.
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# --- auto-ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()
