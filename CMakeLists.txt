cmake_minimum_required(VERSION 3.16)
project(clinvar_tool LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(clinvar_tool
  src/clinvar_tool.cpp
)

find_package(ZLIB REQUIRED)
target_link_libraries(clinvar_tool PRIVATE ZLIB::ZLIB)

# Drobna higiena kompilacji
if (MSVC)
  target_compile_options(clinvar_tool PRIVATE /W4)
else()
  target_compile_options(clinvar_tool PRIVATE -Wall -Wextra -Wpedantic)
endif()
