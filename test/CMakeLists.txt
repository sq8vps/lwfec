cmake_minimum_required(VERSION 3.11.0)
project(test VERSION 0.1.0 LANGUAGES C)

add_executable(test main.c ../gf.c ../rs.c)
target_compile_definitions(test PUBLIC LWFEC_TEST=1)
#target_compile_definitions(test PUBLIC RS_USE_BUFFERS_PER_INSTANCE=1)
target_include_directories(test PRIVATE ..)

