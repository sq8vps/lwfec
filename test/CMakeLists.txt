cmake_minimum_required(VERSION 3.5.0)
project(test VERSION 0.1.0 LANGUAGES C)

add_executable(test main.c ../gf.c ../rs.c)
target_include_directories(test PRIVATE ..)

