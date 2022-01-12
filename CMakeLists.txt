cmake_minimum_required(VERSION 2.8)
#set(CMAKE_C_COMPILER D:/MinGW/bin/gcc.exe)
#set(CMAKE_CXX_COMPILER D:/MinGW/bin/g++.exe)
#set(CMAKE_MAKE_PROGRAM D:/MinGW/bin/mingw32-make.exe)
project(ARCLENGTH)
aux_source_directory(./src DIR_SRCS)
add_executable(arclength ${DIR_SRCS})

#cmd执行命令：cmake -G "MinGW Makefiles" ..
#注：".."代表上一级目录