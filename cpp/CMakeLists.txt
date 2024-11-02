cmake_minimum_required(VERSION 3.26)

project(microtrainsim
  VERSION 1.0
  LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS system)

add_executable(microtrainsim src/main.cpp)
target_link_libraries(microtrainsim
  PRIVATE
  Boost::system)


option(TESTING "Enable tests." ON)

if(TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

install(TARGETS microtrainsim DESTINATION bin)
