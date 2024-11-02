include(CTest)

find_package(Catch2 REQUIRED)
include(Catch)

add_executable(microtrainsim-test
  main.cpp
  test.cpp
  )
target_compile_options(microtrainsim-test
  PRIVATE -Wall -Wextra -Werror)
target_link_libraries(microtrainsim-test Catch2::Catch2)

catch_discover_tests(microtrainsim-test)
