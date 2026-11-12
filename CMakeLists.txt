cmake_minimum_required(VERSION 3.21)
project(nsga2 VERSION 0.2.0 LANGUAGES Fortran)

add_library(nsga2 src/nsga2.f90)
set_target_properties(nsga2 PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
target_include_directories(nsga2 PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/mod>)

if(PROJECT_IS_TOP_LEVEL)
  enable_testing()

  add_library(nsga2_test_support OBJECT test/check.f90 test/benchmarks.f90)
  target_link_libraries(nsga2_test_support PUBLIC nsga2)
  set_target_properties(nsga2_test_support PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/test_mod)
  target_include_directories(nsga2_test_support PUBLIC ${PROJECT_BINARY_DIR}/test_mod)

  foreach(t IN ITEMS test_operators test_benchmarks test_errors)
    add_executable(${t} test/${t}.f90)
    target_link_libraries(${t} PRIVATE nsga2_test_support)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()

  # Expected failures: each must exit with its specific error (a crash does not count).
  set(xfail_invalid_input_without_stat_expect "nsga2: population_size must be >= 2")
  set(xfail_no_evolution_expect "FAILED: ZDT1 solved without any generation")
  foreach(t IN ITEMS xfail_invalid_input_without_stat xfail_no_evolution)
    add_executable(${t} test/xfail/${t}.f90)
    target_link_libraries(${t} PRIVATE nsga2_test_support)
    add_test(NAME ${t} COMMAND ${t})
    set_tests_properties(${t} PROPERTIES PASS_REGULAR_EXPRESSION "${${t}_expect}")
  endforeach()
endif()
