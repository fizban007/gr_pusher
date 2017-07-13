add_custom_target(check)

set(tests_src "test.cpp")
add_executable(tests EXCLUDE_FROM_ALL ${tests_src})
target_link_libraries(tests gr)
add_custom_target(run_tests
  COMMAND "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests"
  DEPENDS tests
  WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
add_dependencies(check run_tests)
