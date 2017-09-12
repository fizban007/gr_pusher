set(gr_src
  "grid.cpp" "multi_array.cpp" "metrics.cpp" "fields.cpp"
  "memory.cpp" "Logger.cpp" "timer.cpp" "solve.cpp")
add_library(gr STATIC ${gr_src})

add_executable(pusher-schwarzschild "pusher-schwarzschild.cpp")
target_link_libraries(pusher-schwarzschild gr)

add_executable(pusher-cks "pusher-cks.cpp")
target_link_libraries(pusher-cks gr)

add_executable(pusher-boyer "pusher-boyer.cpp")
target_link_libraries(pusher-boyer gr)

# add_executable(test "test.cpp")
