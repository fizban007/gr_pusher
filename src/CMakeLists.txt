set(gr_src
  "grid.cpp" "multi_array.cpp" "metrics.cpp" "fields.cpp"
  "memory.cpp" "Logger.cpp" "timer.cpp")
add_library(gr STATIC ${gr_src})

add_executable(pusher "main.cpp")
target_link_libraries(pusher gr)
