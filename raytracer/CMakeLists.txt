function(add_target NAME FILE)
  add_catch(${NAME} ${FILE})

  target_include_directories(${NAME} PRIVATE ../raytracer-geom)
  target_include_directories(${NAME} PRIVATE ../raytracer-reader)

  target_link_libraries(${NAME} PRIVATE ${PNG_LIBRARY})
  target_include_directories(${NAME} PRIVATE ${PNG_INCLUDE_DIRS})
endfunction()

add_target(test_raytracer_asan test_asan.cpp)
add_target(test_raytracer_release test_release.cpp)
