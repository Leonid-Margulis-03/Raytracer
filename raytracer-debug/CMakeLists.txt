function(add_target NAME FILE)
  add_catch(${NAME} ${FILE})

  target_include_directories(${NAME} PRIVATE ../raytracer-geom)
  target_include_directories(${NAME} PRIVATE ../raytracer-reader)
  target_include_directories(${NAME} PRIVATE ../raytracer)

  target_link_libraries(${NAME} PRIVATE ${PNG_LIBRARY} ${JPEG_LIBRARIES})
  target_include_directories(${NAME} PRIVATE ${PNG_INCLUDE_DIRS} ${JPEG_INCLUDE_DIRS})
endfunction()

add_target(test_raytracer_debug_asan test_asan.cpp)
add_target(test_raytracer_debug_release test_release.cpp)
