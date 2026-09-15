cmake_minimum_required(VERSION 3.16)
project(addisons_experiments LANGUAGES CXX)

# Consume an already built MFEM; this project is never added by MFEM itself.
find_package(MFEM CONFIG REQUIRED)
set(MFEM_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH
    "MFEM source tree supplying meshes and the unit-test headers")
set(JM_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/johnson_mercier")

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/data")
foreach(mesh ref-triangle.mesh inline-tri.mesh)
  configure_file("${MFEM_SOURCE_DIR}/data/${mesh}"
                 "${CMAKE_CURRENT_BINARY_DIR}/data/${mesh}" COPYONLY)
endforeach()

foreach(driver hx_compare mg_compare)
  add_executable(${driver} "${JM_SOURCE_DIR}/${driver}.cpp")
  target_link_libraries(${driver} PRIVATE mfem)
  target_compile_features(${driver} PRIVATE cxx_std_11)
  target_compile_definitions(${driver} PRIVATE
    JM_DEFAULT_MESH="${CMAKE_CURRENT_BINARY_DIR}/data/ref-triangle.mesh")
endforeach()

include(CTest)
if(BUILD_TESTING)
  add_executable(jm_experiment_tests
    "${JM_SOURCE_DIR}/tests/test_main.cpp"
    "${JM_SOURCE_DIR}/tests/test_jm_hx.cpp")
  target_include_directories(jm_experiment_tests PRIVATE
    "${JM_SOURCE_DIR}" "${MFEM_SOURCE_DIR}/tests/unit")
  target_link_libraries(jm_experiment_tests PRIVATE mfem)
  target_compile_features(jm_experiment_tests PRIVATE cxx_std_11)
  add_test(NAME jm_experiment_tests COMMAND jm_experiment_tests)
  add_test(NAME hx_compare_smoke COMMAND hx_compare -r 1 -repeat 1)
  add_test(NAME mg_compare_smoke COMMAND mg_compare -r 2)
  set_tests_properties(jm_experiment_tests hx_compare_smoke mg_compare_smoke
    PROPERTIES ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
endif()
