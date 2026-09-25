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

# Keep the existing drivers buildable without these optional dependencies.
option(JM_BUILD_SPECTRUM "Build the LAPACK/ARPACK HX spectral experiment" OFF)
if(JM_BUILD_SPECTRUM)
  find_package(LAPACK REQUIRED)
  find_library(ARPACK_LIBRARY NAMES arpack libarpack.so.2)
  if(NOT ARPACK_LIBRARY)
    message(FATAL_ERROR "JM_BUILD_SPECTRUM requires an LP64 ARPACK-ng library")
  endif()
  add_executable(hx_spectrum "${JM_SOURCE_DIR}/hx_spectrum.cpp")
  target_link_libraries(hx_spectrum PRIVATE mfem ${LAPACK_LIBRARIES} ${ARPACK_LIBRARY})
  target_compile_features(hx_spectrum PRIVATE cxx_std_11)
  if(BUILD_TESTING)
    add_test(NAME hx_spectrum_dense_check COMMAND hx_spectrum -r 0 -select base)
    add_test(NAME hx_spectrum_projection_sweep_check COMMAND hx_spectrum -r 0 -select mult_SHCHS)
    add_test(NAME hx_spectrum_middle_check COMMAND hx_spectrum -r 0 -select middle -tol 1e-8 -max-it 100)
    add_test(NAME hx_spectrum_hybrid_weight_check COMMAND hx_spectrum -r 0 -select custom -order "H(C+S)H" -ws 1.3 -wh .7 -wc .6)
    set_tests_properties(hx_spectrum_dense_check hx_spectrum_projection_sweep_check hx_spectrum_middle_check hx_spectrum_hybrid_weight_check PROPERTIES ENVIRONMENT "OPENBLAS_NUM_THREADS=1")
  endif()
endif()
