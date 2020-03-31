if(TEST_COVERAGE)
    if(NOT ${CMAKE_BUILD_TYPE} STREQUAL Debug)
        message("Setting build type to 'Debug' for test coverage!")
        set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    endif()

    find_program( GCOV_PATH gcov )
    if(NOT GCOV_PATH)
        message(FATAL_ERROR "gcov not found! Aborting...")
    endif()

    find_program( GCOVR gcovr )
    if(NOT GCOVR)
        message(FATAL_ERROR "gcovr not found! Aborting... Please run `pip install gcovr`.")
    endif()

    set(HTML_DIR "${PROJECT_BINARY_DIR}/coverage")
    file(MAKE_DIRECTORY ${HTML_DIR})

    add_custom_command(OUTPUT _run_gcovr_parser
      POST_BUILD
      COMMAND ${GCOVR}
              --root=${PROJECT_SOURCE_DIR}
              --object-dir=${PROJECT_BINARY_DIR}
              --filter="${PROJECT_SOURCE_DIR}/core/"
              --filter="${PROJECT_SOURCE_DIR}/fastslam1/"
              --exclude=".+_test.cpp"
              --print-summary

      COMMAND ${GCOVR}
              --root=${PROJECT_SOURCE_DIR}
              --object-dir=${PROJECT_BINARY_DIR}
              --filter="${PROJECT_SOURCE_DIR}/core/"
              --filter="${PROJECT_SOURCE_DIR}/fastslam1/"
              --exclude=".+_test.cpp"
              --html-details=${HTML_DIR}/coverage.html
      COMMENT "Find coverage html output at `${HTML_DIR}`"
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    add_custom_target (coverage DEPENDS _run_gcovr_parser)
endif()
