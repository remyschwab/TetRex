cmake_minimum_required (VERSION 3.8)

add_custom_target (cli_test ALL)

macro (add_cli_test test_filename)
    add_app_test (${test_filename})

    seqan3_test_component (target "${test_filename}" TARGET_NAME)
    add_dependencies (cli_test ${target})
    add_dependencies (${target} "${PROJECT_NAME}")
endmacro ()

add_cli_test (kbioreg_test.cpp)
target_use_datasources (kbioreg_test FILES ibf_idx.ibf file1.fa file2.fa)
