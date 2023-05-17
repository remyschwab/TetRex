cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

declare_datasource (FILE file1.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/file1.fa
                    URL_HASH SHA256=71e7416fe1d7e10c633f253d30e8bbd27e0b8b7dcd4d982e0958c5ddaf19dc27)
declare_datasource (FILE file2.fa
                    URL ${CMAKE_SOURCE_DIR}/test/data/file2.fa
                    URL_HASH SHA256=1fef24d4e2ed643a89e6aaab2e355231569bcbc1002c2ed2715d8becf95ebad3)
declare_datasource (FILE ibf_idx.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/ibf_idx.ibf
                    URL_HASH SHA256=0adaa7b9810c12dbf9b6e88f2eb09cd3cd553a5082480acd3bfda96a865d09dd)
