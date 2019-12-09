set(BASE_SOURCES
    "${PROJECT_SOURCE_DIR}/src/base/messages.cc"
    "${PROJECT_SOURCE_DIR}/src/base/filewrapper.cc"
    "${PROJECT_SOURCE_DIR}/src/base/filepath.cc"
    "${PROJECT_SOURCE_DIR}/src/base/iowrapper.cc"
    "${PROJECT_SOURCE_DIR}/src/base/exceptions.cc"
    "${PROJECT_SOURCE_DIR}/src/base/tictoc.cc"
#    "${PROJECT_SOURCE_DIR}/src/base/node.cc"
    "${PROJECT_SOURCE_DIR}/src/base/node_list.cc"
#    "${PROJECT_SOURCE_DIR}/src/base/inventoried.cc"
    "${PROJECT_SOURCE_DIR}/src/base/inventory.cc"
    "${PROJECT_SOURCE_DIR}/src/base/stream_func.cc"
    "${PROJECT_SOURCE_DIR}/src/base/tokenizer.cc"
    "${PROJECT_SOURCE_DIR}/src/base/glossary.cc"
    "${PROJECT_SOURCE_DIR}/src/base/property.cc"
    "${PROJECT_SOURCE_DIR}/src/base/property_list.cc"
    "${PROJECT_SOURCE_DIR}/src/base/backtrace.cc"
    "${PROJECT_SOURCE_DIR}/src/base/ansi_colors.cc"
)

add_library(${BASE_LIB_TARGET} STATIC ${BASE_SOURCES})

set_source_files_properties(
    PROPERTIES
    EXTERNAL_OBJECT true
    GENERATED true
)

target_include_directories(${BASE_LIB_TARGET} PUBLIC
    "${BASE_SRC_DIR}"
    "${MATH_SRC_DIR}"
)
