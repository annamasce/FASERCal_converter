cmake_minimum_required(VERSION 3.16)
project(FaserRawToRoot)

set(CMAKE_CXX_STANDARD 17)

# ---- ROOT ----
find_package(ROOT REQUIRED COMPONENTS RIO Tree)
include(${ROOT_USE_FILE})

# ---- Dictionary generation ----
set(DICT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/DictTypes.h)
set(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)

ROOT_GENERATE_DICTIONARY(G__NestedVectors
    ${DICT_HEADERS}
    LINKDEF ${LINKDEF}
)

# ---- Add faser-common ----
add_subdirectory(faser-common)

# ---- Executable ----
add_executable(FaserCal_converter src/FaserCal_converter.cpp G__NestedVectors.cxx)

target_include_directories(FaserCal_converter PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link against EventFormats target (propagates includes)
target_link_libraries(FaserCal_converter
    PRIVATE
        EventFormats
        ROOT::RIO
        ROOT::Tree
)