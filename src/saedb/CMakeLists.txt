project(saedb)

add_subdirectory(io)

file(GLOB SAE_SOURCES
    "*.hpp"
    "*.cpp"
    )

add_library(saedb ${SAE_SOURCES})
target_link_libraries(saedb mgraph)
