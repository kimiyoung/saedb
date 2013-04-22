project(approximate_diameter)

file (GLOB APPROXIMATE_DIAMETER_SOURCES
    "*.hpp"
    "*.cpp"
    )

add_executable(approximate_diameter ${APPROXIMATE_DIAMETER_SOURCES})
target_link_libraries(approximate_diameter saedb)
