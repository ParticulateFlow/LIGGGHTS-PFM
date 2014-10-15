find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

FILE(GLOB TEST_SOURCES *.cpp)
add_executable(runtests ${TEST_SOURCES})
target_link_libraries(runtests ${GTEST_BOTH_LIBRARIES} liggghts)


add_custom_command(TARGET runtests POST_BUILD COMMAND cp -R "${PROJECT_SOURCE_DIR}/tests/scripts/" $<TARGET_FILE_DIR:runtests>)

add_test(liggghts_tests runtests)
#add_test(NAME parallel_tests_4
#  COMMAND mpirun -np 4 $<TARGET_FILE:runtests>
#  WORKING_DIRECTORY $<TARGET_FILE_DIR:runtests>)
#add_test(NAME parallel_tests_8
#  COMMAND mpirun -np 8 $<TARGET_FILE:runtests>
#  WORKING_DIRECTORY $<TARGET_FILE_DIR:runtests>)
#add_test(NAME parallel_tests_32
#  COMMAND mpirun -np 32 $<TARGET_FILE:runtests>
#  WORKING_DIRECTORY $<TARGET_FILE_DIR:runtests>)
