cmake_minimum_required(VERSION 3.10)
project(ultra_precision_pendulum VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /fp:fast")
endif()

# Add the executable
add_executable(pendulum src/UltraPrecisionPendulum.cpp)

# Link math library if needed
if(UNIX)
    target_link_libraries(pendulum m)
endif()

# Install target
install(TARGETS pendulum DESTINATION bin)

# Add testing
enable_testing()
