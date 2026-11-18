cmake_minimum_required(VERSION 3.16)
project(orpheic-decoders LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

# The comments in the sources are Greek, in UTF-8: MSVC must be told, or on some Windows code pages it misreads them.
if(MSVC)
  add_compile_options(/utf-8)
else()
  # Files of 2 GB and more on 32-bit systems too (a whole SACD as one .dff is several GB).
  add_compile_definitions(_FILE_OFFSET_BITS=64)
endif()

# Three decoders, each its own library with no dependencies beyond the C++17 standard library: take one, two or all three.
# Each picks its processor-specific path (SSE2, AVX2 with FMA, NEON) at run time; nothing has to be set here.
add_library(orpheic_mp3 STATIC src/mp3decoder.cpp src/mp3huff.cpp src/mp3hybrid.cpp src/mp3synth.cpp)
add_library(orpheic_flac STATIC src/flacdecoder.cpp)
add_library(orpheic_dst STATIC src/dstdecoder.cpp)
foreach(lib orpheic_mp3 orpheic_flac orpheic_dst)
  target_include_directories(${lib} PUBLIC src)
endforeach()
# The same bits on every system: GCC may fuse a multiply and an add into one instruction where the processor allows it, Clang
# and MSVC do not. In the MP3 decoder every operation is done as written (the explicit FMA of the AVX2 code stays as it is).
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(orpheic_mp3 PRIVATE -ffp-contract=off)
endif()

# All three together, for whoever wants one name to link.
add_library(orpheic_decoders INTERFACE)
target_link_libraries(orpheic_decoders INTERFACE orpheic_mp3 orpheic_flac orpheic_dst)

option(ORPHEIC_DECODERS_TOOLS "Build the small command-line tools under tools/" ON)
if(ORPHEIC_DECODERS_TOOLS)
  file(GLOB tool_sources CONFIGURE_DEPENDS tools/*.cpp)
  list(FILTER tool_sources EXCLUDE REGEX "/(speedtest|refbench)\\.cpp$")   # need FFmpeg (and refbench libFLAC, libmpg123): below
  foreach(src ${tool_sources})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE orpheic_decoders)
    if(MSVC)
      target_compile_definitions(${name} PRIVATE _CRT_SECURE_NO_WARNINGS)   # the tools use plain fopen
    endif()
    if(WIN32)
      target_sources(${name} PRIVATE tools/utf8.manifest)   # file names in UTF-8 on Windows
    endif()
  endforeach()
  # speedtest compares the decoders with the FFmpeg installed on this machine, so it is built only when FFmpeg's
  # development files are there (pkg-config: libavformat, libavcodec, libavutil).
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    # FFmpeg 5.1 or newer (the channel layout API speedtest uses).
    pkg_check_modules(FFMPEG IMPORTED_TARGET libavformat>=59.27 libavcodec>=59.37 libavutil>=57.28)
  endif()
  if(FFMPEG_FOUND)
    add_executable(speedtest tools/speedtest.cpp)
    target_link_libraries(speedtest PRIVATE orpheic_decoders PkgConfig::FFMPEG)
    if(WIN32)
      target_sources(speedtest PRIVATE tools/utf8.manifest)
    endif()
  else()
    message(STATUS "FFmpeg 5.1+ development files not found: speedtest is not built (the other tools are)")
  endif()
  # refbench compares the FLAC decoder with libFLAC and the MP3 decoder with libmpg123, through their public APIs only (and
  # reads MP3 packets with libavformat): built only when all three are there (pkg-config: flac, libmpg123, FFmpeg as above).
  if(PKG_CONFIG_FOUND AND FFMPEG_FOUND)
    pkg_check_modules(REFLIBS IMPORTED_TARGET flac libmpg123)
  endif()
  if(REFLIBS_FOUND)
    add_executable(refbench tools/refbench.cpp)
    target_link_libraries(refbench PRIVATE orpheic_decoders PkgConfig::REFLIBS PkgConfig::FFMPEG)
    if(WIN32)
      target_sources(refbench PRIVATE tools/utf8.manifest)
    endif()
  else()
    message(STATUS "libFLAC, libmpg123 or FFmpeg development files not found: refbench is not built")
  endif()
endif()
