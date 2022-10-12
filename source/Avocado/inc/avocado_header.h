#pragma once

#define AVOCADO_MAJOR_VERSION 1
#define AVOCADO_MINOR_VERSION 0
#define AVOCADO_BUILD_VERSION 1
#define AVOCADO_VERSION                                       \
  "Avocado v" + std::to_string(AVOCADO_MAJOR_VERSION) + "." + \
      std::to_string(AVOCADO_MINOR_VERSION) + "." +           \
      std::to_string(AVOCADO_BUILD_VERSION)
#define AVOCADO_LAST_UPDATE "Oct. 12 2022"

#ifdef _WIN32
#define AVOCADO_USE_WINDOWS
#elif __linux__
#define AVOCADO_USE_LINUX
#endif

#ifndef AVOCADO_USE_LINUX
#ifndef AVOCADO_USE_WINDOWS
#define AVOCADO_USE_WINDOWS
#endif
#endif

#ifdef AVOCADO_USE_WINDOWS
#include <Windows.h>
#include <direct.h>
#include <io.h>
#undef GetCurrentTime
#undef min
#undef max
#elif defined(AVOCADO_USE_LINUX)
#include <sys/stat.h>
#include <unistd.h>
#endif