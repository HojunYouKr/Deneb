#pragma once

#define DENEB_MAJOR_VERSION 1
#define DENEB_MINOR_VERSION 0
#define DENEB_BUILD_VERSION 1
#define DENEB_VERSION                                     \
  "Deneb v" + std::to_string(DENEB_MAJOR_VERSION) + "." + \
      std::to_string(DENEB_MINOR_VERSION) + "." +         \
      std::to_string(DENEB_BUILD_VERSION)
#define DENEB_LAST_UPDATE "June 10 2022"
