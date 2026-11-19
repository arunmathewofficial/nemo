///
/// \file nemo-fmt.h
/// \author Jonathan Mackey
/// \date 2026.11.19
///
/// Include fmt's ranges.h (for formatting containers) from the fmt library
/// that spdlog uses.  spdlog may be built with an external fmt library or
/// with its own bundled copy; <spdlog/fmt/ranges.h> picks the right one.
/// spdlog older than 1.10 doesn't have that header, but the packaged
/// versions of those (e.g. Debian 11, Ubuntu 22.04) use an external fmt.
///

#ifndef NEMO_FMT_H
#define NEMO_FMT_H

#if defined(__has_include)
#if __has_include(<spdlog/fmt/ranges.h>)
#include <spdlog/fmt/ranges.h>
#define NEMO_HAVE_SPDLOG_FMT_RANGES
#endif
#endif

#ifndef NEMO_HAVE_SPDLOG_FMT_RANGES
#include <fmt/ranges.h>
#endif

#endif  // NEMO_FMT_H
