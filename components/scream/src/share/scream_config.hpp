#ifndef SCREAM_CONFIG_HPP
#define SCREAM_CONFIG_HPP

#include <string>

// Include this file, not any lower-level configuration file such as that
// generated by CMake from scream_config.h.in. The intent is to funnel all
// configuration decisions through this header.

#ifdef SCREAM_CONFIG_IS_CMAKE
# include "scream_config.h"
#else
// Purposely error out.
"A non-cmake build of scream is currently not supported."
#endif

#ifdef SCREAM_FORCE_BUILD_FAIL
choke_build();
#endif

namespace scream {

std::string scream_config_string();

} // namespace scream

#endif // SCREAM_CONFIG_HPP
