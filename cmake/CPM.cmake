# SPDX-License-Identifier: MIT
#
# SPDX-FileCopyrightText: Copyright (c) 2019-2023 Lars Melchior and contributors

# Option 1: Use tag (preferred if set)
# Option 2: Use version number (fallback)
if(DEFINED CPM_DOWNLOAD_TAG)
  set(CPM_DOWNLOAD_TAG ${CPM_DOWNLOAD_TAG})
  set(CPM_DOWNLOAD_VERSION "")  # Clear version when tag is used
else()
  set(CPM_DOWNLOAD_VERSION 0.42.0)
  set(CPM_HASH_SUM "2020b4fc42dba44817983e06342e682ecfc3d2f484a581f11cc5731fbe4dce8a")
endif()

# Determine the download URL based on whether tag or version is specified
if(DEFINED CPM_DOWNLOAD_TAG AND NOT CPM_DOWNLOAD_TAG STREQUAL "")
  # Use tag - note: tags might not follow the 'v' prefix pattern
  if(CPM_DOWNLOAD_TAG MATCHES "^v?[0-9]")
    # If tag looks like a version (with optional v prefix), construct URL accordingly
    string(REGEX REPLACE "^v" "" CPM_TAG_CLEAN "${CPM_DOWNLOAD_TAG}")
    set(CPM_DOWNLOAD_URL "https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_TAG_CLEAN}/CPM.cmake")
  else()
    # For non-version tags (like 'main', 'develop', etc.), use the tag directly
    set(CPM_DOWNLOAD_URL "https://raw.githubusercontent.com/cpm-cmake/CPM.cmake/${CPM_DOWNLOAD_TAG}/cmake/CPM.cmake")
  endif()
  set(CPM_DOWNLOAD_FILENAME "CPM_${CPM_DOWNLOAD_TAG}.cmake")
else()
  # Use version number (original behavior)
  set(CPM_DOWNLOAD_URL "https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake")
  set(CPM_DOWNLOAD_FILENAME "CPM_${CPM_DOWNLOAD_VERSION}.cmake")
  set(CPM_USE_HASH_CHECK TRUE)
endif()

if(CPM_SOURCE_CACHE)
  set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/${CPM_DOWNLOAD_FILENAME}")
elseif(DEFINED ENV{CPM_SOURCE_CACHE})
  set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/${CPM_DOWNLOAD_FILENAME}")
else()
  set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/${CPM_DOWNLOAD_FILENAME}")
endif()

# Expand relative path. This is important if the provided path contains a tilde (~)
get_filename_component(CPM_DOWNLOAD_LOCATION ${CPM_DOWNLOAD_LOCATION} ABSOLUTE)

# Download with or without hash checking
if(CPM_USE_HASH_CHECK)
  file(DOWNLOAD
       ${CPM_DOWNLOAD_URL}
       ${CPM_DOWNLOAD_LOCATION} EXPECTED_HASH SHA256=${CPM_HASH_SUM}
  )
else()
  file(DOWNLOAD
       ${CPM_DOWNLOAD_URL}
       ${CPM_DOWNLOAD_LOCATION}
  )
endif()

include(${CPM_DOWNLOAD_LOCATION})