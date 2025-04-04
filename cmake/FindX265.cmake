#
# Copyright (C) 2024 The Android Open Source Project
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.
#

#
# Finds the x265 Encoder library. This module defines:
#
#  X265_FOUND            - True if x265 encoder is found, False otherwise
#  X265_LIBRARIES        - x265 encoder library
#  X265_INCLUDE_DIRS     - Include dir
#

find_path(X265_INCLUDE_DIRS x265.h)

find_library(X265_LIBRARIES NAMES libx265 x265)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(X265 DEFAULT_MSG X265_INCLUDE_DIRS X265_LIBRARIES)

mark_as_advanced(X265_INCLUDE_DIRS X265_LIBRARIES)
