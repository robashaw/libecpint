#
#    Copyright (c) 2020 Robert Shaw
#		 This file is a part of Libecpint.
#
#    Permission is hereby granted, free of charge, to any person obtaining
#    a copy of this software and associated documentation files (the
#    "Software"), to deal in the Software without restriction, including
#    without limitation the rights to use, copy, modify, merge, publish,
#    distribute, sublicense, and/or sell copies of the Software, and to
#    permit persons to whom the Software is furnished to do so, subject to
#    the following conditions:
#
#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


find_program(CLANG_FORMAT_EXECUTABLE
             NAMES clang-format
                   clang-format-11
                   clang-format-10
                   clang-format-9
                   clang-format-8
                   clang-format-7
                   clang-format-6.0
                   clang-format-5.0
                   clang-format-4.0
                   clang-format-3.9
                   clang-format-3.8
                   clang-format-3.7
                   clang-format-3.6
                   clang-format-3.5
                   clang-format-3.4
                   clang-format-3.3
             DOC "clang-format executable")
mark_as_advanced(CLANG_FORMAT_EXECUTABLE)

# Extract version from command "clang-format -version"
if(CLANG_FORMAT_EXECUTABLE)
  execute_process(COMMAND ${CLANG_FORMAT_EXECUTABLE} -version
                  OUTPUT_VARIABLE clang_format_version
                  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(clang_format_version MATCHES "^.*clang-format version .*")
    # clang_format_version sample: "clang-format version 3.9.1-4ubuntu3~16.04.1
    # (tags/RELEASE_391/rc2)"
    string(REGEX
           REPLACE "^.*clang-format version ([.0-9]+).*"
                   "\\1"
                   CLANG_FORMAT_VERSION
                   "${clang_format_version}")
    # CLANG_FORMAT_VERSION sample: "3.9.1"
  else()
    set(CLANG_FORMAT_VERSION 0.0)
  endif()
else()
  set(CLANG_FORMAT_VERSION 0.0)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CLANG_FORMAT_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CLANG_FORMAT REQUIRED_VARS CLANG_FORMAT_EXECUTABLE VERSION_VAR CLANG_FORMAT_VERSION)
