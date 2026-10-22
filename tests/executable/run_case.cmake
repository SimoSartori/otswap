# Runs one case of the otswap executable and checks what it expects.
#
#   cmake -DOTSWAP=<executable> -DCASE=<case .ini> -DWORK=<directory>
#         -DEXAMPLES_DATA=<dir> -DTEST_DATA=<dir> -DPROJECT_VERSION=<version>
#         -P run_case.cmake
#
# WORK is emptied and the case copied into it, with @EXAMPLES_DATA@,
# @TEST_DATA@ and @PROJECT_VERSION@ replaced; the executable then runs in
# WORK, on the copy by its name alone, so every message is the same on
# every machine. The expectations are the case's own lines that start
# with "#! ", comments to the parameter file:
#
#   #! args <arguments>    the command line, instead of the case's name
#   #! exit <status>       the exit status, required
#   #! stdout <line>       stdout, one directive per line; none: empty
#   #! stderr <line>       stderr, whole, one directive per line; none: empty
#   #! wrote <line>        stderr's "otswap: wrote" lines only, the others
#                          not compared (instead of "#! stderr")
#   #! prepare <path>      a file holding "stale", made before the run
#   #! listing <dir>       after the run <dir> holds exactly the files
#   #! file <path>         of these directives under it, at any depth
#   #! overwritten <path>  the prepared file has been written over
#   #! untouched <path>    the prepared file still holds "stale"
#   #! contains <path> <text>  the file holds the text

cmake_minimum_required(VERSION 3.16)

foreach (variable OTSWAP CASE WORK EXAMPLES_DATA TEST_DATA PROJECT_VERSION)
  if (NOT DEFINED ${variable})
    message(FATAL_ERROR "run_case.cmake: ${variable} is not set")
  endif()
endforeach()

get_filename_component(case_name ${CASE} NAME)
file(REMOVE_RECURSE ${WORK})
file(MAKE_DIRECTORY ${WORK})
configure_file(${CASE} ${WORK}/${case_name} @ONLY)

# The directives. Semicolons and square brackets would break CMake's
# lists, so they are swapped for control characters while the file is cut
# into lines, and restored in each line.
string(ASCII 1 semicolon)
string(ASCII 2 open_bracket)
string(ASCII 3 close_bracket)
file(READ ${WORK}/${case_name} text)
string(REPLACE ";" "${semicolon}" text "${text}")
string(REPLACE "[" "${open_bracket}" text "${text}")
string(REPLACE "]" "${close_bracket}" text "${text}")
string(REPLACE "\n" ";" lines "${text}")

set(arguments ${case_name})
unset(expected_exit)
set(expected_stdout "")
set(expected_stderr "")
set(expected_wrote "")
set(compare_wrote OFF)
set(prepared "")
set(listings "")
set(listed "")
set(checks "")

foreach (line IN LISTS lines)
  string(REPLACE "${semicolon}" ";" line "${line}")
  string(REPLACE "${open_bracket}" "[" line "${line}")
  string(REPLACE "${close_bracket}" "]" line "${line}")
  if (NOT line MATCHES "^#! ([a-z]+)( (.*))?$")
    continue()
  endif()
  set(directive ${CMAKE_MATCH_1})
  set(value "${CMAKE_MATCH_3}")
  if (directive STREQUAL "args")
    separate_arguments(arguments UNIX_COMMAND "${value}")
  elseif (directive STREQUAL "exit")
    set(expected_exit ${value})
  elseif (directive STREQUAL "stdout")
    string(APPEND expected_stdout "${value}\n")
  elseif (directive STREQUAL "stderr")
    string(APPEND expected_stderr "${value}\n")
  elseif (directive STREQUAL "wrote")
    set(compare_wrote ON)
    string(APPEND expected_wrote "${value}\n")
  elseif (directive STREQUAL "prepare")
    get_filename_component(parent ${WORK}/${value} DIRECTORY)
    file(MAKE_DIRECTORY ${parent})
    file(WRITE ${WORK}/${value} "stale\n")
  elseif (directive STREQUAL "listing")
    list(APPEND listings ${value})
  elseif (directive STREQUAL "file")
    list(APPEND listed ${value})
  elseif (directive MATCHES "^(overwritten|untouched|contains)$")
    # Checked after the run; the text of "contains" may hold anything
    # but a semicolon.
    list(APPEND checks "${directive} ${value}")
  else()
    message(FATAL_ERROR "${case_name}: unknown directive '${directive}'")
  endif()
endforeach()

if (NOT DEFINED expected_exit)
  message(FATAL_ERROR "${case_name}: no '#! exit' directive")
endif()

execute_process(COMMAND ${OTSWAP} ${arguments}
                WORKING_DIRECTORY ${WORK}
                RESULT_VARIABLE status
                OUTPUT_VARIABLE out
                ERROR_VARIABLE err)

set(failures "")

if (NOT status STREQUAL expected_exit)
  string(APPEND failures "exit status ${status}, expected ${expected_exit}\n")
endif()

if (NOT out STREQUAL expected_stdout)
  string(APPEND failures "stdout:\n${out}-- expected:\n${expected_stdout}--\n")
endif()

if (compare_wrote)
  set(wrote "")
  string(REPLACE ";" "${semicolon}" err_lines "${err}")
  string(REPLACE "[" "${open_bracket}" err_lines "${err_lines}")
  string(REPLACE "]" "${close_bracket}" err_lines "${err_lines}")
  string(REPLACE "\n" ";" err_lines "${err_lines}")
  foreach (line IN LISTS err_lines)
    if (line MATCHES "^otswap: wrote ")
      string(APPEND wrote "${line}\n")
    endif()
  endforeach()
  if (NOT wrote STREQUAL expected_wrote)
    string(APPEND failures "'wrote' lines of stderr:\n${wrote}-- expected:\n${expected_wrote}-- stderr:\n${err}--\n")
  endif()
elseif (NOT err STREQUAL expected_stderr)
  string(APPEND failures "stderr:\n${err}-- expected:\n${expected_stderr}--\n")
endif()

foreach (directory IN LISTS listings)
  file(GLOB_RECURSE found RELATIVE ${WORK} ${WORK}/${directory}/*)
  set(expected "")
  foreach (path IN LISTS listed)
    string(FIND "${path}" "${directory}/" position)
    if (position EQUAL 0)
      list(APPEND expected ${path})
    endif()
  endforeach()
  list(SORT found)
  list(SORT expected)
  if (NOT found STREQUAL expected)
    string(REPLACE ";" "\n" found_text "${found}")
    string(REPLACE ";" "\n" expected_text "${expected}")
    string(APPEND failures "files in ${directory}:\n${found_text}\n-- expected:\n${expected_text}\n--\n")
  endif()
endforeach()

foreach (check IN LISTS checks)
  if (check MATCHES "^(overwritten|untouched) (.*)$")
    set(kind ${CMAKE_MATCH_1})
    set(path ${CMAKE_MATCH_2})
    if (NOT EXISTS ${WORK}/${path})
      string(APPEND failures "${path} does not exist\n")
      continue()
    endif()
    file(READ ${WORK}/${path} contents)
    if (kind STREQUAL "overwritten" AND contents STREQUAL "stale\n")
      string(APPEND failures "${path} still holds what was there before the run\n")
    elseif (kind STREQUAL "untouched" AND NOT contents STREQUAL "stale\n")
      string(APPEND failures "${path} was written over\n")
    endif()
  elseif (check MATCHES "^contains ([^ ]+) (.*)$")
    set(path ${CMAKE_MATCH_1})
    set(wanted "${CMAKE_MATCH_2}")
    if (NOT EXISTS ${WORK}/${path})
      string(APPEND failures "${path} does not exist\n")
      continue()
    endif()
    file(READ ${WORK}/${path} contents)
    string(FIND "${contents}" "${wanted}" position)
    if (position EQUAL -1)
      string(APPEND failures "${path} does not hold '${wanted}'\n")
    endif()
  endif()
endforeach()

if (NOT failures STREQUAL "")
  message(FATAL_ERROR "${case_name}:\n${failures}")
endif()
