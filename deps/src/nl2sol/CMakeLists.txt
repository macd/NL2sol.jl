#find_program( GCOV_PATH gcov )

if (CMAKE_BUILD_TYPE MATCHES "COVERAGE")
    set (CMAKE_SHARED_LINKER_FLAGS_COVERAGE
         "-lgcov -shared --coverage"
         CACHE STRING "Flags used by the shared library linker for coverage builds."
         FORCE )
endif ()

add_library (nl2sol SHARED assess.f covclc.f dfault.f dotprd.f dupdat.f gqtstp.f 
             imdcon.f itsmry.f linvrt.f litvmu.f livmul.f lmstep.f lsqrt.f lsvmin.f
             ltsqar.f nl2itr.f nl2sno.f nl2sol.f parchk.f qapply.f qrfact.f
             reldst.f rmdcon.f rptmul.f slupdt.f slvmul.f stopx.f today.f
             v2norm.f vaxpy.f vcopy.f vscopy.f)


# "Publish" this library for cmake
target_include_directories (nl2sol PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
