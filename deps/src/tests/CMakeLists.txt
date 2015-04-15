set (TEST_SRCS nlmain.f ivvset.f madj.f madr.f nltest.f testj.f testr.f xinit.f)

add_executable (nlmain ${TEST_SRCS})
target_link_libraries (nlmain LINK_PUBLIC nl2sol)

# OK, I tried to add the test definition to the top level CMakeLists.txt file
# and it did not work, so be aware of that.
add_test(nlmain nlmain)