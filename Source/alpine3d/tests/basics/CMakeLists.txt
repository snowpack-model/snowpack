
# add the tests
ADD_TEST(basics_simple.smoke run_basics.sh "../simple/output_ref/5_2_dischma.met.bz2")
SET_TESTS_PROPERTIES(basics_simple.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ|fail")
