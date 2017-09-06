
# add the tests
ADD_TEST(res5exp.smoke run_res5exp.sh)
SET_TESTS_PROPERTIES(res5exp.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ")
