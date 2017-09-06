
# add the tests
ADD_TEST(res1exp.smoke run_res1exp.sh)
SET_TESTS_PROPERTIES(res1exp.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ")
