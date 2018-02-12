
# add the tests
ADD_TEST(simple.smoke run_simple.sh)
SET_TESTS_PROPERTIES(simple.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ|fail")
