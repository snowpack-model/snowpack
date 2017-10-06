
# add the tests
ADD_TEST(basics_1exp.smoke run_basics.sh "../res1exp/output_ref/MST96_res.met.bz2")
SET_TESTS_PROPERTIES(basics_1exp.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ|fail")

ADD_TEST(basics_5exp_N.smoke run_basics.sh "../res5exp/output_ref/MST961_res4.met.bz2")
SET_TESTS_PROPERTIES(basics_5exp_N.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ|fail")
ADD_TEST(basics_5exp_S.smoke run_basics.sh "../res5exp/output_ref/MST963_res4.met.bz2")
SET_TESTS_PROPERTIES(basics_5exp_S.smoke
                     PROPERTIES LABELS smoke
                     FAIL_REGULAR_EXPRESSION "error|differ|fail")