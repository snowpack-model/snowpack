
# add the tests
ADD_TEST(massbalance_check_simple.smoke test_mass_balance.sh)
SET_TESTS_PROPERTIES(massbalance_check_simple.smoke 
                                          PROPERTIES LABELS smoke
                                          FAIL_REGULAR_EXPRESSION "ERROR")

ADD_TEST(energybalance_check_simple.smoke test_energy_balance.sh)
SET_TESTS_PROPERTIES(energybalance_check_simple.smoke 
                                          PROPERTIES LABELS smoke
                                          FAIL_REGULAR_EXPRESSION "ERROR")
