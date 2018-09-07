An example to run the solver for sphere plane contact problem with the query:

$ bazel build multibody/boussinesq_solver:example_test
$ bazel-bin/multibody/boussinesq_solver/example_tests

By running the commands above, the solver will run for a set of sphere plane contact model setups. The params of the tests are specified in the function GTEST_TEST(ExampleTest, SpherePlane) in the example_test.cc file.





An example to run the solver for sphere plane contact problem assuming only the plane is compliant. Solve the model without involving queries

$ bazel build multibody/boussinesq_solver:sphere_compliance_test
$ bazel-bin/multibody/boussinesq_solver/sphere_compliance_test

The Yound modulus can be adjust in the function SolveCaTimesF() in the shpere_compliance_test.cc


