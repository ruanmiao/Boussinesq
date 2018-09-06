An example to run the solver:

$ bazel build multibody/boussinesq_solver:example_test
$ bazel-bin/multibody/boussinesq_solver/example_tests

By running the commands above, the solver will run for a set of sphere plane contact model setups. The params of the tests are specified in the function GTEST_TEST(ExampleTest, SpherePlane) in the example_test.cc file.
