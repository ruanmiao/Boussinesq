
#include "Pendulum.h"  // to get some types
#include <iostream>

using namespace std;
using namespace Drake;

struct OutputTest {
    template <typename ScalarType> using OutputVector = EigenVector<2>::type<ScalarType>;

    template <typename ScalarType>
    OutputVector<ScalarType> output(const ScalarType& t) const;
};

struct OutputTestTwo {
    template <typename ScalarType> using OutputVector = EigenVector<2>::type<ScalarType>;

    OutputVector<double> output(const double& t) const;
};

int main(int argc, char* argv[])
{
  try {
  Eigen::Vector2d x;  x << 0.2, 0.4;

  PendulumState<double> state;
  state.theta = 0.2;
  state.thetadot = .3;

  assert(size(state)==2);

  state = x;
  assert(state.thetadot == 0.4);

  state.theta = 0.5;
  x = toEigen(state);
  assert(x(0) = 0.5);

  {
    Eigen::VectorXd y = toEigen(state);
    assert((x - y).isZero());
  }

  PendulumInput<double> input;
  input.tau = 0.2;

  Eigen::Vector3d abc;  abc << 1,2,3;
  {
    Drake::CombinedVector<double, PendulumState, PendulumInput> test(abc);
    test=2*abc;
    assert(test.first().theta == 2);
    assert(test.first().thetadot == 4);
    assert(test.second().tau == 6);
  }
  {
    Drake::CombinedVectorBuilder<PendulumState,PendulumInput>::type<double> test(abc);
    test=2*abc;
    assert(test.first().theta == 2);
    assert(test.first().thetadot == 4);
    assert(test.second().tau == 6);
  }
  {
    // combining a vector with an unused or empty vector should return the original type
    {
      Drake::CombinedVectorBuilder<PendulumState, NullVector>::type<double> test;
      if (!is_same<PendulumState<double>,decltype(test)>::value)
	throw std::runtime_error("combined vector builder returned " + static_cast<string>(typeid(test).name()));
    }
    {
      Drake::CombinedVectorBuilder<NullVector, PendulumState>::type<double> test;
      if (!is_same<PendulumState<double>,decltype(test)>::value)
	throw std::runtime_error("combined vector builder returned " + static_cast<string>(typeid(test).name()));
    }
  }

  static_assert(Eigen::Matrix<double,2,1>::RowsAtCompileTime == 2,"failed to evaluate RowsAtCompileTime");

/*
  { // test for a polynomial-based algorithm
    static_assert(isPolynomial<Pendulum>,"requires polynomial dynamics");

    PendulumState<Polynomial<double>> x;
    PendulumInput<Polynomial<double>> u;
    auto out = p->dynamicsImplementation(x,u);
  }
*/
  } catch (const exception& e) {
    cout << "ERROR: " << e.what() << endl;
    return -1;
  }

  cout << "all tests passed" << endl;
  return 0;
}
