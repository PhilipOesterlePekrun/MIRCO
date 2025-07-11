#include "mirco_linearsolver.h"

#include <Teuchos_RCP.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialSpdDenseSolver.hpp>
#include <Teuchos_SerialSymDenseMatrix.hpp>

#include "Timer.hpp"

Teuchos::SerialDenseVector<int, double> MIRCO::LinearSolver::Solve(
    Teuchos::SerialSymDenseMatrix<int, double>& matrix,
    Teuchos::SerialDenseVector<int, double>& vector_b)
{
  linearSystemSizeSummed += vector_b.length();
  numLinearCallsTotal++;
  ScopedTimer timer("LinearSolver::Solve()");

  StandardTimer timerMat("LinearSolver::Solve()__solver.setMatrix()");
  timerMat.start();
  Teuchos::SerialSpdDenseSolver<int, double> solver;
  int err = solver.setMatrix(Teuchos::rcpFromRef(matrix));
  if (err != 0)
  {
    std::cout << "Error setting matrix for linear solver (1)";
  }
  timerMat.stop();

  StandardTimer timer0("LinearSolver::Solve()__solver.SetVectors()");
  timer0.start();
  Teuchos::SerialDenseVector<int, double> vector_x;
  vector_x.size(vector_b.length());

  err = solver.setVectors(Teuchos::rcpFromRef(vector_x), Teuchos::rcpFromRef(vector_b));
  if (err != 0)
  {
    std::cout << "Error setting vectors for linear solver (2)";
  }
  timer0.stop();

  StandardTimer timer1("LinearSolver::Solve()__solver.factorWithEquilibration()");
  timer1.start();
  solver.factorWithEquilibration(true);
  timer1.stop();

  StandardTimer timer2("LinearSolver::Solve()__solver.solve()");
  timer2.start();
  err = solver.solve();
  timer2.stop();
  if (err != 0)
  {
    std::cout << "Error setting up solver (3)";
  }

  return vector_x;
}
