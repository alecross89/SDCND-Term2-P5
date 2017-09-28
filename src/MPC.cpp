#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

/*
chose N=10 because when the numbers became larger it became more computationally intensive for the system beacuse
it has to go further out in time. Also, predicting too far out in advance is pointless because the prediction would extend
past the waypoints and so predicting about a second or less into the future was optimal.  N=10 also gives the car enough
time to come back toward the center line.

I chose dt = 0.1, because any order of magnitudes smaller and the program would have to run way too many steps which makes 
it computationally intensive, but if it is too large, then you're not processing the data fast enough to make quality 
predictions
*/
size_t N = 10;
double dt = 0.1;


// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;


double ref_cte = 0; // want car to sit on the line
double ref_epsi = 0; // want car to align on the line
double ref_v = 50; // speed
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

    fg[0] = 0;

    // Impose cost based on deviations of reference trajectory.
    for (int i=0; i < N; i++) {

      // Penalize cross track error 
      fg[0] += 500 * CppAD::pow(vars[cte_start  + i] - ref_cte, 2);

      // Penalize deviations from reference heading angle.
      fg[0] += 6000 * CppAD::pow(vars[epsi_start + i] - ref_epsi, 2);

      // Penalize deviations from reference speed.
      fg[0] += CppAD::pow(vars[v_start + i] - ref_v, 2);
    }

    // Penalize control energy expenditure
    for (int i=0; i < N-1; i++) {

      // Penalize non-zero steering angles 
      fg[0] += 15000 * CppAD::pow(vars[delta_start + i], 2);

      // Penalize non-zero longitudinal acceleration commands (braking or acceleration)
      fg[0] += 10* CppAD::pow(vars[a_start + i], 2);
    }

    // Penalize large changes in controller commands across time steps.
    for (int i=0; i < N-2; i++) {
      // Penalize changes in steering angle from one timestep to the next
      fg[0] += 5000 * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);

      // Penalize changes in the acceleration command from one timestep to the next.
      fg[0] += 10000 * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    }

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints
    for (int i = 0; i < N - 1; i++) {

      AD<double> x1 = vars[x_start + i + 1];
      AD<double> y1 = vars[y_start + i + 1];
      AD<double> psi1 = vars[psi_start + i + 1];
      AD<double> v1 = vars[v_start + i + 1];
      AD<double> cte1 = vars[cte_start + i + 1];
      AD<double> epsi1 = vars[epsi_start + i + 1];

      // The state at time t.
      AD<double> x0 = vars[x_start + i];
      AD<double> y0 = vars[y_start + i];
      AD<double> psi0 = vars[psi_start + i];
      AD<double> v0 = vars[v_start + i];
      AD<double> cte0 = vars[cte_start + i];
      AD<double> epsi0 = vars[epsi_start + i];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + i];
      AD<double> a0 = vars[a_start + i];

      // Since coeffs is a cubic polynomial (for the project), this 
      // code needed to be changed to reflect the higher order. 
      AD<double> x2 = x0 * x0;
      AD<double> x3 = x2 * x0;
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x2+ coeffs[3] * x3;

      // atan(Derivative of f0)
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x2);
      //

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // NOTE: The use of `AD<double>` and use of `CppAD`!
      // This is also CppAD can compute derivatives and pass
      // these to the solver.
      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt
      fg[2 + x_start + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[2 + y_start + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[2 + psi_start + i] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[2 + v_start + i] = v1 - (v0 + a0 * dt);
      fg[2 + cte_start + i] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[2 + epsi_start + i] =
          epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }

  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // number of model variables
  size_t n_vars = N * 6 + (N - 1) * 2;
  // number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Initialize the first elements of vars with the initial state.
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];


  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Lower and upper limits for the constraints should be 0 besides initial state.
  for (int i=0; i < delta_start; i++) {
    vars_lowerbound[i] = -1e6;
    vars_upperbound[i] = 1e6;
  }

  // Set bounds for steering
  // 25 degrees converted to radians.
  for (int i=delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436331944;
    vars_upperbound[i] = 0.436331944;
  }

  // Set bounds for acceleration
  for (int i=a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Populate initial constraint states to fix the initial state by the optimizer.
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;
  
  // - Upper bounds.
  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // Load solution trajectory into trajectory_x and trajectory_y for visualization.
  trajectory_x.clear();
  trajectory_y.clear();
  for (int i=0; i<N-1; i++) {
    double x = solution.x[x_start + i + 1];
    double y = solution.x[y_start + i + 1];
    trajectory_x.push_back(x);
    trajectory_y.push_back(y);
  }
  double steering_cmd = -solution.x[delta_start];
  double accel_cmd = solution.x[a_start];
  return {steering_cmd, accel_cmd};
}
