#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double throttle = j[1]["throttle"];
          double steering_angle = j[1]["steering_angle"];

          // Turn waypoints into vehicle coordinates
          for (int i=0; i<ptsx.size(); i++) {
            
            // make waypoints centered at origin of vehicle
            double x_shift = ptsx[i] - px;
            double y_shift = ptsy[i] - py;

            //rotate variables into reference frame of vehicle
            double angle = -psi;

            // caclulate angle beforehand to save computation time
            double cosineAngle = cos(angle);
            double sineAngle = sin(angle);

            ptsx[i] = x_shift * cosineAngle - y_shift * sineAngle;
            ptsy[i] = x_shift * sineAngle + y_shift * cosineAngle;
          }

          // transform waypoints
          Eigen::VectorXd x_glob_wp(ptsx.size());
          Eigen::VectorXd y_glob_wp(ptsy.size());
          for (int i=0; i<ptsx.size(); i++) {
            x_glob_wp[i] = ptsx[i];
            y_glob_wp[i] = ptsy[i];
          }

          // Fit the cubic polynomial to the waypoints to which have been transformed the vehicles frame of reference
          Eigen::VectorXd fit = polyfit(x_glob_wp, y_glob_wp, 3);

          // Now, evaluate the cross track error (cte) relative to this fit cubic polynomial.
          double cte = polyeval(fit, 0);
          std::cout << "Crosstrack Error (CTE) = " << cte << " [m] " << std::endl;

          double epsi = psi - atan(fit[1]);

          // to account for latency, the model is predicted forward by the delay amount so it is updating at the right time
          double latency_dt = 0.1;
          double Lf = 2.67;


          double latency_x = v * latency_dt;
          double latency_y = 0;
          double latency_psi = -(v / Lf) * steering_angle * latency_dt;
          double latency_v = v + throttle * latency_dt;
          double latency_cte = cte + v * sin(epsi) * latency_dt;

          // Compute expected_psi based on the fit.
          double expected_psi = atan(fit[1] + 
                                2.0 * fit[2] * latency_x + 
                                3.0 * fit[3] * latency_x*latency_x);

          // Compute latency error in psi.
          double latency_epsi = psi - expected_psi;
          
          // Compose the state to pass into the solver.
          Eigen::VectorXd state(6);
          state << latency_x, latency_y, latency_psi, latency_v, latency_cte, latency_epsi;

          // Solve for the steering angle and acceleration solution.
          vector<double> solution = mpc.Solve(state, fit);

          // Assign the solution outputs
          double steer_value = solution[0];
          double throttle_value = solution[1];

          // Package the JSON payload for transmission to the simulator.
          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          msgJson["mpc_x"] = mpc.trajectory_x;
          msgJson["mpc_y"] = mpc.trajectory_y;

          //Display the waypoints/reference line
          msgJson["next_x"] = ptsx;
          msgJson["next_y"] = ptsy;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}