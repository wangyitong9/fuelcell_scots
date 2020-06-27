#include "RungeKutta4.hh"
#include <array>
#include <iostream>
#include <cmath>

#define sDIM 2
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double, 2> state_type;
typedef std::array<double, 2> input_type;

const double Np = 3;
const double Ns = 8;
const double Ncells = 3;


/* sampling time */
const double tau = 5;
/* number of intermediate steps in the ode solver */
const int nint = 5;
OdeSolver ode_solver(sDIM, nint, tau);

double growthbound_helper(double i, double Wair, double T, double SOC) {
  double V = 1.256 - 2.26e-4 * T - 2 * 8.3145* T / 96485.0 * std::asinh((i + 0.001) / (16.0 * 10e-5)) - i*0.05 + 0.12*std::log(1 - i / 0.7);
  double P = 3.0 + 10e-3 * 4650.0 * std::pow(Wair, 2) / (60.0 * 125.0 * 3.0) - 150.0 * i*V * 3.0;
  return 1 / std::sqrt(1 - 4 * P*0.18 / 24 / std::pow(1.2*SOC + 3, 2));
}

/* we integrate the system ode by 0.5 sec (the result is stored in x)  */
auto  system_post = [](state_type &x, input_type &u) -> void {
  /* the ode describing the system */
  auto  system_ode = [](state_type &dxdt, const state_type &x, const input_type &u) -> void {
    dxdt[0] = 1 / 100.0 * ((2.26e-4*x[0] + 2 * 8.3145*x[0] / 96485.0 * std::asinh((u[0] + 0.001)*1e5 / 16.0) + 0.05*u[0] - 0.12*std::log(1 - u[0] / 0.7))*u[0] * 150.0 - 1 / 26.0 * x[0] * std::pow(u[1] / Ncells, 1 / 1.35));
    dxdt[1] = (1.2*x[1] + 3 - std::sqrt(std::pow(1.2*x[1] + 3, 2) - 4.0 * 0.18 / Np / Ns * (Ncells + 4.65*std::pow(u[1], 2) / Ncells / 125.0 / 60.0 - u[0] * 150.0 * Ncells * (1.256 - 2.26e-4*x[0] - 2 * 8.3145*x[0] / 96485.0 * std::asinh((u[0] + 0.001) * 100000.0 / 16.0) - 0.05*u[0] + 0.12*std::log(1 - u[0] / 0.7))))) / 3600.0 / 1.875 / 2.0 / 0.18;
    //dxdt[0] = 0.1*x[0] + x[1];
    //dxdt[1] = u[0] + x[1];
    
  };
  ode_solver(system_ode, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) -> void {
  /*
  double a[2][2];
  double x0 = 400, x1 = 0;
  a[0][0] = (2.26e-4 + 2 * 8.3145 / 96485 * std::asinh((u[0] + 0.001) / 16 * 1e5))*u[0] * 150 / 100;// -1 / 26 * pow(u[1] / Ncells, 1 / 1.35) / 100;
  a[0][1] = 0;
  a[1][0] = 1 / 4 / 0.18 / 3600 / 1.875 * pow(pow(1.2*x1 + 3, 2) - 4 * 0.18 / Np / Ns * (Ncells + 4.65e-3*pow(u[1], 2) / Ncells / 125 * 1000 / 60 - u[0] * 150 * Ncells * (1.256 - 2.26e-4*x0 - 2 * 8.3145*x0 / 96485 * std::asinh((u[0] + 0.001)*1e5 / 16) - 0.05*u[0] + 0.12*log(1 - u[0] / 0.7))), -1 / 2) * 4 * 0.18 / Np / Ns * u[0] * 150 * Ncells*(2.26e-4 + 2 * 8.3145 / 96485 * asinh((u[0] + 0.001) / 16 * 1e5));
  //a[1][0] = std::abs(-1 / 4/0.18/3600/1.875 * pow(pow(1.2*x1 + 3, 2) - 4 * 0.18 / Np/Ns * (Ncells + 4.65e-3*pow(u[1], 2) / Ncells / 125 * 1000 / 60 - u[0] * 150 * Ncells * (1.256 - 2.26e-4*x0 - 2 * 8.3145*x0 / 96485 * asinh((u[0] + 0.001)*1e5 / 16) - 0.05*u[0] + 0.12*log(1 - u[0] / 0.7))), -1 / 2) * 4 * 0.18 / Np/Ns * u[0] * 150*Ncells*(-2.26e-4 - 2 * 8.3145 / 96485 * asinh((u[0] + 0.001) / 16 * 1e5)));
  a[1][1] = 1.2 / 3600 / 1.875 / 2 / 0.18;
  r[0] = r[0] + a[0][0] * 0.5 * r[0];
  r[1] = r[1] + a[1][0] * r[0] * 0.5 + a[1][1] * r[1] * 0.5;*/

  auto growth_bound_ode = [](state_type &drdt, const state_type &r, const input_type &u) {
    double a[2][2];
    double x0 = 400, x1 = 0;
    double temp = std::min(growthbound_helper(u[0], u[1], 273, 0), std::min(growthbound_helper(u[0], u[1], 273, 1), std::min(growthbound_helper(u[0], u[1], 400, 0), growthbound_helper(u[0], u[1], 400, 1))));
    a[0][0] = (2.26e-4 + 2 * 8.3145 / 96485.0 * std::asinh((u[0] + 0.001) / 16.0 * 1e5))*u[0] * 150.0 / 100.0 - 1 / 26.0 * std::pow(u[1] / Ncells, 1 / 1.35) / 100.0;
    a[0][1] = 0.0;
    a[1][0] = 1 / 4.0 / 0.18 / 3600.0 / 1.875 * 1 / std::sqrt(std::pow(1.2*x1 + 3, 2) - 4.0 * 0.18 / Np / Ns * (Ncells + 4.65e-3*std::pow(u[1], 2) / Ncells / 125.0 * 1000.0 / 60.0 - u[0] * 150.0 * Ncells * (1.256 - 2.26e-4*x0 - 2 * 8.3145*x0 / 96485.0 * std::asinh((u[0] + 0.001)*1e5 / 16.0) - 0.05*u[0] + 0.12*std::log(1 - u[0] / 0.7)))) * 4.0 * 0.18 / Np / Ns * u[0] * 150.0 * Ncells*(2.26e-4 + 2 * 8.3145 / 96485.0 * std::asinh((u[0] + 0.001) / 16.0 * 1e5));
    //a[1][0] = std::abs(-1 / 4/0.18/3600/1.875 * pow(pow(1.2*x1 + 3, 2) - 4 * 0.18 / Np/Ns * (Ncells + 4.65e-3*pow(u[1], 2) / Ncells / 125 * 1000 / 60 - u[0] * 150 * Ncells * (1.256 - 2.26e-4*x0 - 2 * 8.3145*x0 / 96485 * asinh((u[0] + 0.001)*1e5 / 16) - 0.05*u[0] + 0.12*log(1 - u[0] / 0.7))), -1 / 2) * 4 * 0.18 / Np/Ns * u[0] * 150*Ncells*(-2.26e-4 - 2 * 8.3145 / 96485 * asinh((u[0] + 0.001) / 16 * 1e5)));
    a[1][1] = 1.2 / 3600.0 / 1.875 / 2.0 / 0.18 - temp * 2.0 * 1.2 / (4 * 3600.0 * 1.875*0.18);//- 1 / 2/3600/1.875/2/0.18 * pow(pow(1.2*r[1] + 3, 2) - 4 * 0.18 / Np/Ns * (Ncells + 4.65e-3*pow(u[1], 2) / Ncells / 125 * 1000 / 60 - u[0] * 150 * Ncells * (1.256 - 2.26e-4*r[0] - 2 * 8.3145*r[0] / 96485 * asinh((u[0] + 0.001)*1e5 / 16) - 0.05*u[0] + 0.012*log(1 - u[0] / 0.7))), -1 / 2) * 2 * 1.2*(1.2*r[1] + 3);
    drdt[0] = a[0][0] * r[0] + a[0][1] * r[1];
    drdt[1] = a[1][0] * r[0] + a[1][1] * r[1];
    //drdt[0] = 0.1*r[0] + r[1];
    //drdt[1] = r[1];
  };
  ode_solver(growth_bound_ode, r, u);
};

int main() {
  state_type x, r;
  input_type u;
  x = { 300.00,0.15 };
  u = { 0, 800 };
  r = { 4,0.06 };
  system_post(x, u);
  radius_post(r, u);
  std::cout << x[0] << " " << x[1] << "\n";
  std::cout << r[0] << " " << r[1] << "\n";
}