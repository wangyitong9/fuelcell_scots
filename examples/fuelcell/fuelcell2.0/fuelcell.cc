/*
* fuelcell.cc
*
*
*
*/

#include <array>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicModelGrowthBound.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double, 3> state_type;
typedef std::array<double, 2> input_type;

const double Np = 3;
const double Ns = 8;
const double Ncells = 3;
const double Pload = 40;
const double w2 = std::sqrt(4 * 0.18*Pload / 24.0) / 2.0 / 0.18 / 3600.0 / 1.875;

/* sampling time */
const double tau = 10;
/* number of intermediate steps in the ode solver */
const int nint = 5;
OdeSolver ode_solver(sDIM, nint, tau);

double growthbound_helper(double i, double Wair, double T, double SOC, double Pm) {
  double V = 1.256 - 2.26e-4 * T - 2 * 8.3145* T / 96485.0 * std::asinh((i + 0.001) / (16.0 * 10e-5)) - i*0.05 + 0.12*std::log(1 - i / 0.7);
  double P = Pm + 3.0 + 10e-3 * 4650.0 * std::pow(Wair, 2) / (60.0 * 125.0 * 3.0) - 150.0 * i*V * 3.0;
  return 1 / std::sqrt(1 - 4 * P*0.18 / 24 / std::pow(1.2*SOC + 3, 2));
}

/* we integrate the system ode by 0.5 sec (the result is stored in x)  */
auto  system_post = [](state_type &x, input_type &u) -> void {
  /* the ode describing the system */
  auto  system_ode = [](state_type &dxdt, const state_type &x, const input_type &u) -> void {
    dxdt[0] = 1 / 100.0 * ((2.26e-4*x[0] + 2 * 8.3145*x[0] / 96485.0 * std::asinh((u[0] + 0.001)*1e5 / 16.0) + 0.05*u[0] - 0.12*std::log(1 - u[0] / 0.7))*u[0] * 150.0 - 1 / 26.0 * x[0] * std::pow(u[1] / Ncells, 1 / 1.35));
    dxdt[1] = (1.2*x[1] + 3 - std::sqrt(std::pow(1.2*x[1] + 3, 2) - 4.0 * 0.18 / Np / Ns * (x[2]+Ncells + 4.65*std::pow(u[1], 2) / Ncells / 125.0 / 60.0 - u[0] * 150.0 * Ncells * (1.256 - 2.26e-4*x[0] - 2 * 8.3145*x[0] / 96485.0 * std::asinh((u[0] + 0.001) * 100000.0 / 16.0) - 0.05*u[0] + 0.12*std::log(1 - u[0] / 0.7))))) / 3600.0 / 1.875 / 2.0 / 0.18;
    dxdt[2] = 0;
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
    double a[3][3];
    double x0 = 400, x1 = 0;
    double temp = std::min(growthbound_helper(u[0], u[1], 273, 0, 0), std::min(growthbound_helper(u[0], u[1], 273, 1, 0), std::min(growthbound_helper(u[0], u[1], 300, 0, 0), growthbound_helper(u[0], u[1], 300, 1, 0))));
    double temp2 = std::max(growthbound_helper(u[0], u[1], 273, 0, 285), std::max(growthbound_helper(u[0], u[1], 273, 1, 285), std::max(growthbound_helper(u[0], u[1], 300, 0, 285), growthbound_helper(u[0], u[1], 300, 1, 285))));
    a[0][0] = (2.26e-4 + 2 * 8.3145 / 96485.0 * std::asinh((u[0] + 0.001) / 16.0 * 1e5))*u[0] * 150.0 / 100.0 -1 / 26.0 * std::pow(u[1] / Ncells, 1 / 1.35) / 100.0;
    a[0][1] = 0.0;
    a[0][2] = 0.0;
    a[1][0] = 1 / 4.0 / 0.18 / 3600.0 / 1.875 * 1/std::sqrt(std::pow(1.2*x1 + 3, 2) - 4.0 * 0.18 / Np / Ns * (Ncells + 4.65e-3*std::pow(u[1], 2) / Ncells / 125.0 * 1000.0 / 60.0 - u[0] * 150.0 * Ncells * (1.256 - 2.26e-4*x0 - 2 * 8.3145*x0 / 96485.0 * std::asinh((u[0] + 0.001)*1e5 / 16.0) - 0.05*u[0] + 0.12*std::log(1 - u[0] / 0.7)))) * 4.0 * 0.18 / Np / Ns * u[0] * 150.0 * Ncells*(2.26e-4 + 2 * 8.3145 / 96485.0 * std::asinh((u[0] + 0.001) / 16.0 * 1e5));
    //a[1][0] = std::abs(-1 / 4/0.18/3600/1.875 * pow(pow(1.2*x1 + 3, 2) - 4 * 0.18 / Np/Ns * (Ncells + 4.65e-3*pow(u[1], 2) / Ncells / 125 * 1000 / 60 - u[0] * 150 * Ncells * (1.256 - 2.26e-4*x0 - 2 * 8.3145*x0 / 96485 * asinh((u[0] + 0.001)*1e5 / 16) - 0.05*u[0] + 0.12*log(1 - u[0] / 0.7))), -1 / 2) * 4 * 0.18 / Np/Ns * u[0] * 150*Ncells*(-2.26e-4 - 2 * 8.3145 / 96485 * asinh((u[0] + 0.001) / 16 * 1e5)));
    a[1][1] = 1.2 / 3600.0 / 1.875 / 2.0 / 0.18 - temp * 2.0 * 1.2 / (4 * 3600.0 * 1.875*0.18);//- 1 / 2/3600/1.875/2/0.18 * pow(pow(1.2*r[1] + 3, 2) - 4 * 0.18 / Np/Ns * (Ncells + 4.65e-3*pow(u[1], 2) / Ncells / 125 * 1000 / 60 - u[0] * 150 * Ncells * (1.256 - 2.26e-4*r[0] - 2 * 8.3145*r[0] / 96485 * asinh((u[0] + 0.001)*1e5 / 16) - 0.05*u[0] + 0.012*log(1 - u[0] / 0.7))), -1 / 2) * 2 * 1.2*(1.2*r[1] + 3);
    a[1][2] = temp2 / Np / Ns / 3600 / 1.875;
    drdt[0] = a[0][0] * r[0] + a[0][1] * r[1] + a[0][2] * r[2];
    drdt[1] = a[1][0] * r[0] + a[1][1] * r[1] + a[1][2] * r[2];
    drdt[2] = 0;
    //drdt[0] = 0.1*r[0] + r[1];
    //drdt[1] = r[1];
  };
  ode_solver(growth_bound_ode, r, u);
};


/****************************************************************************/
/* main computation */
/****************************************************************************/
int main() {
  /* to measure time */
  TicToc tt;
  /* there is one unique manager to organize the bdd variables */
  Cudd mgr;

  /****************************************************************************/
  /* construct SymbolicSet for the state space */
  /****************************************************************************/
  /* setup the workspace of the synthesis problem and the uniform grid */
  /* lower bounds of the hyper rectangle */
  double lb[sDIM] = { 273,0,0 };
  /* upper bounds of the hyper rectangle */
  double ub[sDIM] = { 300,1,285 };
  /* grid node distance diameter */
  double eta[sDIM] = { 1,0.0001,10 };
  scots::SymbolicSet ss(mgr, sDIM, lb, ub, eta);
  ss.addGridPoints();
  ss.printInfo(1);
  /****************************************************************************/
  /* construct SymbolicSet for the input space */
  /****************************************************************************/
  double ilb[iDIM] = { 0, 0 };
  double iub[iDIM] = { 0.65, 200 };
  double ieta[iDIM] = { 0.01, 4 };
  scots::SymbolicSet is(mgr, iDIM, ilb, iub, ieta);
  is.addGridPoints();
  is.printInfo(1);
  /****************************************************************************/
  /* setup class for symbolic model computation */
  /****************************************************************************/
  scots::SymbolicSet sspost(ss, 1); /* create state space for post variables */
  scots::SymbolicModelGrowthBound<state_type, input_type> abstraction(&ss, &is, &sspost);
  /* compute the transition relation */
  tt.tic();
  abstraction.computeTransitionRelation(system_post, radius_post);
  std::cout << std::endl;
  tt.toc();
  /* get the number of elements in the transition relation */
  std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;

  /****************************************************************************/
  /* we continue with the controller synthesis for FG (target) */
  /****************************************************************************/
  /* construct SymbolicSet for target (it is a subset of the state space)  */
  scots::SymbolicSet target(ss);
  /* add inner approximation of P={ x | H x<= h } form state space */
  double H[4 * sDIM] = { -1, 0, 0,
    1, 0, 0,
    0,-1, 0,
    0, 1, 0 };
  double h[4] = { -280,290,-0.5, 0.7 };
  target.addPolytope(4, H, h, scots::INNER);
  std::cout << "Target set details:" << std::endl;
  target.writeToFile("system_target.bdd");

  /* we setup a fixed point object to compute the reachable set */
  scots::FixedPoint fp(&abstraction);
  /* the fixed point algorithm operates on the BDD directly */
  BDD T = target.getSymbolicSet();
  /* we implement the nested fixed point algorithm
  *
  * mu X. nu Y. ( pre(Y) & T ) | pre(X)
  *
  */
  size_t i, j;
  /* outer fp*/
  BDD X = mgr.bddOne();
  BDD XX = mgr.bddZero();
  /* inner fp*/
  BDD Y = mgr.bddZero();
  BDD YY = mgr.bddOne();
  /* the controller */
  BDD C = mgr.bddZero();
  BDD U = is.getCube();
  /* as long as not converged */
  for (i = 1; XX != X; i++) {
    X = XX;
    BDD preX = fp.pre(X);
    YY = mgr.bddOne();
    for (j = 1; YY != Y; j++) {
      Y = YY;
      YY = (fp.pre(Y) & T) | preX;
    }
    XX = YY;
    std::cout << "Iterations inner: " << j << std::endl;
    BDD N = XX & (!(C.ExistAbstract(U)));
    C = C | N;
    //std::cout << C.CountMinterm(17) << std::endl;
  }
  std::cout << "Iterations outer: " << i << std::endl;
  scots::SymbolicSet controller0(ss, is);
  controller0.setSymbolicSet(C);
  std::cout << "Domain size: " << controller0.getSize() << std::endl;
  controller0.writeToFile("reachandstay_set6.bdd");





  C = fp.reach(T, 1);
  /****************************************************************************/
  /* last we store the controller as a SymbolicSet
  * the underlying uniform grid is given by the Cartesian product of
  * the uniform gird of the space and uniform gird of the input space */
  /****************************************************************************/
  scots::SymbolicSet controller1(ss, is);
  controller1.setSymbolicSet(C);
  std::cout << "Domain size: " << controller1.getSize() << std::endl;
  controller1.writeToFile("reach_set6.bdd");






  C = fp.safe(T, 1);
  scots::SymbolicSet controller2(ss, is);
  controller2.setSymbolicSet(C);
  std::cout << "Domain size: " << controller2.getSize() << std::endl;
  controller2.writeToFile("safe_set6.bdd");
  return 1;
}

