#include <iostream>
#include <fstream>
#include <math.h> 
#include <vector>

using namespace std;

// constant definition 
const double e = 1.602176634e-19; // C 
const double k = 8.987551792e9; // m**3 * kg / (s**2 * C**2)

// Convert energy (Mev) to velocity (m/s) for alpha particle
const double alpha_mass = 6.644657230e-27; // kg
const double alpha_charge = 2*e; // C
const double alpha_energy = 8.0108831e-13; // Joule ( ig * um**2 * ns**−2)

const double init_alpha_x = 0.0;
const double init_alpha_y = 0.0;

const double init_alpha_vx = sqrt(2.0 * alpha_energy / alpha_mass);
const double init_alpha_vy = 0.0;

// Set up the initial position of Gold nucleus
const double impact_parm = 1e-15; // m
const double gold_charge = 79.0 * e;

const double gold_x = 0.01;  //m
const double gold_y = -impact_parm;

vector<double> dxxdtt (double a_x , double a_y){
    double position_x = gold_x - a_x;
    double position_y = gold_y - a_y;
    double r = position_x * position_x + position_y * position_y;
    double force = -k * gold_charge * alpha_charge / ( r * alpha_mass );
    double F_x = force * position_x / sqrt(r);
    double F_y = force * position_y / sqrt(r);
    return {F_x, F_y};
}

int main (){
    // setup simulation parameters

    double t = 0.0; // s
    double t_max = 1; // s
    int iteration = 1e9;

    double dt = t_max/iteration;

    vector<double> temp;
    vector<double> alpha_t   = {t};
    vector<double> alpha_x   = {init_alpha_x};
    vector<double> alpha_y   = {init_alpha_y};
    vector<double> alpha_v_x = {init_alpha_vx};
    vector<double> alpha_v_y = {init_alpha_vy};

    double k_v_x , k_v_y;
    double k1_v_x , k1_v_y;
    double k2_v_x , k2_v_y;
    double k3_v_x , k3_v_y;
    double k4_v_x , k4_v_y;

    double k_x  , k_y;
    double k1_x , k1_y;
    double k2_x , k2_y;
    double k3_x , k3_y;
    double k4_x , k4_y;

    double half_dt = dt / 2.0;

    alpha_t.reserve(iteration);
    alpha_x.reserve(iteration);
    alpha_y.reserve(iteration);
    alpha_v_x.reserve(iteration);
    alpha_v_y.reserve(iteration);
    
    for (int i=0;i<iteration;i++){
        // calculate velocity RK4
        half_dt= dt/2.0;

        temp = dxxdtt(alpha_x[i],alpha_y[i]);
        k1_v_x  = temp[0];
        k1_v_y  = temp[1];
        k1_x    = alpha_v_x[i];
        k1_y    = alpha_v_y[i];

        temp = dxxdtt(alpha_x[i]+half_dt*k1_x,alpha_y[i]+half_dt*k1_y);
        k2_v_x = temp[0];
        k2_v_y = temp[1];
        k2_x   = alpha_v_x[i]+half_dt*k1_v_x; 
        k2_y   = alpha_v_y[i]+half_dt*k1_v_y;

        temp = dxxdtt(alpha_x[i]+half_dt*k2_x,alpha_y[i]+half_dt*k2_y);
        k3_v_x = temp[0];
        k3_v_y = temp[1];
        k3_x   = alpha_v_x[i]+half_dt*k2_v_x;
        k3_y   = alpha_v_y[i]+half_dt*k2_v_y;
        
        temp = dxxdtt(alpha_x[i]+dt*k3_x,alpha_y[i]+dt*k3_y);
        k4_v_x = temp[0];
        k4_v_y = temp[1];
        k4_x   = alpha_v_x[i]+dt*k3_v_x;
        k4_y   = alpha_v_y[i]+dt*k3_v_y; 

        k_v_x = dt/6.0*(k1_v_x+2*k2_v_x+2*k3_v_x+k4_v_x);
        k_v_y = dt/6.0*(k1_v_y+2*k2_v_y+2*k3_v_y+k4_v_y);
        k_x   = dt/6.0*(k1_x+2*k2_x+2*k3_x+k4_x);
        k_y   = dt/6.0*(k1_y+2*k2_y+2*k3_y+k4_y);

        alpha_v_x.push_back(alpha_v_x[i]+k_v_x);
        alpha_v_y.push_back(alpha_v_y[i]+k_v_y);
        alpha_x.push_back(alpha_x[i]+k_x);
        alpha_y.push_back(alpha_y[i]+k_y);
        alpha_t.push_back(alpha_t[i]+dt);
    }

    // Output the position and velocity of alpha particle to output.csv
    ofstream output("data.csv");
    output << "t" << "," << "x" << "," << "y" << "," << "v_x" << "," << "v_y" << endl;
    for (int i=0;i<iteration;i++){
        output << alpha_t[i] << "," << alpha_x[i] << "," << alpha_y[i] << "," << alpha_v_x[i] << "," << alpha_v_y[i] << endl;
    }
    output.close();


    return 0;
}