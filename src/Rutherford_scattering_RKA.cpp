#include <iostream>
#include <fstream>
#include <cmath> 
#include <vector>
#include <stdexcept>
#include <iomanip>

using namespace std;

// Constants define  
const double _e = 1/1.602176634e-19 ; // normalize charge
const double _m = 1/9.1093837e-31   ; // normalize mass
const double _l = 1/10e-15          ; // normalize length
const double _v = 1.0/299792458.0   ; // normalize speed
const double _s = _l / _v           ; // normalize time


// Convert energy (Mev) to velocity (m/s) for alpha particle
const double alpha_mass_SI = 6.644657230e-27 ; // Unit : kg
const double alpha_charge = 2                ; // Unit : _e

const double k = 8.987551792e9*_l*_l*_l*_m/(_s*_s*_e*_e); // Unit : m**3 * kg / (s**2 * C**2)
const double alpha_energy = 8.0108865e-13               ; // Unit : Joule ( kg * m**2 * s**−2)

const double init_alpha_vx = sqrt(2.0 * alpha_energy / alpha_mass_SI) * _v ; //Unit : _v
const double init_alpha_vy = 0.0                                           ; //Unit : _v

const double alpha_mass = _m * alpha_mass_SI ; //Unit : _m

// Set up the initial position of Gold nucleus
const double impact_parm = 1e-15 *_l *10 ; //Unit : _l
const double gold_charge = 79.0          ; //Unit : _e

const double gold_x = 100                ; //Unit : _l
const double gold_y = 0                  ; //Unit : _l


vector<double> dxxdtt (double a_x, double a_y){
    // Function dxxdtt calculates the acceleration of alpha particle.
    /*
    Input Arguments:
        a_x : x position of alpha particle
        a_y : y position of alpha particle
    return :
        F_x : x-axis acceleration of alpha particle
        F_y : y-axis acceleration of alpha particle
    */

    double position_x = gold_x - a_x;
    double position_y = gold_y - a_y;
    double r     = position_x * position_x + position_y * position_y   ;
    double force = -k * gold_charge * alpha_charge / ( r * alpha_mass );
    double F_x   = force * position_x / sqrt(r);
    double F_y   = force * position_y / sqrt(r);
    return {F_x, F_y};
}


vector<double> forwarder(double x, double y, double vx, double vy, double dt){
    // Function Forwarder performs a single RK4 step.
    /*
    Input Arguments:
        x   : x position
        y   : y position
        vx  : x velocity
        vy  : y velocity
        dt  : time step
    
    return:
        vector<double> contains four elements in order:
            1. step_x   : x position step
            2. step_y   : y position step
            3. step_v_x : x velocity step
            4. step_v_y : y velocity step
    */

    double half_dt = dt / 2.0;
    vector<double> temp      ;
    
    temp = dxxdtt(x,y)       ;
    double k1_v_x   = temp[0];
    double k1_v_y   = temp[1];
    double k1_x     = vx     ; 
    double k1_y     = vy     ;

    temp = dxxdtt(x + half_dt * k1_x, y + half_dt * k1_y);
    double k2_v_x = temp[0];
    double k2_v_y = temp[1];
    double k2_x   = vx + half_dt * k1_v_x;
    double k2_y   = vy + half_dt * k1_v_y;   

    temp = dxxdtt(x + half_dt * k2_x, y + half_dt * k2_y);
    double k3_v_x = temp[0];
    double k3_v_y = temp[1];
    double k3_x   = vx + half_dt * k2_v_x;
    double k3_y   = vy + half_dt * k2_v_y;

    temp = dxxdtt(x + dt * k3_x, y + dt * k3_y);
    double k4_v_x = temp[0];
    double k4_v_y = temp[1];
    double k4_x   = vx + dt * k3_v_x;
    double k4_y   = vy + dt * k3_v_y ;

    double step_v_x = dt/6.0*(k1_v_x + 2*k2_v_x + 2*k3_v_x + k4_v_x );
    double step_v_y = dt/6.0*(k1_v_y + 2*k2_v_y + 2*k3_v_y + k4_v_y );
    double step_x   = dt/6.0*(k1_x   + 2*k2_x   + 2*k3_x   + k4_x   );
    double step_y   = dt/6.0*(k1_y   + 2*k2_y   + 2*k3_y   + k4_y   );
    
    return {step_x , step_y , step_v_x , step_v_y};
}

vector<double> update_dt (vector<double> system_state, double dt, double ideal_error, double S1, double S2, double max_attempts){
    // Perform a big time step with dt = dt using rk4 routine
    vector<double> big_step  = forwarder(system_state[0], system_state[1], system_state[2], system_state[3], dt);
    vector<double> big_state = {big_step[0]+system_state[0], big_step[1]+system_state[1], big_step[2]+system_state[2], big_step[3]+system_state[3]};

    // Perform two small time steps with dt = dt / 2 using rk4 routine
    vector<double> small_step  = forwarder(system_state[0], system_state[1], system_state[2], system_state[3], dt/2.0);
    vector<double> small_state = {small_step[0]+system_state[0], small_step[1]+system_state[1], small_step[2]+system_state[2], small_step[3]+system_state[3]};
    small_step  = forwarder(small_state[0], small_state[1], small_state[2], small_state[3], dt/2.0);
    small_state = {small_step[0]+small_state[0], small_step[1]+small_state[1], small_step[2]+small_state[2], small_step[3]+small_state[3]};
    
    // Estimate truncation error using Error between two small steps and one big step
    double trunc_err = abs(big_state[0] - small_state[0])+abs(big_state[1] - small_state[1]);
    
    // Calculate score to determine whether dt is suitable
    double score = ideal_error / trunc_err;

    // Start iteration, if score is not in the range of 0.9 ~ 1.1 then update dt
    // Storage the number of attempts, if attempts > max_attempts then throw error
    double attempts = 0;
    while ( trunc_err > ideal_error) {  // Tranditional RKA 
    //while ( score > 1.1 || score < 0.9 ) { // Enables stronger restriction with this while condition
        if (attempts > max_attempts) {
            throw invalid_argument( "max attempts" );
        }

        double err_factor = pow((ideal_error/trunc_err) , 1.0/5.0);
        double dt_est = dt * err_factor;
        // dt_new is the new dt, which is the smaller one between dt_est and dt/S2
        double dt_new = S1*dt_est > S2*dt ? S2*dt : S1*dt_est < dt/S2 ? dt/S2 : S1*dt_est;
        
        dt = dt_new ;
        attempts ++ ;

        // Perform a big time step and two small steps again using rk4 routine to update new truncation error and score
        big_step    = forwarder(system_state[0], system_state[1], system_state[2], system_state[3], dt);
        big_state   = {big_step[0]+system_state[0], big_step[1]+system_state[1], big_step[2]+system_state[2], big_step[3]+system_state[3]};
        small_step  = forwarder(system_state[0], system_state[1], system_state[2], system_state[3], dt/2.0);
        small_state = {small_step[0]+system_state[0], small_step[1]+system_state[1], small_step[2]+system_state[2], small_step[3]+system_state[3]};
        small_step  = forwarder(small_state[0], small_state[1], small_state[2], small_state[3], dt/2.0);
        small_state = {small_step[0]+small_state[0], small_step[1]+small_state[1], small_step[2]+small_state[2], small_step[3]+small_state[3]};
        
        // Update Truncatiom Error and score
        trunc_err = abs(big_state[0] - small_state[0])+abs(big_state[1] - small_state[1]);
        score = ideal_error / trunc_err;
    }

    return {dt, trunc_err};
}


vector<vector<double>> rka(double dt, double t_max, vector<double> initial_condition ,double S1, double S2, double ideal_error, double max_attempts){

    vector<double> temp;
    vector<double> alpha_t   = {0};
    vector<double> alpha_x   = {initial_condition[0]};
    vector<double> alpha_y   = {initial_condition[1]};
    vector<double> alpha_v_x = {initial_condition[2]};
    vector<double> alpha_v_y = {initial_condition[3]};
    
    double i = 0  ;
    double t = 0.0;
    double p = 0.0;

    while(t<t_max){
        // current state
        vector<double> system_state = {alpha_x[i], alpha_y[i], alpha_v_x[i], alpha_v_y[i]};
        
        // update dt
        temp = update_dt(system_state, dt, ideal_error, S1, S2, max_attempts);
        dt  = temp[0];
        double err = temp[1];

        
        // Evolution for next step
        vector<double> evelution_step = forwarder(system_state[0], system_state[1], system_state[2], system_state[3], dt);
        double k_x   = evelution_step[0];
        double k_y   = evelution_step[1];
        double k_v_x = evelution_step[2];
        double k_v_y = evelution_step[3];

        alpha_v_x.push_back(alpha_v_x[i] + k_v_x);
        alpha_v_y.push_back(alpha_v_y[i] + k_v_y);
        alpha_x  .push_back(alpha_x  [i] + k_x  );
        alpha_y  .push_back(alpha_y  [i] + k_y  );
        alpha_t  .push_back(alpha_t  [i] + dt   );
        
        i += 1  ;
        t += dt ;
        
        cout << setprecision(20);
        if (round(t/t_max*100) > p){
            p = round(t/t_max*100);
            cout << p << "%, " << t << " s_bar, " << alpha_x[i] << " m_bar, " << alpha_y[i] << " m_bar, " << err << " m_bar" << endl;
        }
    }
    return {alpha_t, alpha_x, alpha_y, alpha_v_x, alpha_v_y};
}

int main (){
    vector<int> impact_parm_fac = {-10,-8,-6,-4,-2,0,2,4,6,8,10};
    for (int i=0;i<impact_parm_fac.size();i++){
        // Set up initial condition
        double init_alpha_x = 0.0                            ; //Unit : _l
        double init_alpha_y = impact_parm_fac[i]*impact_parm ; //Unit : _l
        vector<double> initial_condition = {init_alpha_x, init_alpha_y, init_alpha_vx, init_alpha_vy};

        // Adaptive Runge-Kutta
        vector<vector<double>> result = rka(10, 4000, initial_condition, 1-1E-2,1.1, 1E-13, 1e6);
        string name;
        name =  "./output/data/" + to_string(impact_parm_fac[i]) + ".csv";

        // Output the position and velocity of alpha particle to output.csv
        ofstream output(name);

        // Output csv header
        output << "t" << "," << "x" << "," << "y" << "," << "v_x" << "," << "v_y" << endl;

        // Get RKA iteration steps
        int iteration = result[0].size();

        // Output csv data
        for (int i=0;i<iteration;i++){
            output << result[0][i] << "," << result[1][i] << "," << result[2][i] << "," << result[3][i] << "," << result[4][i] << endl;
        }
        // Close file stream
        output.close();
    }
    

    
    
    // Calculate the final angle, theoretical angle and error
    // double final_angle = atan((result[2][iteration-1]-result[2][iteration-2])/(result[1][iteration-1]-result[1][iteration-2]))*180/M_PI;
    // if (final_angle < 0){
    //     final_angle = final_angle + 180;
    // }

    // double D = 8.987551792e9 * (gold_charge / _e) * (alpha_charge / _e) / alpha_energy;
    // double theoretical_angle = atan(D / 2 / (impact_parm / _l)) * 2 * 180 / M_PI;
    // double error             = abs(theoretical_angle - final_angle) / theoretical_angle * 100;
    // cout << "==========================================" << endl;
    // cout << "[Done!] data saved to "<< name << endl;
    // cout << "==========================================" << endl;
    // cout << "Final  angle : " << final_angle << " degree" << endl;
    // cout << "Theory angle : " << theoretical_angle << " degree" << endl;
    // cout << "==========================================" << endl;
    // cout << "Error : " << error << " %" << endl;
    // cout << "==========================================" << endl;

    
}

