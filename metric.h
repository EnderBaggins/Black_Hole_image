#pragma once // I believe to make class file only called once
#include <iostream>
#include <vector>
#include <cmath>

class boyer_lindquist_metric { 
public:
    boyer_lindquist_metric(double a0, double M0) {
        //Don't need to find a in this section (that is done in the main file) which inputs a here
        //Initialize the parameter a an M
        a = a0;
        M = M0;
        //Potential other initializations
    }
    //This function computes the metric and stores it in a vector
    //The metric is stored in the following order: dr/dt, dtheta/dt, dphi/dt, dmur/dt, dmutheta/dt, dmuphi/dt where mu_i is the covariant momentum
    void compute_metric(double r, double theta) {
        //Compute the metric
        //Starting by computing the intermediate variable and their derivatives
        sintheta = std::sin(theta); //Computing sin(theta) once to save computation time
        costheta = std::cos(theta); //Computing cos(theta) once to save computation time

        rho2 = r * r + a * a * costheta * costheta;//checked2
        d_rho2_dr = 2.0 * r; //checked2
        d_rho2_dtheta = -2.0 * a * a * costheta * sintheta; //checked and fixed one thing2

        Delta = r * r - 2.0 * M * r + a * a;//checked2
        d_Delta_dr = 2.0 * r - 2.0 * M;//checked2
        d_Delta_dtheta = 0.0; //checked2

        Sigma = (r*r+a*a)*(r*r+a*a) - a*a*Delta*sintheta*sintheta;//checked2
        d_Sigma_dr = 4.0 *r *(r*r+a*a)-a*a*sintheta*sintheta*d_Delta_dr;//fixed something2
        d_Sigma_dtheta = -2.0 * a * a * Delta * sintheta * costheta; //checked2

        //Now compute the kerr metric components
        g_00 = -1+2.0 * M * r/rho2; // g_(t,t) //checked
        g_03 = -2.0*M*a*r*sintheta*sintheta/rho2; // g_(t,phi)//checked
        g_11 = rho2/Delta; // g_(r,r)//checked
        g_22 = rho2; // g_(theta,theta)//checked
        g_33 = Sigma*sintheta*sintheta/rho2; // g_(phi,phi)//checked
        //Now computing the lapse function and shift vector
        alpha = std::sqrt(rho2*Delta/Sigma); //checked2
        beta3 = -2*M*a*r/Sigma; //changed to equivalent expression //checked2
        //Now computing the upper gamma^ij (i.e the inverse of the spatial metric)
        gamma11 = Delta/rho2; //changed to equivalent expression //checked2
        gamma22 = 1.0/rho2;//changed to equivalent expression //checked2
        gamma33 = rho2/(Sigma*sintheta*sintheta);//changed to equivalent ex2pression //checked2
        //Computing the derivative values with respect to r
        // I analytically solved each one then plugged in
        d_alpha_dr = 1.0/(2.0*alpha) * (d_rho2_dr*Delta/Sigma + rho2*d_Delta_dr/Sigma - rho2*Delta*d_Sigma_dr/(Sigma*Sigma));//checked why is his differenet?
        d_beta3_dr = -2.0 * M * a *(1/Sigma - r/(Sigma*Sigma) * d_Sigma_dr); // checked
        d_gamma11_dr = 1.0/rho2 * d_Delta_dr - Delta/(rho2*rho2) * d_rho2_dr; //checked2
        d_gamma22_dr = -1/(rho2*rho2) * d_rho2_dr;//checked2
        d_gamma33_dr = d_rho2_dr/(Sigma * sintheta * sintheta) - rho2/(Sigma * Sigma * sintheta*sintheta) * d_Sigma_dr;//checked2
        //Computing the derivative values with respect to theta
        // I analytically solved each one then plugged in
        d_alpha_dtheta = 1.0/(2.0*alpha) * (d_rho2_dtheta*Delta/Sigma + rho2*d_Delta_dtheta/Sigma - rho2*Delta*d_Sigma_dtheta/(Sigma*Sigma));//checked2
        d_beta3_dtheta = (2.0 * M * a * r)/(Sigma*Sigma)*d_Sigma_dtheta; // changed to equivalent expression //checked2
        d_gamma11_dtheta = d_Delta_dtheta/rho2-Delta/(rho2*rho2)*d_rho2_dtheta; //checked2
        d_gamma22_dtheta = -1/(rho2*rho2)*d_rho2_dtheta; //checked
        d_gamma33_dtheta = d_rho2_dtheta/(Sigma*sintheta*sintheta) - rho2/(Sigma*Sigma*sintheta*sintheta)*d_Sigma_dtheta-2.0*rho2/(Sigma*sintheta*sintheta*sintheta)*costheta; // This was my issue just forgot some factors
        
        }
    double a, M;
    double rho2, Delta, Sigma; //intermediate variables
    double sintheta, costheta; //intermediate variables of sin and cos that will shorten computation time
    double d_rho2_dr, d_Delta_dr, d_Sigma_dr; //partial derivatives with respect to r
    double d_rho2_dtheta, d_Delta_dtheta, d_Sigma_dtheta; //partial derivatives with respect to theta
    double alpha, beta3; // only beta3 since the only the phi component of beta is non-zero Note this is the beta^phi not beta_phi
    double gamma11, gamma22, gamma33; //Upper gamma^ij, 0 = t, 1 = r, 2 = theta, 3 = phi
    double g_00, g_03, g_11, g_22, g_33; //components of the kerr metric
    double d_alpha_dr, d_beta3_dr, d_gamma11_dr, d_gamma22_dr, d_gamma33_dr; //partial derivatives with respect to r
    double d_alpha_dtheta, d_beta3_dtheta, d_gamma11_dtheta, d_gamma22_dtheta, d_gamma33_dtheta; //partial derivatives with respect to theta
    
    
};