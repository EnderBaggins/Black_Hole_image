#include "metric.h"
#include "dorman.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h>
int main() {
    clock_t tStart = clock();
    double a = 0.0;                             // Will be my variable a
    double M = 1.0;                             //unit mass
    double theta0 = 10*2.0*M_PI/360;            // THIS IS THE ANGLE OF THE SCREEN IN DEGREES
    
    double rin = 5.0*M;                         //The radii of the accretion disk
    double rout = 20.0*M;
    double radius = M + std::sqrt(M*M-a*a);     //The radius of the black hole

    double D = 100.0*M;                         //The distance to the screen

    boyer_lindquist_metric metric(a, M);//Initializes the metric object
    
    rk45_dormand_prince integrator(6, 1.0e-12, 1.0e-12);//Initializes the integrator object

    auto dydt = [&metric] (double t, const std::vector<double> &y) {
        //Compute the metric

        metric.compute_metric(y[0], y[1]); // y[0] = r, y[1] = theta
        //finding the upper mu_t
        double mu_t_upper = std::sqrt(metric.gamma11 * y[3] * y[3] + metric.gamma22 * y[4] * y[4] + metric.gamma33 * y[5] * y[5])/metric.alpha;
        
        double dr_dt = metric.gamma11 * y[3] / mu_t_upper; //Checked
        double dtheta_dt = metric.gamma22 * y[4] / mu_t_upper; //Checked
        double dphi_dt = metric.gamma33 * y[5] / mu_t_upper - metric.beta3; //checked
        double dmur_dt = -metric.alpha*mu_t_upper*metric.d_alpha_dr + y[5]*metric.d_beta3_dr - 1/(2.0*mu_t_upper)*(y[3]*y[3]*metric.d_gamma11_dr + y[4]*y[4]*metric.d_gamma22_dr + y[5]*y[5]*metric.d_gamma33_dr); //checked
        double dmutheta_dt = -metric.alpha*mu_t_upper*metric.d_alpha_dtheta + y[5]*metric.d_beta3_dtheta - 1/(2.0*mu_t_upper)*(y[3]*y[3]*metric.d_gamma11_dtheta + y[4]*y[4]*metric.d_gamma22_dtheta + y[5]*y[5]*metric.d_gamma33_dtheta);//checked

        return std::vector<double> {dr_dt, dtheta_dt, dphi_dt, dmur_dt, dmutheta_dt,0.0}; // last element is zero since d muphi/dt = 0
    };

    double BH_rad_mod = 1.001; // I needed this otherwise would take forever to run (I think that the step size as nearing horizon becomes infinitely small)
    double accretion_angle_mod = 0.001; // additive modulator (aka how thick the disk is)
    auto stop = [&D,&radius,&rout,&rin,&BH_rad_mod,&accretion_angle_mod] (double t, const std::vector<double> &y) {
        bool toofar = y[0] >= D+1; // when distance bewteen photon and BH is greater than screen
        bool tooclose = y[0]<= radius*BH_rad_mod;//when photon is inside the black hole
        bool hitaccresion = y[0] <= rout && y[0] >= rin && y[1] <= M_PI/2.0+accretion_angle_mod && y[1] >= M_PI/2.0-accretion_angle_mod;
        //when photon hits the accretion disk
        return toofar or tooclose or hitaccresion; 
    };
    
    

    double screenxlim = 30.0;                       //The size of the screen in M distance units
    double screenylim = 20;
    
    int screenxres = 200;//1280;                    //Num of pixels along the x axis
    int screenyres = 100;//800;                     //Num of pixels along the y axis
    double pixels[screenxres][screenyres];          // STAR to let it be allocated dynamically
    
    for (int i =0;i<screenxres;i++){
        double x = 2.0*i*screenxlim/screenxres-screenxlim;
        for (int j =0;j<screenyres;j++){
            double y = 2.0*j*screenylim/screenyres-screenylim;
            // THis is run for each pixel in this for loop
            double alpha = y/D;
            double beta = x/D;
            double r = std::sqrt(D*D + x*x + y*y);
            double theta = theta0-alpha;
            double phi = beta;
            metric.compute_metric(r, theta); // y[0] = r, y[1] = theta
            double mur = -std::cos(beta)*std::cos(alpha)*std::sqrt(metric.g_11);
            double mutheta = std::sin(alpha)*std::sqrt(metric.g_22);
            double muphi = std::sin(beta)*std::cos(alpha)*std::sqrt(metric.g_33);
            std::vector<double> y0 = {r, theta, phi, mur, mutheta, muphi};
            
            integrator.integrate(dydt, stop, 0.01, 0.0, y0);
            auto finalpos = integrator.result.back();// NNEED TO SEE IF THIS IS APPROPRIATE
            
            if (finalpos[0] >= D){
                pixels[i][j] = 0.05; 
            }
            else if (finalpos[0] <= radius*BH_rad_mod){
                pixels[i][j] = 0.0; 
            }
            
            else if (finalpos[0] <= rout && finalpos[0] >= rin && finalpos[1] <= M_PI/2.0 + accretion_angle_mod && finalpos[1] >= M_PI/2.0 - accretion_angle_mod){
                // THis is where I implement brightness information
                double omega = 1.0/(a+std::pow(finalpos[0],1.5)/std::sqrt(M));
                metric.compute_metric(finalpos[0], finalpos[1]);
                double mu0 =-metric.alpha*std::sqrt(metric.gamma11*finalpos[3]*finalpos[3] + metric.gamma22*finalpos[4]*finalpos[4] + metric.gamma33*finalpos[5]*finalpos[5])+(-metric.beta3*finalpos[5]);
                double onez = (1 - omega*finalpos[5]/mu0)/std::sqrt(-metric.g_00-omega*omega*metric.g_33-2.0*omega*metric.g_03);
                pixels[i][j] = 1/(onez*onez*onez);
            }
            else{
                pixels[i][j] = 1000.0;// Really just to troubleshoot in case something else breaks
            }

        };
    };

    //std::ofstream output_file("BH_Image.csv");
    //std::ofstream output_file("close"+std::to_string(screenxres)+"_"+std::to_string(screenyres)+"a="+std::to_string(std::lround(100*a)/100)+"Th="+std::to_string(std::lround(theta0*360/(2*M_PI)))+".csv");
    std::ofstream output_file("test1.csv");
    for (int i = 0; i < screenxres; i++) {// I might have weird transpose things (but fixed in python code)
        output_file << pixels[i][0];
        for (int j = 1; j < screenyres; j++){
            output_file <<"," << pixels[i][j];
        };
        output_file << "\n";
    };
    output_file.close();
   tStart = clock() - tStart;
   double time_taken = ((double)tStart)/CLOCKS_PER_SEC; // in seconds
    std::cout << "Time taken: " << time_taken / 60 << " minutes" << std::endl;
     return 0;
}