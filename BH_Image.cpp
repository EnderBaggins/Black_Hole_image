#include "metric.h"
#include "dorman.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <omp.h> 
#include <nlohmann/json.hpp>

int main() {
    clock_t tStart = clock();                //To time the code

    //Loading In the JSON config file
    std::ifstream i("config.json");
    nlohmann::json config;
    i >> config;

    //Setting the parameters from the JSON file
    std::string savename = config["savename"];       //Name of the file to save the image
    //Visual Parameters
    double theta0 = config["theta0"];         // inclination angle of the observer in Degrees
    int xres = config["xres"];       //Num of pixels along the x axis
    int yres = config["yres"];       //Num of pixels along the y axis

    //Black Hole Parameters
    double a = config["a"];                         // Will be my variable a
    double M = config["M"];                          //unit mass
    double rin = config["rin"];                      //The inner radii of the accretion disk
    double rout = config["rout"];                  //The outer radii of the accretion disk
    double thetathick = config["thetathick"];          //The angular thickness of the disk
    double radius = M + std::sqrt(M*M-a*a);        //The radius of the black hole

    //Screen parameters
    double D = config["D"];         //The distance to the screen
    double xsize = config["xsize"];   //The size of the screen in M distance units
    double ysize = config["ysize"];




    /*
    std::string savename = "alexBH_Image";       //Name of the file to save the image
    //Black Hole Parameters
    double a = 0.99;                         // Will be my variable a
    double M = 1.0;                          //unit mass
    double rin = 5.0*M;                      //The radii of the accretion disk
    double rout = 20.0*M;
    double radius = M + std::sqrt(M*M-a*a);  //The radius of the black hole
    
    //Screen parameters
    double theta0 = 80;         // inclination angle of the observer in Degrees
    double D = 500.0*M;         //The distance to the screen
    int xres = 1280;       //Num of pixels along the x axis
    int yres = 800;       //Num of pixels along the y axis
    double xsize = 30.0;   //The size of the screen in M distance units
    double ysize = 20;
*/
    //Construct the solver for each thread
    int nthreads = omp_get_max_threads();
    std::cout << "number of threads: " << nthreads << std::endl;

    std::vector<boyer_lindquist_metric> metric_v;
    std::vector<rk45_dormand_prince> integrator_v;
    //Initializing the integrator and metric objects for each thread
    for (int i = 0; i < nthreads; i++) {
        integrator_v.emplace_back(6, 1.0e-12, 1.0e-12);
        metric_v.emplace_back(a, M);
    }

    double BH_rad_mod = 1.01;              // I needed this otherwise would take forever to run (I think that the step size as nearing horizon becomes infinitely small)
    //double accretion_angle_mod = 0.01;     // additive modulator (aka how thick the disk is)
    
    //Stop condition of the photon path
    auto stop = [&D,&radius,&rout,&rin,
                 &BH_rad_mod,&thetathick]
                (double t, const std::vector<double> &y) {
        bool toofar = y[0] >= D+1;                          // when distance bewteen photon and BH is greater than screen
        bool tooclose = y[0]<= radius*BH_rad_mod;           //when photon is inside the black hole
        bool hitaccresion = y[0] <= rout && y[0] >= rin     //when photon hits the accretion disk
                         && y[1] <= M_PI/2.0+thetathick 
                         && y[1] >= M_PI/2.0-thetathick;
        
        return toofar or tooclose or hitaccresion; 
    };//end of stop
    
    

    
    //Creating the image
    double pixels[xres][yres];
    theta0 = theta0*2.0*M_PI/360; // Converts theta0 to radians

   //Looping over the pixels
#pragma omp parallel for  
    for (int i =0;i<xres;i++){
        //get the id of the thread
        int thread_id = omp_get_thread_num();
    
        double x = 2.0*i*xsize/xres-xsize;    // This is the x position of the pixel
        for (int j =0;j<yres;j++){
            double y = 2.0*j*ysize/yres-ysize;// This is the y position of the pixel

            // THis is run for each pixel in this for loop
            double alpha = y/D;
            double beta = x/D;
            double r = std::sqrt(D*D + x*x + y*y);
            double theta = theta0-alpha;
            double phi = beta;

            auto &metric = metric_v[thread_id];
            metric.compute_metric(r, theta); // y[0] = r, y[1] = theta
            
            double mur = -std::cos(beta)*std::cos(alpha)*std::sqrt(metric.g_11);
            double mutheta = std::sin(alpha)*std::sqrt(metric.g_22);
            double muphi = std::sin(beta)*std::cos(alpha)*std::sqrt(metric.g_33);
            std::vector<double> y0 = {r, theta, phi, mur, mutheta, muphi};
            
            auto dydt = make_dydt(metric);
            
            //integrating the pixel until the stop condition
            integrator_v[thread_id].integrate(dydt, stop, 0.01, 0.0, y0);
            auto finalpos = integrator_v[thread_id].result.back();
            
            if (finalpos[0] >= D){
                pixels[i][j] = 0.0; // if the photon escapes
            }
            else if (finalpos[0] <= radius*BH_rad_mod){
                pixels[i][j] = 0.0; // if the photon gets trapped
            }
                                    // if the photon hits the accretion disk
            else if (finalpos[0] <= rout && finalpos[0] >= rin 
                  && finalpos[1] <= M_PI/2.0 + thetathick 
                  && finalpos[1] >= M_PI/2.0 - thetathick){
                
                
                // Calculating the brightness of the pixel
                // This is the doppler factor
                // and maybe a gravitational effect
                double omega = 1.0/(a+std::pow(finalpos[0],1.5)/std::sqrt(M));
                metric.compute_metric(finalpos[0], finalpos[1]);

                double mu0 =-metric.alpha*std::sqrt(metric.gamma11*finalpos[3]*finalpos[3] 
                                                  + metric.gamma22*finalpos[4]*finalpos[4] 
                                                  + metric.gamma33*finalpos[5]*finalpos[5])
                                          +(-metric.beta3*finalpos[5]);

                double onez = (1 - omega*finalpos[5]/mu0)
                             /std::sqrt(-metric.g_00-omega*omega*metric.g_33
                                        -2.0*omega*metric.g_03);

                pixels[i][j] = 1/(onez*onez*onez);}  // The brightness of the pixel
            else{
                pixels[i][j] = 1000.0;// Troubleshooter only if non of the above conditions are met
            }

        };
    };

    std::ostringstream filepath;
    filepath << "Output_files/" << savename << ".csv";
    std::ofstream output_file(filepath.str());
    // Save configuration information
output_file << "Configuration:\n";
output_file << "savename: " << savename << "\n";
output_file << "theta0: " << theta0 << "\n";
output_file << "xres: " << xres << "\n";
output_file << "yres: " << yres << "\n";
output_file << "a: " << a << "\n";
output_file << "M: " << M << "\n";
output_file << "rin: " << rin << "\n";
output_file << "rout: " << rout << "\n";
output_file << "thetathick: " << thetathick << "\n";
output_file << "D: " << D << "\n";
output_file << "xsize: " << xsize << "\n";
output_file << "ysize: " << ysize << "\n\n"; // Two newlines to separate config from data
    // Save the pixel data
    for (int i = 0; i < xres; i++) {// Saves as transposed matrix (just data.transpose() in python)
        output_file << pixels[i][0];
        for (int j = 1; j < yres; j++){
            output_file <<"," << pixels[i][j];
        };
        output_file << "\n";
    };
    output_file.close();
   
   //Time taken
   tStart = clock() - tStart;
   double time_taken = ((double)tStart)/CLOCKS_PER_SEC; // in seconds
    std::cout << "Time taken: " << time_taken / 60 << " minutes" << std::endl;
     return 0;



}