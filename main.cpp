#include <iostream>
#include <fstream>
#include <armadillo>
#include <time.h>
#include <string>

using namespace std;
using namespace arma;

void writestring2file(string arg_filename, string arg_outstring);
double general_tridiag(vec &arg_a, vec &arg_b, vec &arg_c, vec &arg_u, vec &arg_y, int arg_n);
double specific_tridiag(vec &arg_u, vec &arg_y, int arg_n);
double LU_decomp(vec &arg_y, int n);

int main(int argc, char *argv[]){
    cout << "Armadillo version: "
	 << arma_version::as_string()
     << endl;
    //find cmd-line args
    int n; bool LU;
    if (argc == 1){
      n = 5;
      LU = false;
      cout << "No cmd-line args: n="
	   << n << " LU=" << LU
	   << endl;
    }
    else if (argc == 2) {
      n = atoi(argv[1]);
      LU = false;
      cout << "1 cmd-line args: n="
	   << n << " LU=" << LU
	   << endl;
    }
    else if (argc == 3) {
      n = atoi(argv[1]);
      string cmd_arg2(argv[2]);
      if (cmd_arg2 == "1") {
          LU = true;
      };
      cout << "2 cmd-line args: n="
	   << n << " LU=" << LU
       << endl;
    }
    else {
        cout << "ERROR you gave to many cmd-line exercise"
            << endl;
    }
    
    double t_gen, t_spec, t_LU; //CPU time in seconds it takes to calculate the diff. methods
    
    //generate x,f,y -vectors
    double x0 = 0.0;
    double x1 = 1.0;
    double h = (x1 - x0)/(n-1.0);
    vec x = linspace<vec>(x0,x1,n);
    vec y1 = zeros<vec>(n); // y = h^2 * f = h^2*100*e^(-10*x)
    vec y2 = zeros<vec>(n);
    y1 = h*h*100.0*exp(-10.0*x);
    y2 = h*h*100.0*exp(-10.0*x); //need to arrays since they will be overridden

    //generate diagonal vector elements and matrix
    vec a = ones<vec>(n); a *= -1.0;
    vec b = ones<vec>(n); b *= 2.0;
    vec c = ones<vec>(n); c *= -1.0;

    //generate vectors for u
    vec u_gen = zeros<vec>(n);
    vec u_spec = zeros<vec>(n);
    vec u_LU = y1; //starts out as y1, will be overridden

    //allocate parameters for storing data
    string time_filename = "../project1/data/dderiv_time_c++.dat";
    string u_filename = "../project1/data/dderiv_u_c++_n" + to_string(n) + "_";

    //start exercises
    if (LU) {
        /*calculate u using LU-decomposition*/
        t_LU = LU_decomp(u_LU, n); // turns empty array U_LU into solution

        //write timing results to file
        string string_time_data = "method: n=" + to_string(n) + " time=";
        writestring2file(time_filename, "LU " + string_time_data + to_string(t_LU));

    } else {
        /*calculate u using the general tridiagonal method*/
        t_gen = general_tridiag(a, b, c, u_gen, y1, n); //turns empty array u_gen into solution

        /*calculate u using the specific tridiagonal method*/
        t_spec = specific_tridiag(u_spec, y2, n); //turns empty array u_spec into solution

        //write timing-results to file
        string string_time_data = "tridiagonal method: n=" + to_string(n) + " time=";
        writestring2file(time_filename, "general " + string_time_data + to_string(t_gen));
        writestring2file(time_filename, "specific " + string_time_data + to_string(t_spec));

        /*write u(x)-results to file IF n is not too large */
        if (n <= 1e+4) {
            //make new data-file
            u_filename += "tridiag.dat";
            ofstream outfile; outfile.open(u_filename.c_str()); outfile.close();
            cout << "Started new file to /data/-directory" << endl;

            //add first line
            string string_u_data = "x, u_gen, u_spec";
            writestring2file(u_filename, string_u_data);
            //add all x_i, u_gen_i, u_spec_i to data-file
            for (int i=0; i < n; i++) {
                string_u_data = to_string(x(i)) + ", " + to_string(u_gen(i)) + ", " + to_string(u_spec(i));
                writestring2file(u_filename, string_u_data);
            }
        }
    }
}

double general_tridiag(vec &arg_a, vec &arg_b, vec &arg_c, vec &arg_u, vec &arg_y, int arg_n){
    double k;
    clock_t t0, t1;

    t0 = clock();
    for (int i=1; i<=arg_n-3; i++){
        k = arg_a(i+1)/( (double) arg_b(i) );
        arg_b(i+1) -= k*arg_c(i);
        arg_y(i+1) -= k*arg_y(i);
    }
    for (int i=arg_n-2; i>=1; i--){
        arg_u(i) = (arg_y(i) - arg_u(i+1)*arg_c(i))/( (double) arg_b(i) );
    }
    t1 = clock();

    return (t1 - t0)/((double) CLOCKS_PER_SEC); //measure time of forward and backward substitution
}

double specific_tridiag(vec &arg_u, vec &arg_y, int arg_n){
    clock_t t0,t1;
    vec d = zeros<vec>(arg_n);

    t0 = clock();
    for (int i=1; i<=arg_n-1; i++){
        d(i) = (i+1)/( (double) i );
    }

    for (int i=2; i<=arg_n-2; i++){
        arg_y(i) += arg_y(i-1)/d(i-1);
    }

    for (int i=arg_n-2; i>=1; i--){
        arg_u(i) = (arg_y(i) + arg_u(i+1))/d(i);
    }
    t1 = clock();

    return (t1 - t0)/((double) CLOCKS_PER_SEC); //measure time of forward and backward substitution
}

double LU_decomp(vec &arg_y, int arg_n){
    /*Solve the equation A*u = y were the matrix A
     * is fetched as 'arg_a', and y is fetched as 'arg_y'.
    */
    clock_t t0, t1;

    mat A (arg_n,arg_n, fill::eye); A *= 2.0;
    for (int i=0; i<arg_n-1; i++){
      A(i,i+1) = -1.0;
      A(i+1,i) = -1.0;
    } //filling A elementwise

    t0 = clock();
    arg_y = solve(A, arg_y); //solve A*x = L*U*x = L*w = y for x
    //array of argument 'arg_y' has now become the solution u of 'A*u = y'
    t1 = clock();

    return (t1 - t0)/((double) CLOCKS_PER_SEC); //measure time of solve()-function two times
}

void writestring2file (string arg_filename, string arg_outstring) {
  /* Take one single string and add it to the 
     last line of [filename] (defined locally),
     then add endline and close file.
     This is slow and unefficient, but not more is needed.
  */
  ofstream outfile;
  outfile.open(arg_filename.c_str(), ofstream::app);
  outfile << arg_outstring << endl;
  outfile.close();
  //cout << "wrote string to file" << endl;
} // writestring2file
