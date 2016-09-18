//libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <cstdlib>
#include <sstream>

//declaring functions
int writestring2file(char filename[], char outstring[],
		     bool delete_file);
int write3vars2file (char filename[], double a,
		     double b, double c);
double func(double x);
void dfdx_2c(double& dfdx, double x, double h);
void dfdx_3c(double& dfdx, double x, double h);

//start main
int main (int argn, char* argv[]) {
  writestring2file(data_loc, "", 1); //make sure file is deleted
  writestring2file(data_loc, title, 0); //make new file with new 'title' as first line"
  writestring2file(data_loc, "h, f2c, f3c", 0); 

  for (int i=0; i<h_size; i++) {
    h[i+1] = h[i]/h_const; //making h-values
    dfdx_2c(df2c[i], x, h[i]);
    dfdx_3c(df3c[i], x, h[i]);

    write3vars2file(data_loc, h[i], df2c[i], df3c[i]); 
  } //for-loop over h

  delete h, df2c, df3c;
  return 0;
} // main

//define functions
int writestring2file (char filename[], char outstring[], bool delete_file) {
  /* Take one single string and add it to the 
     last line of [filename] (defined locally),
     then add endline and close file.
     This is slow and unefficient, but not more is needed.
  */
  if (delete_file) {
    /*delete preexisting file and exit*/
    if( std::remove( filename ) == 0 ) {
      std::cout << "one file successfully deleted:" << std::endl
		<< "\t" << filename << std::endl;
      return 0;
    } //if removed succesfully
    else {
      std::cout << "could not delete file:" << std::endl
		<< "\t" << filename << std::endl;
      return 0;
    } //else not removed succesfully
  } //if boolean 'delete_file'
  
  std::ofstream outfile;
  outfile.open(filename, std::ios::app); 
  outfile << outstring << std::endl;
  outfile.close();
  std::cout << "wrote string to: " << std::endl
	    << "\t" << filename << std::endl;
  return 0;
} // write2file

int write3vars2file (char filename[], double a, double b, double c) {
  /* Open a file and append three variables to it,
     using the csv-format.*/
  std::ofstream outfile;
  outfile.open(filename, std::ios::app); 
  outfile << std::scientific << std::setprecision(20)
	  << a << ", "
	  << b << ", "
	  << c << std::endl;
  outfile.close();
  /*std::cout << "writing skalars to file: " << std::endl
    << "\t" << filename << std::endl;*/
  return 0;
} // write3vars2file
