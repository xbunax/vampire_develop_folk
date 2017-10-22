//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
void extract( std::string arg_string, std::vector<double> arg_vector );

//------------------------------------------------------------------------------
// Command line parsing function
//------------------------------------------------------------------------------
void command( int argc, char* argv[] ){
// Command line options for utility to be implemented:
//    --xyz - generate xyz file
//    --povray - generate povray files
//    --vector - generate raw xyz vector data
//    --vtk - generate vtk files
//    --cells - collate data into cells
//    --cell-size = x - define cell size
//    --povray-cells - generate renderable cell positions for povray
//    --verbose [= true, false] - set verbose output to screen
//    --colours = default, rwb [red-white-blue], oyb [orange-yellow-blue], jet, gs [grey-scale], cw [colour-wheel], <filename>
//    --objects = spins, cones, spheres, cubes
//    --slice = x,x,y,y,z,z
//    --multiscale = gradient, material, region

   // check for axis initialisations
   bool x_initialised = false;
   bool z_initialised = false;

   // temporary string for storing command line argument
   std::string temp_str;

   for (int arg = 1; arg < argc; arg++){

      // read prefix
      std::string sw=argv[arg]
      if (sw == "--vector-z"){

         // check number of args not exceeded
         if (arg+1 < argc){
            arg++;
            temp_str = string(argv[arg]);
         }
         else {
            terminaltextcolor(RED);
            std::cerr << "Error - expected 3 comma separated variables in brackets."
                      << "\n" << "Check for spaces in command-line arguments"
                      << std:endl;
            terminaltextcolor(WHITE);
            return EXIT_FAILURE;
         }

         // work through vector and extract values
         extract(temp_str,vdc::vector_z);

         // confirm initialisation of z-axis
         z_initialised = true;
      }
      else if (sw == "--vector-x"){

         // check number of args not exceeded
         if (arg+1 < argc){
            arg++;
            temp_str = string(argv[arg]);
         }
         else {
            terminaltextcolor(RED);
            std::cerr << "Error - expected 3 comma separated variables in brackets."
                      << "\n" << "Check for spaces in command-line arguments"
                      << std::endl;
            terminaltextcolor(WHITE);
            return EXIT_FAILURE;
         }

         // work through vector and extract values
         extract(temp_str,vdc::vector_x);

         // confirm initialisation of x-axis
         x_initialised = true;
      }

   }

   // check for valid axis initialisations
   if ( z_initialised && !x_initialised ){

      // find a line on the plane going through origin with normal vector_z
      vdc::vector_x = {1.0, 0.0, -1.0*vdc::vector_z[0]/vdc::vector_z[2]};
   }
   else if ( !z_initialised && x_initialised ){

      // x-axis cannot be initialised alone
      terminaltextcolor(RED);
      std::cerr << "Error - x-axis cannot be initialised alone."
                << "\n" << "To use 1D colour scheme, initialise z-axis instead"
                << std::endl;
      terminaltextcolor(WHITE);
      return EXIT_FAILURE;
   }
   else if ( z_initialised && x_initialised ){

      // check if input axes are orthogonal
      double zdotx;
      zdotx = vdc::vector_z[0]*vdc::vector_x[0] + vdc::vector_z[1]*vdc::vector_x[1] + vdc::vector_z[2]*vdc::vector_x[2];

      if ( (zdotx > 0.000001) || (zdotx < -0.000001) ){
         terminaltextcolor(RED);
         std::cerr << "Error - input axes are not orthogonal." << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }

   }

}

//------------------------------------------------------------------------------
// Extracts 3D vector coordinates from string: {x,y,z} or (x,y,z)
// where x,y and z are doubles
//------------------------------------------------------------------------------
void extract( std::string arg_string, std::vector<double> arg_vector ){
   int marker = 0; // position in the vector string

   // check for opening brackets
   if ( arg_string[marker] != ( "(" || "{" ) ){
      terminaltextcolor(RED);
      std::cerr << "Error - brackets required around 3 comma separated values"
                << std::endl;
      terminaltextcolor(WHITE);
      return EXIT_FAILURE;
   }

   // move to next character
   marker++;

   // read coordinates
   for ( int i = 0; i < 3; i++){
      std::string tmp_string;

      // read coordinate-value
      int j = 0;
      while ( arg_string[marker] != "," ){
         tmp_string.resize( tmp_string.size() +1 );
         tmp_string[j] = arg_string[marker];

         // move through number
         marker++;
         j++;
      }

      arg_vector[i] = std::stod(tmp_string);

      // skip comma, check for closing brackets
      if ( arg_string[marker] == "," ){
         marker++;
      }
      else if ( (arg_string[marker] != ( ")" || "}" )) && ( i == 2 ){
         terminaltextcolor(RED);
         std::cerr << "Error - brackets required around 3 comma separated values"
                   << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }
      else if ( (arg_string[marker] == ( ")" || "}" )) && ( i != 2 ){
         terminaltextcolor(RED)
         std::cerr << "Error - three coordinates required"
                   << std::endl;
         terminaltextcolor(WHITE);
         return EXIT_FAILURE;
      }

   }

   // normalise arg_vector
   double length;
   length = std:sqrt( arg_vector[0]*arg_vector[0] + arg_vector[1]*arg_vector[1] + arg_vector[2]*arg_vector[2] );
   arg_vector[0] = arg_vector[0]/length;
   arg_vector[1] = arg_vector[1]/length;
   arg_vector[2] = arg_vector[2]/length;

}

} // end of namespace vdc
