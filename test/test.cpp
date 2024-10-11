/*--------------------------------------------------------------------------*/
/*----------------------------- File main.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for testing MCFBlock. It just creates one and loads it from a
 * file, then prints it back to a file. It also de-serializes it,
 * serialize a copy out of it, and finally de-serialize the copy.
 *
 * Little more than a compilation check.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>

#include "BundleSolver.h"
#include <fstream>
#include <sstream>

#include "QPPnltMP.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;

using namespace SMSpp_di_unipi_it;

#if( OPT_USE_NAMESPACES )
 using namespace NDO_di_unipi_it;
#endif

const char *const ParF = "ParValue.qp";

/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
 if( argc < 2 ) {
  cerr << "Usage: " << argv[ 0 ]
       << " DIMACS_in [ DIMACS_out netCDF_out netCDF_out_2 ]" << endl;
  return( 1 );
  }

 ifstream ParFile( ParF );
 if( ! ParFile.is_open() )
  cerr << "Warning: cannot open parameters file """ << ParF << """" << endl;
 QPPenaltyMP *MP = new QPPenaltyMP( &ParFile );

 ParFile.close();

 auto Slv = Solver::new_Solver( "BundleSolver" );

 // pass the MPSolver to the Bundle  - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 BundleSolver * BndSlv = static_cast< BundleSolver * >( Slv );
 if( BndSlv == nullptr )
  throw( std::logic_error( "the solver is not of the Bundle type" ) );
                    
 // all done
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File main.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/

