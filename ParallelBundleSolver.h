/*--------------------------------------------------------------------------*/
/*-------------------- File ParallelBundleSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ParallelBundleSolver class, which extends BundleSolver
 * to implement the computation of the C0Function in parallel, with a simple
 * master/slave scheme.
 *
 * Apart from that ParallelBundleSolver is completely equivalent to
 * BundleSolver, in particular regarding the class of :Block it can solve.
 * However, having ParallelBundleSolver solve a :Block with a single
 * C0Function hardly makes sense (although it also makes little difference,
 * in that for a single C0Function ParallelBundleSolver behaves exactly as
 * BundleSolver with an extremely limited overhead).
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ParallelBundleSolver
 #define __ParallelBundleSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BundleSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ParallelBundleSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a BundleSolver where multiple C05Function computations are parallelized
/** The ParallelBundleSolver class extends BundleSolver to implement the
 * computation of the C0Function in parallel, with a basic master/slave
 * scheme.
 *
 * Apart from that ParallelBundleSolver is completely equivalent to
 * BundleSolver, in particular regarding the class of :Block it can solve.
 * However, having ParallelBundleSolver solve a :Block with a single
 * C0Function hardly makes sense (although it also makes little difference,
 * in that for a single C0Function ParallelBundleSolver behaves exactly as
 * BundleSolver with an extremely limited overhead). */

class ParallelBundleSolver : public BundleSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 *  @{ */

/*----------------------------- CONSTANTS ----------------------------------*/

/*--------------------------------------------------------------------------*/
 /// public enum for the int algorithmic parameters
 /** Public enum describing the different algorithmic parameters of int type
  * that ParallelBundleSolver has in addition to these of BundleSolver (and
  * CDASolver, Solver, ThinComputeInterface). The value intLastPBndSlvPar is
  * provided so that the list can be easily further extended by derived
  * classes.
  *
  * In fact ParallelBundleSolver currently has no extra parameters w.r.t.
  * BundleSolver, except the fact that it does really react to the
  * intMaxThread one of ThinComputeInterface which BundleSolver instead
  * ignores (since it's completely sequential implementation).
 */

 enum int_par_type_PBndSlv {

 intLastPBndSlvPar = intLastBndSlvPar ,
 ///< first allowed new int parameter for derived classes
 /**< Convenience value for easily allow derived classes
  * to extend the set of int algorithmic parameters. */

 };  // end( int_par_type_PBndSlv )

/*--------------------------------------------------------------------------*/
 /// public enum for the double algorithmic parameters
 /** Public enum describing the different algorithmic parameters of double
  * type that ParallelBundleSolver has in addition to these of BundleSolver
  * (and CDASolver, Solver, ThinComputeInterface). The value dblLastPBndSlvPar
  * is provided so that the list can be easily further extended by derived
  * classes. */

 enum dbl_par_type_PBndSlv {
  dblPoolingInt = dblLastBndSlvPar ,
  ///< parameter for declaring the frequency of pooling of results

  dblLastPBndSlvPar ///< first new double parameter for derived classes
                    /**< Convenience value for easily allow derived classes
		     * to extend the set of double algorithmic parameters. */

  };  // end( dbl_par_type_BndSlv )

/*@} -----------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING ParallelBundleSolver ----------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ParallelBundleSolver
 *  @{ */

 /// constructor: ensure every field is initialized

 ParallelBundleSolver( void ) : BundleSolver() {
  // ensure all parameters are properly given their default value
  MaxThread = ThinComputeInterface::get_dflt_int_par( intMaxIter );
  PoolingInt = 1e-4;
  }

/*--------------------------------------------------------------------------*/
 /// destructor: does nothing special (explicitly)

 virtual ~ParallelBundleSolver() {}

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

 /// set the int parameters of ParallelBundleSolver
 /** Set the int parameters specific of ParallelBundleSolver, which actually
  * is a parameter of ThinComputeInterface that ParallelBundleSolver
  * "listens to" while BundleSolver does not:
  *
  * - intMaxThread [0]: maximum number of threads that compute() can spawn.
  *                     Actually this is not "threads" but "tasks", as it
  *   is implemented by calling intMaxThread std::asynch and therefore the
  *   actual number of threads may be smaller (down to actually none)
  *   depending on the C++ scheduler; yet, clearly this gives an upper bound
  *   on the total number of extra threads spawned when compute() is called
  *   (which are actually spawned each time that InnerLoop() is called within
  *   compute(), i.e., at every function iteration round, and then reined in
  *   when InnerLoop() ends). */

 void set_par( idx_type par , int value ) override {
  if( par == intMaxThread ) {
   MaxThread = value;
   }
  else
   BundleSolver::set_par( par , value );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// set the double parameters of ParallelBundleSolver
 /** Set the double parameters specific of ParallelBundleSolver, or calls the
  * BundleSolver version to deal with the rest:
  *
  * - dblPoolingInt [1e-4]: waiting time, in seconds, between each round of 
  *                         InnerLoop() checking for any of the std::asynch
  *   having finished and produced a result. That is, ParallelBundleSolver
  *   uses a dampened active wait strategy whereby there is no std::mutex or
  *   suchlike to automatically wake up the main thread when one of the tasks
  *   have finished, but the main thread waits the given amount between each
  *   test. Clearly, too long a waiting time means the main thread being slow
  *   to react and therefore wasted time, while too short a waiting time means
  *   too much time spent in active wait. The trade-off is dependent on the
  *   expected duration of each (or, better, of the faster among the)
  *   C05Function computations, and therefore left to the user's choice by
  *   means of this parameter. */

 void set_par( idx_type par , double value ) override {
  if( par == dblPoolingInt ) {
   if( value < 0 )
    throw( std::invalid_argument( "PoolingInt must be >= 0" ) );
   PoolingInt = value;
   }
  else
   BundleSolver::set_par( par , value );
  }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the BundleSolver
 *
 *  @{ */

 /*!!
 idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastBndSlvPar ) );
  }
  !!*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastPBndSlvPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dflt_dbl_par( idx_type par ) const override {
  if( par == dblPoolingInt )
   return( 1e-4 );
  else
   return( BundleSolver::get_dflt_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/
 
 int get_int_par( idx_type par ) const override {
  if( par == intMaxThread )
   return( MaxThread );
  else
   return( BundleSolver::get_int_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 
 double get_dbl_par( idx_type par ) const override {
  if( par == dblPoolingInt )
   return( PoolingInt );
  else
   return( BundleSolver::get_dbl_par( par ) );
  }

/*--------------------------------------------------------------------------*/

 /*!!
 idx_type int_par_str2idx( const std::string & name ) const override {
  const auto it = int_pars_map.find( name );
  if( it != int_pars_map.end() )
   return( it->second );
  else
   return( CDASolver::int_par_str2idx( name ) );
  }
  !!*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 idx_type dbl_par_str2idx( const std::string & name ) const override {
  if( name == "dblPoolingInt" )
   return( dblPoolingInt );
 
  return( BundleSolver::dbl_par_str2idx( name ) );
  }

/*--------------------------------------------------------------------------*/

 /*!!
 const std::string & int_par_idx2str( idx_type idx ) const override {
  if( ( idx >= intLastParCDAS ) && ( idx < intLastBndSlvPar ) )
   return( int_pars_str[ idx - intBPar1 ] );

  return( CDASolver::int_par_idx2str( idx ) );
  }
  !!*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::string & dbl_par_idx2str( idx_type idx ) const override {
  static const std::string __psname = "dblPoolingInt";
  if( idx == dblPoolingInt )
   return( __psname );

  return( BundleSolver::dbl_par_idx2str( idx ) );
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
 /* Performs the parallel inner loop: runs (at most) MaxThread std::asynch,
  * each one compute()-ing a different component, up until the conditions to
  * stop are satisfied or there no longer are available components to
  * evaluate. */

 Index InnerLoop( bool extrastep = false ) override;

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index MaxThread;    ///< maximum number of different threads (tasks)

 double PoolingInt;  ///< waiting time between each pooling round

 // generic fields- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*------------------------------ PRIVATE FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class ParallelBundleSolver )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ParallelBundleSolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File ParallelBundleSolver.h --------------------*/
/*--------------------------------------------------------------------------*/
