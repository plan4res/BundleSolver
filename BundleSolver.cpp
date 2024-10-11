/*--------------------------------------------------------------------------*/
/*------------------------ File BundleSolver.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BundleSolver class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Dipartimento di Matematica ed Informatica \n
 *         Universita' di Cagliari \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BundleSolver.h"

#include "LagBFunction.h"

#include "QPPnltMP.h"

#include "OSIMPSolver.h"

#include "OsiClpSolverInterface.hpp"

#include <iomanip>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define BLOG( l , x ) if( f_log && ( LogVerb > l ) ) *f_log << x

#define BLOG2( l , c , x ) if( f_log && ( LogVerb > l ) && c ) *f_log << x

/*--------------------------------------------------------------------------*/

#define USE_MPTESTER 0

// if USE_MPTESTER is nonzero, the MPSolver is a MPTester. in particular, if
// USE_MPTESTER == 1 then the master of the MPTester is an OSIMPSolver and
// the slave is a QPPenaltyMP, while if USE_MPTESTER != 1 then the master is
// a QPPenaltyMP and the slave is an OSIMPSolver

#if USE_MPTESTER
 #include "MPTester.h"
#endif

/*--------------------------------------------------------------------------*/
/* If WHICH_OSI_QP is nonzero, the OsiSolverInterface used in OSIMPSolver (if
 * any. i.e., unless a QPPenaltyMP is used) will be one of the two specific
 * OsiXXXSolverInterface that support a QP Master problem, i.e.,
 *
 * - 1: a OsiCpxSolverInterface
 *
 * - 2: a OsiGrbSolverInterface
 *
 * This requires some includes that can be avoided otherwise, especially
 * since both Cplex and Gurobi are commercial and therefore the corresponding
 * OsiXXXSolverInterface are in general not be available to all users. */

#ifndef WHICH_OSI_QP
 #define WHICH_OSI_QP 2
#endif

#if WHICH_OSI_QP == 1
 #include "OsiCpxSolverInterface.hpp"
#elif WHICH_OSI_QP == 2
 #include "OsiGrbSolverInterface.hpp"
 #include "gurobi_c++.h"
#endif

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Perform long and costly checks on the data structures, coded bit-wise:
  *
  * - CHECK_DS & 1 == checks the data structures representing the bundle and
  *                   the global pools against the MPSolver and the
  *                   C05Function(s)
  *
  * - CHECK_DS & 2 == checks that the aggregated linearization produced by
  *                   the C05Function agrees with that produced by the
  *                   MPSolver
  *
  * - CHECK_DS & 4 == checks that the aggregated linearization errors
  *                   directly computed with fresh data (linearization +
  *                   constant) out of the C05Function agree with these
  *                   stored in the MPSolver
  *
  * - CHECK_DS & 8 == checks that the lower bounds out of the C05Function
  *                   agree with these stored in the MPSolver
  */
#else
 #define CHECK_DS 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr HpNum Nearly  = 1.01;
static constexpr HpNum Nearly2 = 1.02;

static constexpr char LogBnd = 16;        // log Bundle changes
static constexpr char LogVar = 32;        // log variables changes

static constexpr Index tSP1Msk = 12;  // mask for tSPar1: the long-term t-s
static constexpr Index kSLTTS =  4;   // "soft" long-term t-strategy
static constexpr Index kHLTTS =  8;   // "hard" long-term t-strategy
static constexpr Index kBLTTS = 12;   // "balancing" long-term t-strategy
static constexpr Index kEGTTS = 16;   // "endgame" long-term t-strategy
static constexpr Index tSPHMsk1 = 192;  // mask for heuristics: bits 6 and 7
static constexpr Index tSPHMsk2 = 768;  // mask for heuristics: bits 7 and 8

static constexpr unsigned char RstAlg = 1;  // don't reset algorithmic params
static constexpr unsigned char RstCrr = 2;  // don't reset current point to
                                            // all-0, use Variable value()

static constexpr auto InINF = SMSpp_di_unipi_it::Inf< Index >();

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/
// set precision for long floats (10 digits) in scientific notation

static inline std::ostream & def( std::ostream & os ) {
 os.setf( std::ios::scientific , std::ios::floatfield );
 os << setprecision( 10 );
 return( os );
 }

/*--------------------------------------------------------------------------*/
// set precision for short floats (2 digits) in scientific notation

static inline std::ostream & shrt( std::ostream & os ) {
 os.setf( std::ios::scientific , std::ios::floatfield );
 os << setprecision( 2 );
 return( os );
 }

/*--------------------------------------------------------------------------*/
// set precision for short floats (4 digits) in fixed notation

static inline std::ostream & fixd( std::ostream & os ) {
 os.setf( std::ios::fixed , std::ios::floatfield );
 os << setprecision( 4 );
 return( os );
 }

/*--------------------------------------------------------------------------*/
// cleanly print +/-INF

static inline void pval( std::ostream & os , double val ) {
 if( val == BundleSolver::INFshift )
  os << "INF";
 else
  if( val == -BundleSolver::INFshift )
   os << "-INF";
  else
   os << val;
 }

/*--------------------------------------------------------------------------*/

static void Compact( BundleSolver::Vec_VarValue & g ,
		     BundleSolver::c_Subset & B )
{
 // takes a "dense" n-vector g and "compacts" it deleting the elements whose
 // indices are in B; all elements of B must be in the range 0 .. n, B must
 // be ordered in increasing sense
 // the remaining entries in g are shifted left of the minimum possible
 // amount in order to fill the holes left by the deleted ones
 // g is *not* resized in here

 auto Bit = B.begin();
 auto i = *(Bit++);
 auto git = g.begin() + (i++);

 for( ; Bit != B.end() ; ++i ) {
  auto h = *(Bit++);
  while( i < h )
   *(git++) = g[ i++ ];
  }

 std::copy( g.begin() + i , g.end() , git );

 }  // end( Compact )

/*--------------------------------------------------------------------------*/

static void set_difference_in_place( BundleSolver::Subset & S1 ,
				     BundleSolver::c_Subset & S2 )
{
 // removes from S1 all elements in S2, resizing it accordingly
 // both S1 and S2 are assumed to be ordered and with unique elements

 if( S1.empty() )  // nothing to delete from
  return;          // nothing to do

 auto S1it = S1.begin();
 auto S2it = S2.begin();

 // first phase: find the first element present in both S1 and S2

 for( ; ; ) {
  while( ( S1it != S1.end() ) && ( *S1it < *S2it ) )
   ++S1it;
  if( S1it == S1.end() )
   break;
  while( ( S2it != S2.end() ) && ( *S1it > *S2it ) )
   ++S2it;
  if( S2it == S2.end() )
   break;
  if( *S1it == *S2it )
   break;
  }

 if( ( S1it == S1.end() ) || ( S2it == S2.end() ) ) // if there are none
  return;                                           // all done

 // now S1it points to the first element in S1 == than the first in S2
 // elements in S1 after the common one(s) will have to be moved
 auto S1wit = S1it++;  // skip the first equal element
 S2it++;

 for( ; ( S1it != S1.end() ) && ( S2it != S2.end() ) ; ) {
  while( ( S1it != S1.end() ) && ( *S1it < *S2it ) )
   *(S1wit++) = *(S1it++);
  if( S1it == S1.end() )
   break;
  while( ( S2it != S2.end() ) && ( *S1it > *S2it ) )
   ++S2it;
  if( S2it == S2.end() )
   break;
  if( *S1it == *S2it ) { ++S1it; ++S2it; }
  }

 while( S1it != S1.end() )  // copy the part remaining after the end of S2
  *(S1wit++) = *(S1it++);

 S1.resize( std::distance( S1.begin() , S1wit ) );

 }  // end( set_difference_in_place )

/*--------------------------------------------------------------------------*/

static void set_union_in_place( BundleSolver::Subset & S1 ,
				BundleSolver::c_Subset & S2 )
{
 // make S1 to be the union of S1 and S2
 if( S2.empty() )
  return;

 if( S1.empty() )
  S1 = S2;
 else {
  BundleSolver::Subset tmp;
  std::set_union( S1.begin() , S1.end() , S2.begin() , S2.end() ,
		  std::back_inserter( tmp ) );
  S1 = std::move( tmp );
  }
 }  // end( set_union_in_place )

/*--------------------------------------------------------------------------*/

static void set_union_in_place( BundleSolver::Subset & S1 ,
				BundleSolver::Subset && S2 )
{
 // make S1 to be the union of S1 and S2, if useful destroy S2 in the process
 if( S2.empty() )
  return;

 if( S1.empty() )
  S1 = std::move( S2 );
 else {
  BundleSolver::Subset tmp;
  std::set_union( S1.begin() , S1.end() , S2.begin() , S2.end() ,
		  std::back_inserter( tmp ) );
  S1 = std::move( tmp );
  }
 }  // end( set_union_in_place )

/*--------------------------------------------------------------------------*/

static double norm( const BundleSolver::Vec_VarValue & v , char t )
{
 double res = 0;
 if( t == 0 ) {    // INF-norm
  for( auto el : v )
   if( std::abs( el ) > res )
    res = std::abs( el );
  }
 else
  if( t == 1 )    // 1-norm
   for( auto el : v )
    res += std::abs( el );
  else {          // 2-norm
   for( auto el : v )
    res += el * el;

   res = sqrt(  res );
   }

 return( res );
 }

/*--------------------------------------------------------------------------*/

static void vect_sum( BundleSolver::Vec_VarValue & v1 , double * v2 )
{
 for( auto & el : v1 )
  el += *(v2++);
 }

/*--------------------------------------------------------------------------*/

static void chgsign( double * v , Index n )
{
 for( const auto ev = v + n ; v < ev ; ++v )
  *v = - *v;
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BundleSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0( BundleSolver );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register BundleSolverState to the State factory

SMSpp_insert_in_factory_cpp_0( BundleSolverState );

/*--------------------------------------------------------------------------*/
// define and initialize here the vector of int parameters names
const std::vector< std::string > BundleSolver::int_pars_str = {
 "intBPar1" ,
 "intBPar2" ,
 "intBPar3" ,
 "intBPar4" ,
 "intBPar6" ,
 "intBPar7" ,
 "intMnSSC" ,
 "intMnNSC" ,
 "inttSPar1" ,
 "intMaxNrEvls" ,
 "intDoEasy" ,
 "intWZNorm" ,
 "intFrcLstSS" ,
 "intTrgtMng" ,
 "intMPName" ,
 "intMPlvl" ,
 "intQPmp1" ,
 "intQPmp2",
 "OSImp1" ,
 "OSImp2" ,
 "OSImp3" ,
 "intRstAlg"
 };

// define and initialize here the vector of double parameters names
const std::vector< std::string > BundleSolver::dbl_pars_str = {
 "dblNZEps" ,
 "dbltStar" ,
 "dblMinNrEvls" ,
 "dblBPar5" ,
 "dblm1" ,
 "dblm2" ,
 "dblm3" ,
 "dblmxIncr" ,
 "dblmnIncr" ,
 "dblmxDecr" ,
 "dblmnDecr" ,
 "dbltMaior" ,
 "dbltMinor" ,
 "dbltInit" ,
 "dbltSPar2" ,
 "dbltSPar3" ,
 "dblCtOff"
 };

// define and initialize here the vector of string parameters names
const std::vector< std::string > BundleSolver::str_pars_str = {
 "strEasyCfg" ,
 "strHardCfg"
 };

// define and initialize here the map for int parameters names
const std::map< std::string , BundleSolver::idx_type >
 BundleSolver::int_pars_map = {
 { "intBPar1" , BundleSolver::intBPar1  } ,
 { "intBPar2" , BundleSolver::intBPar2  } ,
 { "intBPar3" , BundleSolver::intBPar3 } ,
 { "intBPar4" , BundleSolver::intBPar4 } ,
 { "intBPar6" , BundleSolver::intBPar6 } ,
 { "intBPar7" , BundleSolver::intBPar7 } ,
 { "intMnSSC" , BundleSolver::intMnSSC } ,
 { "intMnNSC" , BundleSolver::intMnNSC } ,
 { "inttSPar1" , BundleSolver::inttSPar1 } ,
 { "intMaxNrEvls" , BundleSolver::intMaxNrEvls } ,
 { "intDoEasy" , BundleSolver::intDoEasy } ,
 { "intWZNorm" , BundleSolver::intWZNorm } ,
 { "intFrcLstSS" , BundleSolver::intFrcLstSS } ,
 { "intTrgtMng" , BundleSolver::intTrgtMng } ,
 { "intMPName" , BundleSolver::intMPName } ,
 { "intMPlvl" , BundleSolver::intMPlvl } ,
 { "intQPmp1" , BundleSolver::intQPmp1 } ,
 { "intQPmp2" , BundleSolver::intQPmp2 } ,
 { "intOSImp1" , BundleSolver::intOSImp1 } ,
 { "intOSImp2" , BundleSolver::intOSImp2 } ,
 { "intOSImp3" , BundleSolver::intOSImp3 } ,
 { "intRstAlg" , BundleSolver::intRstAlg } ,
 };

// define and initialize here the map for double parameters names
const std::map< std::string , BundleSolver::idx_type >
 BundleSolver::dbl_pars_map = {
 { "dblNZEps" , BundleSolver::dblNZEps } ,
 { "dbltStar" , BundleSolver::dbltStar } ,
 { "dblMinNrEvls" , BundleSolver::dblMinNrEvls } ,
 { "dblBPar5" , BundleSolver::dblBPar5 } ,
 { "dblm1" , BundleSolver::dblm1 } ,
 { "dblm2" , BundleSolver::dblm2 } ,
 { "dblm3" , BundleSolver::dblm3 } ,
 { "dblmxIncr" , BundleSolver::dblmxIncr } ,
 { "dblmnIncr" , BundleSolver::dblmnIncr } ,
 { "dblmxDecr" , BundleSolver::dblmxDecr } ,
 { "dblmnDecr" , BundleSolver::dblmnDecr } ,
 { "dbltMaior" , BundleSolver::dbltMaior } ,
 { "dbltMinor" , BundleSolver::dbltMinor } ,
 { "dbltInit" , BundleSolver::dbltInit } ,
 { "dbltSPar2" , BundleSolver::dbltSPar2 } ,
 { "dbltSPar3" , BundleSolver::dbltSPar3 } ,
 { "dblCtOff" , BundleSolver::dblCtOff }
 };

// define and initialize here the default int parameters
const std::vector< int > BundleSolver::dflt_int_par = {
 10 ,  // intBPar1
100 ,  // intBPar2
  1 ,  // intBPar3
  1 ,  // intBPar4
  0 ,  // intBPar6
  3 ,  // intBPar7
  0 ,  // intMnSSC
  3 ,  // intMnNSC
 12 ,  // inttSPar1
  2 ,  // intMaxNrEvls
  1 ,  // intDoEasy
  2 ,  // intWZNorm
  0 ,  // intFrcLstSS
  0 ,  // intTrgtMng
  0 ,  // intMPName
  0 ,  // intMPlvl
  0 ,  // intQPmp1
  0 ,  // intQPmp2
  4 ,  // intOSImp1
  0 ,  // intOSImp2
  1 ,  // intOSImp3
  2    // intRstAlg, default value:
       // RstAlg = 0  -  reset algorithmic parameters
       // RstCrr = 1  -  set current point to using values of the Variable
 };

// define and initialize here the default double parameters
const std::vector< double > BundleSolver::dflt_dbl_par = {
 0 ,      // dblNZEps
 1e+2 ,   // dbltStar
 0 ,      // dblMinNrEvls
 30 ,     // dblBPar5
 0.01 ,   // dblm1
 0.99 ,   // dblm2
 0.99 ,   // dblm3
 10 ,     // dblmxIncr
 1.5 ,    // dblmnIncr
 0.1 ,    // dblmxDecr
 0.66 ,   // dblmmDecr
 1e+6 ,   // dbltMaior
 1e-6,    // dbltMinor
 1 ,      // dbltInit
 1e-3 ,   // dbltSPar2
 0 ,      // dbltSPar3
 1e-1     // dblCtOff
 };

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF BundleSolver --------------------------*/
/*--------------------------------------------------------------------------*/

int BundleSolver::compute( bool changedvars )
{
 // ensure no concurrent accesses - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 lock();                        // ... either from other threads

 if( Result == kStillRunning )  // ... or from the same
  throw( std::logic_error( "BundleSolver::compute() called within itself" )
	 );

 // basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! f_Block ) {    // no Block is there to compute
  Result = kBlockLocked;

  BundleSolver_error_return:
  f_mutex.unlock();  // unlock the mutex
  return( Result );
  }

 if( MaxIter == 0 ) {  // no iteration must be performed
  Result = kStopIter;
  goto BundleSolver_error_return;
  }

 Result = kStillRunning;    // still working

 // start timer now (so that processing Modification is included) - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 c_start = std::chrono::system_clock::now();

 double tot_time = 0;  // independently keep function evaluation time
 long tot_NrEvls = 0;  // total number of C05Function evaluations

 // first, process any outstanding Modification - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // v_mod is atomically copied in a temporary data structure to be processed,
 // but while the latter happens new Modification may come in; hence,
 // process_outstanding_Modification() may be called more than once

 while( num_outstanding_Modification() ) {
  bool owned = f_Block->is_owned_by( f_id );       // check if already locked
  if( ( ! owned ) && ( ! f_Block->read_lock() ) ) {  // if not read_lock now
   Result = kBlockLocked;                           // return error on failure
   goto BundleSolver_error_return;
   }

  // if there are "easy" components and changing them is supported, process
  // the corresponding Modification (stored it v_FakeSolver)
  if( NrEasy && ( DoEasy & ~1 ) )
   process_outstanding_easy_Modification();

  // process any other Modification
  process_outstanding_Modification();

  if( ! owned )             // if the Block was actually read_locked
   f_Block->read_unlock();  // read_unlock it
  }

 // initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there is a 0-th component and the value had not been computed before
 // (maybe because the 0-th component has changed), do it now
 if( f_lf && ( Fi0Lmb == INFshift ) ) {
  f_lf->compute( true );
  Fi0Lmb = rs( f_lf->get_upper_estimate() );
  if( UpFiLmbdef == NrFi ) {  // ready to compute the total upper bound
   ++UpFiLmbdef;              // do so
   UpFiLmb.back() = std::accumulate( UpFiLmb.begin() , --(UpFiLmb.end()) ,
				     Fi0Lmb );
   }
  if( LwFiLmbdef == NrFi ) {  // ready to compute the total lower bound
   ++LwFiLmbdef;              // do so
   LwFiLmb.back() = std::accumulate( LwFiLmb.begin() , --(LwFiLmb.end()) ,
				     Fi0Lmb );
   }
  }

 f_wFi = NrFi - 1;  // since not all components are necessarily evaluated
                    // at all iterations, the order in which they are seen
 // may be important; keep track of the last evaluated component so as to
 // proceed round-robin-like across multiple iterations
 double lastETT = 0;  // last "time" eEveryTTime events have been called
 ParIter = 0;         // number of iterations in this call
 ++SCalls;            // one more call
 RifeqFi = ( UpRifFi == UpFiLmb );  // true if the reference values are right

 if( NeedsG1() )
  G1.resize( NrFi );
 else
  G1.clear();

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( ; ; ) {
  // check if time is over- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( MaxTime < INFshift ) && ( get_elapsed_time() > MaxTime ) ) {
   BLOG( 1 , " ~ stop due to max time" << std::endl );
   Result = kStopTime;
   break;
   }

  // run time-periodic events - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( EveryTTm && ( get_elapsed_time() >= lastETT + EveryTTm ) ) {
   for( auto & ev : v_events[ eEveryTTime ] ) {
    auto res = ev();
    if( res == eStopOK ) { Result = kOK; break; }
    if( res == eStopError ) { Result = kError; break; }
    }

   lastETT = get_elapsed_time();  // reset counter
   }

  if( Result != kStillRunning ) {
   BLOG( 1 , " ~ stop due to time-periodic event" << std::endl );
   break;
   }

  // construct the direction d- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //!! PrintBundle();

  FormD();

  // some log - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // another iteration (master problem solution)- - - - - - - - - - - - - - -

  ++ParIter;

  // check for "bad" termination- - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( Result == kStopTime )  // time is up
   continue;                 // return at the start and stop

  if( Result == kInfeasible ) {  // the Master Problem is infeasible
   BLOG( 1 , " ~ stop (infeasible)" << std::endl );
   break;
   }

  if( Result >= kError ) {  // problems in the Master Problem solver
   BLOG( 1 , " ~ error in the MPSolver" << std::endl );
   break;
   }

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( IsOptimal() ) {  // if optimality is detected
   // run optimality events - - - - - - - - - - - - - - - - - - - - - - - - -
   int res = eContinue;
   for( auto & ev : v_events[ eBeforeTermination ] )
    if( ( res = ev() ) != eContinue )
     break;

   if( res == eForceContinue ) {
    BLOG( 1 , " ~ optimal stop aborted by optimality event" << std::endl );
    continue;  // go back to master problem solution
    }

   if( res == eStopError ) {
    BLOG( 1 , " ~ stop (error) by optimality event" << std::endl );
    Result = kError;
    break;
    }

   BLOG( 1 , " ~ stop (optimal)" << std::endl );
   Result = kOK;
   break;
   }

  // run iteration-periodic events- - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( EverykIt && ( ! ( ParIter % EverykIt ) ) )
   for( auto & ev : v_events[ eEverykIteration ] ) {
    auto res = ev();
    if( res == eStopOK ) { Result = kOK; break; }
    if( res == eStopError ) { Result = kError; break; }
    }

  if( Result != kStillRunning ) {
   BLOG( 1 , " ~ stop due to time-periodic event" << std::endl );
   break;
   }

  // check if "ex-ante" Noise Reduction is needed - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ensure that the Sigma* is "not too negative", if it is increase t (if
  // possible) and re-solve the MP; note that this kind of NR only happens if
  // the oracle is "unfaithful", i.e., it pretends to provide information with
  // the required accuracy but in fact it does not
  //
  // however, avoid doing any of this if the linearization errors are not
  // computed w.r.t. the "true" value of UpFiLmb but w.r.t. a "random"
  // reference value, since then the fact that linearization errors are
  // negative is not meaningful

  if( RifeqFi && ( Sigma < - max_error( UpRifFi.back() , RelAcc ) ) &&
      ( Sigma <= - m3 * DST ) ) {
   if( t >= tMaior ) {
    BLOG( 1 , " ~ stop: NR required but t maximum" << std::endl );
    Result = kLowPrecision;
    break;
    }

   t = std::min( t * mxIncr , tMaior );
   BLOG( 2 , " ~ NR: t increased to " << shrt << t << std::endl );
   tHasChgd = true;
   continue;
   }

  // update out-of-base counters- - - - - - - - - - - - - - - - - - - - - - -

  UpdtCntrs();

  // Hard Long-Term t-strategy for quadratic stabilization- - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the hard long-term t-strategy requires t to increase if the step is too
  // small, and therefore has to be checked before the others
  // however, it is only viable under a quadratic stabilization
  //
  // however, avoid doing any of this if the linearization errors are not
  // computed w.r.t. the "true" value of UpFiLmb

  if( ( ( ! ( MPName & 1 ) ) || ( MPName & 4 ) ) &&
      ( tStar > 0 ) && ( ( tSPar1 & tSP1Msk ) == kHLTTS ) && RifeqFi ) {

   double AFL = std::abs( UpFiLmb.back() );
   if( AFL < 1 )
    AFL = 1;

   if( abs( vStar.back() ) <= tSPar2 * EpsU * AFL ) {
    BLOG( 1 , "small v => increase t" << std::endl << "           " );

    // collect two numbers vc and vl such that v( tNew ) >= vc + tNew * vl
    // we require that v( tNew ) >= vc + tNew * vl = tSPar2 * EpsU * AFL
    // ==> tNew = ( tSPar2 * EpsU * AFL - vc ) / vl

    double vl , vc;
    Master->SensitAnals( vl , vc );

    double tt;
    if( - vl < Eps< HpNum >() )  // v( t ) is [almost] constant ==> D*_t [~]= 0
     tt = tStar;                 // ==> the CP model is ~bounded
    else
     tt = std::min( tStar , ( tSPar2 * EpsU * AFL * Nearly + vc ) /
		            ( - vl ) );

    if( ( tHasChgd = ( tt != t ) ) ) {
     t = tt;
     continue;         // loop only if t changes
     }
    }
   }  // end if( Hard t-strategy )  - - - - - - - - - - - - - - - - - - - - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // compute Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormLambda1( t );

  // update the number of items to be fetched from the oracle - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UpdtaBP3();

  // eliminate outdated info- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // This is done *after* the call to Master->SensitAnals() in the Hard
  // Long-Term t-strategy and to FormLambda1(), because elimination of items
  // from the bundle may make the current solution of the master problem
  // invalid, and therefore all solution information may be lost. In theory
  // this should not happen, since only items "out of base" are eliminated,
  // and therefore the solution remains optimal; however, not all MPSolvers
  // may behave in this respect.

  SimpleBStrat();

  // run the inner loop - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // first some initializations - - - - - - - - - - - - - - - - - - - - - - -
  // all stuff that must be computed/changed inside InnerLoop()

  Alfa1 = 0;
  ScPr1 = NeedsScPr1() ? Master->ReadGid() : 0;
  if( NeedsG1() ) {
   G1Norm = INFshift;
   G1.assign( NrFi , double( 0 ) );
   }

  CurrNrEvls.assign( NrFi , Index( 0 ) );
  MPchgs = 0;  // != 0 if the MP is guaranteed to change enough after the
               // insertion of new information to ensure convergence

  auto start = std::chrono::system_clock::now();

  auto cnt = InnerLoop();

  auto end = std::chrono::system_clock::now();
  std::chrono::duration< double > elapsed = end - start;

  tot_time += elapsed.count();
  tot_NrEvls += std::accumulate( CurrNrEvls.begin() , CurrNrEvls.end() , 0 );

  CmptdinL = false;

  // compute DeltaFi- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1.back() == INFshift )
   DeltaFi = INFshift;
  else
   if( UpFiLmb1.back() == -INFshift )
    DeltaFi = -INFshift;
   else
    DeltaFi = UpFiLmb1.back() - UpRifFi.back();

  // update FiBest- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1.back() < UpFiBest ) {
   UpFiBest = UpFiLmb1.back();
   if( MaxSol > 1 )
    LmbdBst = Lambda1;
   }

  // update the "aggregated" Alfa1 and ScPr1- - - - - - - - - - - - - - - - -

  UpdateHeuristicInfo();

  // some log about the newly obtained information- - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log2( tot_time );

  // check whether either any error has occurred or time has expired- - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( UpFiLmb1.back() == -INFshift ) {
   BLOG2( 1 , f_convex , " ~ stop (Fi = -INF)" << std::endl );
   BLOG2( 1 , ! f_convex , " ~ stop (Fi = INF)" << std::endl );
   Result = kUnbounded;
   break;
   }

  if( Result == kError ) {
   BLOG( 1 , " ~ stop (error)" << std::endl );
   break;
   }

  if( Result == kStopTime )  // time has ran up inside InnerLoop()
   continue;                 // go back at the beginning to stop

  // check for the conditional lower bound- - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ! TrueLB ) &&
      ( UpFiBest <= LowerBound.back() *
	            ( 1 - ( LowerBound.back() > 0 ? RelAcc : - RelAcc ) ) )
      ) {
   BLOG( 1 , "            FiBest " );
   BLOG2( 1 , f_convex , "< conditional LB" );
   BLOG2( 1 , ! f_convex , "> conditional UB" );
   BLOG( 1 , ": unbounded " << std::endl );
   if( UpFiLmb1.back() < UpFiLmb.back() )  // Lambda1 is better than Lambda
    GotoLambda1();                         // go to Lambda1
   else                                    // if not
    if( ! RifeqFi )    // and the alfas are not computed w.r.t. UpFiLmb
     GotoLambda();     // ensure they are so

   Result = kUnbounded;
   break;
   }

  // avoid the t-changing phase if a vertical linearization has been found- -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is because the vertical linearization making Lambda1 unfeasible,
  // surely "change enough the master problem already"
  // yet, one possible t-strategy would be to set t to the largest value
  // that would have produced a feasible point: t := Alfa1 / ( - ScPr1 )
  // (with Alfa1 and ScPr1 of that particular constraint, though, not the
  // "global" ones)

  if( MPchgs > 1 ) {
   if( ParIter >= MaxIter ) {  // if we have done too many iterations
    BLOG( 1 , " ~ stop due to max iter" << std::endl );
    Result = kStopIter;        // stop already
    break;
    }
   else                        // otherwise
    continue;                  // go to the next one
   }

  // avoid the t-changing phase if the linearization errors are not reliable-
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is because we have firmly established one feasible (finite) upper
  // estimate in Lambda1, which ends the "phase 0" in which the linearization
  // errors were computed against an arbitrary value and starts the "phase 1"
  // in which the real optimization takes place

  if( ( ! RifeqFi ) && ( UpFiLmb1.back() < INFshift ) ) {
   // if we are still in "phase 0", and we just found a point where the
   // function value is finite, end the "phase 0" by immediately jumping
   // there. note that one may expect the thing on the function value to be
   // redundant since any component evaluating to +INF should generate a
   // vertical linearization and therefore set MPchgs = 2, which is acted
   // upon right above, but this may not happen. which is a problem if
   // MPchgs == 0 (but this is acted upon right below) but not otherwise,
   // since a "normal" NS will be done which is the right thing to do
   BLOG( 1 , "            Fi1 defined ==> SS " << std::endl );
   GotoLambda1();         // go to the feasible point
   continue;              // and start the actual minimization of Fi()
   }

  // check if noise reduction has to be done- - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ! MPchgs ) {
   if( t >= tMaior ) {
    BLOG( 1 , "            stop: NR required but t maximum" << std::endl );
    Result = kLowPrecision;
    break;
    }
   t = std::min( t * mxIncr , tMaior );
   BLOG( 1 , "            NR: t increased to " << shrt << t << std::endl );
   tHasChgd = true;
   continue;
   }

  // the NS / SS decision - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note again the "<" in the SS condition below (which means this is ever
  // so slightly stronger than it should), which is there to avoid the
  // condition to fire when UpFiLmb1.back() == INF == UpTrgt

  SSDone = ( UpFiLmb1.back() < UpTrgt ) ? true : false;

  VarValue tt = t , tm = t , tp = t;  // setup for the heuristic t

  if( SSDone ) {  // SS - - - - - - - - - - - - - - - - - - - - - - - - - - -

   BLOG( 1 , std::endl << " SS[" << CSSCntr << "]: DFi = " << shrt );
   if( f_convex ) {
    BLOG( 1 , DeltaFi << def << " ~ Up1(" << UpFiLmb1.back()
	      << ") <= UpTrgt(" << UpTrgt << ")" );
    }
   else
    BLOG( 1 , - DeltaFi << def << " ~ Lw1(" << - UpFiLmb1.back()
	      << ") >= LwTrgt(" << - UpTrgt << ")" );

   if( tSPar1 & 1 ) {
    tt = Heuristic( tSPar1 >> 6 );
    BLOG( 1 , " ~ Ht = " << shrt << tt );
    }

   if( tSPar3 ) {
    tp *= std::abs( tSPar3 );
    if( tSPar3 > 0 )
     tm /= tSPar3;
    }

   if( ++CSSCntr > MnSSC ) {  // increasing t is possible: note the ">"
    // due to the fact that the counter has just been increased
    if( ( ( tSPar1 & tSP1Msk ) == kBLTTS )  &&
	( DSTS <= tSPar2 * Sigma ) && ( CSSCntr < 10 ) ) {  //!! 10!
     // if the "balancing" long-term t-strategy is active and D*_t( 1 )
     // is small already, inhibit t increases (but not small heuristic
     // decreases, if active) unless "too many SS happened"
     BLOG( 1 , " ~ small D*_t( 1 )" );
     tp = t;
     }
    else {
     tm = t * mnIncr;  // minimum significant increase
     tp = t * mxIncr;  // maximum significant increase
     CSSCntr = 0;      // a significant increase happened, reset counter
     }
    }

   BLOG( 1 , std::endl );

   GotoLambda1();
   CNSCntr = 0;
   CmptdinL = ( cnt == NrFi - NrEasy );
   }
  else {        // NS - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   BLOG( 1 , std::endl << " NS[" << CNSCntr << "]: " );
   BLOG2( 1 , DeltaFi < INFshift , "DFi = " << shrt << rs( DeltaFi )
	      <<  " ~ " << def );
   if( f_convex ) {
    BLOG( 1 , "Lw1(" << def << LwFiLmb1.back() << ") >= LwTrgt(" << LwTrgt
	      << ")" );
    }
   else
    BLOG( 1 , "Up1(" << - LwFiLmb1.back() << ") <= UpTrgt(" << - LwTrgt
	      << ")" );

   if( tSPar1 & 2 ) {
    tt = Heuristic( tSPar1 >> 8 );
    BLOG( 1 , " ~ Ht = " << shrt << tt );
    }

   if( tSPar3 ) {
    tm /= std::abs( tSPar3 );
    if( tSPar3 > 0 )
     tp *= tSPar3;
    }

   if( ++CNSCntr > MnNSC ) {  // decreasing t is possible: note the ">"
    // due to the fact that the counter has just been increased
    if( ( ( ( tSPar1 & tSP1Msk ) == kSLTTS ) ||
	  ( ( tSPar1 & tSP1Msk ) == kHLTTS ) ) &&
	( abs( vStar.back() ) <= tSPar2 * EpsU * max_error() ) ) {
     // if either the "hard" or the "soft" long-term t-strategy is active
     // and v* is small already, inhibit t decreases (but not small
     // heuristic increases, if active)
     BLOG( 1 , " small v" );
     tm = t;
     }
    else
     if( ( ( tSPar1 & tSP1Msk ) == kBLTTS ) &&
	 ( tSPar2 * DSTS >= abs( Sigma ) ) ) {
      // if the "balancing" long-term t-strategy is active and D*_t( 1 )
      // is large already, inhibit t decreases (but not small heuristic
      // increases, if active); note that one may add the clause "unless
      // too many NS happened", i.e., "&& ( CNSCntr < 20 )": this version
      // avoids problems which may occur with ill-set tStar or tSPar2, but
      // it may give worse performances with "difficult" problems
      // also note the "abs( Sigma )": Sigma should be positive, but in
      // case it is not the control would always be true irrespectively of
      // the magnitude of tSPar2 and tStar just because of the sign
      BLOG( 1 , " ~ large D*_t( 1 )" );
      tm = t;
      }
     else {
      tm = t * mxDecr;  // maximum significant decrease
      tp = t * mnDecr;  // minimum significant decrease
      CNSCntr = 0;      // a significant decrease happened, reset counter
      }
    }


   BLOG( 1 , std::endl );
   CSSCntr = 0;

   }   // end else( NS )- - - - - - - - - - - - - - - - - - - - - - - - - - -

  // actually update t- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // if the endgame t-strategy fires (note the "/ 10"!!), the regular
  // t-updating mechanism is superseeded
  if( ( tSPar1 & kEGTTS ) && ( UpFiLmb.back() < INFshift ) &&
      ( DSTS < max_error() / 10 ) ) {
    tt = std::max( t * ( mxDecr + mnDecr ) / 2 , tMinor );
    BLOG( 1 , " ~ endgame, t = " << shrt << tt );
    //!! the reverse should also be done: if sigma is small and D*( t* ) is
    //!! large, t should be increased --> but this would happen surely at
    //!! the beginning, it should be done only near the end
    }
  else             // regular update mechanism
   if( tm != tp )  // if t can change, select it in [ tm , tp ]
    tt = std::min( std::min( tMaior , tp ) ,
		   std::max( std::max( tMinor , tm ) , tt ) );
   else            // else
    tt = t;        // keep it as it is

  if( ( tHasChgd = ( t != tt ) ) )
   t = tt;

  // check max number of iterations - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ParIter >= MaxIter ) {
   BLOG( 1 , " ~ stop due to max iter" << std::endl );
   Result = kStopIter;
   break;
   }
  }  // end( main loop )

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if necessary, force one last SS to the stability center - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FrcLstSS && ( ! CmptdinL ) &&
     ( ( Result == kOK ) || ( Result == kStopIter ) ||
       ( Result == kLowPrecision ) ) ) {
  BLOG( 1 , "            Recomputing the current point" << std::endl );

  UpFiLmb1 = UpFiLmb;
  LwFiLmb1 = LwFiLmb;
  Fi0Lmb1 = Fi0Lmb;
  UpTrgt = UpFiLmb1.back();
  LwTrgt = LwFiLmb1.back();

  FiStatus.assign( NrFi , kUnEval );
  for( Index i = 0 ; i < NumVar ; i++ )
   LamVcblr[ i ]->set_value( Lambda[ i ] );

  // note that Alfa1, ScPr1, G1 are computed inside GetGi() that is not
  // called inside this call to InnerLoop(), so they are not initialised
  CurrNrEvls.assign( NrFi , Index( 0 ) );
  MPchgs = 0;  // != 0 if the MP is guaranteed to change enough after the
               // insertion of new information to ensure convergence

  auto start = std::chrono::system_clock::now();

  auto cnt = InnerLoop( true );

  auto end = std::chrono::system_clock::now();
  std::chrono::duration< double > elapsed = end - start;

  tot_time += elapsed.count();

  CmptdinL = ( cnt == NrFi - NrEasy );

  // not being able to compute all non-easy components is an error
  if( ( ! CmptdinL ) && ( Result != kStopTime ) )
   Result = kError;
  }

 // final printouts - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_log && ( LogVerb >= 1 ) ) {
  *f_log << std::endl << "Call " << SCalls << ": "  << fixd << ParIter
	 << " ~ " << tot_NrEvls  << " ~ " << get_elapsed_time() << " ~ "
	 << tot_time << " -> ";
  switch( Result ) {
   case( kOK ):           *f_log << "optimal"; break;
   case( kStopTime ):     *f_log << "max time"; break;
   case( kStopIter ):     *f_log << "max iter"; break;
   case( kInfeasible ):   *f_log << "infeasible"; break;
   case( kUnbounded ):    *f_log << "unbounded"; break;
   case( kLowPrecision ): *f_log << "inexact oracle"; break;
   default:               *f_log << "error";
   }
  if( ( Result != kInfeasible ) && ( Result != kUnbounded ) ) {
   *f_log << " ~ Fi* = " << def;
   pval( *f_log , rs( UpRifFi.back() ) );
   }
  *f_log << std::endl;
  }

 unlock();  // unlock the mutex

 return( Result );

 }  // end( BundleSolver::compute )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::set_Block( Block * block )
{
 if( f_Block == block )  // registering to the same Block
  return;                // cowardly and silently return

 if( f_Block ) {  // changing from a previous oracle - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  guts_of_destructor();   // deallocate memory
  }

 Solver::set_Block( block );  // attach to the new Block

 if( ! f_Block )  // that was actually clearing the Block
  return;         // all done

 // lock the Block - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
                       "LagrangianDualSolver: unable to lock the Block" ) );

 // generate the abstract representation
 f_Block->generate_abstract_variables();
 f_Block->generate_abstract_constraints();
 f_Block->generate_objective();

 /* Two types of block can be handled by the BundleSolver:

     1. Only one single non-smooth function
     2. a sum of some non-smooth functions

    The algorithm here developed aims at solving non-constrained non-smooth
    optimization. The block can have box constraints at the most.
    It is expected the block to have in the first case a FRealObjective
    whose the function is a C05Function one and having no children, while in
    the second case a FRealObjective one whose the function is a
    LinearFunction and having as many sub-blocks as the number of components.
    In the latter case, each sub-block must not contain any Variable or
    Constraint.cVariable may have a lower and upper bound. If the lower
    bound  has a finite value, it must be 0. */

 const auto & sb = f_Block->get_nested_Blocks();

 if( sb.empty() ) {  // no sub-Block
  // the objective function of the Block must be a C05Function  - - - - - - -

  auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
  if( ! obj )
   throw( std::invalid_argument( "objective is not a FRealObjective" ) );

  auto c05f = dynamic_cast< C05Function * >( obj->get_function() );
  if( ! c05f )
   throw( std::invalid_argument( "the objective is not a C05Function" ) );

  f_convex = c05f->is_convex();
  if( ( ! f_convex ) && ( ! c05f->is_concave() ) )
   throw( std::invalid_argument(
			     "only convex or concave objectives allowed " ) );

  if( ( f_convex && ( obj->get_sense() == Objective::eMax ) ) ||
      ( ( ! f_convex ) && ( obj->get_sense() == Objective::eMin ) ) )
   throw( std::invalid_argument( "can only minimize convex / maximize concave"
				 ) );
  v_c05f.push_back( c05f );
  f_lf = nullptr;
  }
 else {  // there are sub-Block
  // the objective function of the block must be a LinearFunction- - - - - - -

  if( ! f_Block->get_objective() )  // there is no Objective
   f_lf = nullptr;
  else {
   auto obj = dynamic_cast< FRealObjective * >( f_Block->get_objective() );
   if( ! obj )
    throw( std::logic_error( "the objective is not a real function" ) );

   if( ! obj->get_function() ) // the FRealObjective has no Function
    f_lf = nullptr;
   else {
    f_lf = dynamic_cast< LinearFunction * >( obj->get_function() );
    if( ! f_lf )
     throw( std::logic_error( "the objective is not a LinearFunction" ) );

    if( ! f_lf->get_num_active_var() )  // the LinearFunction has no Variable
     f_lf = nullptr;
    }
   }

  v_c05f.resize( sb.size() );

  for( Index i = 0 ; i < sb.size() ; ++i ) {  // for each sub-block
   // the objective function of each sub-block must be a C05Function - - - - -

   auto obj = dynamic_cast< FRealObjective * >( sb[ i ]->get_objective() );
   if( ! obj )
    throw( std::logic_error( "the objective is not a real function" ) );

   auto c05f = dynamic_cast< C05Function * >( obj->get_function() );
   if( ! c05f )
    throw( std::logic_error( "the objective is not a C05Function" ) );

   // all have to be the same convexity - - - - - - - - - - - - - - - - - - -
   v_c05f[ i ] = c05f;
   if( ! i )
    f_convex = c05f->is_convex();
   else
    if( c05f->is_convex() != f_convex )
     throw( std::invalid_argument(
			  "objectives must be all convex or all concave" ) );

   // all have to be max/min in the right way- - - - - - - - - - - - - - -
   if( ( ! f_convex ) && ( ! c05f->is_concave() ) )
   throw( std::invalid_argument(
			     "only convex or concave objectives allowed " ) );

   if( ( f_convex && ( obj->get_sense() == Objective::eMax ) ) ||
       ( ( ! f_convex ) && ( obj->get_sense() == Objective::eMin ) ) )
    throw( std::invalid_argument( "can only minimize convex/maximize concave"
				  ) );

   // nephews are not allowed- - - - - - - - - - - - - - - - - - - - - - - - -
   if( sb[ i ]->get_nested_Blocks().size() )
    throw( std::logic_error( "nephew are not allowed" ) );

   // Variable not allowed - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( sb[ i ]->get_static_variables().size() )
    throw( std::logic_error( "static Variable are not allowed" ) );

   if( sb[ i ]->get_dynamic_variables().size() )
    throw( std::logic_error( "dynamic Variable are not allowed" ) );

   // neither are Constraint - - - - - - - - - - - - - - - - - - - - - - - - -
   if( sb[ i ]->get_static_constraints().size() )
    throw( std::logic_error( "static Constraint are not allowed" ) );

   if( sb[ i ]->get_dynamic_constraints().size() )
    throw( std::logic_error( "dynamic Constraint are not allowed" ) );
   }  // end( for each sub-Block )
  }  // end( there are sub-Block )

 // the set of "active" Variable in all Function must be the same- - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NumVar = v_c05f[ 0 ]->get_num_active_var();
 LamVcblr.resize( NumVar );
 auto vi = std::as_const( v_c05f[ 0 ] )->begin();
 for( Index i = 0 ; i < NumVar ; ++vi )
  LamVcblr[ i++ ] = static_cast< ColVariable * >( & (*vi) );

 if( f_lf ) {
  if( f_lf->get_num_active_var() != NumVar )
   throw( std::logic_error( "the list of active Variable do not match" ) );

  auto v = f_lf->begin();
  for( auto vi = LamVcblr.begin() ; vi != LamVcblr.end() ; ++v , ++vi )
   if( static_cast< ColVariable * >( & (*v) ) != *vi )
    throw( std::logic_error( "the list of active Variable do not match" ) );
  }

 for( Index i = 1 ; i < sb.size() ; ++i ) {
  if( v_c05f[ i ]->get_num_active_var() != NumVar )
   throw( std::logic_error( "the list of active Variable do not match" ) );

  auto v = v_c05f[ i ]->begin();
  for( auto vi = LamVcblr.begin() ; vi != LamVcblr.end() ; ++v , ++vi )
   if( static_cast< ColVariable * >( & (*v) ) != *vi )
    throw( std::logic_error( "the list of active Variable do not match" ) );
  }

 // if some Variable are present, they are of the ColVariable type - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index NumBVar = 0;  // count the number of Variable in the Block
 auto v_s_Variable = f_Block->get_static_variables();
 for( auto & el : v_s_Variable ) {
  if( un_any_thing_0( ColVariable , el , ++NumBVar ) )
   continue;
  if( un_any_thing_1( ColVariable , el , NumBVar += var.size() ) )
   continue;
  if( un_any_thing_K( ColVariable , el , NumBVar += var.num_elements() ) )
   continue;
  throw( std::logic_error( "some static Variable is not a ColVariable" ) );
  }

 auto v_d_Variable = f_Block->get_dynamic_variables();
 for( auto & el : v_d_Variable ) {
  if( un_any_thing_0( std::list< ColVariable > , el , ++NumBVar ) )
   continue;
  if( un_any_thing_1( std::list< ColVariable > , el ,
		      NumBVar += var.size() ) )
   continue;
  if( un_any_thing_K( std::list< ColVariable > , el ,
		      NumBVar += var.num_elements() ) )
   continue;
  throw( std::logic_error( "some dynamic Variable is not a ColVariable" ) );
  }

 if( NumBVar < NumVar )
  throw( std::logic_error( "too few ColVariable in the Block" ) );

 // check that the Variable in the Block agree with that in the C05Function- -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector< ColVariable * > LamBVcblr( NumBVar );

 Index cnt = 0;
 for( auto & el : v_s_Variable )
  un_any_static( el , [ & ]( ColVariable & sv ) { LamBVcblr[ cnt++ ] = & sv;
                  } , un_any_type< ColVariable >() );

 for( auto & el : v_d_Variable )
  un_any_dynamic( el , [ & ]( ColVariable & sv ) { LamBVcblr[ cnt++ ] = & sv;
                  } , un_any_type< ColVariable >() );

 std::sort( LamBVcblr.begin() , LamBVcblr.end() );

 std::vector< ColVariable * > LamVcblrO( LamVcblr );
 std::sort( LamVcblrO.begin() , LamVcblrO.end() );

 if( ! std::includes( LamBVcblr.begin() , LamBVcblr.end() ,
		      LamVcblrO.begin() , LamVcblrO.end() ) )
 throw( std::logic_error(
		   "some ColVariable in C05Function are not in the Block" ) );

 LamVcblrO.clear();
 LamBVcblr.clear();

 // if some Constraint are present, their can only be either BoxConstraint
 // (with LHS == 0), LB0Constraint or NNConstraint - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // one day general linear constraints will be allowed
 //
 // note that un_any_thing() only serves to verify that the stuff is of the
 // right type, and therefore it has to do nothing; this is obtained by
 // passing it as, the "function" argument, a void --> void lambda doing
 // nothing immediately applied to nothing, cue the curios list of
 // parentheses "[](){}()"

 for( auto & el : f_Block->get_static_constraints() ) {
  if( un_any_thing( BoxConstraint , el , [](){}() ) )
   continue;
  if( un_any_thing( LB0Constraint , el , [](){}() ) )
   continue;
  if( un_any_thing( NNConstraint , el , [](){}() ) )
   continue;
  if( un_any_const_static( el , []( BoxConstraint & b ){} ,
                           un_any_type< BoxConstraint >() ) )
   continue;
  throw( std::logic_error( "unsupported type of static Constraint" ) );
  }

 for( auto & el : f_Block->get_dynamic_constraints() ) {
  if( un_any_thing( std::list< BoxConstraint > , el , [](){}() ) )
   continue;
  if( un_any_thing( std::list< LB0Constraint > , el , [](){}() ) )
   continue;
  if( un_any_thing( std::list< NNConstraint > , el , [](){}() ) )
   continue;
  throw( std::logic_error( "unsupported type of dynamic Constraint" ) );
  }

 // read information about the C05Function - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrFi = v_c05f.size();

 // check if there are "easy" components and deal with them- - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NrEasy = 0;
 if( DoEasy & 1 ) {
  // retrieve the ComputeConfig for the "easy" components, if any
  ComputeConfig * eCC = nullptr;
  if( ! EasyCfg.empty() ) {
   auto cfg = Configuration::deserialize( EasyCfg );
   if( ! ( eCC = dynamic_cast< ComputeConfig * >( cfg ) ) )
    delete cfg;
   }

  IsEasy.resize( NrFi , nullptr );
  auto NEit = NoEasy.begin();
  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( ( NEit != NoEasy.end() ) && ( Index( *NEit ) == k ) ) {
    ++NEit;
    continue;
    }

   auto LagB = dynamic_cast< LagBFunction * >( v_c05f[ k ] );
   if( LagB ) {
    auto MILPs = new MILPSolver();
    try {  // check if the inner Block of the LagBFunction is all-linear
     // do this by trying to register the MILPSolver to the inner Block; if
     // the operation succeeds than the component may be easy (provided that
     // also all variables are continuous), otherwise it surely is not,
     // which is captured by the fact that exception is thrown; note that
     // [MILP]Solver::set_Block() does *not* call Block::register_Solver(),
     // which therefore may have to be done later
     MILPs->set_Block( LagB->get_inner_block() );
     // the component is easy only if all variables are continuous
     if( ! MILPs->get_num_integer_vars() ) {
      IsEasy[ k ] = MILPs;
      ++NrEasy;

      // this is done since OSIMPSolver does not deal with constant term
      constant_value += LagB->get_constant_term();

      // if dynamic updating of the easy component is allowed for any piece
      // of data, register the MILPSolver with the inner Block, as this is not
      // done by [MILP]Solver::set_Block(); note that this calls set_Block()
      // again, which is why it is important that [MILP]Solver::set_Block()
      // check that the Block is the same and ignores it
      if( DoEasy & ~1 )
       LagB->get_inner_block()->register_Solver( MILPs );
      }
     else  // everything is linear, but there are integer variables
      delete MILPs;
     }
    catch( ... ) {  // exception means that something nonlinear is there
     delete MILPs;
     }
    }
   }

  if( ! NrEasy )
   IsEasy.clear();
  else {
   // ComputeConfig-ure the easy components
   if( eCC || ( ! CmpCfg.empty() ) )
    for( Index k = 0 ; k < NrFi ; ++k ) {
     if( ! IsEasy[ k ] )
      continue;

     ComputeConfig * cfg = nullptr;
     if( ( k < Index( CmpCfg.size() ) ) && ( ! CmpCfg[ k ].empty() ) ) {
      auto tcfg = Configuration::deserialize( CmpCfg[ k ] );
      if( ! ( cfg = dynamic_cast< ComputeConfig * >( tcfg ) ) )
       delete tcfg;
      }

     if( ( ! cfg ) && eCC )
      cfg = eCC->clone();

     if( cfg )
      v_c05f[ k ]->set_ComputeConfig( cfg );
     }

   // if easy components can be dynamically changed, attach a FakeSolver to
   // the inner Block of each LagBFunction to record the changes
   if( DoEasy & ~1 ) {
    v_FakeSolver.resize( NrEasy );
    auto FSit = v_FakeSolver.begin();
    for( Index k = 0 ; k < NrFi ; ++k )
     if( IsEasy[ k ] ) {
      auto LagB = static_cast< LagBFunction * >( v_c05f[ k ] );
      *FSit = new FakeSolver();
      LagB->get_inner_block()->register_Solver( *(FSit++) );
      }
    }
   }

  delete eCC;  // TODO: do not clone() the last time

  }  // end( if( DoEasy ) )

 // configure all non-easy components- - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if the non-easy ComputeConfig is provided, apply it.
 //
 // in all cases, set the global pool size, *after* having configured them.
 // this is necessary in that with BPar2 == 0, BundleSolver just takes
 // whatever size of the global pool it finds in the C05Function (that may
 // have been just set by the non-easy ComputeConfig). with BPar2 > 0,
 // instead, BundleSolver ensures that the size of the global pool is *at
 // least* BPar2 by increasing it if it is below. this means that:
 // - BundleSolver never *decreases* the size of the global pool
 // - BundleSolver only uses the first BPar2 linearizations in each global
 //   pool; if there are more, the other ones are ignored
 //
 // meanwhile, also set the accuracy of multipliers
 // MinQuad requires a "high" accuracy to work (1e-12) while standard solvers
 // do not, and in fact may complain if such a tight accuracy is set, but
 // since we don't really trust the accuracy of MPSolver, we give the
 // C05Function more slack
 const double eps = ( MPName & 1 ) ? 1e-7 : 1e-10;

 vBPar2.resize( NrFi + 1, 0 );
 InvItemVcblr.resize( NrFi );

 // retrieve the ComputeConfig for the non-easy components, if any
 ComputeConfig * hCC = nullptr;
 if( ! HardCfg.empty() ) {
  auto cfg = Configuration::deserialize( HardCfg );
  if( ! ( hCC = dynamic_cast< ComputeConfig * >( cfg ) ) )
   delete cfg;
  }

 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( NrEasy && IsEasy[ k ] )
   continue;

  // ComputeConfig-ure the non-easy component
  ComputeConfig * cfg = nullptr;
  if( ( k < Index( CmpCfg.size() ) ) && ( ! CmpCfg[ k ].empty() ) ) {
   auto tcfg = Configuration::deserialize( CmpCfg[ k ] );
   if( ! ( cfg = dynamic_cast< ComputeConfig * >( tcfg ) ) )
    delete tcfg;
   }

  if( ( ! cfg ) && hCC )
   cfg = hCC->clone();

  if( cfg )
   v_c05f[ k ]->set_ComputeConfig( cfg );

  // ensure that the accuracy of multipliers is at least eps
  if( v_c05f[ k ]->get_dbl_par( C05Function::dblAAccMlt ) > eps )
   v_c05f[ k ]->set_par( C05Function::dblAAccMlt , eps );

  // manage the global pool size
  Index gps = v_c05f[ k ]->get_int_par( C05Function::intGPMaxSz );
  if( BPar2 == 0 ) {  // use the current global pool size
   if( gps < 2 )
    throw( std::logic_error( "BPar2 == 0 but too small global pool" ) );
   vBPar2.back() += gps;
   vBPar2[ k ] = gps;
   }
  else {              // force the global pool size to be *at least* BPar2
   if( gps < BPar2 )
    v_c05f[ k ]->set_par( C05Function::intGPMaxSz , int( BPar2 ) );
   vBPar2.back() += BPar2;
   vBPar2[ k ] = BPar2;
   }

  InvItemVcblr[ k ].resize( vBPar2[ k ] , InINF );
  }

 delete hCC;  // TODO: do not clone() the last time

 // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 t = tInit;
 Prevt = INFshift;

 Lambda.resize( NumVar );    // the default starting point
 Lambda1.resize( NumVar );   // the tentative point

 if( MaxSol > 1 )  // best point found so far
  LmbdBst.resize( NumVar );

 OOBase.resize( vBPar2.back() , Inf< SIndex >() );
 // counter for eliminating outdated items: Inf< SIndex >() means empty

 ItemVcblr.resize( vBPar2.back() , make_pair( InINF , InINF ) );

 NrItems.resize( NrFi + 1 , 0 );
 FrFItem.resize( NrFi , 0 );
 MaxItem.resize( NrFi , 0 );

 FreList = {};
 whisZ.resize( NrFi , InINF );
 Zvalid.resize( NrFi , false );

 CurrNrEvls.resize( NrFi , 0 );

 FiStatus.resize( NrFi , kUnEval );
 TrueLB = false;

 UpFiBest = INFshift;      // best, ...
 UpRifFi.resize( NrFi + 1 , 0 );  // and reference Fi() values
 RifeqFi = false;                 // reference values != UpFiLmb
 UpFiLmb1.resize( NrFi + 1 );     // upper and lower function value
 LwFiLmb1.resize( NrFi + 1 );     // ... at the tentative point
 UpFiLmb.resize( NrFi + 1 ,  INFshift );  // upper
 LwFiLmb.resize( NrFi + 1 , -INFshift );  // ... and lower Fi-value
 UpFiLmbdef = LwFiLmbdef = 0;             // ... at the current point
 LowerBound.resize( NrFi + 1 , -INFshift );  // global lower bounds
 f_global_LB = -INFshift;         // algorithmic global LB

 vStar.resize( NrFi + 1 , 0 );
 whisG1.resize( NrFi , InINF );  // no representative yet

 Result = kError;
 SSDone = false;

 // initialize the MP Solver - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Master )        // a MPSolver is set already
  Master->SetDim();  // clear all its internal state
 else {              // a MPSolver is not set yet, create it now
  #if( ! USE_MPTESTER )
   if( MPName & 1 ) {  // the MPSolver is a OSIMPSolver
  #endif
    auto osi_mps = new OSIMPSolver();
    Master = osi_mps;

    if( MPName & 2 ) {
     #if WHICH_OSI_QP == 1
      osi_mps->SetOsi( new OsiCpxSolverInterface() );
     #elif WHICH_OSI_QP == 2
      // externally create the Gurobi environment in order to be able to
      // set it in "silent mode" before it is started
      GRBenv * envP;
      auto err = GRBemptyenv( & envP );
      if( err )
       throw( std::logic_error( "cannot create Gurobi environment" ) );
      GRBsetintparam( envP , "OutputFlag" , 0 );
      err = GRBstartenv( envP );
      if( err )
       throw( std::logic_error( "cannot initialize Gurobi environment" ) );
      osi_mps->SetOsi( new OsiGrbSolverInterface( envP ) );
     #else
      throw( std::invalid_argument( "MPName & 2 not supported" ) );
     #endif
     }
    else
     osi_mps->SetOsi( new OsiClpSolverInterface() );

    osi_mps->SetStabType( MPName & 4 ? OSIMPSolver::quadratic
			             : OSIMPSolver::boxstep );

    osi_mps->SetAlgo( OSIMPSolver::OsiAlg( algo ) ,
		      OSIMPSolver::OsiRed( reduction ) );

    osi_mps->SetThreads( threads );
   #if( ! USE_MPTESTER )
    }
   else {  // the MPSolver is a QPPenaltyMP
  #endif
    QPPenaltyMP *qp = new QPPenaltyMP();
    qp->SetPricing( CtOff );
    qp->SetMaxVarAdd( MxAdd );
    qp->SetMaxVarRmv( MxRmv );
    #if( USE_MPTESTER )
     #if( USE_MPTESTER == 1 )
      Master = new MPTester( Master , qp );
     #else
      Master = new MPTester( qp , Master );
     #endif
    #else
     Master = qp;
    }
    #endif
  }

 InitMP();

 // cleanup MILPSolver data- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now that the MPSolver has read all the information it needs out of the
 // MILPSolver, cleanup all un-necessary data; note that clear_problem() tells
 // which data to delete with exactly the same bit mapping as DoEasy tells
 // which to keep, hence an xor works; the first bit gets zeroed (hence
 // ultimately set to 1) so that the coefficient matrix is always discarded,
 // as changing it is not managed by the MPSolver (and neither really is by
 // MILPSolver, currently)
 //
 // note that if "easy" components are fully static one may want to get rid
 // of the MILPSolver entirely; however, if new variables are added then
 // GetADesc() can be called, which relies on the MILPSolver. getting rid of
 // the MILPSolver would require knowing that the variables set is also
 // static, which would require one parameter; doable, but not now

 if( NrEasy ) {
  const char which = ( DoEasy & ~1 ) ^ 15;
  for( Index k = 0 ; k < NrFi ; ++k )
   if( IsEasy[ k ] )
    IsEasy[ k ]->clear_problem( which );
  }

 // reset algorithm  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note: this has to be done before the next step since it sets Lambda, and
 //       if linearizations are added to the bundle their linearization error
 //       depends on Lambda

 ReSetAlg( RstAlgPrm );  // Lambda is reset inside

 // deal with existing linearizations- - - - - - - - - - - - - - - - - - - - -
 // - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if BPar7 & 8, read from all the C05Function and immediately add to the
 // master problem each and every linearization found in their global pools
 //
 // otherwise read from all the C05Function and mark into InvItemVcblr each
 // and every linearization found in their global pools
 //
 // note that if ( BPar7 & 3 ) >= 2, BundleSolver will happily delete from
 // the global pool any linearization it deletes from the bundle; yet we do
 // not immediately delete existing linearizations from the global pools here.
 // these will likely be overwritten during the optimization, and if memory
 // is a problem they can be cleaned up by the user before set_Block() is
 // called. besides, in many scenarios there will be no linearizations anyway
 // 
 // however, if ( BPar7 & 3 ) == 3 then BundleSolver does not care at all
 // about what linearizations are there in the global pool because it will
 // treat any position in the global pool as available for it regardless to
 // if there is anything there; hence, in this case we do not bother to even
 // look

 if( BPar7 & 8 ) {
  for( Index k = 0 ; k < NrFi ; ++k )
   for( Index i = 0 ; i < vBPar2[ k ] ; ++i )
    if( v_c05f[ k ]->is_linearization_there( i ) )
     add_to_bundle( k , i );
  }
 else
  if( ( BPar7 & 3 ) < 3 ) {
   for( Index k = 0 ; k < NrFi ; ++k )
    for( Index i = 0 ; i < vBPar2[ k ] ; ++i )
     if( v_c05f[ k ]->is_linearization_there( i ) )
      add_to_global_pool( k , i );
   }

 //!! PrintBundle();
 #if CHECK_DS & 1
  CheckBundle();
 #endif
 #if CHECK_DS & 4
  CheckAlpha();
 #endif

 // finally, release the Block - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( BundleSolver::set_Block )

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intMaxIter ):
   if( value < 0 )
    throw( std::invalid_argument( "MaxIter must be >= 0" ) );
   MaxIter = value;
   break;
  case( intMaxSol ):
   if( value < 1 )
    throw( std::invalid_argument( "MaxSol must be >= 1" ) );
   MaxSol = value;
   break;
  case( intEverykIt ):
   EverykIt = value;
   break;
  case( intLogVerb ):
   LogVerb = value;
   break;
  case( intBPar1 ):
   if( value < 1 )
    throw( std::invalid_argument( "BPar1 must be >= 1" ) );
   BPar1 = value;
   break;
  case( intBPar2 ):
   if( value < 2 )
    throw( std::invalid_argument( "BPar2 must be >= 2" ) );
   if( BPar2 == Index( value ) )
    break;
   if( f_Block )
    throw( std::invalid_argument( "changing BPar2 not supported yet" ) );
   BPar2 = value;
   break;
  case( intBPar3 ):
   if( Index( value ) < BPar4 )
    throw( std::invalid_argument( "BPar3 must be >= BPar4" ) );
   BPar3 = value;
   break;
  case( intBPar4 ):
   if( value < 1 )
    throw( std::invalid_argument( "BPar4 must be >= 1" ) );
   BPar4 = value;
   if( BPar4 > BPar3 )
    BPar3 = BPar4;
   break;
  case( intBPar6 ):
   BPar6 = value;
   break;
  case( intBPar7 ):
   BPar7 = value;
   break;
  case( intMnSSC ):
   MnSSC = value;
   break;
  case( intMnNSC ):
   MnNSC = value;
   break;
  case( inttSPar1 ):
   tSPar1 = value;
   break;
  case( intMaxNrEvls ):
   MaxNrEvls = value;
   break;
  case( intDoEasy ):
   DoEasy = char( value & 15 );
   break;
  case( intWZNorm ):
   if( WZNorm != char( value ) ) {
    WZNorm = char( value );
    if( ! ( WZNorm & ~3 ) )  // the easy case, constant
     NrmZFctr = 1;           // factor is known
    else
     NrmZFctr = INFshift;    // factor to be computed
    }
   break;
  case( intFrcLstSS ):
   FrcLstSS = bool( value );
   break;
  case( intTrgtMng ):
   TrgtMng = Index( value );
   break;
  case( intMPName ):
   if( ( value < 0 ) || ( value > 15 ) )
    throw( std::invalid_argument( "MPName must be in [0, 15]" ) );
   MPName = value;
   break;
  case( intMPlvl ):
   MPlvl = value;
   break;
  case( intQPmp1 ):
   MxAdd = value;
   break;
  case( intQPmp2 ):
   MxRmv = value;
   break;
  case( intOSImp1 ):
   if( algo != Index( value ) ) {
    algo = value;
    if( auto osi_mps = dynamic_cast< OSIMPSolver * >( Master ) )
     osi_mps->SetAlgo( OSIMPSolver::OsiAlg( algo ) ,
		       OSIMPSolver::OsiRed( reduction ) );
    }
   break;
  case( intOSImp2 ):
   if( reduction != Index( value ) ) {
    reduction = value;
    if( auto osi_mps = dynamic_cast< OSIMPSolver * >( Master ) )
     osi_mps->SetAlgo( OSIMPSolver::OsiAlg( algo ) ,
		       OSIMPSolver::OsiRed( reduction ) );
    }
   break;
  case( intOSImp3 ):
   if( threads != Index( value ) ) {
    threads = value;
    if( auto osi_mps = dynamic_cast< OSIMPSolver * >( Master ) )
     osi_mps->SetThreads( threads );
    }
   break;
  case( intRstAlg ):
   RstAlgPrm = value;
   break;
  default:
   CDASolver::set_par( par , value );
  }
 }  // end( BundleSolver::set_par( int ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ):
   if( value <= 0 )
    throw( std::invalid_argument( "dblMaxTime must be > 0" ) );
   MaxTime = value;
   break;
  case( dblRelAcc ):
   if( ( value <= 0 ) || ( value >= INFshift ) )
    throw( std::invalid_argument( "RelAcc must be > 0 and finite" ) );
   RelAcc = value;
   break;
  case( dblAbsAcc ):
   if( value <= 0 )
    throw( std::invalid_argument( "AbsAcc must be > 0" ) );
   AbsAcc = value;
   break;
  case( dblEveryTTm ):
   if( value < 0 )
    throw( std::invalid_argument( "EveryTTm must be >= 0" ) );
   EveryTTm = value;
   break;
  case( dblNZEps ):
   if( value < 0 )
    throw( std::invalid_argument( "NZEps must be >= 0" ) );
   NZEps = value;
   break;
  case( dbltStar ):
   tStar = value;
   break;
  case( dblMinNrEvls ):
   MinNrEvls = std::max( double( -1 ) , value );
   break;
  case( dblBPar5 ):
   BPar5 = value;
   break;
  case( dblm1 ):
   if( std::abs( value ) >= 1 )
    throw( std::invalid_argument( "| m1 | must be in (0, 1)" ) );
   m1 = value;
   break;
  case( dblm2 ):
   if( ( value < std::abs( m1 ) ) || ( value >= 1 ) )
    throw( std::invalid_argument( "m2 must be in [ | m1 |, 1)" ) );
   m2 = value;
   break;
  case( dblm3 ):
  if( ( value <= 0 ) || ( value >= 1 ) )
    throw( std::invalid_argument( "m3 must be in (0, 1)" ) );
   m3 = value;
   break;
  case( dblmxIncr ):
  if( value <= 1 )
    throw( std::invalid_argument( "mxIncr must be > 1" ) );
   mxIncr = value;
   break;
  case( dblmnIncr ):
   if( value <= 1 )
    throw( std::invalid_argument( "mnIncr must be > 1" ) );
   mnIncr = std::min( value , mxIncr );
   break;
  case( dblmxDecr ):
   if( ( value <= 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "mxDecr must be in (0, 1)" ) );
   mxDecr = value;
   break;
  case( dblmnDecr ):
   if( ( value <= 0 ) || ( value > 1 ) )
    throw( std::invalid_argument( "mnDecr must be in (0, 1)" ) );
   mnDecr = std::max( value , mxDecr );
   break;
  case( dbltMaior ):
   if( value <= 0 )
    throw( std::invalid_argument( "tMaior must be > 0" ) );
   tMaior = value;
   break;
  case( dbltMinor ):
   if( ( value <= 0 ) || ( value > tMaior ) )
    throw( std::invalid_argument( "tMinor must be in (0, tMaior]" ) );
   tMinor = value;
   break;
  case( dbltInit ):
   if( ( value < tMinor ) || ( value > tMaior ) )
    throw( std::invalid_argument( "tInit must be in [tMinor, tMaior]" ) );
   tInit = value;
   break;
  case( dbltSPar2 ):
   if( value <= 0 )
    throw( std::invalid_argument( "tSPar2 must be > 0" ) );
   tSPar2 = value;
   break;
  case( dbltSPar3 ):
   tSPar3 = std::abs( value ) > 1 ? value : 0;
   break;
  case( dblCtOff ):
   if( value < 0 )
    throw( std::invalid_argument( "CtOff must be >= 0" ) );
   CtOff = value;
   break;
  default:
   CDASolver::set_par( par , value );
  }
 }  // end( BundleSolver::set_par( double ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par , std::string && value )
{
 switch( par ) {
  case( strEasyCfg ):
   EasyCfg = std::move( value );
   break;
  case( strHardCfg ):
   HardCfg = std::move( value );
   break;
  default:
   CDASolver::set_par( par , value );
  }
 }  // end( BundleSolver::set_par( std::string && ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par , std::vector< int > && value )
{
 if( par == vintNoEasy )
  NoEasy = std::move( value );
 else
  CDASolver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::set_par( idx_type par ,
			    std::vector< std::string > && value )
{
 if( par == vstrCmpCfg )
  CmpCfg = std::move( value );
 else
  CDASolver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::set_log( std::ostream * log_stream )
{
 f_log = log_stream;
 if( Master )
  Master->SetMPLog( f_log , MPlvl );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::get_dual_solution( Configuration * solc )
{
 // construct the important linearization for each non-easy component (unless
 // it is already there, and signal to the C05Functions which one it is

 for( Index k = 0 ; k < NrFi ; ++k )         // for all components
  if( ( ! NrEasy ) || ( ! IsEasy[ k ] ) ) {  // but skip the easy ones
   C05Function::LinearCombination lc;

   if( Zvalid[ k ] ) {
    // the optimal aggregated linearization for component k is in the
    // bundle (and, therefore, global pool) already: the optimal
    // coefficients are very simple, it's just that one
    lc.resize( 1 );
    lc[ 0 ].first = ItemVcblr[ whisZ[ k ] ].second;
    lc[ 0 ].second = 1;
    }
   else {
    // retrieve optimal multipliers from the Master
    Index MBDm;
    cIndex_Set MBse;
    cHpRow Mlt = Master->ReadMult( MBse , MBDm , k + 1 , false );

    // copy them in the LinearCombination
    lc.resize( MBDm );
    auto lcit = lc.begin();

    if( MBse )
     for( Index h ; ( h = *(MBse++) ) < InINF ; ) {
      lcit->first = ItemVcblr[ h ].second;
      (lcit++)->second = *(Mlt++);
      }
    else
     for( Index h = 0 ; h < MBDm ; ) {
      lcit->first = ItemVcblr[ h++ ].second;
      (lcit++)->second = *(Mlt++);
      }
    }

   v_c05f[ k ]->set_important_linearization( std::move( lc ) );

   }  // end( if( not easy ) )

 }  // end( BundleSolver::get_dual_solution() )

/*--------------------------------------------------------------------------*/

int BundleSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intMaxIter ):   return( MaxIter );
  case( intMaxSol ):    return( MaxSol );
  case( intEverykIt ):  return( EverykIt );
  case( intLogVerb ):   return( LogVerb );
  case( intBPar1 ):     return( BPar1 );
  case( intBPar2 ):     return( BPar2 );
  case( intBPar3 ):     return( BPar3 );
  case( intBPar4 ):     return( BPar4 );
  case( intBPar6 ):     return( BPar6 );
  case( intBPar7 ):     return( BPar7 );
  case( intMnSSC ):     return( MnSSC );
  case( intMnNSC ):     return( MnNSC );
  case( inttSPar1 ):    return( tSPar1 );
  case( intMaxNrEvls ): return( MaxNrEvls );
  case( intDoEasy ):    return( DoEasy );
  case( intWZNorm ):    return( WZNorm );
  case( intFrcLstSS ):  return( FrcLstSS );
  case( intTrgtMng ):   return( TrgtMng );
  case( intMPName ):    return( MPName );
  case( intMPlvl ):     return( MPlvl );
  case( intQPmp1 ):     return( CtOff );
  case( intQPmp2 ):     return( MxRmv );
  case( intOSImp1 ):    return( algo );
  case( intOSImp2  ):   return( reduction );
  case( intOSImp3 ):    return( threads  );
  case( intRstAlg ):    return( RstAlgPrm  );
  default:              return( CDASolver::get_int_par( par ) );
  }
 }  // end( BundleSolver::get_int_par )

/*--------------------------------------------------------------------------*/

double BundleSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):   return( MaxTime );
  case( dblRelAcc ):    return( RelAcc );
  case( dblAbsAcc ):    return( AbsAcc );
  case( dblEveryTTm ):  return( EveryTTm );
  case( dblNZEps ):     return( NZEps );
  case( dbltStar ):     return( tStar );
  case( dblMinNrEvls ): return( MinNrEvls );
  case( dblBPar5 ):     return( BPar5 );
  case( dblm1 ):        return( m1 );
  case( dblm2 ):        return( m2 );
  case( dblm3 ):        return( m3 );
  case( dblmxIncr ):    return( mxIncr );
  case( dblmnIncr ):    return( mnIncr );
  case( dblmxDecr ):    return( mxDecr );
  case( dblmnDecr ):    return( mnDecr );
  case( dbltMaior ):    return( tMaior );
  case( dbltMinor ):    return( tMinor );
  case( dbltInit ):     return( tInit );
  case( dbltSPar2 ):    return( tSPar2 );
  case( dbltSPar3 ):    return( tSPar3 );
  case( dblCtOff ):     return( CtOff );
  default:              return( CDASolver::get_dbl_par( par ) );
  }
 }  // end( BundleSolver::get_dbl_par )

/*--------------------------------------------------------------------------*/

const std::string & BundleSolver::get_str_par( idx_type par ) const
{
 switch( par ) {
  case( strEasyCfg ):   return( EasyCfg );
  case( strHardCfg ):   return( HardCfg );
  default:              return( CDASolver::get_str_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

const std::vector< int > & BundleSolver::get_vint_par( idx_type par ) const
{
 if( par == vintNoEasy )
  return( NoEasy );

 return( CDASolver::get_vint_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > & BundleSolver::get_vstr_par( idx_type par )
 const
{
 if( par == vstrCmpCfg )
  return( CmpCfg );

 return( CDASolver::get_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/
/*----------- METHODS FOR HANDLING THE State OF THE BundleSolver -----------*/
/*--------------------------------------------------------------------------*/

State * BundleSolver::get_State( void ) const {
  return( new BundleSolverState( this ) );
  }

/*--------------------------------------------------------------------------*/

void BundleSolver::put_State( const State & state )
{
 // if state is not a BundleSolverState &, exception will be thrown
 const auto & s = dynamic_cast< const BundleSolverState & >( state );

 guts_of_put_State( s );

 Lambda = s.Lambda;

 for( Index i = 0 ; i < NrFi ; ++i )
  if( s.v_comp_State[ i ] )
   v_c05f[ i ]->put_State( *(s.v_comp_State[ i ]) );

 }  // end( BundleSolver::put_State( const & ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::put_State( State && state )
{
 // if state is not a BundleSolverState &&, exception will be thrown
 auto && s = dynamic_cast< BundleSolverState && >( state );

 guts_of_put_State( s );

 Lambda = std::move( s.Lambda );

 for( Index i = 0 ; i < NrFi ; ++i )
  if( s.v_comp_State[ i ] ) {
   v_c05f[ i ]->put_State( std::move( *(s.v_comp_State[ i ]) ) );
   delete s.v_comp_State[ i ];
   }

 s.v_comp_State.clear();

 }  // end( BundleSolver::put_State( && ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::serialize_State( netCDF::NcGroup & group ,
				    const std::string & sub_group_name ) const
{
 if( ! sub_group_name.empty() ) {
  auto gr = group.addGroup( sub_group_name );
  serialize_State( gr );
  return;
  }

 // do it "by hand" since there is no BundleSolverState available to call
 // State::serialize() from
 group.putAtt( "type", "BundleSolverState" );

 auto nv = group.addDim( "BundleSolver_NumVar" , NumVar );

 ( group.addVar( "BundleSolver_Lambda" , netCDF::NcDouble() , nv ) ).putVar(
			               { 0 } , {  NumVar } , Lambda.data() );

 ( group.addVar( "BundleSolver_t" , netCDF::NcDouble() ) ).putVar( &t );

 auto nfi = group.addDim( "BundleSolver_NrFi" , NrFi + 1 );

 if( UpFiLmbdef ) {
  group.addDim( "BundleSolver_UpFiLmbdef" , UpFiLmbdef );
  ( group.addVar( "BundleSolver_UpFiLmb" , netCDF::NcDouble() , nfi )
    ).putVar( { 0 } , { NrFi + 1 } , UpFiLmb.data() );
  }

 if( LwFiLmbdef ) {
  group.addDim( "BundleSolver_LwFiLmbdef" , LwFiLmbdef );
  ( group.addVar( "BundleSolver_LwFiLmb" , netCDF::NcDouble() , nfi )
    ).putVar( { 0 } , { NrFi + 1 } , LwFiLmb.data() );
  }

 if( Fi0Lmb != 0 )
  ( group.addVar( "BundleSolver_Fi0Lmb" , netCDF::NcDouble() )
    ).putVar( & Fi0Lmb );

 if( f_global_LB > - INFshift )
  ( group.addVar( "BundleSolver_global_LB" , netCDF::NcDouble() )
    ).putVar( & f_global_LB );

 for( Index i = 0 ; i < NrFi ; ++i ) {
  if( NrEasy && IsEasy[ i ] )
   continue;

  v_c05f[ i ]->serialize_State( group ,
				"Component_State_" + std::to_string( i ) );
  }
 }  // end( BundleSolver::serialize_State )

/*--------------------------------------------------------------------------*/
/*----------------------- OTHER PROTECTED METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::guts_of_put_State( const BundleSolverState & state )
{
 if( Result == kStillRunning )
  throw( std::logic_error(
	"BundleSolver::put_State() called within BundleSolver::compute()" ) );

 if( ( NrFi != state.NrFi ) || ( NumVar != state.NumVar ) )
  throw( std::invalid_argument(
			  "BundleSolver::put_State(): inconsistent State" ) );

 if( t != state.t ) {
  t = state.t;
  tHasChgd = true;
  }

 if( state.UpFiLmbdef ) {
  UpFiLmbdef = state.UpFiLmbdef;
  UpFiLmb = state.UpFiLmb;
  }
 else {
  UpFiLmbdef = 0;
  std::fill( UpFiLmb.begin() , UpFiLmb.end() ,  INFshift );
  }

 if( state.LwFiLmbdef ) {
  LwFiLmbdef = state.LwFiLmbdef;
  LwFiLmb = state.LwFiLmb;
  }
 else {
  LwFiLmbdef = 0;
  std::fill( LwFiLmb.begin() , LwFiLmb.end() ,  -INFshift );
  }

 // if Lambda has changed, the master problem need be informed
 if( Lambda != state.Lambda ) {
  Vec_VarValue foo( NrFi + 1 , 0 );

  if( UpFiLmbdef == NrFi + 1 ) {
   // the function value is known, so it has to be passed to the master:
   // see the comments inside GotoLambda1() for the rationale of the
   // curios definition
   foo.front() = UpFiLmb.back() - UpRifFi.back();
   std::transform( UpFiLmb.begin() , --(UpFiLmb.end()) , UpRifFi.begin() ,
		   ++(foo.begin()) , std::minus< double >() );
   UpRifFi = UpFiLmb;
   RifeqFi = true;
   }
  // else no change in the (unknown) f-values, so all-0 is OK

  Master->ChangeCurrPoint( state.Lambda.data() , foo.data() );
  }

 Fi0Lmb = state.Fi0Lmb;
 f_global_LB = state.global_LB;
 
 }  // end( BundleSolver::guts_of_put_State )

/*--------------------------------------------------------------------------*/

void BundleSolver::FormD( void )
{
 // initialize the Master Problem Solver- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // change/set t as required- - - - - - - - - - - - - - - - - - - - - - - - -

 // Special treatment of the "empty Master Problem" case: there are no
 // subgradients, so this is only a feasibility problem, which should give
 // a feasible point as close as possible to the starting one. This is
 // "free" with some stabilizing terms (e.g. the quadratic one), but not
 // necessarily so with others (e.g. the trust region). In order to "stay
 // as close as possible", t is temporarily decreased to its minimum value.
 // As soon as there is something in the bundle, the current value of t is
 // restored (Prevt is used to hold it).

 if( Master->BCSize() >= Master->BSize() ) {
  if( ( t > tMinor ) && ( Prevt == INFshift ) ) {
   Prevt = t;
   t = tMinor;
   tHasChgd = true;
   }
  }
 else
  if( Prevt < INFshift ) {
   if( t != Prevt ) {
    t = Prevt;
    tHasChgd = true;
    }
   Prevt = INFshift;
   }

 if( tHasChgd ) {
  Master->Sett( t );
  tHasChgd = false;
  }

 // collect and set individual lower bounds, if any - - - - - - - - - - - - -
 // if the MPSolver accepts them, collect and if necessary set the individual
 // lower bounds. note that if all of them are finite and the 0-th component
 // is not there their sum would give an alternative valid global lower bound.
 // however, the same information is already encoded in the individual bounds,
 // hence it's of no use. the exception is that QPPenaltyMP does not allow
 // individual lower bounds, but this is not a permanent issue and it'll go
 // away when we'll get rid of MPSolver; besides, it very unlikely to ever
 // really happen

 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( NrEasy && IsEasy[ k ] )  // skip easy components
   continue;

  // get the lower bound out of the C05Function
  auto LwrBndk = f_convex ?   v_c05f[ k ]->get_global_lower_bound()
                          : - v_c05f[ k ]->get_global_upper_bound();

  // use it to update the lower estimate in Lambda
  update_LwFiLambd( k , LwrBndk );

  if( LwrBndk != LowerBound[ k ] ) {  // if it has changed
   LowerBound[ k ] = LwrBndk;         // record the new value

   #if( ! USE_MPTESTER )
    // QPPenaltyMP does not allow individual lower bounds, and if a MPTester
    // is used then a QPPenaltyMP is involved anyway
    if( MPName & 1 ) {
     if( LwrBndk > -INFshift )  // translate it w.r.t. UpRifFi
      LwrBndk -= UpRifFi[ k ];

     // pass it to the master problem
     Master->SetLowerBound( LwrBndk , k + 1 );
     }
   #endif
   }
  }

 // collect and set the global lower bound, if any- - - - - - - - - - - - - -
 // note: this used to be done elsewhere, in particular after each function
 //       evaluation. the rationale was that in the Lagrangian case one would
 //       do heuristics as a part of the computation, and these could produce
 //       a better lower bound that one may immediately want to check. but
 //       the truth is that one does not: the real way in which a lower bound
 //       is useful is when it enters the master problem and helps in getting
 //       the optimality conditions. thus, the right place to check and
 //       update the lower bound(s) is right before the master problem is
 //       solved. the conditional lower bound ( TrueLB == false ) is indeed
 //       useful right after the function computation to prove unboundedness,
 //       but it is available then, and anyway it typically does not change
 //       when the function is computed, unlike the "hard" one

 auto LwrBnd = f_convex ?   f_Block->get_valid_lower_bound( false )
                        : - f_Block->get_valid_upper_bound( false );

 if( LwrBnd != LowerBound.back() ) {
  if( TrueLB ) {  // if the global lower bound was a "true" one, its value
   // is "baked in" the total lower and upper estimate of the function value
   // in Lambda: ensure it is recomputed. this works both if the bound was
   // there and it is changed and if it is reset 
   if( UpFiLmbdef > NrFi ) {
    --UpFiLmbdef;
    UpFiLmb.back() = INFshift;
    }
   if( LwFiLmbdef > NrFi ) {
    --LwFiLmbdef;
    LwFiLmb.back() = -INFshift;
    }
   }

  TrueLB = ( LwrBnd > -INFshift );  // if it's finite it's a true LB

  // if a true global lower bound was not there and it is set, then ensure
  // it will be properly "baked in" the total lower and upper estimate
  if( TrueLB ) {
   if( UpFiLmbdef > NrFi ) {
    --UpFiLmbdef;
    UpFiLmb.back() = INFshift;
   }
   if( LwFiLmbdef > NrFi ) {
    --LwFiLmbdef;
    LwFiLmb.back() = -INFshift;
    }
   }

  // if the total upper estimate in Lambda needs be recomputed, do it now
  // the only role of the total upper estimate in Lambda in the master
  // problem is to translate the global lower bound (if any); this will
  // be done during the final call to SetLowerBound()
  if( UpFiLmbdef == NrFi ) {
   ++UpFiLmbdef;  // all components + the sum computed
   UpFiLmb.back() = std::accumulate( UpFiLmb.begin() , --(UpFiLmb.end()) ,
				     Fi0Lmb );
   // note that here the global lower bound (if any) is "baked in" the total
   // upper function estimate
   if( TrueLB && ( UpFiLmb.back() < LwrBnd ) )
    UpFiLmb.back() = LwrBnd;
   }

  // if the total lower estimate in Lambda needs be recomputed, do it now
  if( LwFiLmbdef == NrFi ) {
   ++LwFiLmbdef;  // all components + the sum computed
   LwFiLmb.back() = std::accumulate( LwFiLmb.begin() , --(LwFiLmb.end()) ,
				     Fi0Lmb );
   // note that here the global lower bound (if any) is "baked in" the total
   // lower function estimate
   if( TrueLB && ( LwFiLmb.back() < LwrBnd ) )
    LwFiLmb.back() = LwrBnd;
   }

  LowerBound.back() = LwrBnd;        // in all cases, record it
  if( TrueLB ) {   // if the bound value is finite
   // translate it using the reference value of the hard components. the
   // "easy" components are not translated, but the non-easy ones are,
   // hence the global lower bound has to be translated by the contribution
   // of the non-easy components (and the linear one, that somewhat
   // un-intuitively is treated in the same way). this *would* be the sum of
   // their reference value *if* the (previous) global lower bound had not
   // impacted the total reference value. to avoid checking it (and because
   // it likely is faster) we compute the correction term by subtracting the
   // value of the "easy" components by the total reference value
   auto rf = UpRifFi.back();
   if( NrEasy )
    for( Index k = 0 ; k < NrFi ; ++k )
     if( IsEasy[ k ] )
      rf -= UpRifFi[ k ];

   LwrBnd -= rf;
   }

  // set it in the master problem, translated if it's finite
  // note that the bound has to be set in the master problem even if it is
  // -INF, because before it was not so, hence it has to be reset
  Master->SetLowerBound( LwrBnd );

  }  // end( if( the global lower bound has changed ) )

 if( ! TrueLB )  // if no true LB, see if "conditional" one is there
  LowerBound.back() = f_convex
                      ?   f_Block->get_valid_lower_bound( true )
                      : - f_Block->get_valid_upper_bound( true );

 for( ; ; )  // error-handling loop - - - - - - - - - - - - - - - - - - - - - -
 {        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // ensure the MPSolver will not take too much time
  if( MaxTime < INFshift ) {
   Master->SetMPTime();
   Master->SetPar( MPSolver::kMaxTme , MaxTime - get_elapsed_time() );
   }

  MPSolver::MPStatus mps = Master->SolveMP();  // solve the MP

  if( mps == MPSolver::kOK )        // everything's alright
   break;

  if( mps == MPSolver::kUnfsbl ) {  // the MP is empty
   if( ! Master->BCSize() )         // but no vertical linearizations
    mps = MPSolver::kError;         // it must be a numerical error
   else {                           // there are vertical linearizations
    Result = kInfeasible;           // the MP can really be infeasible
    return;                         // nothing else to do
    }
   }

  if( mps == MPSolver::kUnbndd )  // the MP is unbounded
   /* With a "loosely stabilised" MPSolver (one with, say, a 1-norm or
    * INF-norm primal penalty), this may happen and it could be mended
    * by decreasing t, which would lead to the following piece of code:
    *
    * if( ( t <= tMinor ) || ( Master->BCSize() >= Master->BSize() ) ) {
    *  // ... but t must always be >= tMinor, and it is already == tMinor
    *  // in the "empty" case of the initial iteration with empty bundle
    *  BLOG( 1 , std::endl << "Bundle::FormD: failure in MPSolver." );
    *  Result = kError;
    *  return;
    *  }
    *
    * BLOG( 1 , std::endl << "Bundle::FormD: MP unbounded, decreasing t" );
    * Master->Sett( t = std::max( t / 2 , tMinor ) );
    * continue;
    *
    * However, currently BundleSolver only admits QPPenaltyMP or OSiMPSolver
    * with either Quadratic or BoxStep stabilisation, which can never be
    * unbounded: hence, the MPSolver returning kUnbndd can only be a
    * numerical error in disguise. */
   mps = MPSolver::kError;

  if( mps == MPSolver::kStppd ) {  // stopped by time limit
   //!! so far, the time limit in the MPSolver is only due to the global time
   //!! limit in the NDOSolver, but one day we may want to set it
   //!! independently; then, some checks will have to be done if the
   //!! solution is feasible and it can still be used and v is sufficiently
   //!! < 0: in this case we can use the direction as well, otherwise we
   //!! have to give the MPSolver more time
   Result = kStopTime;
   break;
   }

  // mps == MPSolver::kError, i.e., there has been a numerical problem in- -
  // the MP Solver; it's not yet time to despair, as by eliminating items- -
  // it may be possible to solve the problem - - - - - - - - - - - - - - - -

  BLOG( 2 , std::endl << "Bundle::FormD: error in MP, emergency delete" );

  Index MBDm;
  cIndex_Set MBse;
  cHpRow Mlt = Master->ReadMult( MBse , MBDm );
  Index i = InINF;

  // the last *removable* item in Base is eliminated - - - - - - - - - - - -

  if( MBse ) {
   for( ; MBDm-- ; )
    if( ( OOBase[ MBse[ MBDm ] ] >= 0 ) &&
        ( Mlt[ MBDm ] >= Eps< HpNum >() ) ) {
     i = MBse[ MBDm ];
     break;
     }
   }
  else
   for( ; MBDm-- ; )
    if( ( OOBase[ MBDm ] >= 0 ) && ( Mlt[ MBDm ] >= Eps< HpNum >() ) ) {
     i = MBDm;
     break;
     }

  if( i == InINF )  // there are no *removable* items in Base - - - - - - - -
  for( Index j = Master->MaxName() ; j-- ; )  // pick any removable item
    if( ( OOBase[ j ] >= 0 ) && ( OOBase[ j ] < Inf< SIndex >() ) ) {
     i = j;
     break;
     }

  if( i == InINF ) {  // there are no removable items at all - - - - - - - -
   BLOG( 1 , std::endl << "Bundle::FormD: unrecoverable MP failure." );
   Result = kError;
   return;
   }

  // i can be deleted, do it and try again hoping this will solve the
  // numerical issue in the master problem

  if( ( BPar7 & 3 ) == 3 ) {
   // if BPar7 tells that the removal must be "hard", actually remove the
   // linearization from the global pool, as Delete() does not do that
   inhibit_Modification( true );
   v_c05f[ ItemVcblr[ i ].first ]->delete_linearization(
						   ItemVcblr[ i ].second );
   inhibit_Modification( false );
   }
  Delete( i );

  }  // end ( error-handling loop )- - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // MinQuad (in QPPenaltyMP) has the habit to increase eR when things go
 // wrong, but never reset it to the original value. over long sequences
 // of calls this may make eR to become "big", which in turn ultimately
 // reduces accuracy and created problems. ensure that eR is properly reset
 // after every successful call
 if( ! ( MPName & 1 ) )
  Master->SetPar( MPSolver::kOptEps , 1e-12 );

 Sigma = Master->ReadSigma();              // read Sigma*
 vStar.back() = Master->ReadFiBLambda();   // read the total v*
 // v* is the predicted decrease in  Lambda1 w.r.t. the value in Lambda;
 // however, if any (non-easy) component does not have any subgradient
 // in the bundle this value is not well-defined (the master problem is
 // added an artificial constraint to make v[ k ] bounded) and +INF is
 // returned for that component, and therefore for the total v* (the sum)

 // now retrieve vStar[ k ] for each component: for easy ones, the master
 // problem produces the *exact* Fi-value (up = lw) at Lambda1, which we
 // store unmodified in vStar[ k ] (that is not used anyway) for it to be
 // retrieved later by FormLambda1()
 for( Index k = 0 ; k < NrFi ; ++k )
  vStar[ k ] = Master->ReadFiBLambda( k + 1 );  // read model value

 if( NrEasy ) {  // there are easy components
  // adjust the contribution of the easy components to v* and Sigma*
  // easy components are treated differently from hard ones in that their
  // value in the model is *not* translated by their (reference) Fi-value in
  // Lambda; however, v* has to measure the total decrease predicted by the
  // model in Lambda1 w.r.t. the (reference) Fi-value in Lambda, and Sigma*
  // the error, again, w.r.t. the (reference) Fi-value in Lambda. since the
  // MPSolver provides non-translated Fi-values, the correction has to be
  // made here.
  VarValue EasyRifFi = 0;
  for( Index k = 0 ; k < NrFi ; ++k )
   if( IsEasy[ k ] )
    EasyRifFi += UpRifFi[ k ];

  // note that v* and Sigma* have "opposite signs" (the former is negative,
  // the latter positive) which justifies why one finds a "-=" and a "+="
  if( vStar.back() < INFshift )
   vStar.back() -= EasyRifFi;
  Sigma += EasyRifFi;
  }

 DSTS = Master->ReadDStart( std::abs( tStar ) );  // D_{t*}( z* )

 // Sigma* + D*_{t*}( -z* ) is the "maximum expected increase" used in
 // the stopping criterion, EpsU is that relative to Fi( Lambda )
 if( UpFiLmb.back() < INFshift )
  EpsU = ( DSTS + Sigma ) / std::max( std::abs( UpFiLmb.back() ) ,
				      double( 1 ) );
 else
  EpsU = 1;  // ensure EpsU is initialized somehow

 Zvalid.assign( NrFi , false );    // the z[ i ] are no longer valid

 DST = Master->ReadDStart( t );  // D_t( z* )

 // Delta* = D_t( z* ) + Sigma* is <= - v*, and a weaker requirement about
 // how much the (total) function must increase for a NS to be declared;
 // however, this only holds if v* is "true", which means that all the
 // components have some diagonal linearization in their bundle; otherwise
 // v* is "fake" and so is Delta* (but anyway, this means that any finite
 // lower bound is much better than what we currently have)

 // compute (very easy if the stabilization is quadratic) || d* ||_2
 if( ( ! ( MPName & 1 ) ) || ( MPName & 4 ) )  // quadratic stabilization
  // ReadDStart( t ) == t || z* ||_2^2 / 2   and   d = - t z*   ==>
  // || d* ||_2 == t || z* ||_2 == t * sqrt( 2 * DST / t )
  NrmD = t * sqrt( 2 * DST / t );
 else {                                        // boxstep stabilization
  auto tdir = Master->Readd();
  NrmD = 0;                                    // || d* ||_2 need be computed
  for( Index i = 0 ; i < NumVar ; ++i )
   NrmD += tdir[ i ] * tdir[ i ];
  NrmD = sqrt(  NrmD );
  }

 // in the easy case NrmD also gives the (2-)norm of z*
 // important note: the relationship d* = - t z* upon which the following is
 // based is actually NOT VALID WHEN THERE ARE CONSTRAINTS, if z* is to be
 // interpreted as the aggregate subgradient of the objective. however, it
 // IS VALID IF z* IS TO BE INTERPRETED AS THE AGGREGATE SUBGRADIENT OF THE
 // ESSENTIAL OBJECTIVE (f + i_X), WHICH IS EXACTLY WHAT IS NEEDED HERE.
 // in fact, consider the case where the constraints are just Lambda >= 0;
 // what we have is that d* = t proj_{>= 0}( - z* ); but it is exactly the
 // condition proj_{>= 0}( - z* ) == 0 that gives the optimality condition.
 // in the Lagrangian case, \Lambda >= 0 corresponds to having relaxed
 // inequalities A u <= b, and z* = b - A u*; hence, it is only required that
 // z* == 0 (that is, proj_{>= 0}( - z* ) == 0) to ensure that u* is feasible
 // and therefore Sigma*-optimal, as opposed to requiring z* == 0
 if( ( ( ! ( MPName & 1 ) ) || ( MPName & 4 ) ) && ( ( WZNorm & 3 ) == 2 ) )
  NrmZ = NrmD / t;
 else {  // otherwise it has to be computed the hard way
  // NOTE: THIS CODE IS BOTH HORRIBLY INEFFICIENT DUE TO A CRAP IMPLEMENTATION
  // OF OsiMPSolver::ReadZ AND INCORRECT WHEN THERE ARE CONSTRAINTS, AS THE
  // z* COMPUTED BY ReadZ() IS THAT OF THE OBJECTIVE BUT NOT OF THE ESSENTIAL
  // OBJECTIVE. the right vector should be easy to compute since it's basically
  // the slack s that is explicit in the dual formulation of the master problem,
  // adjusted with the slacks, but this is no time to dawdle with this
  Index dim;
  const Index * nms;
  std::vector< double > tZ( NumVar );
  Master->ReadZ( tZ.data() , nms , dim );
  tZ.resize( dim );
  NrmZ = ::norm( tZ , WZNorm & 3 );
  }

 // if still needed, compute the scaling factor for z*
 if( NrmZFctr == INFshift )
  compute_NrmZFctr();

 // if the scaling factor could be computed one can check if z* == 0
 // has happened and declare a globally valid LB
 if( ( UpFiLmb.back() < INFshift ) && ( vStar.back() < INFshift ) &&
     ( NrmZFctr < INFshift ) && ( NrmZ <= NrmZFctr * NZEps ) ) {
  if( f_global_LB < UpFiLmb.back() + vStar.back() )
   f_global_LB = UpFiLmb.back() + vStar.back();
  }
 }  // end( BundleSolver::FormD )

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtCntrs( void )
{
 // increase all the OOBase[] counters but those == +/-Inf< SIndex >() - - - - -
 // items whose OOBase[] becomes 0 (e.g. the newly entered items, which have
 // OOBase[] == -1) are set to +1, in such a way that only the items in the
 // optimal base have OOBase[] == 0; note that the converse is not true, as
 // items in the optimal base may have OOBase[] < 0 instead

 for( auto OOit = OOBase.begin() ;
      OOit != OOBase.begin() + Master->MaxName() ; ++OOit )
  if( ( *OOit < Inf< SIndex >() ) && ( *OOit > -Inf< SIndex >() ) ) {
   ++(*OOit);
   if( ! *OOit )
    ++(*OOit);
   }

 // set to 0 the OOBase[] counter for items in base (if not < 0)- - - - - - -
 // note that chechking if the multiplier is strictly positive should be
 // redundant, if one was trusting the MPSolver
 // MinQuad requires a "high" accuracy to work (1e-12) while standard solvers
 // do not, and in fact may complain if such a tight accuracy is set
 double eps = ( MPName & 1 ) ? 1e-9 : 1e-12;

 Index MBDim;
 const Index * MBse;
 const double * Mlt = Master->ReadMult( MBse , MBDim );
 if( MBse ) {
  for( Index i ; ( i = *(MBse++) ) < InINF ; ++Mlt )
   if( ( *Mlt >= eps ) && ( OOBase[ i ] > 0 ) )
    OOBase[ i ] = 0;
  }
 else
  for( Index i = 0 ; i < MBDim ; ++i , ++Mlt )
   if( ( *Mlt >= eps ) && ( OOBase[ i ] > 0 ) )
    OOBase[ i ] = 0;

 /*!!
 // note that there is a case in which a component wFi has Z[ wFi ] "for free"
 // in the bundle: this is when wFi only has *one* subgradient in base (or, in
 // practice, a subgradient with multiplier very close to one). This is
 // checked here (it is basically for free), and in case whisZ[] is properly
 // set so as to avoid pointless aggregations and OOBase[] is set to -1,
 // because under no circumstances such a subgradient can ever be removed
 // from the bundle
 //
 // it now seems to me that this is stupid, since if there is only one
 // subgradient in base no aggregation is ever performed; the only issue
 // is if the base is all (but one) taken by constraints, but then even
 // aggregating does not help

 if( MBse ) {
  for( Index i ; ( i = *(MBse++) ) < InINF ; Mlt++ )
   if( *Mlt >= Eps< HpNum >() ) {
    if( ( *Mlt >= 1 - RMPAccSol ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) - 1 ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }
  }
 else
  for( Index i = 0 ; i < MBDim ; i++ , Mlt++ )
   if( *Mlt >= Eps< double >() ) {
    if( ( *Mlt >= 1 - RAccSol ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) - 1 ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }
    !!*/

 }  // end( UpdtCntrs ) - - - - - - - - - - - - - - - - - - - - - - - - - - -

/*--------------------------------------------------------------------------*/

void BundleSolver::FormLambda1( HpNum Tau )
{
 Master->MakeLambda1( Lambda.data() , Lambda1.data() , Tau );

 if( Master->NumBxdVars() ) {
  // as the relative precision required to the MPSolver is not enough to
  // ensure that the bounds on the variables will be satisfied with the
  // precision required by the FiOracle, the (upper and lower) bounds are
  // strictly enforced here
  //
  //!! this part either to be updated with the bounds from the
  //   OneVarConstraint or, more likely, to be completely removed

  std::vector< VarValue > tL1 = Lambda1;

  if( Master->NumNNVars() )             // there are NN vars and UB vars
   if( Master->NumNNVars() == NumVar )  // actually, all variables are NN
    for( Index i = 0 ; i < NumVar ; ++i ) {
     if( tL1[ i ] < 0 )
      tL1[ i ] = 0;

     const double UBh = LamVcblr[ i ]->get_ub();
     if( tL1[ i ] > UBh )
      tL1[ i ] = UBh;
     }
   else                                 // not all variables are NN
    for( Index i = 0 ; i < NumVar ; ++i ) {
     if( Master->IsNN( i ) && ( tL1[ i ] < 0 ) )
      tL1[ i ] = 0;

     const double UBh = LamVcblr[ i ]->get_ub();
     if( tL1[ i ] > UBh )
      tL1[ i ] = UBh;
     }
  else  // there are only UB vars
   for( Index i = 0 ; i < NumVar ; ++i ) {
    const double UBh = LamVcblr[ i ]->get_ub();
    if( tL1[ i ] > UBh )
     tL1[ i ] = UBh;
    }

  Lambda1 = tL1;

  }  // end( if( the bounds have to be enforced ) )

 // move the value from Lambda1 to the ColVariable - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiStatus.assign( NrFi , kUnEval );
 whisG1.assign( NrFi , InINF );

 for( Index i = 0 ; i < NumVar ; i++ )
  LamVcblr[ i ]->set_value( Lambda1[ i ] );

 // compute the upper and lower model at the tentative point - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 UpFiLmb1def = LwFiLmb1def = 0;
 for( Index k = 0 ; k < NrFi ; ++k )
  if( NrEasy && IsEasy[ k ] ) {  // k is an easy component
   UpFiLmb1[ k ] = LwFiLmb1[ k ] = vStar[ k ];  // we know the exact value
   ++UpFiLmb1def;                               // thus both UpFi1 and LwFi1
   ++LwFiLmb1def;                               // are known (and equal)
   }
  else {                         // k is a hard component
   // compute upper and lower bound for k (possibly +/- INF)

   UpFiLmb1[ k ] = INFshift;

   // computing the upper bound requires the Lipschitz constant and is
   // only done if bit 4 of TrgtMng == 1
   if( ( TrgtMng & 16 ) && ( UpFiLmb[ k ] < INFshift ) ) {
    c_VarValue Lk = v_c05f[ k ]->get_Lipschitz_constant();
    if( Lk < INFshift ) {
     UpFiLmb1[ k ] = UpFiLmb[ k ] + Lk * NrmD;
     ++UpFiLmb1def;
     }
    }

   if( vStar[ k ] < INFshift ) {
    LwFiLmb1[ k ] = std::max( UpRifFi[ k ] + vStar[ k ] , LowerBound[ k ] );
    ++LwFiLmb1def;
    }
   else
    LwFiLmb1[ k ] = -INFshift;
   }

 // now compute total upper and lower bound (possibly +/- INF)
 // this requires the value of the linear function, so ensure it is computed
 if( f_lf ) {
  f_lf->compute( true );
  Fi0Lmb1 = rs( f_lf->get_upper_estimate() );
  }
 else
  Fi0Lmb1 = 0;

 if( UpFiLmb1def == NrFi ) {
  ++UpFiLmb1def;  // all components + the sum computed
  // can now compute the total value: 
  UpFiLmb1.back() = std::accumulate( UpFiLmb1.begin() , --(UpFiLmb1.end()) ,
				     Fi0Lmb1 );
  // note that this is the point where the lower bound (if any) is taken
  // into account when computing the total upper function estimate
  if( TrueLB && ( UpFiLmb1.back() < LowerBound.back() ) )
   UpFiLmb1.back() = LowerBound.back();
  }
 else
  UpFiLmb1.back() = INFshift;

 // note that this computation does not make explicit use of the global
 // lower bound; this is not necessary since if it is a "true" lower bound
 // then it is inserted in the master problem, and therefore it determines
 // the values of the vStar[ k ] in such a way that their sum (with all the
 // required corrections) is never < than the global bound
 if( LwFiLmb1def == NrFi ) {
  ++LwFiLmb1def;  // all components + the sum computed
  LwFiLmb1.back() = std::accumulate( LwFiLmb1.begin() , --(LwFiLmb1.end()) ,
				     Fi0Lmb1 );
  }
 else
  LwFiLmb1.back() = -INFshift;


 // update the upper and lower targets - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that vStar.back() == INFshift if the bundle of any of the hard
 // components is empty, in which case there are no targets: any "finite"
 // information about them is better than what we currently have

 if( vStar.back() < INFshift ) {
  UpTrgt = UpRifFi.back() + ( 1.0 - m2 ) * vStar.back();

  if( m1 > 0 )
   LwTrgt = UpRifFi.back() + vStar.back() + m1 * ( DST + Sigma );
  else
   LwTrgt = UpRifFi.back() + ( 1.0 + m1 ) * vStar.back();
  }
 else {
  UpTrgt = INFshift;
  LwTrgt = -INFshift;
  }

 // print the new information - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_log && ( LogVerb > 4 ) ) {
  *f_log << std::endl << "    Lambda1 = [ ";
  for( auto el : Lambda1 )
   *f_log << el << " ";
  *f_log << "]";

  if( LogVerb > 5 )
   for( Index k = 0 ; k < NrFi ; ++k ) {
    *f_log << std::endl << "    UB[ " << k << " ] = " << def;
    if( f_convex )
     pval( *f_log , UpFiLmb1[ k ] );
    else
     pval( *f_log , - LwFiLmb1[ k ] );
    *f_log << ", LB[ " << k << " ] = ";
    if( f_convex )
     pval( *f_log , LwFiLmb1[ k ] );
    else
     pval( *f_log , - UpFiLmb1[ k ] );
    }
  }
 }  // end( BundleSolver::FormLambda1 )

/*--------------------------------------------------------------------------*/

BundleSolver::Index BundleSolver::InnerLoop( bool extrastep )
{
 // here one might change the value of wFi, corresponding to the first
 // component to be evaluated, if a non-strictly-round-robin order is
 // sought for

 // compute the minimum number of components to evaluate
 Index minceval = ( MinNrEvls >= 0 ? Index( MinNrEvls )
		                   : ( NrFi - NrEasy ) * ( - MinNrEvls ) );
 Index ceval = 0;  // how many components have been evaluated so far

 for( bool insrtd = false ; ; ) {
  // round-robin-like loop between the different components
  if( ! FindNext() ) {  // find next component
   if( ! ceval )        // no component found to evaluate
    Result = kError;    // this is an error (or is it?)
   break;               // anyway, nothing else to do but stop
   }

  if( FiAndGi( f_wFi , ! extrastep ) )
   insrtd = true;

  // return if an unrecoverable error happens
  if( ( FiStatus[ f_wFi ] <= kUnEval ) || ( FiStatus[ f_wFi ] >= kError ) ) {
   Result = kError;
   break;
   }

  // if any component evaluates to -INF, then the whole problem evaluates to
  // -INF and it is therefore unbounded below; due to convexity, it is "very
  // seriously unbounded" in that a convex function being -INF anywhere is
  // -INF everywhere, hence if this ever happens it will do it "very soon"
  // (the very first time the offending component is evaluated)
  if( UpFiLmb1[ f_wFi ] == -INFshift )  {
   UpFiLmb1.back() = -INFshift;
   break;
   }

  if( ! CurrNrEvls[ f_wFi ] )  // not evaluated before
   ++ceval;                    // one more evaluated
  ++CurrNrEvls[ f_wFi ];       // evaluated once more

  // a SS can be performed: note the "<" in the SS condition below (which
  // means it is ever so slightly stronger than it should), which is there
  // to avoid the condition to work when UpFiLmb1.back() == INF == UpTrgt
  if( ( ! MPchgs ) && ( UpFiLmb1.back() < UpTrgt ) )
   MPchgs = 1;

  if( ( ! MPchgs ) && insrtd && RifeqFi && ( LwFiLmb1.back() > LwTrgt ) )
   // doing a NS without possibly evaluating all the components is inhibited
   // if the linearization errors are not computed w.r.t. the "true" value
   // of (every component of); this corresponds to the assumption in the
   // theory that a finite upper bound is known for every component. this
   // implies that eventually all components will be evaluated, which will
   // typically yield a SS
   //
   // for a NS to be performed, LwFiLmb1 must be > than the lower target;
   // again, note the ">" instead of the ">=" (which means this is ever so
   // slightly stronger than it should), which is there to avoid the
   // condition to work when LwFiLmb1.back() == -INF == LwTrgt
   //
   // however, for a NS to guarantee no cycling, at least something must
   // have been inserted (on top of all the other conditions)
   MPchgs = 1;

  if( ( MaxTime < INFshift ) && ( get_elapsed_time() > MaxTime ) ) {
   Result = kStopTime;     // time has ran up
   break;                  // nothing else to do but stop
   }

  if( ( ceval < minceval ) || extrastep )  // too few components compute()-d
   continue;               // do not stop regardless of MPchgs

  if( MPchgs )             // the MP is guaranteed to change
   break;                  // happily stop

  }  // end( for( functions evaluation loop ) )

 return( ceval );

 }  // end( BundleSolver::InnerLoop )

/*--------------------------------------------------------------------------*/

bool BundleSolver::FiAndGi( Index wFi , bool getgi )
{
 // compute and set upper and lower cutoffs and the accuracy- - - - - - - - -

 if( f_log && ( LogVerb > 3 ) )
  *f_log << std::endl << "            Fi[ " << wFi;

 if( getgi )
  SetupFiLambda1( wFi );
 else
  SetupFiLambda( wFi );

 // compute the C05Function and retrieve upper and lower estimates- - - - - -

 auto fwFi = v_c05f[ wFi ];

 auto start = std::chrono::system_clock::now();

 FiStatus[ wFi ] = fwFi->compute( ( FiStatus[ wFi ] == kUnEval ) );

 auto end = std::chrono::system_clock::now();
 std::chrono::duration< double > elapsed = end - start;

 if( ( FiStatus[ wFi ] <= kUnEval ) || ( FiStatus[ wFi ] >= kError ) ) {
  if( f_log && ( LogVerb > 3 ) )
   *f_log << " ] = Error #" <<  FiStatus[ wFi ] << ", stop";
  return( false );
  }

 auto ue = fwFi->get_upper_estimate();
 auto le = fwFi->get_lower_estimate();

 if( f_log && ( LogVerb > 3 ) ) {
  *f_log << " ]: UB = "<< def;
  pval( *f_log , ue );
  *f_log << ", LB = ";
  pval( *f_log , le );
  *f_log << " [" << fixd << elapsed.count() << "] " << def;
  }

 if( f_convex ) {
  if( ue == -INFshift )  // very special case: value == -INF
   return( true );       // immediately terminate
  }
 else
  if( le == INFshift )   // very special case: value == +INF
   return( true );       // immediately terminate

 // if getgi == false the method is actually being called on Lambda, hence
 // it is Lambda's estimates that need be updated, not Lambda1's
 if( ! getgi ) {
  update_UpFiLambd( wFi , f_convex ? ue : - le );
  update_LwFiLambd( wFi , f_convex ? le : - ue );

  // furthermore one immediately returns before getting any linearization
  return( false );
  }

 // update UpFiLambd1[ wFi ] (and possibly UpFiLambd1[ NrFi ])
 update_UpFiLambd1( wFi , f_convex ? ue : - le );

  // if bit 4 of TrgtMng == 1, then compute the upper bound in Lambda
  // provided by the upper bound in Lambda1 and try to update UpFiLmb[ wFi ]
  // (and possibly UpFiLambd1[ NrFi ])
  // note that, even if this suceeds and therefore decreases UpFiLmb[ wFi ]
  // (and possibly UpFiLambd[ NrFi ], which would be a "rather big" decrease
  // from +INF to something finite), as the theory requires the upper target
  // is *not* changed
  if( ( TrgtMng & 16 ) && ( UpFiLmb1[ wFi ] < INFshift ) ) {
   c_VarValue LwFi = fwFi->get_Lipschitz_constant();
   if( LwFi < INFshift )
    update_UpFiLambd( wFi , UpFiLmb1[ wFi ] + LwFi * NrmD );
   }

 // update LwFiLambd1[ wFi ] (and possibly LwFiLambd1[ NrFi ])
 update_LwFiLambd1( wFi , f_convex ? le : - ue );

 // get new linearizations- - - - - - - - - - - - - - - - - - - - - - - - - -

 return( GetGi( wFi ) );

 }  // end( BundleSolver::FiAndGi )

/*--------------------------------------------------------------------------*/

void BundleSolver::SetupFiLambda1( Index wFi )
{
 auto fwFi = v_c05f[ wFi ];

 // start by setting a "time cutoff" with the remaining total time
 if( MaxTime < INFshift )
  fwFi->set_par( dblMaxTime , MaxTime - get_elapsed_time() );

 if( ! ( TrgtMng & 15 ) )
  return;

 // compute upper and lower cutoffs and the accuracy
 auto UpCutOff = INFshift;
 auto LwCutOff = -INFshift;

 if( ( TrgtMng & 7 ) &&
     ( UpFiLmb1.back() > UpTrgt ) && ( LwFiLmb1.back() < LwTrgt ) ) {
  // finite upper and lower cutoffs can (and make sense to) be set only if the
  // upper and lower targets are finite (they either both are or none of them
  // are), which depends on the fact that v* = vStar.back() (and its "cousin"
  // Delta* = D_t( z* ) + Sigma*) is finite, which means that *all* non-easy
  // components have diagonal linearizations in their bundle
  //
  // note that this implies that vStar[ wFi ] (for this particular component)
  // must also be well-defined (< INF)
  //
  // if some component does not have any linearization yet, any "finite"
  // information about them is better than what we currently have, and hence
  // the cutoffs should be "as weak as they can be"
  //
  // also, the upper cutoff only makes sense to be set if the SS condition
  //
  //  \bar{f}_{tot}( x ) <= \bar{tau}  \equiv  UpFiLmb1.back() <= UpTrgt
  //
  // it is not satisfied already, i.e., if \bar{f}_{tot}( x ) > \bar{tau}
  // (note that this can never hold if \bar{tau} = +INF, due to the ">",
  // which ensures that the upper target is finite when the if() is entered).
  // indeed, the second part of the UpCutOff computation corresponds to the
  // following argument: we want the SS condition
  //
  //  \bar{f}_{tot}( x ) = \sum_i \bar{f}_i( x ) <= \bar{tau}
  //
  // to hold, which means
  //
  //  \bar{f}_k( x ) <= \bar{tau} - \sum_{i \neq k} i \bar{f}_i( x )
  //
  // but
  //
  //  \sum_{i \neq k} i \bar{f}_i( x ) = \bar{f}_{tot}( x ) - \bar{f}_k( x )
  //
  // and hence we will be content if
  //
  //  \bar{f}_k( x ) <= \bar{tau} - ( \bar{f}_{tot}( x ) - \bar{f}_k( x ) )
  //
  // this is of course conditional to the fact that the function value is
  // available, i.e., \bar{f}_{tot}( x ) < INF ==> \bar{f}_k( x ) < INF, but
  // it is also conditional to the fact that the SS condition does not hold
  // already. in fact, in the above condition \bar{f}_k( x ) on the left means
  // "after having computed the function value at x", whereas \bar{f}_k( x )
  // on the right means "the current value it has". if the SS condition
  //
  //  \bar{f}_{tot}( x ) <= \bar{tau}
  //
  // holds already when the method is called, then
  //
  //  \bar{f}_k( x ) <= \bar{tau} - ( \bar{f}_{tot}( x ) - \bar{f}_k( x ) )
  //
  // also holds with \bar{f}_k( x ) meaning "the current value it has" in
  // both places, and therefore *a fortiori* after that the function is
  // computed. thus, if the SS condition holds already, then there is no
  // reason to put a finite upper cutoff, since any value would do.
  //
  // of course the symmetric argument holds for the lower cutoff and the
  // NS condition
  //
  //  \underline{f}_{tot}( x ) >= \underline{tau}
  //
  // and note that if either one of the two conditions holds already, then
  // no cutoff (be it upper or lower) need be set, and neither does the
  // accuracy, because the fate of the step is already sealed whatever the
  // information that the oracle provides

  // note: the reason why upper and lower cutoffs are not entirely separated
  // is the computation of BetaK(), which is currently very easy but it may
  // one day become more sophisticated and therefore costly
  const auto bk = BetaK( wFi );
  const auto LwFiK = UpRifFi[ wFi ] + vStar[ wFi ];

  UpCutOff = LwFiK - m2 * bk * vStar.back();
  LwCutOff = LwFiK + std::abs( m1 ) * bk * ( DST + Sigma );

  if( UpCutOff < LwCutOff )  // this should never happen, but account for
   LwCutOff = UpCutOff;      // numerical issues

  // note that here UpCutOff >= LwCutOff and the next step can only increase
  // UpCutOff and decrease LwCutOff, so the relationship will keep holding

  // note that UpFiLmb1.back() < INFshift ==> UpFiLmb1[ wFi ] < INFshift
  if( UpFiLmb1.back() < INFshift )
   UpCutOff = std::max( UpTrgt - ( UpFiLmb1.back() - UpFiLmb1[ wFi ] ) ,
			UpCutOff );

  // note that LwFiLmb1.back() > -INFshift ==> LwFiLmb1[ wFi ] > -INFshift
  if( LwFiLmb1.back() > -INFshift )
   LwCutOff = std::min( LwTrgt - ( LwFiLmb1.back() - LwFiLmb1[ wFi ] ) ,
			LwCutOff );
  }

 // assign the cutoff values to the C05Function - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ( TrgtMng & 1 ) && f_convex ) ||
     ( ( TrgtMng & 2 ) && ( ! f_convex ) ) ) {
  auto lt = f_convex ? LwCutOff : - UpCutOff;
  if( f_log && ( LogVerb > 3 ) )
   *f_log << " ~ lt = " << def << lt;
  fwFi->set_par( dblLwCutOff , lt );
  }

 if( ( ( TrgtMng & 2 ) && f_convex ) ||
     ( ( TrgtMng & 1 ) && ( ! f_convex ) ) ) {
  auto ut = f_convex ? UpCutOff : - LwCutOff;
  if( f_log && ( LogVerb > 3 ) )
   *f_log << " ~ ut = " << def << ut;
  fwFi->set_par( dblUpCutOff , ut );
  }

 if( TrgtMng & 12 ) {
  double EpsCurr = 100;  // 1e+2 relative error is "a finite INF"

 // set EpsCurr to the difference between upper and lower cutoffs if they
 // are finite, and leave it to "INF" == "anything goes" otherwise
  if( ( TrgtMng & 4 ) &&
      ( UpCutOff < INFshift ) && ( LwCutOff > -INFshift ) )
   EpsCurr = ( UpCutOff - LwCutOff ) /
                           std::max( 1.0 , std::abs( UpRifFi[ wFi ] ) );

  // set EpsCurr to (almost) the current achieved accuracy EpsU as
  // computed via tStar if this is larger than the one before (or no
  // value has been set before)
  if( ( TrgtMng & 8 ) &&
      ( ( EpsCurr < EpsU / Nearly ) || ( EpsCurr == 100 ) ) )
   EpsCurr = EpsU / Nearly;

  // EpsCurr never needs be (much) smaller than the required relative
  // accuracy for the overall computation,
  EpsCurr = std::max( EpsCurr , RelAcc / Nearly );

  if( f_log && ( LogVerb > 3 ) )
   *f_log << " ~ eps = " << shrt << EpsCurr;
  fwFi->set_par( dblRelAcc , EpsCurr );
  }
 }  // end( BundleSolver::SetupFiLambda1 )

/*--------------------------------------------------------------------------*/

void BundleSolver::SetupFiLambda( Index wFi )
{
 auto fwFi = v_c05f[ wFi ];

 // start by setting a "time cutoff" with the remaining total time
 if( MaxTime < INFshift )
  fwFi->set_par( dblMaxTime , MaxTime - get_elapsed_time() );

 if( ! ( TrgtMng & 15 ) )
  return;

 // compute upper and lower cutoffs and the accuracy- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this is just to recompute an already computed value, hence getting back
 // what one already had is OK

 double UpCutOff = UpFiLmb1[ wFi ];
 double LwCutOff = LwFiLmb1[ wFi ];

 double EpsCurr = ( UpCutOff - LwCutOff ) /
                                std::max( 1.0 , std::abs( UpRifFi[ wFi ] ) );

 // assign the cutoff values to the C05Function - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_convex ) {
  if( TrgtMng & 1 )
   fwFi->set_par( dblLwCutOff , LwCutOff );

  if( TrgtMng & 2 )
   fwFi->set_par( dblUpCutOff , UpCutOff );
  }
 else {
  if( TrgtMng & 1 )
   fwFi->set_par( dblUpCutOff , - LwCutOff );

  if( TrgtMng & 2 )
   fwFi->set_par( dblLwCutOff , - UpCutOff );
  }

 if( TrgtMng & 12 )
  fwFi->set_par( dblRelAcc , EpsCurr );

 }  // end( BundleSolver::SetupFiLambda )

/*--------------------------------------------------------------------------*/

bool BundleSolver::GetGi( Index wFi )
{
 bool insrtd = false;  // keep track if anything new at all was inserted
 auto fwFi = v_c05f[ wFi ];

 for( Index Ftchd = 0 ; Ftchd < aBP3 ; ++Ftchd ) {
  bool diagonal = true;
  bool HasLinearization;

  if( ! Ftchd ) {  // first look for a constraint then for a sub-gradient
   if( UpFiLmb1[ wFi ] == INFshift ) {
    HasLinearization = fwFi->has_linearization( diagonal = false );
    if( ! HasLinearization )
     HasLinearization = fwFi->has_linearization( diagonal = true );
    }
   else
    HasLinearization = fwFi->has_linearization( diagonal );
   }
  else {
   if( UpFiLmb1[ wFi ] == INFshift ) {
    HasLinearization = fwFi->compute_new_linearization( diagonal = false );
    if( ! HasLinearization )
     HasLinearization = fwFi->compute_new_linearization( diagonal = true );
    }
   else
    HasLinearization = fwFi->compute_new_linearization( diagonal );
   }

  if( ! HasLinearization )  // no new linearization of either type available
   break;                   // nothing else to do

  if( ! diagonal )          // a vertical linearization changes the MP
   MPchgs = 2;              // no matter what else happens

  // check if aggregation has to be performed - - - - - - - - - - - - - - - -
  // doing this now could occasionally result in useless aggregations, but it
  // is necessary due to limitations in the MPSolver interface (there can be
  // only one "un-named item being inserted", so inserting Z[ wFi ] while
  // inserting the new item is complicated

  auto wh = BStrategy( wFi );

  // get the space for the item from the MPSolver - - - - - - - - - - - - - -

  auto G1k = Master->GetItem( wFi + 1 );

  // fetch the item from the Oracle - - - - - - - - - - - - - - - - - - - - -

  fwFi->get_linearization_coefficients( G1k );
  if( ! f_convex )
   chgsign( G1k , NumVar );

  auto Alfa1k = rs( fwFi->get_linearization_constant() );
  HpNum eps;

  // pass the base to the MP Solver - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set SGBse = nullptr;
  Master->SetItemBse( SGBse , NumVar );

  // compute ScPr1k and Alfa1k- - - - - - - - - - - - - - - - - - - - - - - -

  Index cp;
  HpNum ScPr1k;
  bool is_rep = diagonal && ( ! Ftchd ) && ( ! CurrNrEvls[ wFi ] );
  // if it is the first subgradient of the first call to GetGi() for this
  // component, it may be the "representative subgradient"

  if( diagonal ) {  // it is a subgradient
   // compute the lower bound in Lambda provided by the subgradient
   auto FikLmb = Alfa1k +
    std::inner_product( Lambda.begin() , Lambda.end() , G1k , double( 0 ) );

   // try to update LwFiLmb[ wFi ] (and possibly LwFiLambd.back())
   // note that, even if this suceeds and therefore increases LwFiLmb[ wFi ]
   // (and possibly LwFiLambd.back(), which would be a "rather big" increase
   // from -INF to something finite), as the theory requires the lower target
   // is *not* changed
   update_LwFiLambd( wFi , FikLmb );

   // note that given FikLmb, the linearization error in Lambda is obtained
   // for free as UpFiLmb[ wFi ] - FikLmb; however, the MPSolver interface
   // requires CheckSubG() to be called and the computation to be done there
   // (possibly with a larger error, but one day MPSolver will go ...)

   // compute the linearization error in Lambda1
   Alfa1k = UpFiLmb1[ wFi ] - Alfa1k -
    std::inner_product( Lambda1.begin() , Lambda1.end() , G1k , double( 0 ) );

   // this is the eps so that G1 is an eps-subgradient in Lambda1
   eps = Alfa1k;

   // CheckSubG changes Alfa1k so that G1 is an Alfa1k-subgradient in Lambda
   cp = Master->CheckSubG( UpFiLmb1[ wFi ] - UpRifFi[ wFi ] , t , Alfa1k ,
			   ScPr1k );
   }
  else {             // it is a constraint
   // the definition of constraint in FiOracle (hence MPSolver) is
   //
   //      SubG * Lambda <= GetVal()
   //
   // i.e., SubG * Lambda - GetVal() <= 0, whereas in C05Function it is
   //
   //       ( 0 , - g ) ( v , x ) >= \alpha
   //
   // i.e., g x + \alpha <= 0; this means that the \alpha produced by
   // get_linearization_constant() is the opposite than that of GetVal()
   // in fact, the standard form of the constraints in the master problem is
   // [ v / 0 ] >= g d - \alpha while the diagonal linearizations are
   //
   //       ( 1 , - g ) ( v , x ) >= \alpha
   //
   // and indeed Alfa1k is also "changed in sign" in the diagonal
   // linearization (subgradient) case

   Alfa1k = - Alfa1k;
   cp = Master->CheckCnst( Alfa1k , ScPr1k , Lambda.data() );
   }

  Index gpp = InINF;  // position in the global pool where to put it

  if( f_log && ( LogVerb > 2 ) ) {
   *f_log << std::endl << "            New " << shrt;
   if( diagonal ) {
    if( eps >= std::max( std::abs( UpRifFi[ wFi ] ) , double( 1 ) )
	       * RelAcc / 10 )
     *f_log << "eps-subgradient with eps = " << eps;
    else
     *f_log << "subgradient";
    *f_log << " for Fi[ " << wFi << " ] ~ Alfa1 = " << Alfa1k
	   << " ~ gd = " << rs( ScPr1k );
    }
   else
    *f_log << "constraint " << wh << " ~ rhs = " << Alfa1k;
   }

  bool to_insert = true;  // if it has to be inserted

  if( cp < InINF ) {  // the item is a copy - - - - - - - - - - - - - - - - -
   BLOG( 2 , " is copy of " << cp << " (" << ItemVcblr[ cp ].second << ")" );

   wh = cp;  // we have it already

   auto OldA1k = (Master->ReadLinErr())[ cp ];

   assert( ( ItemVcblr[ cp ].first == wFi ) &&
           ( ItemVcblr[ cp ].second < vBPar2[ wFi ] ) &&
	   ( InvItemVcblr[ wFi ][ ItemVcblr[ cp ].second ] == cp ) );

   if( OldA1k >= Alfa1k + std::max( std::abs( Alfa1k ) , double( 1 ) )
                          * RelAcc / 10 ) {
    // if the copy has a *substantially* smaller Alfa than the original,
    // replace the original with the copy; in principle relative differences
    // smaller than RelAcc could be ignored, but we use RelAcc / 10 for safety

    BLOG( 2 , " with smaller Alfa" );

    gpp = ItemVcblr[ cp ].second;
    if( ( BPar7 & 3 ) < 3 ) {
     // BundleSolver does not immediately replace the copy unless necessary,
     // but clearly if one linearization in the global pool has to be
     // sacrificed, it'll be the copy
     auto ngpp = find_place_in_global_pool( wFi );
     if( ngpp < InINF ) {       // a free place has been found
      // although the old linearization is kept, it is removed from the
      // bundle: the position cp is now associated with ngpp, which
      // means that position gpp is now free
      remove_from_global_pool( wFi , gpp , false );
      gpp = ngpp;                      // store the copy there
      BLOG( 2 , " (" << gpp << ")" );  // print the chosen place
      }
     }

    // if it is the "representative subgradient", add its contribution to
    // the required ones of Alfa1, ScPr1 and G1 (if any); do this before
    // the call to SubstItem() because the state of the G1k memory after
    // the call is unclear
    if( is_rep ) {
     whisG1[ wFi ] = cp;
     if( NeedsAlfa1() )
      Alfa1 += Alfa1k;
     if( NeedsScPr1() )
      ScPr1 += ScPr1k;
     if( NeedsG1() )
      vect_sum( G1 , G1k );
     }

    Master->SubstItem( cp );  // substitute it in the master problem
    // note that the number of items of component wFi in the master problem
    // is unchanged
    }
   else                 // the item is a copy, not better than the original
    to_insert = false;  // do nothing
   }
  else {           // the item is not a copy- - - - - - - - - - - - - - - - -
   // insert the item, if there is space

   if( wh == InINF )  // the position has not been selected in BStrategy()
    wh = FindAPlace( wFi );  // find a free spot in the bundle

   if( wh == InINF ) {  // no space found ...
    if( ! Ftchd ) {     // ... and this was the first item
     BLOG( 1 , std::endl << " ERROR: No space in the bundle" << std::endl );
     Result = kError;   // signal an error to end the outer Fi-cycle
     }
    else
     BLOG( 1 , std::endl << " WARNING: No space in the bundle" << std::endl );
    break;              // the cycle ends
    }

   if( ItemVcblr[ wh ].second < vBPar2[ wFi ] )
    // the place is occupied already: this happens if the bundle was full
    // (and, possibly aggregation has been performed for safety)
    Master->RmvItem( wh );  // the old item has to be removed first
   else {                   // the place is unoccupied
    ++NrItems[ wFi ];       // one more item in the bundle (otherwise the
    ++NrItems[ NrFi ];      // number remains the same as one is replaced)
    }

   // if it is the "representative subgradient", add its contribution to
   // the required ones of Alfa1, ScPr1 and G1 (if any); do this before
   // the call to SetItem() because the state of the G1k memory after
   // the call is unclear
   if( is_rep ) {
    whisG1[ wFi ] = wh;
    if( NeedsAlfa1() )
     Alfa1 += Alfa1k;
    if( NeedsScPr1() )
     ScPr1 += ScPr1k;
    if( NeedsG1() )
     vect_sum( G1 , G1k );
    }

   Master->SetItem( wh );   // insert the new item in the MP Solver

   // now find a position in the global pool of component wFi where to store
   // the new linearization
   gpp = find_place_in_global_pool( wFi );

   if( gpp == InINF ) {     // there is none
    // this means that not only the global pool is full, but also the bundle
    // is full: one can therefore put it in the very place of the item it
    // replaces, which must be an item of the same component because
    // BStrategy() ensures this
    assert( ItemVcblr[ wh ].first == wFi );
    gpp = ItemVcblr[ wh ].second;
    assert( gpp < vBPar2.back() );
    }

   BLOG( 2 , " stored in " << wh << " (" << gpp << ")"  );
   }

  // if something was inserted, bookkeeping is needed - - - - - - - - - - - -

  if( to_insert ) {
   inhibit_Modification( true );
   v_c05f[ wFi ]->store_linearization( gpp );
   inhibit_Modification( false );

   insrtd = true;
   add_to_global_pool( wFi , gpp , wh );

   if( diagonal )       // it is a subgradient
    OOBase[ wh ] = -1;  // ensure it won't be touched again this round
   else                 // it is a constraint
    // mark it as permanently fixed: this may be a bad choice in practice,
    // although it is required by the theory (we'll see ...)
    OOBase[ wh ] = -Inf< SIndex >();
   }

  #if CHECK_DS & 1
   CheckBundle();
  #endif

  }  // end( items-collecting loop )- - - - - - - - - - - - - - - - - - - - -

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( insrtd );  // returns true if at least one item was inserted

 }  // end( BundleSolver::GetGi )

/*--------------------------------------------------------------------------*/

void BundleSolver::GotoLambda1( void )
{
 // compute DeltaFi - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector< VarValue > DF( NrFi + 1 );
 /* DeltaFi = UpFiLmb1 - UpRifFi; note that the code one may expect
  *
  * std::transform( UpFiLmb1.begin() , UpFiLmb1.end() , UpRifFi.begin() ,
  *                 DF.begin() , std::minus< double >() );
  *
  * is wrong since the format of DeltaFi expected by ChangeCurrPoint() is
  * different from the one used in BundleSolver; in particular, the total
  * value need be in DF.front() rather than in DF.back(), and the value for
  * component i need be in DF[ i + 1 ] rather than in DF[ i ]. */

 DF.front() = UpFiLmb1.back() - UpRifFi.back();
 std::transform( UpFiLmb1.begin() , --(UpFiLmb1.end()) , UpRifFi.begin() ,
 		 ++(DF.begin()) , std::minus< double >() );

 // do the move - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // Lambda = Lambda1 and all associated data structures

 Lambda.swap( Lambda1 );
 UpFiLmb.swap( UpFiLmb1 );
 LwFiLmb.swap( LwFiLmb1 );
 UpRifFi = UpFiLmb;
 RifeqFi = true;
 UpFiLmbdef = UpFiLmb1def;
 LwFiLmbdef = LwFiLmb1def;
 Fi0Lmb = Fi0Lmb1;

 // change the current point in the MP Solver - - - - - - - - - - - - - - - -

 Master->ChangeCurrPoint( t , DF.data() );

 #if CHECK_DS & 4
  CheckAlpha();
 #endif
 #if CHECK_DS & 8
  CheckLBs();
 #endif

 }  // end( BundleSolver::GotoLambda1 )

/*--------------------------------------------------------------------------*/

void BundleSolver::GotoLambda( void )
{
 // compute DeltaFi - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector< VarValue > DF( NrFi + 1 );  // DF = UpFiLmb - UpRifFi

 DF.front() = UpFiLmb.back() - UpRifFi.back();
 std::transform( UpFiLmb.begin() , --(UpFiLmb.end()) , UpRifFi.begin() ,
 		 ++(DF.begin()) , std::minus< double >() );

 UpRifFi = UpFiLmb;        // set UpFiLmb as the reference values
 RifeqFi = true;

 // "change the current point in the MP Solver" - - - - - - - - - - - - - - -
 // use the special case of ChangeCurrPoint() with Tau == 0, whereby only
 // the DF component of the change is effective: the current point does not
 // really change, but the linearization errors (and bounds) do

 Master->ChangeCurrPoint( double( 0 ) , DF.data() );

 #if CHECK_DS & 4
  CheckAlpha();
 #endif
 #if CHECK_DS & 8
  CheckLBs();
 #endif

 }  // end( BundleSolver::GotoLambda )

/*--------------------------------------------------------------------------*/

void BundleSolver::ResetAlfa( Index k )
{
 std::vector< VarValue > Gi( NumVar );
 std::vector< VarValue > Alfa( Master->MaxName( k == NrFi ? InINF : k + 1 ) );
 if( k == NrFi ) {  // all components need be reset
  for( Index i = 0 ; i < Master->MaxName() ; ++i )
   if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
    auto Ai = rs( v_c05f[ ItemVcblr[ i ].first
			  ]->get_linearization_constant(
						   ItemVcblr[ i ].second ) );
    #ifndef NDEBUG
     if( std::isnan( Ai ) )  // linearization no longer valid
      throw( std::logic_error( "inconsistent ItemVcblr" ) );
    #endif

    // compute the linearization error in Lambda
    v_c05f[ ItemVcblr[ i ].first ]->get_linearization_coefficients(
		   Gi.data() , Range( 0 , NumVar ) , ItemVcblr[ i ].second );
    if( ! f_convex )
     chgsign( Gi.data() , NumVar );

    Alfa[ i ] = UpRifFi[ ItemVcblr[ i ].first ] - Ai -
              std::inner_product( Lambda.begin() , Lambda.end() , Gi.data() ,
				  double( 0 ) );
    }
  }
 else {             // only that specific component need be reset
  for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
   if( InvItemVcblr[ k ][ i ] < vBPar2.back() ) {
    auto Ai = rs( v_c05f[ k ]->get_linearization_constant( i ) );

    #ifndef NDEBUG
     if( std::isnan( Ai ) )  // linearization no longer valid
      throw( std::logic_error( "inconsistent ItemVcblr" ) );
    #endif

    // compute the linearization error in Lambda
    v_c05f[ k ]->get_linearization_coefficients( Gi.data() ,
						 Range( 0 , NumVar ) , i );
    if( ! f_convex )
     chgsign( Gi.data() , NumVar );

    Alfa[ InvItemVcblr[ k ][ i ] ] = UpRifFi[ k ] - Ai -
	       std::inner_product( Lambda.begin() , Lambda.end() , Gi.data() ,
				   double( 0 ) );
    }
  }

 Master->ChgAlfa( Alfa.data() , k + 1 );

 }  // end( BundleSolver::ResetAlfa )

/*--------------------------------------------------------------------------*/

void BundleSolver::SimpleBStrat( void )
{
 if( ( BPar7 & 3 ) == 3 ) {  // "eager" deletion
  std::vector< Subset > tbdltd( NrFi );
  for( Index i = 0 ; i < Master->MaxName() ; ++i )
   if( ( OOBase[ i ] < Inf< SIndex >() ) &&
       ( OOBase[ i ] > SIndex( BPar1 ) ) ) {
    tbdltd[ ItemVcblr[ i ].first ].push_back( ItemVcblr[ i ].second );
    Delete( i );
    }

  inhibit_Modification( true );
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ! tbdltd[ k ].empty() )
    v_c05f[ k ]->delete_linearizations( std::move( tbdltd[ k ] ) , false );
  inhibit_Modification( false );
  }
 else                        // "lazy" deletion
  for( Index i = 0 ; i < Master->MaxName() ; ++i )
   if( ( OOBase[ i ] < Inf< SIndex >() ) && ( OOBase[ i ] > SIndex( BPar1 ) ) )
    Delete( i );

 #if CHECK_DS & 1
  CheckBundle();
 #endif

 }  // end( BundleSolver::SimpleBStrat )

/*--------------------------------------------------------------------------*/

double BundleSolver::BetaK( Index wFi ) {
 return( 1.0 / double( NrFi - NrEasy ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::Log1( void )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "{" << SCalls << "-" << ParIter << "-"
	<< NrItems.back() << "-" << fixd << get_elapsed_time() << "} t = "
	<< shrt << t << " ~ D*_1( z* ) = " << Master->ReadDStart( 1 )
	<< " ~ Sigma = " << Sigma << std::endl << "           ";

 if( UpFiLmb.back() == INFshift )
  *f_log << " Fi undefined";
 else
  *f_log << " Fi = " << def << rs( UpFiLmb.back() ) << " ~ eU = "
	 << shrt << EpsU;

 if( BPar6 )
  *f_log << " ~ BP3 = " << aBP3;

 }  // end( BundleSolver::Log1 )

/*--------------------------------------------------------------------------*/

void BundleSolver::Log2( double ft )
{
 if( ( ! f_log ) || ( LogVerb <= 1 ) )
  return;

 *f_log << std::endl << "            [" << fixd << ft << "] " << def;

 if( LowerBound.back() > -INFshift ) {
  if( f_convex )
   *f_log << "LB = " << LowerBound.back() << " ~ ";
  else
   *f_log << "UB = " << - LowerBound.back() << " ~ ";
  }

 *f_log << "Fi1 = ";

 if( f_convex ) {
  if( UpFiLmb1.back() == -INFshift ) {
   *f_log << "-INF => STOP." << std::endl;
   return;
   }
  else
   if( UpFiLmb1.back() == INFshift ) {
    *f_log << "+INF" << std::endl;
    return;
    }
   else
    *f_log << UpFiLmb1.back() << shrt;
  }
 else
  if( UpFiLmb1.back() == -INFshift ) {
   *f_log << "+INF => STOP." << std::endl;
   return;
   }
  else
   if( UpFiLmb1.back() == INFshift ) {
    *f_log << "-INF" << std::endl;
    return;
    }
   else
    *f_log << - UpFiLmb1.back() << shrt;

 if( NeedsAlfa1() )
  *f_log << " ~ Alfa1 = " << Alfa1;

 if( NeedsScPr1() )
  *f_log << " ~ Gi1xd = " << ScPr1;

 *f_log << std::endl;

 }  // end( BundleSolver::Log2 )

/*--------------------------------------------------------------------------*/

void BundleSolver::compute_NrmZFctr( void )
{
 auto wf = ( WZNorm << 2 );
 // if we need to sum but some component has no linearization, return:
 // NrmZFctr remains undefined
 if( wf > 1 )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ( ! NrEasy ) || ( ! IsEasy[ k ] ) ) && ( NrItems[ k ] == 0 ) )
    return;

 // now we sum: note that if ! f_convex one should change the sign, but
 // the norm is invariant w.r.t. the sign
 std::vector< VarValue > tg( NumVar , 0 );
 if( f_lf ) {      // the linear 0-th component is there
  auto & cf = f_lf->get_v_var();
  for( Index i = 0 ; i < NumVar ; ++i )
   tg[ i ] = cf[ i ].second;
  }
 else              // there is no 0-th component
  if( wf <= 1 ) {  // and we just wanted is subgradient
   NrmZFctr = 1;   // ... which is all-0, so use 1
   return;
   }

 if( wf > 1 ) {  // also need to sum a subgradient for each component
  std::vector< VarValue > tg1( NumVar );

  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( NrEasy && IsEasy[ k ] )
    continue;

   Index i = 0;
   while( InvItemVcblr[ k ][ i ] == InINF )
    ++i;

   v_c05f[ k ]->get_linearization_coefficients( tg1.data() ,
						Range( 0 , NumVar ) , i );

   std::transform( tg.begin() , tg.end() , tg1.begin() , tg.begin() ,
		   std::plus< VarValue >() );
   }
  }

 NrmZFctr = ::norm( tg , WZNorm & 3 );

 }  // end( compute_NrmZFctr )

/*--------------------------------------------------------------------------*/

bool BundleSolver::FindNext( void )
{
 Index InitwFi = f_wFi;
 do {
  f_wFi = ( f_wFi + 1 ) % NrFi;    // next patient, please
  if( NrEasy && IsEasy[ f_wFi ] )  // skip easy components
   continue;
  if( ( FiStatus[ f_wFi ] == kUnEval ) ||
      ( ( FiStatus[ f_wFi ] < kError ) && ( FiStatus[ f_wFi ] > kOK ) &&
	( CurrNrEvls[ f_wFi ] < MaxNrEvls ) ) )
   return( true );

  } while( f_wFi != InitwFi );

 return( false );

 }  // end( BundleSolver::FindNext )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolver::InitMP( void )
{
 // set the size- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Master->SetDim( vBPar2.back() , &FakeFi , false );

 // MinQuad requires a "high" accuracy to work (1e-12) while standard solvers
 // do not, and in fact may complain if such a tight accuracy is set
 double eps = ( MPName & 1 ) ? 1e-9 : 1e-12;

 Master->SetPar( MPSolver::kOptEps , eps );
 Master->SetPar( MPSolver::kFsbEps , eps );

 // insert the constant subgradient of the 0-th component - - - - - - - - - -
 // TODO: check if the 0-th component is sparse, if so pass a proper base
 //       to the MPSolver

 if( f_lf ) {
  auto G0 = Master->GetItem( 0 );
  f_lf->get_linearization_coefficients( G0 );
  if( ! f_convex )
   chgsign( G0 , NumVar );
  Master->SetItemBse( nullptr , NumVar );
  Master->SetItem( InINF );
  }

 tHasChgd = true;

 if( MPName & 8 )
  Master->CheckIdentical();

 // set the log file
 Master->SetMPLog( f_log , MPlvl );

 }  // end( BundleSolver::InitMP )

/*--------------------------------------------------------------------------*/

Index BundleSolver::BStrategy( Index wFi )
{
 // this method implements the B-strategies of the code, i.e., which "old"
 // items are discarded if the bundle is full and a new item belonging to
 // component wFi has to be inserted
 //
 // this is called *before* that we know if the place will actually be
 // required, so it has a "loose attitude"; in particular, it will return
 // InINF under two opposite set of conditions:
 //
 // - there is plenty of space left in the bundle, so that no B-strategy (no
 //   removal or aggregation) is required;
 //
 // - there is no way in which space can be found, i.e., the bundle (for this
 //   component) is full, and removal/aggregation are not successful (a very
 //   strange occurrence due to an extremely small bundle or it being chock
 //   full of constraints)
 //
 // Picking a specific spot in the free space is the task of FindAPlace(),
 // which however is not called right away because the place may end up not
 // being needed. If BStrategy() return InINF because there is plenty of
 // space then FindAPlace() will suceed, if BStrategy() return InINF because
 // there is no way space can be found then FindAPlace() will fail and it
 // will be clear that disaster looms
 //
 // important note: if the returned wFi is not InINF, *and* it is the name
 // of an item still in the bundle, then:
 //
 // - the item belongs to the component wFi, so as to keep the number of
 //   items of that component constant
 //
 // - the item is *not* deleted, since it is not 100% sure this will need
 //   to be done, so deletion will be responsibility of the caller

 // there are "free" items in the global pool: there is "plenty of space"
 if( FrFItem[ wFi ] < vBPar2[ wFi ] )
  return( InINF );

 // there is not plenty of space, take 1- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // there are no "free" items in the gloal pool, but there are items in the
 // global pool that are not in the bundle; these will go first

 if( NrItems[ wFi ] < vBPar2[ wFi ] )
  return( InINF );

 // there is not plenty of space, take 2- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the bundle for component wFi is full, which implies that the global pool
 // is also full: among the items in the bundle for component wFi, find the
 // removable one with largest OOBase[], snd among these with the largest
 // OOBase[] select that with largest Alfa[]
 // note: the Z[ wFi ] in the bundle (if any) is not removable and
 //       therefore cannot be selected, which in particular happens if
 //       wFi has only *one* subgradient in base

 Index wh;
 SIndex OOwh = -Inf< SIndex >();
 HpNum Awh = -Inf< HpNum >();
 cHpRow tA = Master->ReadLinErr();
 for( auto i : InvItemVcblr[ wFi ] ) {
  assert( i < vBPar2.back() );
  if( ( OOBase[ i ] > OOwh ) ||
      ( ( OOBase[ i ] == OOwh ) && ( tA[ i ] > Awh ) ) ) {
   wh = i;
   OOwh = OOBase[ i ];
   Awh = tA[ i ];
   }
  }

 if( OOBase[ wh ] < 0 )    // all items are non-removable: nothing else to
  return( InINF );         // do (except maybe complaining very loudly)
 else
  if( OOBase[ wh ] > 0 )   // a place is found
   return( wh );           // there are no problems, all done

 // wh is a basic item- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this means that *all* items of component wFi are either in base or not
 // removable, for otherwise we would have selected an item with OOBase > 0;
 // we cannot discard anything before having performed aggregation, but in
 // order to do so we also need to free some space for the Z[ wFi ]
 //
 // note: there is an easy case where z is "naturally" in the base without
 //       any aggregation: there is only one subgradient in base (for the
 //       component wFi), and therefore, its Mlt[] is == 1. however this
 //       can never happen here because it would mean vBPar2[ wFi ] == 1
 //       which is not permitted. the other case in which this could happen
 //       is that all the items in the bundle for component wF are
 //       constraints, but then aggregation would not be useful (in fact
 //       this is checked and reported as a failure)
 //
 // note: this also means that wh is the item in base with largest Alpha
 //
 // note: that since *all* items of this component are either in base or
 //       not removable, we can scan MBse[] for the items to be removed,
 //       possibly ignoring items with Mlt[] == 0 -- but in fact not doing
 //       it because there is not any

 if( Zvalid[ wFi ] )  // a valid Z[ wFi ] is in the bundle: it is safe to
  return( wh );       // replace wh with anything the oracle provides us

 // a valid Z[ wFi ] is not already in: aggregation has to be performed

 Index MBDm;
 cIndex_Set MBse;
 cHpRow Mlt = Master->ReadMult( MBse , MBDm , wFi + 1 , false );

 Index whZ = InINF;  // the position where Z[ wFi ] has to go

 if( ( whisZ[ wFi ] < InINF ) && Master->IsSubG( whisZ[ wFi ] ) ) {
  whZ = whisZ[ wFi ];  // preferably re-use the last position

  if( whZ == wh ) {  // this was the slot selected for the new item
   // re-select wh as the one with min Mult among all the removable ones
   // different from whZ

   wh = InINF;
   cHpRow tMlt = Mlt;
   HpNum tMin = Inf< HpNum >();
   if( MBse ) {
    cIndex_Set tMBse = MBse;
    for( Index h ; ( h = *(tMBse++) ) < InINF ; ++tMlt )
     if( ( h != whZ ) && ( *tMlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
      wh = h;
      tMin = *tMlt;
      }
    }
   else
    for( Index h = 0 ; h < MBDm ; ++h , ++tMlt )
     if( ( h != whZ ) && ( *tMlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
      wh = h;
      tMin = *tMlt;
      }
    }

  if( wh == InINF )  // nothing valid found (??)
   return( wh );
  }
 else {
  // there is no last position for Z[ wFi ], choose the one with min Mult
  // among all the removable ones different from wh

  cHpRow tMlt = Mlt;
  HpNum tMin = Inf< HpNum >();
  if( MBse ) {
   cIndex_Set tMBse = MBse;
   for( Index h ; ( h = *(tMBse++) ) < InINF ; ++tMlt )
    if( ( h != wh ) && ( *tMlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
     whZ = h;
     tMin = *tMlt;
     }
   }
  else
   for( Index h = 0 ; h < MBDm ; ++h , ++tMlt )
    if( ( h != wh ) && ( *tMlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
     whZ = h;
     tMin = *tMlt;
     }
  }

 if( whZ == InINF )  // there is no removable item apart from wh
  return( InINF );   // nothing else to do except complaining very loudly

 // tell the C05Function what is going to happen- - - - - - - - - - - - - - -
 // note that this only happens when the bundle (for component wFi) is "very
 // full", and therefore also the global pool (for component wFi) is such.
 // hence, whZ is an item already in the bundle, and therefore in the global
 // pool. the natural choice is to put the new aggregate linearization in the
 // same position in the global pool where whZ was, i.e.,
 // ItemVcblr[ whZ ].second. hence, ItemVcblr, InvItemVcblr and so on need
 // not be changed

 LinearCombination coeff( MBDm );
 if( MBse )
  for( Index i = 0 ; i < MBDm ; ++i ) {
   coeff[ i ].first = ItemVcblr[ MBse[ i ] ].second;
   coeff[ i ].second = Mlt[ i ];
   }
 else
  for( Index i = 0 ; i < MBDm ; ++i ) {
   coeff[ i ].first = ItemVcblr[ i ].second;
   coeff[ i ].second = Mlt[ i ];
   }

 inhibit_Modification( true );
 v_c05f[ wFi ]->store_combination_of_linearizations( coeff ,
						   ItemVcblr[ whZ ].second );
 inhibit_Modification( false );

 Master->RmvItem( whZ );  // remove the old item in position whZ

 // ask the MPSolver the memory for keeping Z[ wFi ]- - - - - - - - - - - - -
 // note: Mlt and MBse could very well be "temporary" memory belonging to the
 // MPSolver, and any call to a method of the MPSolver may invalidate it;
 // the calls start now, and in fact MBse and Mlt are no longer used

 SgRow tZ = Master->GetItem( wFi + 1 );

 // read Z[ wFi ] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index ZBDm;
 cIndex_Set ZBse;
 Master->ReadZ( tZ , ZBse , ZBDm , wFi + 1 );

 #if CHECK_DS & 2
  // MinQuad requires a "high" accuracy to work (1e-12) while standard solvers
  // do not, and in fact may complain if such a tight accuracy is set
  double eps = ( MPName & 1 ) ? 1e-9 : 1e-12;

  std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                            : f_log;
  *wlog << def;
  std::vector< VarValue > Z( NumVar );
  v_c05f[ wFi ]->get_linearization_coefficients( Z.data() ,
						 Range( 0 , NumVar ) ,
						 ItemVcblr[ whZ ].second );
  if( ! f_convex )
   chgsign( Z.data() , NumVar );
  if( ZBse ) {
   Index j = 0;
   for( Index i = 0 ; i < NumVar ; ++i )
    if( ( j < ZBDm ) && ( ZBse[ j ] == i ) ) {
     if( std::abs( Z[ i ] - tZ[ j ] ) >=
	 eps * std::max( Z[ i ] , double( 1 ) ) )
      *wlog << std::endl << "Z[ " << i << " ]: F = " << Z[ i ]
	    << " ~ M = " << tZ[ j ];
     ++j;
     }
    else
     if( std::abs( Z[ i ] ) >= eps )
      *wlog << std::endl << "Z[ " << i << " ]: F = " << Z[ i ]
	    << " ~ M = 0";
   }
  else
   for( Index i = 0 ; i < NumVar ; ++i )
    if( std::abs( Z[ i ] - tZ[ i ] ) >=
	eps * std::max( Z[ i ] , double( 1 ) ) )
     *wlog << std::endl << "Z[ " << i << " ]: F = " << Z[ i ]
	   << " ~ M = " << tZ[ i ];
 #endif

 // now pass Z[ wFi ] back to the MP Solver - - - - - - - - - - - - - - - - -

 Master->SetItemBse( ZBse , ZBDm );

 HpNum ScPri;
 HpNum Ai = Master->ReadSigma( wFi + 1 );      // its alfa is Sigma[ wFi ]

 #if CHECK_DS & 2
  HpNum tAi = rs( v_c05f[ wFi ]->get_linearization_constant(
						ItemVcblr[ whZ ].second ) );
  tAi = UpRifFi[ wFi ] - tAi -
                         std::inner_product( Lambda.begin() , Lambda.end() ,
					     Z.begin() , double( 0 ) );

  if( std::abs( tAi - Ai ) >=
      eps * std::max( std::max( Ai , UpRifFi[ wFi ] ) , double( 1 ) ) )
   *wlog << std::endl << "Sigma: F = " << tAi << " ~ M = " << Ai;
 #endif

 // note that Tau == -1, meaning that Ai need not be changed since
 // Ai is already the linearization error in Lambda, but still the
 // ScPri need be computed
 Master->CheckSubG( 0 , -1 , Ai , ScPri );

 Master->SetItem( whZ );  // set Z[ wFi ] in position whZ

 whisZ[ wFi ] = whZ;      // Z[ wFi ] is in the bundle in position whZ
 Zvalid[ wFi ] = true;    // ... and it is valid
 OOBase[ whZ ] = -1;      // ... and it won't be removed in this iteration

 BLOG( 2 , std::endl << "Aggregation performed into " << whZ );

 // at this point, Z[ wFi ] is in the bundle, hence it is safe to replace wh
 // with anything the oracle provides us

 return( wh );

 }  // end( BundleSolver::BStrategy )

/*--------------------------------------------------------------------------*/

Index BundleSolver::FindAPlace( Index wFi )
{
 // this method is used to return the index of an available position in the
 // bundle where to store a new item belonging to "component" wFi; if there
 // are no possible positions left, then InINF is returned
 //
 // note that positions in the bundle are not "reserved to components", so
 // the task of FindAPlace() is easy. it is the business of other parts of
 // the code to ensure that a component does not use more than its fair
 // share of positions in the bundle

 Index wh = InINF;

 if( ! FreList.empty() ) {       // there are deleted items
  wh = FreList.top();            // pick the one with smaller name
  if( wh >= Master->MaxName() )  // FreList has all items with "large" names
   FreList = {};                 // clear FreList
  else {                         // wh is a "small" name
   FreList.pop();                // take it away
   return( wh );                 // all done
   }
  }

 // if there are no deleted items (with suitably small names)
 if( Master->MaxName() < vBPar2.back() )  // ... but there is still space
  wh = Master->MaxName();                 // next name

 return( wh );

 }  // end( BundleSolver::FindAPlace )

/*--------------------------------------------------------------------------*/

bool BundleSolver::NeedsAlfa1( void )
{
 // Alfa1 is only used in Heuristic2
 return( ( ( ( tSPar1 & 1 ) && ( ( tSPar1 & tSPHMsk1 ) == 64 ) ) ) ||
	 ( ( ( tSPar1 & 2 ) && ( ( tSPar1 & tSPHMsk2 ) == 256 ) ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

bool BundleSolver::NeedsScPr1( void )
{
 // ScPr1 is used by everyone save for Heuristic1
 return( ( ( ( tSPar1 & 1 ) && ( tSPar1 & tSPHMsk1 ) ) ) ||
	 ( ( ( tSPar1 & 2 ) && ( tSPar1 & tSPHMsk2 ) ) ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

bool BundleSolver::NeedsG1( void )
{
 // G1 is only used in Heuristic4
 return( ( ( ( tSPar1 & 1 ) && ( ( tSPar1 & tSPHMsk1 ) == 172 ) ) ) ||
	 ( ( ( tSPar1 & 2 ) && ( ( tSPar1 & tSPHMsk2 ) == 768 ) ) ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdateHeuristicInfo( void )
{
 // if required, update the "aggregated" Alfa1 and ScPr1, that are used in
 // the t heuristics, using the "representatives" of all components that
 // have not been compute()-d in the latest InnerLoop()

 if( NeedsAlfa1() && NeedsG1() ) {
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ! CurrNrEvls[ k ] ) && ( whisG1[ k ] < InINF ) ) {
    Alfa1 += (Master->ReadLinErr())[ whisG1[ k ] ];
    ScPr1 += Master->ReadGid( whisG1[ k ] );
    }

  return;
  }

 if( NeedsAlfa1() ) {
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ! CurrNrEvls[ k ] ) && ( whisG1[ k ] < InINF ) )
    Alfa1 += (Master->ReadLinErr())[ whisG1[ k ] ];

  return;
  }

 if( NeedsG1() )
  for( Index k = 0 ; k < NrFi ; ++k )
   if( ( ! CurrNrEvls[ k ] ) && ( whisG1[ k ] < InINF ) )
    ScPr1 += Master->ReadGid( whisG1[ k ] );

 }  // end( UpdateHeuristicInfo )

/*--------------------------------------------------------------------------*/
/* Front-end for four different heuristics for short-term t management.
 *
 * Three of the heuristics (1, 2, and 3) are based on three slightly different
 * variants of the same idea: considering the objective f( x ) from the
 * current stability center \bar{x} along direction d* seen as a function of
 * t, assuming that d* = - t z* (which is, strictly speaking, only true in
 * the pure quadratic proximal case and therefore does not cleanly generalise
 * to the generalised one; yet these are heuristics).
 *
 * That is, we consider the translated function along z*
 *
 *    q( v ) = f( \bar{x} - v z* ) - f( \bar{x} )
 *
 * Assuming (only for notational simplicity) differentiability, we thus have
 *
 *    q'( v ) = < - z* , f'( \bar{x} - v z* ) >
 *
 * After that f( \bar{x} + d* ) = f( \bar{x} - t z* ) = q( t ) has been
 * computed, we know:
 *
 * - the aggregated subgradient z*, which is a Sigma*-subgradient in \bar{x}
 *
 * - the newly obtained subgradient g, which is an Alfa1-subgradient in 0 and
 *   eps-subgradient in t, where
 *
 *     eps = DeltaFi - ( Alfa1 + < g , d* > )
 *
 * We thus assume:
 *
 * - q( 0 ) = 0
 *
 * - q'( 0 ) = < - z* , z* > = - NrmZ^2
 *
 * Note, however, that z* is a Sigma*-subgradient in \bar{x}, and therefore
 * the value of the linearization there is rather -Sigma*; thus, we could
 * alternatively assume q( 0 ) = - Sigma*.
 *
 * - q( t ) = f( \bar{x} - t z* ) - f( \bar{x} ) = DeltaFi
 *
 * - q'( t ) = < - z* , g >; since we have ScPr1 = < d* , g > =
 *   < - t z* , g > (note again that this only holds in the quadratic case),
 *   we conclude q'( t ) = ScPr1 / t
 *
 * Note, however, that g* is a Alfa1-subgradient of \bar{x}, and therefore we
 * could alternatively take the value in t as that of the corresponding
 * linearization, i.e., q( t ) = ScPr1 - Alfa1.
 *
 * We can then consider the quadratic function
 *
 *    m( v ) = a v^2 + b v + c
 *
 * and construct different forms of it corresponding to different choices of
 * three of the four information we have, then use its minimum
 *
 *   v* = - b / ( 2 a )
 *
 * as the suggested new value for t. Since v* only makes sense if a > 0,
 * when a <= 0 we use the best possible convex approximation of a concave
 * function by setting a = 0, in which case the minimum is the extreme of
 * the interval [ tMinor , tMaior ] dictated by the sign of b.
 *
 * The fourth heuristic is based on an entirely different idea related to the
 * Moreau-Yoshida regularization, called "reversal form of the poorman's
 * quasi-Newton update". */

HpNum BundleSolver::Heuristic( Index whch )
{
 switch( whch & 3 ) {
  case( 0 ): return( Heuristic1() );
  case( 1 ): return( Heuristic2() );
  case( 2 ): return( Heuristic3() );
  }

 return( Heuristic4() );
 }

/*--------------------------------------------------------------------------*/
/* With the notation above, in the first case we impose
 *
 *    m( 0 )  = c = - Sigma*
 *    m( t )  = a t^2 + b t + c = DeltaFi
 *    m'( 0 ) = [ 2 a 0 ] + b = - NrmZ^2
 *
 * which yields c = - Sigma*, b = - NrmZ^2,
 * a = ( DeltaFi + NrmZ^2 t + Sigma*  ) / t^2.
 * Note that z* is a Sigma*-subgradient in \bar{x}, and therefore
 *
 *    f( \bar{x} + d* ) >= f( \bar{x} ) + < d* , z* > - Sigma*
 *                       = f( \bar{x} ) - t < z* , z* > - Sigma*
 *    ==> DeltaFi = f( \bar{x} + d* ) - f( \bar{x} ) >= - NrmZ^2 t - Sigma*
 *    ==> a = DeltaFi + NrmZ^2 t + Sigma* >= 0
 *
 * which guarantees that m() is convex and therefore the minimum of m() is
 *
 *   v* = - ( - NrmZ^2 ) / ( 2 ( DeltaFi + NrmZ^2 t + Sigma* ) / t^2 )
 *      =   NrmZ^2 / ( 2 ( DeltaFi + NrmZ^2 t + Sigma* ) / t^2 )
 *      =   t^2 NrmZ^2 / ( 2 ( DeltaFi + NrmZ^2 t + Sigma* ) )
 *
 * Note that, conveniently, v* >= 0 always holds. This corresponds to the
 * fact that m'( 0 ) = b < 0, i.e., m() is surely decreasing in 0.
 * 
 * However, this formula has a serious issue: we need to know DeltaFi,
 * which may well not be defined when a NS is performed and multiple
 * components are present since the incremental approach may stop the
 * inner loop before having computed them all. The obvious solution is
 * to replace DeltaFi with
 *
 *    \underline{f}( \bar{x} + d* ) - \bar{f}( \bar{x} ) =
 *    LwFiLmb1.back() - UpRifFi.back()
 *
 * which is always well-defined since a finite lower bound is always
 * available (unless some component evaluates to -INF, in which case
 * the algorithm stops and this method is not invoked). */

HpNum BundleSolver::Heuristic1( void )
{
 auto DF = DeltaFi < INFshift ? DeltaFi : LwFiLmb1.back() - UpRifFi.back();
 auto NZ2 = NrmZ * NrmZ;
 if( DF + NZ2 * t + Sigma > 1e-16 )  // this should always be >=
  return( t * t * NZ2 / ( 2 * ( DF + NZ2 * t + Sigma ) ) );
 else                 // a == 0, all depends on the sign of b
  /* there is no "if" here, NrmZ >= 0 by definition, a fortiori NZ2
  if( - NZ2 <= 0 )    // b < 0  */
   return( tMaior );  // ==> tMaior
  /* there is no "else" here, - NZ2 > 0 cannot happen
  else                // b > 0
   return( tMinor );  // ==> tMinor */
 }

/*--------------------------------------------------------------------------*/
/* With the notation above, in the second case we rather impose
 *
 *    m( 0 )  = c = 0
 *    m( t )  = a t^2 + b t [ + 0 ] = ScPr1 - Alfa1
 *    m'( t ) = 2 a t + b = ScPr1 / t
 *
 * which yields c = 0, a = Alfa1 / t^2, b = ( ScPr1 - 2 Alfa1 ) / t
 *
 * Since Alfa1 >= 0, a >= 0 which implies that m() is surely convex and the
 * minimum is
 * 
 *   v* = - [ ( ScPr1 - 2 Alfa1 ) / t ] / [ 2 Alfa1 / t^2 ]
 *      = t ( 2 Alfa1 - ScPr1 ) / ( 2 Alfa1 )
 *
 * Note that, unlike in the first case, there is no guarantee that v* >= 0,
 * because we fix the derivative in t and therefore m'( 0 ) = b may turn up
 * to be positive (m() in increasing in 0).
 *
 * Since this formula does not really use DeltaFi, it being undefined is not
 * an issue here. However, this formula has a somewhat similar issue with NS
 * (and SS alike) in that not all components may have been evaluated, and
 * therefore only a "partial" g may be available. ScPr1 and Alfa1 are
 * computed for all non-"easy" components using the "representative
 * subgradients" out of the previous iterations, *provided they have not by
 * chance been deleted* (which should not happen unless the bundle is very
 * very small). Yet, "easy" components are left out. There may be some way
 * put of this, e.g. by using z*_i in place of g_i for the "easy" components,
 * but this is nontrivial and therefore avoided for now. */

HpNum BundleSolver::Heuristic2( void )
{
 if( Alfa1 > 1e-16 )            // it is always >= 0, but it may be ==
  return( t * ( 2 * Alfa1 - ScPr1 ) / ( 2 * Alfa1 ) );
 else                           // a == 0,  all depends on the sign of b
  if( ScPr1 - 2 * Alfa1 <= 0 )  // b < 0
   return( tMaior );            // ==> tMaior
  else                          // b > 0
   return( tMinor );            // ==> tMinor
 }

/*--------------------------------------------------------------------------*/
/* With the notation above, in the third case we rather impose
 *
 *    m'( 0 ) = [ 2 a 0 ] + b = - NrmZ^2
 *    m( t )  = a t^2 + b t + c = < something >
 *    m'( t ) = 2 a t + b = ScPr1 / t
 *
 * which yields b = - NrmZ^2, a  = ( ScPr1 / t + NrmZ^2 ) / ( 2 t ), and c
 * ... something depending on which value we choose for m( t ), but we
 * don't care about because c does not appear in the computation of v*.
 * If m() is convex, i.e.
 *
 *    ScPr1 / t + NrmZ^2 > 0
 *
 * yields
 * 
 *   v* = - [ - NrmZ^2 ] / [ 2 ( ScPr1 / t + NrmZ^2 ) / ( 2 t ) ]
 *      = t NrmZ^2 / ( ScPr1 / t + NrmZ^2 )
 *
 * Note that, if m() is convex, then v* >= 0 holds because, as usual, we have
 * fixed m'( 0 ) = b < 0 and therefore m() is decreasing in 0. 
 *
 * See above for the "issue" about ScPr1 having been computed with a
 * "partial" g; however, since this formula does not really use DeltaFi,
 * it being undefined is not an issue here.
 *
 * Note also that the possible fourth case
 *
 *    m( 0 )  = c = 0 [ or - Sigma* ]
 *    m'( 0 ) = [ 2 a 0 ] + b = - NrmZ^2
 *    m'( t ) = 2 a t + b = ScPr1 / t
 *
 * only changes c w.r.t. the current one, hence it does not change v*, and
 * therefore need not be separately considered. */

HpNum BundleSolver::Heuristic3( void )
{
 auto NZ2 = NrmZ * NrmZ;
 if( ScPr1 / t + NZ2 > 1e-16 )
  return( t * NZ2 / ( ScPr1 / t + NZ2 ) );
 else                 // a == 0, all depends on the sign of b
  /* there is no "if" here, NrmZ >= 0 by definition, a fortiori NZ2
  if( - NZ2 <= 0 )    // b < 0  */
   return( tMaior );  // ==> tMaior
  /* there is no "else" here, - NZ2 > 0 cannot happen
  else                // b > 0
   return( tMinor );  // ==> tMinor */
 }

/*--------------------------------------------------------------------------*/
/* This heuristic is instead based on a completely different approach. It is
 * called "reversal form of the poorman's quasi-Newton update" and its
 * nontrivial rationale is described in details in
 *
 *  C. Lemarechal and C. Sagastizabal. Variable metric bundle methods: from
 *  conceptual to implementable forms. Mathematical Programming,
 *  76(3):393-410, 1997
 *
 * A more refined version of the same is proposed in
 *
 *  P.A. Rey and C. Sagastizabal. Dynamical adjustment of the prox-parameter
 *  in variable metric bundle methods. Optimization, 51(2):423-447, 2002
 *
 * It should be noted that this heuristic is explicitly developed for being
 * used at SS only.
 *
 * The proposed new value is
 *
 *   t = < v , u > / || v ||^2
 *
 * where
 *
 *   v = g - z*
 *
 *   u = ( \bar{x} - v z* ) - \bar{x} + t v = d* + t v
 *
 * although v would in general be g_{i+1} - g_i, hence the choice of z* as
 * g_i is somewhat arbitrary; but in general z* is considered to be "the best
 * (approximate) subgradient we have at \bar{x}".
 *
 * Hence
 *
 *   v = < v , d* + t v > / || v ||^2
 *     = [ < v , d* > + t < v , v > ] / || v ||^2
 *     = < g - z* , d* > / || g - z* ||^2 + t
 *     = t + [ < g , d* > + t || z* ||^2 ] /
 *           [ || g ||^2  - 2 < g , z* > + || z* ||^2 ]
 *     = t + [ < g , d* > + t || z* ||^2 ] /
 *           [ || g ||^2  + 2 < g , d* > / t + || z* ||^2 ]
 *
 * The issue with this formula is the || g || term. This is the same issue as
 * with ScPr1 = < g , d* >, i.e., what to do with the easy components which
 * do not explicitly compute a subgradient. For the scalar product we could
 * use < z*_i , d* > that should be is available "for free" out of the Master
 * Problem but it currently isn't; in theory z*_i is also available, and we
 * could use it to compute < z*_i , d* >, but in practice due to the current
 * implementation of OSIMPSolver it is too costly to compute. */

HpNum BundleSolver::Heuristic4( void )
{
 auto NZ2 = NrmZ * NrmZ;
 if( G1Norm == INFshift )
  G1Norm = norm( G1 , 2 );
 return( t + ( ScPr1 + t * NZ2 ) / ( G1Norm + 2 * ScPr1 / t + NZ2 ) );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::guts_of_destructor( void )
{
 if( Master ) {
  Master->SetDim();
  #if( USE_MPTESTER )
   auto t = dynamic_cast< MPTester * >( Master );
   assert( t );
   #if( USE_MPTESTER == 1 )
    auto o = dynamic_cast< OSIMPSolver * >( t->get_master() );
   #else
    auto o = dynamic_cast< OSIMPSolver * >( t->get_slave() );
   #endif
   assert( o );
   o->SetOsi();
  #else
   if( MPName & 1 ) {
    auto o = dynamic_cast< OSIMPSolver * >( Master );
    assert( o );
    o->SetOsi();
    }
  #endif

  delete Master;
  Master = nullptr;
  }

 whisG1.clear();
 vStar.clear();

 LowerBound.clear();
 LwFiLmb.clear();
 UpFiLmb.clear();
 LwFiLmb1.clear();
 UpFiLmb1.clear();
 UpRifFi.clear();

 FiStatus.clear();

 CurrNrEvls.clear();

 Zvalid.clear();
 whisZ.clear();
 FreList = {};

 MaxItem.clear();
 FrFItem.clear();
 NrItems.clear();

 ItemVcblr.clear();

 OOBase.clear();

 LmbdBst.clear();
 Lambda1.clear();
 Lambda.clear();

 InvItemVcblr.clear();
 vBPar2.clear();

 if( NrEasy ) {  // if there are "easy" components, delete the MILPSolver
  if( DoEasy & ~1 ) {
   // if easy components can be changed, before doing this unregister the
   // MILPSolver from the inner Block (since it is registered there);
   // meanwhile unregister the FakeSolver that is also registered there
   auto FSit = v_FakeSolver.begin();
   for( Index k = 0 ; k < NrFi ; ++k )
    if( IsEasy[ k ] ) {
     auto iB = static_cast< LagBFunction * >(
					   v_c05f[ k ] )->get_inner_block();
     iB->unregister_Solver( *(FSit++) , true );
     iB->unregister_Solver( IsEasy[ k ] , true );
     }

   v_FakeSolver.clear();
   }
  else  // "easy" components are static, just delete the MILPSolver
   for( Index k = 0 ; k < NrFi ; ++k )
    delete IsEasy[ k ];

  IsEasy.clear();
  NrEasy = 0;
  }

 LamVcblr.clear();

 v_c05f.clear();

 }  // end( BundleSolver:guts_of_destructor )

/*--------------------------------------------------------------------------*/

void BundleSolver::ReSetAlg( unsigned char RstLvl )
{
 if( ! ( RstLvl & RstAlg ) ) {  // reset algorithmic parameters - - - - - - -
  ParIter = 0;             // reset iterations count
  CSSCntr = CNSCntr = 0;   // ... comprised consecutive NS/SS count

  if( t != tInit ) {       // reset t
   t = tInit;
   tHasChgd = true;
   }

  // reset the dynamic number of fetched items
  if( BPar6 && ( BPar5 > 0 ) )
   aBP3 = BPar4;
  else
   aBP3 = BPar3;
  }

 // reset function values in Lambda, since they are changing
 std::fill( UpFiLmb.begin() , UpFiLmb.end() ,  INFshift );  // upper
 std::fill( LwFiLmb.begin() , LwFiLmb.end() , -INFshift );  // lower
 UpFiLmbdef = LwFiLmbdef = 0;             // ... at the current point
 f_global_LB = -INFshift;         // algorithmic global LB

 if( RstLvl & RstCrr ) {  // get an initial point - - - - - - - - - - - - - -
  // note that the thusly constructed Master Problem assumes the stability
  // centre to be all-0, in particular when constructing the Lagrangian costs
  // of the easy components, so if that's not happening the right value has
  // to be explicitly passed

  for( Index i = 0 ; i < NumVar ; ++i )
   Lambda[ i ] = LamVcblr[ i ]->get_value();

  bool nonzero = false;
  for( Index i = 0 ; i < NumVar ; ++i )
   if( Lambda[ i ] ) { nonzero = true; break; }

  if( nonzero ) {
   // call ChangeCurrPoint( DLambda ) with DLambda == Lambda - oldLambda ==
   // Lambda since oldLambda == 0 by construction
   Vec_VarValue foo( NrFi + 1 , 0 );  // no change in the (unknown) f-values
   Master->ChangeCurrPoint( Lambda.data() , foo.data() );
   Fi0Lmb = INFshift;  // the value of the linear part must be computed
   }
  else  // Lambda was all-0 anyway
   Fi0Lmb = 0;  // then the value of the linear part is quite obvious ...
  }
 else {                   // reset the current point to all-0 - - - - - - - -
  // note that the thusly constructed Master Problem precisely assumes
  // the stability centre to be all-0, in particular when constructing the
  // Lagrangian costs of the easy components, so that's fine as it is
  Lambda.assign( NumVar , 0 );
  // "tell" this to the ColVariable of the C05Function(s)
  for( Index i = 0 ; i < NumVar ; ++i )
   LamVcblr[ i++ ]->set_value( 0 );
  Fi0Lmb = 0;  // then the value of the linear part is quite obvious ...
  }
 }  // end( BundleSolver::ReSetAlg )

/*--------------------------------------------------------------------------*/

void BundleSolver::Delete( cIndex i , bool ModDelete )
{
 // deletes from the bundle the item in position i
 //
 // ModDelete == true means that this is called in response of a Modification
 // where the linearization has been removed from the global pool
 //
 // whatever BPar7 says, no linearization is physically deleted here inside,
 // this has to be done by the caller (if needed)

 cIndex k = ItemVcblr[ i ].first;

 // check if this item was the "representative" for its component - - - - - -

 if( whisG1[ k ] == i )  // it is the representative of k
  whisG1[ k ] = InINF;   // a new representative is needed

 // check if this item was the z* for its component - - - - - - - - - - - - -

 if( whisZ[ k ] == i ) {  // it is the aggregate subgradient of k
  whisZ[ k ] = InINF;     // no aggregate subgradient is in the bundle
  Zvalid[ k ] = false;    // a fortiori, no valid one
  }

 // delete the item from the MP - - - - - - - - - - - - - - - - - - - - - - -

 Master->RmvItem( i );

 BLOG( 2 , std::endl << "Item " << i << " removed" );

 // bookkeeping of internal data structures - - - - - - - - - - - - - - - - -
 // note that any item whose name is >= Master->MaxName() is surely not in
 // the bundle (master problem), and therefore it need not be in FreList

 cIndex MxNm = Master->MaxName();
 if( i < MxNm )
  FreList.push( i );

 OOBase[ i ] = Inf< SIndex >();
 --NrItems[ k ];
 --NrItems[ NrFi ];

 // remove from the global pool: the removal is "hard" if either BPar7 says
 // so, or the linearization had been deleted anyway

 remove_from_global_pool( k , ItemVcblr[ i ].second ,
			  ( ( BPar7 & 3 ) == 3 ) || ModDelete );
 ItemVcblr[ i ].second = InINF;

 // check if compacting FreList is appropriate- - - - - - - - - - - - - - - -
 // the issue with having indices of "free" position in the bundle stored in
 // a priority_queue is the following: if the bundle gets "full", but then is
 // "emptied", FreList may end up containing "many" elements, and in
 // particular elements that are >= Master->MaxName(), which therefore are
 // useless since they are obviously not in the bundle. the check above tries
 // to avoid that, but it may clearly fail (say, if small items are deleted
 // before large ones). checking if there are items with name >=
 // Master->MaxName() in FreList and deleting them is not cheap. the only
 // easy-to-check case is the one where FreList.size() > Master->MaxName():
 // if this happens, FreList is cleared and re-initialized

 if( FreList.size() > MxNm ) {
  FreList = {};
  for( Index h = 0 ; h < MxNm ; ++h )
   if( ItemVcblr[ h ].second == InINF )
    FreList.push( h );
  }
 }  // end( BundleSolver::Delete )

/*--------------------------------------------------------------------------*/

void BundleSolver::UpdtaBP3( void )
{
 if( BPar5 == 0 )
  return;

 switch( BPar6 ) {
  case( 4 ):
   if( UpFiLmb[ NrFi ] > -INFshift )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 / std::log10( EpsU / RelAcc ) );
   break;
  case( 3 ):
   if( UpFiLmb[ NrFi ] > -INFshift )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 / std::sqrt( EpsU / RelAcc ) );
   break;
  case( 2 ):
   if( UpFiLmb[ NrFi ] > -INFshift )
    aBP3 = ( BPar5 > 0 ? BPar4 : BPar3 ) +
           Index( BPar5 * ( RelAcc / EpsU ) );
   break;
  case( 1 ):
   if( ! ( ParIter % Index( std::abs( BPar5 ) ) ) ) {
    if( BPar5 > 0 )
     aBP3++;
    else
     aBP3--;
    }
  }

 if( aBP3 > Index( BPar3 ) )
  aBP3 = BPar3;
 else
  if( aBP3 < BPar4 )
   aBP3 = BPar4;

 }  // end( BundleSolver::UpdtaBP3 )

/*--------------------------------------------------------------------------*/

bool BundleSolver::IsOptimal( double eps ) const
{
 if( ! RifeqFi )    // the linearization errors are not "properly computed"
  return( false );  // no way one can detect optimality

 if( eps <= 0 )
  eps = RelAcc;

 if( vStar.back() >= INFshift )  // some components have no subgradients
  return( false );               // no way one can detect optimality

 c_VarValue err = max_error( eps );
 if( err >= INFshift )
  return( false );

 if( ( tStar > 0 ) && ( DSTS + Sigma <= err ) )
  return( true );

 return( ( Sigma <= err ) &&
	 ( NrmZFctr < INFshift ) && ( NrmZ <= NrmZFctr * NZEps ) );

 }  // end( BundleSolver::IsOptimal )

/*--------------------------------------------------------------------------*/

void BundleSolver::FModChg( VarValue shift , Index wFi )
{
 if( ( ! std::isnan( shift ) ) && ( ! f_convex ) )
  shift = - shift;

 if( shift == INFshift ) {      // function changed monotonically up
  if( UpFiLmb[ wFi ] < INFshift ) {
   UpFiLmb[ wFi ] = INFshift;   // reset upper function value for component
   --UpFiLmbdef;                // one less known
   }
  if( UpFiLmb.back() < INFshift ) {
   UpFiLmb.back() = INFshift;   // reset total upper function value
   --UpFiLmbdef;                // one less known
   }
  UpFiBest = INFshift;          // comprised best one
  return;
  }

 if( shift == -INFshift ) {     // function changed monotonically dn
  if( LwFiLmb[ wFi ] > -INFshift ) {
   LwFiLmb[ wFi ] = -INFshift;  // reset lower function value for component
   --LwFiLmbdef;                // one less known
   }
  if( LwFiLmb.back() > -INFshift ) {
   LwFiLmb.back() = -INFshift;  // reset total lower function value
   --LwFiLmbdef;                // one less known
   }
  f_global_LB = -INFshift;      // global LB no longer valid
  return;
  }

 if( std::isnan( shift ) ) {    // function changed unpredictably
  if( UpFiLmb[ wFi ] < INFshift ) {
   UpFiLmb[ wFi ] = INFshift;   // reset upper function value for component
   --UpFiLmbdef;                // one less known
   }
  if( UpFiLmb.back() < INFshift ) {
   UpFiLmb.back() = INFshift;   // reset total upper function value
   --UpFiLmbdef;                // one less known
   }
  UpFiBest = INFshift;          // and of course best one
  if( LwFiLmb[ wFi ] > -INFshift ) {
   LwFiLmb[ wFi ] = -INFshift;  // reset lower function value for component
   --LwFiLmbdef;                // one less known
   }
  if( LwFiLmb.back() > -INFshift ) {
   LwFiLmb.back() = -INFshift;  // reset total lower function value
   --LwFiLmbdef;                // one less known
   }
  f_global_LB = -INFshift;      // global LB no longer valid
  return;
  }

 // function changed by shift: just update everything

 if( UpFiLmb[ wFi ] < INFshift )
  UpFiLmb[ wFi ] += shift;

 if( UpFiLmb.back() < INFshift )
  UpFiLmb.back() += shift;

 if( UpRifFi[ wFi ] < INFshift )
  UpRifFi[ wFi ] += shift;

 if( UpRifFi.back() < INFshift )
  UpRifFi.back() += shift;

 if( UpFiBest < INFshift )
  UpFiBest += shift;

 if( LwFiLmb[ wFi ] > -INFshift )
  LwFiLmb[ wFi ] += shift;

 if( LwFiLmb.back() > -INFshift )
  LwFiLmb.back() += shift;

 f_global_LB += shift;

 }  // end( BundleSolver::FModChg )

/*--------------------------------------------------------------------------*/

void BundleSolver::remove_from_global_pool( Index k , Index i , bool hard )
{
 // update InvItemVcblr and all the associated fields to the fact that the
 // linearization currently in position i of the global pool of component k
 // is removed; actually, there may as well not be any linearization in
 // position i already
 //
 // the removal is "hard" (meaning the linearization is actually deleted) if
 // hard == true, and "soft" (the linearization is kept in the global pool
 // but deleted from the bundle) otherwise
 //
 // the actual removal from the global pool is not handled here

 InvItemVcblr[ k ][ i ] = hard ? InINF : vBPar2.back();
 while( MaxItem[ k ] && ( InvItemVcblr[ k ][ MaxItem[ k ] - 1 ] == InINF ) )
  --MaxItem[ k ];
 if( i < FrFItem[ k ] )   // creating a new "hole" before the FrFItem
  FrFItem[ k ] = i;       // this is the new FrFItem
 else                     // deleting something that may be FrFItem
  while( FrFItem[ k ] &&
	 ( InvItemVcblr[ k ][ FrFItem[ k ] - 1 ] >=
	                        ( ( BPar7 & 3 ) ? vBPar2.back() : InINF ) ) )
   --FrFItem[ k ];

 }  // end( BundleSolver::remove_from_global_pool )

/*--------------------------------------------------------------------------*/

Index BundleSolver::find_place_in_global_pool( Index k )
{
 // returns a suitable position in the global pool of component k, or InINF
 // if there is no free space

 if( FrFItem[ k ] < vBPar2[ k ] )  // there are free positions
  return( FrFItem[ k ] );          // return the first of them

 // if there are no free positions but there are less items in the bundle
 // (for component k) than the size of the global pool, there must be a
 // some linearization in the global pool that is not in the bundle: find
 // and return the first one (smallest position)
 Index gpp = InINF;
 if( NrItems[ k ] < vBPar2[ k ] )
  for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
   if( InvItemVcblr[ k ][ i ] >= vBPar2.back() ) {
    gpp = i;
    break;
    }

 return( gpp );

 }  // end( BundleSolver::find_place_in_global_pool )

/*--------------------------------------------------------------------------*/

void BundleSolver::add_to_global_pool( Index k , Index i , Index wh )
{
 // update ItemVcblr, InvItemVcblr and all the associated fields the to fact
 // that the linearization in position i in the global global_pool of
 // component k will be kept in the bundle at position wh; if wh == InINF,
 // this means the linearization is in the global pool but not in the bundle
 //
 // the actual addition to the global pool is not handled here

 if( wh < InINF ) {
  ItemVcblr[ wh ].first = k;
  ItemVcblr[ wh ].second = i;
  InvItemVcblr[ k ][ i ] = wh;
  }
 else
  InvItemVcblr[ k ][ i ] = vBPar2.back();

 if( i >= MaxItem[ k ] )
  MaxItem[ k ] = i + 1;

 while( ( FrFItem[ k ] < MaxItem[ k ] ) &&
	( InvItemVcblr[ k ][ FrFItem[ k ] ] <
	                        ( ( BPar7 & 3 ) ? vBPar2.back() : InINF ) ) )
  ++FrFItem[ k ];

 }  // end( BundleSolver::add_to_global_pool( k , i , wh ) )

/*--------------------------------------------------------------------------*/

void BundleSolver::add_to_bundle( Index k , Index i )
{
 // add to the bundle (master problem) the item corresponding to the
 // linearization to be found at position i in the global pool of component
 // k; this assumes that the linearization is already there in the global
 // pool. if InvItemVcblr[ k ][ i ] < vBPar2.back(), i.e., the item is
 // already in the bundle, then it is replaced, otherwise it is added
 //
 // note that CheckSubG() or CheckCnst() need be called, but even if the
 // item is identical to some in the bundle already this information is
 // ignored and the item is inserted anyway; hence, if the check is active,
 // it is temporarily deactivated (and then re-activated)

 auto wh = InvItemVcblr[ k ][ i ];
 if( wh >= vBPar2.back() ) {  // the item is not there already
  wh = FindAPlace( k );       // find a "free" spot in the bundle
  if( wh == InINF )           // one must be there
   throw( std::logic_error( "no space found in the bundle" ) );

  ++NrItems[ k ];              // keep count
  ++NrItems[ NrFi ];
  add_to_global_pool( k , i , wh );  // update dictionaries
  }
 else                          // the item is there already
  Master->RmvItem( wh );       // remove it so that it can be replaced

 // ask the MPSolver for the space to write the item to
 auto G1 = Master->GetItem( k + 1 );

 // recover the linearization from the C05Function
 v_c05f[ k ]->get_linearization_coefficients( G1 , Range( 0 , NumVar ) , i );
 if( ! f_convex )
  chgsign( G1 , NumVar );

 // recover the constant and "translate" it w.r.t. Lambda
 auto Ai = rs( v_c05f[ k ]->get_linearization_constant( i ) );

 Master->SetItemBse( nullptr , NumVar );

 if( MPName & 8 )                   // if checking for copies is active
  Master->CheckIdentical( false );  // temporarily de-activate it now

 double ScPri;
 if( v_c05f[ k ]->is_linearization_vertical( i ) )
  Master->CheckCnst( Ai , ScPri , Lambda.data() );
 else {
  Ai = UpRifFi[ k ] - Ai -
   std::inner_product( Lambda.begin() , Lambda.end() , G1 , double( 0 ) );
  Master->CheckSubG( 0 , 0 , Ai , ScPri );
  }

 if( MPName & 8 )                   // if checking for copies is active
  Master->CheckIdentical();         // re-activate it now

 Master->SetItem( wh );  // add the item to the master problem

 }  // end( BundleSolver::add_to_bundle )

/*--------------------------------------------------------------------------*/

void BundleSolver::reset_bundle( void )
{
 // completely resets the bundle, because a (bunch of) Modification(s) saying
 // so has(ve) been received. this only affects the BundleSolver data
 // structures and the MPSolver, not the C05Function(s)

 OOBase.assign( vBPar2.back() , Inf< SIndex >() );

 ItemVcblr.assign( vBPar2.back() , make_pair( InINF , InINF ) );

 for( Index k = 0 ; k < NrFi ; ++k )
  InvItemVcblr[ k ].assign( vBPar2[ k ] , InINF );

 NrItems.assign( NrFi + 1 , 0 );
 FrFItem.assign( NrFi , 0 );
 MaxItem.assign( NrFi , 0 );

 FreList = {};
 whisZ.assign( NrFi , InINF );
 Zvalid.assign( NrFi , false );

 whisG1.assign( NrFi , InINF );

 Master->RmvItems();

 }  // end( BundleSolver::reset_bundle )

/*--------------------------------------------------------------------------*/

Lst_sp_Mod::size_type BundleSolver::num_outstanding_Modification( void )
{
 auto res = v_mod.size();

 if( NrEasy && ( DoEasy & ~1 ) ) {
  auto FSit = v_FakeSolver.begin();
  for( Index k = 0 ; k < NrFi ; )
   if( IsEasy[ k++ ] )
    res += (*FSit)->get_Modification_list().size();
  }

 return( res );
 }

/*--------------------------------------------------------------------------*/

bool BundleSolver::is_special_GroupMod( GroupModification & gmod )
{
 // recognise "special" GroupModification for changing the set of "active"
 // Variable of all the Objective at the same time; note that these
 // contain FunctionModVars* not necessarily C05FunctionModVars* because
 // the Modification may not be strongly quasi-additive

 if( gmod.sub_Modifications().size() != NrFi + ( f_lf ? 1 : 0 ) )
  return( false );

 auto smi = gmod.sub_Modifications().begin();
 auto sm0 = *(smi++);
 for( ; smi !=  gmod.sub_Modifications().end() ; ++smi )
  if( typeid( sm0 ) != typeid( *smi ) )
   return( false );

 smi = gmod.sub_Modifications().begin();
 ++smi;

 // check FunctionModVarsAddd
 if( const auto mod0 =
     std::dynamic_pointer_cast< FunctionModVarsAddd >( sm0 ) ) {
  for( ; smi != gmod.sub_Modifications().end() ; ++smi ) {
   auto modi = std::static_pointer_cast< FunctionModVarsAddd >( *smi );
   if( ( mod0->first() != modi->first() ) ||
       ( mod0->vars() != modi->vars() ) )
    throw( std::logic_error( "different Variable change in components" ) );
   }

  return( true );
  }

 // check FunctionModVarsRngd
 if( const auto mod0 =
     std::dynamic_pointer_cast< FunctionModVarsRngd >( sm0 ) ) {
  for( ; smi != gmod.sub_Modifications().end() ; ++smi ) {
   auto modi = std::static_pointer_cast< FunctionModVarsRngd >( *smi );
   if( mod0->range() != modi->range() )
    throw( std::logic_error( "different Variable change in components" ) );
   }

  return( true );
  }

 // check FunctionModVarsSbst
 if( const auto mod0 =
     std::dynamic_pointer_cast< FunctionModVarsSbst >( sm0 ) ) {
  for( ; smi != gmod.sub_Modifications().end() ; ++smi ) {
   auto modi = std::static_pointer_cast< FunctionModVarsSbst >( *smi );
   if( mod0->subset() != modi->subset() )
    throw( std::logic_error( "different Variable change in components" ) );
   }

  return( true );
  }

 return( false );

 }  // end( BundleSolver::is_special_GroupMod )

/*--------------------------------------------------------------------------*/

void BundleSolver::flatten_Modification_list( Lst_sp_Mod & vmt , sp_Mod mod )
{
 const auto tmod = std::dynamic_pointer_cast< GroupModification >( mod );
 if( tmod && ( ! is_special_GroupMod( *tmod ) ) )
  for( auto submod : tmod->sub_Modifications() )
   flatten_Modification_list( vmt , submod );
 else
  vmt.push_back( mod );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::flatten_easy_Modification_list( Lst_sp_Mod & vmt ,
						   sp_Mod mod )
{
 if( const auto tmod = std::dynamic_pointer_cast< GroupModification >( mod ) )
  for( auto submod : tmod->sub_Modifications() )
   flatten_Modification_list( vmt , submod );
 else
  vmt.push_back( mod );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::process_outstanding_easy_Modification( void )
{
 // look at all easy components; for each of these scan the Modification in
 // the corresponding FakeSolver and divine which changes need be done to
 // the Master Problem

 auto FSit = v_FakeSolver.begin();
 for( Index k = 0 ; k < NrFi ; ++k ) {
  if( ! IsEasy[ k ] )
   continue;

  // ensure that the MILPSolver has "digested" the Modification - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // since some Modification in the inner Block cause LagBFunction to react
  // by doing other changes to it, ensure that this is done before the list
  // of Modification is scanned; however, note that if the list is empty
  // already (in which case this method may not even be called, and therefore
  // neither apply_obj_Modification() is) then no new Modification can be
  // added by LagBFunction: LagBFunction does not "create Modification out
  // of thin air", only reacts to those issued from outside

  static_cast< LagBFunction * >( v_c05f[ k ] )->apply_obj_Modification();

  // MILPSolver::compute() has the only role of scanning the list of
  // Modification and updating its internal data structures accordingly

  IsEasy[ k ]->MILPSolver::compute( false );

  // construct the flattened list of Modification in the FakeSolver - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Lst_sp_Mod v_mod_tmp;

  (*FSit)->lock_Modification_list();

  for( auto mod : (*FSit)->get_Modification_list() )
   flatten_easy_Modification_list( v_mod_tmp , mod );

  (*FSit)->get_Modification_list().clear();

  (*(FSit++))->unlock_Modification_list();

  if( v_mod_tmp.empty() )  // nothing here to see
   continue;               // move over

  // scan all the Modification in the FakeSolver- - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  char whch = 0;

  for( ; ! v_mod_tmp.empty() ; v_mod_tmp.pop_front() ) {
   auto mod = v_mod_tmp.front().get();

   // FunctionMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( const auto tmod = dynamic_cast< const FunctionMod * >( mod ) ) {
    auto obs = tmod->function()->get_Observer();
    if( dynamic_cast< const Objective * >( obs ) ) {

     // C05FunctionModLin - - - - - - - - - - - - - - - - - - - - - - - - - -
     if( auto modl = dynamic_cast< const C05FunctionModLin * >( tmod ) ) {

      whch |= 1;
      continue;
      }

     // C05FunctionMod- - - - - - - - - - - - - - - - - - - - - - - - - - - -
     if( auto modl = dynamic_cast< const C05FunctionMod * >( tmod ) ) {

      if( modl->type() == C05FunctionMod::NothingChanged ) {

       const auto shift = modl->shift();

       if( ( shift ==  INFshift ) ||
           ( shift == -INFshift ) ||
           ( std::isnan( shift ) ) )
        throw( std::logic_error(
         "unexpected *C05FunctionMod* from Objective Function" ) );

       constant_value += shift;
       }
      else {

       whch |= 1;
       continue;
       }
      }
     }

    throw( std::logic_error( "unsupported Modification in easy component" ) );
    }

   // FunctionModVars- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( const auto tmod = dynamic_cast< const FunctionModVars * >( mod ) ) {
    auto obs = tmod->function()->get_Observer();
    if( dynamic_cast< const Objective * >( obs ) ) {
     whch |= 1;
     continue;
     }

    throw( std::logic_error( "unsupported Modification in easy component" ) );
    }

   // VariableMod- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( const auto tmod = dynamic_cast< const VariableMod * >( mod ) ) {
    const auto xj = dynamic_cast< const ColVariable * >( tmod->variable() );

    if( ! xj )     // unknown variable type
     throw( std::logic_error( "unsupported Modification in easy component" ) );

    // fixing variables or changing their type is not supported; all the rest
    // boils down to a change of bounds, and therefore is

    if( ( xj->is_fixed() != xj->is_fixed( tmod->old_state() ) ) ||
        ( xj->is_integer() != xj->is_integer( tmod->old_state() ) ) )
     throw( std::logic_error( "unsupported Modification in easy component" ) );

    whch |= 4;
    continue;

    }  // end( VariableMod )

   // RowConstraintMod - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( const auto tmod = dynamic_cast< const RowConstraintMod * >( mod ) ) {

    if( ( tmod->type() != RowConstraintMod::eChgLHS ) &&
        ( tmod->type() != RowConstraintMod::eChgRHS ) &&
        ( tmod->type() != RowConstraintMod::eChgBTS ) )
     throw( std::logic_error( "unsupported Modification in easy component" ) );

    if( dynamic_cast< const OneVarConstraintMod * >( tmod ) )
     whch |= 4;
    else
     whch |= 2;

    continue;

    }  // end( RowConstraintMod )

   // any other Modification is ignored; we assume it is a "physical"
   // Modification whose corresponding "abstract" one has already been dealt
   // with or it is still in the queue waiting to be dealt with

   }  // end( for( all Modification ) )

  // finally, act upon the detected Modification- - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( whch & 1 )         // changes in costs
   Master->ChgCosts( k + 1 , Lambda.data() );

  if( whch & 2 )         // changes in LHS/RHS
   Master->ChgRLHS( k + 1 );

  if( whch & 4 )         // changes in LBD/UBD
   Master->ChgLUBD( k + 1 );

  }  // end( for( k ) )
 }  // end( process_outstanding_easy_Modification )

/*--------------------------------------------------------------------------*/

void BundleSolver::process_outstanding_Modification( void )
{
 // multiple-loop version, where several passes are done in order to gather
 // which kind of Modification have occurred and avoid doing costly work
 // more than once

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 0-th loop: "atomically flatten" v_mod into a temporary list to better
 // handle it, then clear it

 Lst_sp_Mod v_mod_tmp;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;  // try to acquire lock, spin on failure

 for( auto mod : v_mod )
  flatten_Modification_list( v_mod_tmp , mod );

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the 1st loop is made in reverse, from the latest Modification to the
 // earlies, and does the following:
 // - reset/change the upper/lower bounds that need to (check all shift() of
 //   all *FunctionMod*)
 // - check if some global pool has been "hard" reset, i.e., all the
 //   linearization in there have been deleted; this is brought about by
 //   a C05FunctionMod with type() == GlobalPoolRemoved and which().empty()
 //   or by any FunctionMod that does *not* imply strong quasi-additivity
 //   (i.e., that it is not a *C05*FunctionMod*)
 // - check the LagBFunctionMod coming out the easy components and store in
 //   reset[] the or of their what(), so that this can be acted upon later
 // - if a *FunctionMod* changing linearizations happens before a reset of the
 //   global pool (meaning it is found afterwards in the reverse order) it
 //   is deleted since it is useless (after having checked if it also impacts
 //   the upper/lower bounds)
 // - if there is more than one component, check that no "naked"
 //   *FunctionModVars* is there
 // - check that no ConstraintMod or VariableMod are there, since they are
 //   not handled (yet)
 //
 // TODO: during the 1st loop we could compute the set of components that have
 //       been modified anyhow and use this information to avoid constructing
 //       the numerous data structures like reset[] that are indexed over
 //       NrFi. this might be important if, say, NrFi is 10000 but only a
 //       smattering of the components (say, one) change

 std::vector< bool > reset( NrFi , false );

 // how many non-easy components are reset (easy never are)
 Index cntreset = 0;

 bool Fi0Chgd = false;  // true if the 0-th component changes

 bool to_delete;  // should have been defined inside, but there is not
                  // visible by the lambda

 for( auto rimod = v_mod_tmp.rbegin() ; rimod != v_mod_tmp.rend() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to rimod
      [ & to_delete , & v_mod_tmp ]( decltype( rimod ) & ri ) {
       if( to_delete )
	ri = std::reverse_iterator( v_mod_tmp.erase( std::next( ri ).base()
						     ) );
       else
	++ri;
       }( rimod ) ) {
  to_delete = false;
  auto mod = *rimod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // first check if it is any kind of FunctionMod, since this gives immediate
  // access to the component, and any FunctionMod pertaining to an already
  // reset component can be almost immediately deleted
  if( const auto tmod = std::dynamic_pointer_cast< FunctionMod >( mod ) ) {
   if( tmod->function() == f_lf ) {
    const auto shift = tmod->shift();
    // special immediate treatment of the 0-th component, which is simple
    if( std::isnan( shift ) ) {  // is a C05FunctionModLin*
     Fi0Chgd = true;             // changing the coefficients
     Fi0Lmb = INFshift;          // the value is no longer known
     if( UpFiLmbdef == NrFi + 1 ) {  // the total upper bound was known
      UpFiLmb.back() = INFshift;     // it is no longer so
      UpFiLmbdef = NrFi;             // it will have to be recomputed
      }
     if( LwFiLmbdef == NrFi + 1 ) {  // the total lower bound was known
      LwFiLmb.back() = -INFshift;    // it is no longer so
      LwFiLmbdef = NrFi;             // it will have to be recomputed
      }
     UpFiBest = INFshift;            // and the best value as well
     }
    else {                       // a FunctionMod changing the constant
     #ifndef NDEBUG
      if( ( shift == INFshift ) || ( shift == -INFshift ) )
       throw( std::logic_error( "unexpected *FunctionMod* from LinearFunction" ) );
     #endif
     if( Fi0Lmb < INFshift )
      Fi0Lmb += shift;
     if( UpFiLmb.back() < INFshift )
      UpFiLmb.back() += shift;
     if( UpRifFi.back() < INFshift )
      UpRifFi.back() += shift;
     if( UpFiBest < INFshift )
      UpFiBest += shift;
     if( LwFiLmb.back() > -INFshift )
      LwFiLmb.back() += shift;
     f_global_LB += shift;
     }

    to_delete = true;  // in either case, all that had to be done
    continue;          // has been done
    }

   auto wFi = get_index_of_component( tmod->function() );

   // adjust or reset upper/lower values as needed
   // note that the list is scanned in reverse, hence these changes are
   // applied in reverse order. however, if the upper/lower values are reset
   // at any point in the list they stay reset forever. indeed, even if a
   // function has a finite shift after a reset, this says nothing because
   // there are no known values to shift. if, rather, the values are only 
   // shifted by finite amounts, the total shift is the sum of the shift,
   // and the order of additions does not change the result
   if( wFi < NrFi )
    FModChg( tmod->shift() , wFi );
   else {
    // this is a FunctionMod coming from some unknown Function, not any
    // business of BundleSolver
    to_delete = true;
    continue;
    }

   if( NrEasy && IsEasy[ wFi ] ) {  // coming from an easy component
    if( DoEasy & ~1 ) {
     // some changes in easy components are separately supported: ignoring
     // the Modification is entirely possible
     to_delete = true;
     continue;
     }

    // changes in easy components not separately supported
    throw( std::logic_error( "unsupported change in easy component" ) );
    }
   else                             // coming from a non-easy component
    if( reset[ wFi ] ) {
     // any FunctionMod after (before) one that completely resets the
     // component is useless, delete it and move forward (backward)
     to_delete = true;
     continue;
     }

   // if the component is not reset (yet), one must look in details what
   // exact type the *FunctionMod* is and react accordingly

   // a C05FunctionModRngd only changes existing linearizations, and
   // therefore is never a "hard" reset
   if( const auto ttmod =
       std::dynamic_pointer_cast< C05FunctionModRngd >( tmod ) ) {
    switch( ttmod->type() ) {
     case( C05FunctionMod::AllLinearizationChanged ):
     case( C05FunctionMod::AllEntriesChanged ):
      continue;
     default:
      throw( std::invalid_argument( "wrong type in C05FunctionModRngd" ) );
     }  // end( switch( ttmod->f_type ) )
    }  // end( if( ttmod == C05FunctionModRngd ) )

   // a C05FunctionModSbst only changes existing linearizations, and
   // therefore is never a "hard" reset; the only easy case is
   // NothingChanged, which by definition does nothing save for the
   // shift(), that has been dealt with already
   if( const auto ttmod =
       std::dynamic_pointer_cast< C05FunctionModSbst >( tmod ) ) {
    switch( ttmod->type() ) {
     case( C05FunctionMod::AllLinearizationChanged ):
     case( C05FunctionMod::AllEntriesChanged ):
      continue;
     default:
      throw( std::invalid_argument( "wrong type in C05FunctionModSbst" ) );
     }  // end( switch( ttmod->f_type ) )
    }  // end( if( ttmod == C05FunctionModSbst ) )

   // a C05FunctionMod of type GlobalPoolRemoved with which.empty() resets
   // all the component. NothingChanged by definition does nothing (save
   // for the shift(), that has been dealt with already) ... except
   // possibly doing a lot, i.e., signalling the switch from convex to
   // concave, or vice-versa, which is not allowed; all other cases
   // will have to be dealt with later
   if( const auto ttmod =
       std::dynamic_pointer_cast< C05FunctionMod >( tmod ) ) {
    switch( ttmod->type() ) {
     case( C05FunctionMod::GlobalPoolRemoved ):
      if( ttmod->which().empty() ) {
       reset[ wFi ] = true;
       ++cntreset;
       to_delete = true;
       }
      continue;
     case( C05FunctionMod::NothingChanged ):
      if( v_c05f[ wFi ]->is_convex() != f_convex )
       throw( std::logic_error( "convex/concave switch not allowed" ) );
      to_delete = true;
     case( C05FunctionMod::AllLinearizationChanged ):
     case( C05FunctionMod::AllEntriesChanged ):
     case( C05FunctionMod::AlphaChanged ):
     case( C05FunctionMod::GlobalPoolAdded ):
      continue;
     default:
      throw( std::invalid_argument( "wrong type in C05FunctionMod" ) );
     }  // end( switch( ttmod->f_type ) )
    }  // end( if( ttmod == C05FunctionMod ) )

   // a C05FunctionModLin* only changes existing linearizations, and
   // therefore is never a "hard" reset
   if( std::dynamic_pointer_cast< C05FunctionModLinRngd >( tmod ) )
    continue;

   // if control reaches this point, mod is (indistinguishable from) a base
   // FunctionMod, and in particular it is not a C05FunctionMod*. hence the
   // change in the Function is not quasi-additive, and therefore a fortiori
   // not strongly quasi-additive. as a result, this is a "hard" reset

   reset[ wFi ] = to_delete = true;
   ++cntreset;

   }  // end( if( tmod == FunctionMod ) )

  // a "naked" FunctionModVars is only allowed if there is only one
  // component (comprised the linear one). if it is allowed, it is
  // of no consequence here, except for the possible effect on the
  // function values, if it is a C05FunctionModVars*, meaning that it
  // represents a strongly quasi-additive variable change. if not, the
  // variable change also implies a reset
  // in no case, however, the Modification is removed from the list
  if( const auto tmod = std::dynamic_pointer_cast< FunctionModVars >( mod ) ) {
   if( ( NrFi > 1 ) || f_lf )
    throw( std::invalid_argument( "naked FunctionModVars not allowed" ) );

   auto wFi = get_index_of_component( tmod->function() );

   FModChg( tmod->shift() , wFi );  // change/reset upper/lower values

   if( std::dynamic_pointer_cast< C05FunctionModVarsAddd >( tmod ) )
    continue;

   if( std::dynamic_pointer_cast< C05FunctionModVarsRngd >( tmod ) )
    continue;

   if( std::dynamic_pointer_cast< C05FunctionModVarsSbst >( tmod ) )
    continue;

   // if control reaches here, this is a FunctionModVars* that is not a
   // C05FunctionModVars*, i.e., a non strongly quasi-additive variable
   // change, which implies a "hard" reset for the component
   reset[ wFi ] = true;
   ++cntreset;
   continue;

   }  // end( if( tmod == FunctionModVars ) )

  // a GroupModification here can only be a bunch of identical
  // *FunctionModVar*: pick the first one and act on it
  if( const auto tmod =
      std::dynamic_pointer_cast< GroupModification >( mod ) ) {
   auto fmod = tmod->sub_Modifications().front();

   if( std::dynamic_pointer_cast< C05FunctionModVarsAddd >( fmod ) )
    continue;

   if( std::dynamic_pointer_cast< C05FunctionModVarsRngd >( fmod ) )
    continue;

   if( std::dynamic_pointer_cast< C05FunctionModVarsSbst >( fmod ) )
    continue;

   // if control reaches here, this is a FunctionModVars* that is not a
   // C05FunctionModVars*, i.e., a non strongly quasi-additive variable
   // change, which implies a "hard" reset for *all* components
   if( NrEasy ) {
    for( Index k = 0 ; k < NrFi ; ++k )
     if( ! IsEasy[ k ] )
      reset[ k ] = true;
    }
   else
    reset.assign( NrFi , true );
   cntreset = NrFi - NrEasy;
   continue;

   }  // end( if( tmod == GroupModification ) )

  if( std::dynamic_pointer_cast< ConstraintMod >( mod ) )
   throw( std::invalid_argument( "ConstraintMod not handled (yet)" ) );

  if( std::dynamic_pointer_cast< VariableMod >( mod ) )
   throw( std::invalid_argument( "VariableMod not handled (yet)" ) );

  if( std::dynamic_pointer_cast< BlockMod >( mod ) )
   throw( std::invalid_argument( "BlockMod not handled (yet)" ) );

  if( std::dynamic_pointer_cast< BlockModAD >( mod ) )
   throw( std::invalid_argument( "BlockModAD not handled (yet)" ) );

  // if control reaches here, the Modification is "unknown", probably a
  // "physical" Modification that BundleSolver does not care about

  to_delete = true;

  }  // end( 1st loop, in reverse )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if any global pool has been reset/all of them have reset, then delete all
 // items in the corresponding bundle/reset all the bundle

 if( cntreset == NrFi - NrEasy )  // all (non-easy) components have been reset
  reset_bundle();
 else
  if( cntreset ) {           // at least a (non-easy) component has been reset
   for( Index k = 0 ; k < NrFi ; ++k )
    if( reset[ k ] ) {       // reset[ k ] ==> ! IsEasy[ k ]
     // if BundleSolver "plays nice" with other Solvers, it keeps track of
     // linearizations in the global pool even if they are not in the bundle
     // to avoid overriding them; but there is no longer anything in the
     // global pool, so reset any such information; do this first because
     // MaxItem[ k ] is zeroed during the next loop
     if( ( BPar7 & 3 ) < 3 )
      std::fill( std::next( InvItemVcblr[ k ].begin() , MaxItem[ k ] ) ,
		 InvItemVcblr[ k ].end() , InINF );

     // delete all linearizations in the bundle for this component (in
     // the sense of updating the BundleSolver data structures, since they
     // have been deleted already from the global pool)
     for( Index i = 0 ; i < MaxItem[ k ] ; ++i )
      if( InvItemVcblr[ k ][ i ] < vBPar2.back() )
       Delete( InvItemVcblr[ k ][ i ] , true );
     }
   }

 // After this point, all the Modification adding, deleting or modifying
 // linearizations are significant: they either pertain to components that
 // have never been reset, or are the remaining ones after the (last) one
 // resetting the component

 if( v_mod_tmp.empty() ) {  // no more Modification to process
  // ordinarily, changes in the 0-th component would be dealt with at the
  // end, but if this is the only thing that happened the end is now, so
  // they have to be dealt with immediately
  if( Fi0Chgd )
   Master->ChgSubG( 0 , NumVar , 0 );

  return;                   // all done
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 2nd loop, again in reverse: check for "soft" reset of components, i.e.,
 // when all existing linearization changes. any Modification that changes
 // the linearizations happening before a "soft" reset of the global pool
 // (meaning it is found afterwards in the reverse order) is deleted since it
 // is useless.
 //
 // note that the linearization error of a linearization depends on both the
 // initial constant (\alpha), the linearization itself (g) and the current
 // stability centre (Lambda); thus, if any of those changes, the
 // linearization error need be recomputed. hence reset[ k ] == true
 // implies AlphaC[ k ] == true; the reverse implication does not hold, i.e.,
 // AlphaC[ k ] == true with reset[ k ] == false is possible and it means
 // that only the constants need be changed, but not all the rest. in this
 // loop we consider changes of \alpha and g, while in the 4th loop we will
 // consider changes of Lambda due to the removal of variables; additions
 // never create problems since new variables are always initialized to 0, and
 // therefore they never change the existing linearization error. in fact, not
 // all changes of g necessarily change the linearization error: if a
 // component g_i changes such that Lambda_i == 0, this has no impact. however,
 // in this reverse loop the map between the Lambda[] vector and the indices
 // in the Modification is nontrivial (if additions/removals happened), which
 // makes checking this too complicated. anyway, the issue will go away in the
 // version of BundleSolver that does not use MPSolver since there the
 // linearizations will (likely) be represented by means of their "naked"
 // constant \alpha rather than by their linearization error
 //
 // note that Modification changing the linearizations happening *after* a
 // "soft" reset of the global pool (meaning it is found *before* in the
 // reverse order) is also useless, since a reset forces the re-reading of all
 // linearizations, which by definition happens at their current (final)
 // state. yet, this is not done immediately
 //
 // note that we make no serious attempt at keeping track of the combined
 // effect of all changes, in order to detect if a large set of small
 // changes actually implies a reset. this is complicated for "horizontal"
 // changes (for all linearizations, a range/subset of entries) because the
 // names of the changed Variable may not be current (additions/deletions may
 // happen in the meantime), and keeping track is too burdensome. similarly
 // for "vertical" changes (a set of specific linearizations). some steps
 // in this direction will perhaps be done in later stages of development

 if( cntreset ) {                // if there was any hard reset
  reset.assign( NrFi , false );  // reset reset (couldn't resist)
  cntreset = 0;
  }

 std::vector< bool > AlphaC( NrFi , false );

 for( auto rimod = v_mod_tmp.rbegin() ; rimod != v_mod_tmp.rend() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to rimod
      [ & to_delete , & v_mod_tmp ]( decltype( rimod ) & ri ) {
       if( to_delete )
	ri = std::reverse_iterator( v_mod_tmp.erase( std::next( ri ).base()
						     ) );
       else
	++ri;
       }( rimod ) ) {
  to_delete = false;
  auto mod = *rimod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // a C05FunctionModRngd only changes a range of the linearizations,
  // and therefore is not considered a "soft" reset even if which().empty()
  // in fact the range could be so large as to be (almost) all the
  // variables, which would count as a reset, but so far we don't attempt
  // at detecting this. the only easy case would be NothingChanged, but
  // any such C05FunctionModRngd has been deleted already. however, if the
  // component is "soft" reset already, it can be deleted
  //
  // AllEntriesChanged and AllLinearizationChanged are equivalent, since
  // even if only g changes, also the linearization error does (unless it
  // only changes in places where Lambda == 0, but this cannot be checked
  // efficiently)
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModRngd >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   switch( tmod->type() ) {
    case( C05FunctionMod::AllEntriesChanged ):
    case( C05FunctionMod::AllLinearizationChanged ):
     if( tmod->which().empty() ) {  // reset of the constants
      if( reset[ wFi ] )            // component reset already
       to_delete = true;            // nothing else to do
      else                          // component not reset
       AlphaC[ wFi ] = true;        // reset the constants
      }
     continue;
    default:  // this must not happen
     throw( std::invalid_argument( "wrong type() in C05FunctionModRngd" ) );
    }
   }  // end( if( tmod == C05FunctionModRngd ) )

  // a C05FunctionModSbst only changes a subset of the linearizations,
  // and therefore is not considered a "soft" reset even if which().empty()
  // in fact the subset could be so large as to be (almost) all the
  // variables, which would count as a reset, but so far we don't attempt
  // at detecting this. the only easy case would be NothingChanged, but
  // any such C05FunctionModSbst has been deleted already. however, if the
  // component is "soft" reset already, it can be deleted
  //
  // AllEntriesChanged and AllLinearizationChanged are equivalent, since
  // even if only g changes, also the linearization error does (unless it
  // only changes in places where Lambda == 0, but this cannot be checked
  // efficiently)
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModSbst >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   switch( tmod->type() ) {
    case( C05FunctionMod::AllEntriesChanged ):
    case( C05FunctionMod::AllLinearizationChanged ):
     if( tmod->which().empty() ) {  // reset of the constants
      if( reset[ wFi ] )            // component reset already
       to_delete = true;            // nothing else to do
      else                          // component not reset
       AlphaC[ wFi ] = true;        // reset the constants
      }
     continue;
    default:  // this must not happen
     throw( std::invalid_argument( "wrong type() in C05FunctionModSbst" ) );
    }
   }  // end( if( tmod == C05FunctionModSbst ) )

  // a C05FunctionMod of type AllLinearizationChanged or AllEntriesChanged
  // with which.empty() "soft" resets all the component and Alpha;
  // AlphaChanged only changes the constants (obviously)
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionMod >( mod ) ) {
   if( ! tmod->which().empty() )   // we only react to which().empty()
    continue;

   auto wFi = get_index_of_component( tmod->function() );

   switch( tmod->type() ) {
    case( C05FunctionMod::AlphaChanged ):
     AlphaC[ wFi ] = true;
     break;
    case( C05FunctionMod::AllEntriesChanged ):
    case( C05FunctionMod::AllLinearizationChanged ):
     if( ! reset[ wFi ] )
      ++cntreset;
     reset[ wFi ] = AlphaC[ wFi ] = true;
     break;
    default:  // this must not happen, as GlobalPoolRemoved with
              // which.empty() has been dealt with and deleted before
     throw( std::invalid_argument(
		        "wrong type in C05FunctionMod with empty which()" ) );
    }

   to_delete = true;
   continue;

   }  // end( if( tmod == C05FunctionMod ) )

  // a C05FunctionModLinRngd implies that a specific range in all the
  // linearizations must be changed by adding; this is never considered a
  // "soft" reset even if in fact the range could be so large as to be
  // (almost) all the variables, which would count as a reset, but so far
  // we don't attempt at detecting this. however, if the component is "soft"
  // reset already, it can be deleted
  //
  // again, on the grounds that changing the linearization (presumably)
  // changes the linearization errors as well, this also resets all the
  // Alpha
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModLinRngd >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   if( reset[ wFi ] )            // component reset already
    to_delete = true;            // nothing else to do
   else                          // component not reset
    AlphaC[ wFi ] = true;        // reset the constants

   continue;

   }  // end( if( tmod == C05FunctionModLinRngd ) )

  // a C05FunctionModLinSbst implies that a specific subset in all the
  // linearizations must be changed by adding; this is never considered a
  // "soft" reset even if in fact the subset could be so large as to be
  // (almost) all the variables, which would count as a reset, but so far
  // we don't attempt at detecting this. however, if the component is "soft"
  // reset already, it can be deleted
  //
  // again, on the grounds that changing the linearization (presumably)
  // changes the linearization errors as well, this also resets all the
  // Alpha
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModLinSbst >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   if( reset[ wFi ] )            // component reset already
    to_delete = true;            // nothing else to do
   else                          // component not reset
    AlphaC[ wFi ] = true;        // reset the constants

   continue;

   }  // end( if( tmod == C05FunctionModLinSbst ) )

  // a C05FunctionModLin implies that *all* the linearizations must be
  // changed by adding them \delta; this may in principle be handled in
  // a specialised way by BundleSolver, but is currently not, and
  // therefore it is a "full" reset
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModLin >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );
   if( ! reset[ wFi ] )
    ++cntreset;
   AlphaC[ wFi ] = reset[ wFi ] = to_delete = true;

   }  // end( if( tmod == C05FunctionModLin ) )
  }  // end( 2nd loop, again in reverse )

 // note that even if there were no more Modification to process we could not
 // stop because this means that reset[ k ] == true and/or AlphaC[ k ] == true
 // for some k. in fact v_mod_tmp was not empty(), and elements can be removed
 // from it only if some component is "soft" reset. 

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 3rd loop, forward: prepare for addition/removal/changes of individual
 // linearization for each component by computing the four sets
 // - linearizations that need be removed
 // - linearizations that need be added
 // - linearizations that need be changed
 // - constants that need be changed
 // Note that:
 // - due to limitations in the MPSolver interface, changing a linearization
 //   implies changing its constant; therefore, Cchg[] contains changes in
 //   the constants only and nothing else, which means that has empty
 //   intersection with all other three sets
 // - if a linearization that is added/changed is later removed, it is no
 //   longer added/changed
 // - if a linearization that is removed/changed is later added it is no
 //   longer removed/changed (adding over an existing linearization changes
 //   it anyway, no reason to remove it)
 // - thus, Addd[], Rmvd[], Chgd[] and Cchg[] all have empty intersection
 // - linearization changes to a reset component can be ignored
 // - constant changes when all constants change can be ignored
 // - due to limitations in the MPSolver interface, "horizontal" changes (of
 //   a given range/subset of entries) to a subset of linearizations are
 //   not supported. these must be either mapped in "horizontal" changes to
 //   *all* linearizations (of a given component), or to "vertical" changes
 //   of all components to a subset of linearization. somewhat arbitrarily,
 //   the second option is chosen here. as a consequence, C05FunctionModRngd
 //   and C05FunctionModSbst are considered C05FunctionMod. note that
 //   C05FunctionMod* with which().empty() still remain untreated, as well
 //   as C05FunctionModLinRngd and C05FunctionModLinSbst (that by definition
 //   always concern all the linearizations), while C05FunctionModLin have
 //   been dealt with already

 std::vector< Subset > Addd( NrFi );  // added linearizations
 std::vector< Subset > Rmvd( NrFi );  // removed linearizations
 std::vector< Subset > Chgd( NrFi );  // changed linearizations
 std::vector< Subset > Cchg( NrFi );  // changed constants

 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to imod
      [ & to_delete , & v_mod_tmp ]( decltype( imod ) & it ) {
       if( to_delete )
	it =  v_mod_tmp.erase( it );
       else
	++it;
       }( imod ) ) {
  to_delete = false;
  auto mod = *imod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly
  //
  // actually, only C05FunctionMod need be treated here. we do not need to
  // distinguish between the base and the derived classes because we take
  // the change of a range/subset of the entries for a change of the whole
  // linearization. however, we do in fact distinguish them because we
  // ignore the Modification if ttmod->which().empty() and type() ==
  // AllEntriesChanged or == AllLinearizationChanged. any such C05FunctionMod
  // has been deleted (at worst) in the 2nd loop, and therefore here it must
  // be a C05FunctionModRngd or C05FunctionModSbst
  //
  // also, note that NothingChanged and AlphaChanged cannot happen, since
  // again these cases have been dealt with already

  if( const auto tmod = std::dynamic_pointer_cast< C05FunctionMod >( mod ) ) {
   auto wFi = get_index_of_component( tmod->function() );

   switch( tmod->type() ) {
    case( C05FunctionMod::AlphaChanged ):
     assert( ! tmod->which().empty() );
     // this cannot happen, the case has been dealt with and deleted
     if( AlphaC[ wFi ] )  // constants are reset already
      break;              // nothing else to do

     // note that reset[ wFi ] ==> AlphaC[ wFi ], hence surely reset[ wFi ]
     // is false here
     // add to Cchg[ wFi ] the names in tmod->which(), save those that are
     // in either Addd[ wFi ], or Rmvd[ wFi ], or Chgd[ wFi ]
     if( Addd[ wFi ].empty() && Rmvd[ wFi ].empty() && Chgd[ wFi ].empty() )
      set_union_in_place( Cchg[ wFi ] , tmod->which() );
     else {
      // only those constants corresponding to linearizations that are not
      // added or removed or changed need be changed
      Subset tmp( tmod->which() );

      if( ! Addd[ wFi ].empty() )
       set_difference_in_place( tmp , Addd[ wFi ] );

      if( ( ! Rmvd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Rmvd[ wFi ] );

      if( ( ! Chgd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Chgd[ wFi ] );

      set_union_in_place( Cchg[ wFi ] , std::move( tmp ) );
      }
     break;
    case( C05FunctionMod::AllLinearizationChanged ):
    case( C05FunctionMod::AllEntriesChanged ):
     // AllEntriesChanged is completely equivalent to AllLinearizationChanged
     // since the change in the linearization implies that of the constant
     // if tmod->which().empty(), this must actually be either a
     // C05FunctionModRngd or a C05FunctionModSbst: save it, for it
     // will be dealt with in the next loop
     if( tmod->which().empty() )
      continue;

     if( reset[ wFi ] )  // changes in reset components
      break;             // are ignored

     // add to Chgd[ wFi ] the names in tmod->which(), save those that are
     // in either Addd[ wFi ] or Rmvd[ wFi ]
     if( Addd[ wFi ].empty() && Rmvd[ wFi ].empty() )
      set_union_in_place( Chgd[ wFi ] , tmod->which() );
     else {
      // only those items that are not added and/or removed need be changed
      Subset tmp( tmod->which() );

      if( ! Addd[ wFi ].empty() )
       set_difference_in_place( tmp , Addd[ wFi ] );

      if( ( ! Rmvd[ wFi ].empty() ) && ( ! tmp.empty() ) )
       set_difference_in_place( tmp , Rmvd[ wFi ] );

      set_union_in_place( Chgd[ wFi ] , std::move( tmp ) );
      }

     // changed linearizations imply changed constants
     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , Chgd[ wFi ] );

     break;

    case( C05FunctionMod::GlobalPoolAdded ):
     // add to Addd[ wFi ] the names in tmod->which(), and remove them
     // from Rmvd[ wFi ], and from Chgd[ wFi ] if the component is not
     // reset, and from Cchg[ wFi ] is constants are not reset

     set_union_in_place( Addd[ wFi ] , tmod->which() );

     set_difference_in_place( Rmvd[ wFi ] , tmod->which() );

     if( ! reset[ wFi ] )
      set_difference_in_place( Chgd[ wFi ] , tmod->which() );

     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , tmod->which() );

     break;

    case( C05FunctionMod::GlobalPoolRemoved ):
     // add to Rmvd[ wFi ] the names in tmod->which(), and remove them
     // from Addd[ wFi ], and from Chgd[ wFi ] if the component is not
     // reset, and from Cchg[ wFi ] is constants are not reset

     set_union_in_place( Rmvd[ wFi ] , tmod->which() );

     set_difference_in_place( Addd[ wFi ] , tmod->which() );

     if( ! reset[ wFi ] )
      set_difference_in_place( Chgd[ wFi ] , tmod->which() );

     if( ! AlphaC[ wFi ] )
      set_difference_in_place( Cchg[ wFi ] , tmod->which() );

    }  // end( switch( tmod->type() ) )

   to_delete = true;   // delete it

   }  // end( if( tmod == C05FunctionMod ) )
  }  // end( 3rd loop, forward )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // now delete all linearization that need to, if any

 if( std::find_if( Rmvd.begin() , Rmvd.end() ,
		   []( Subset & Rk ) { return( ! Rk.empty() ); }
		   ) != Rmvd.end() ) {
  // at least a component has had linearizations removed, but note that not
  // all linearizations need be in the bundle; if they are not they are
  // just removed from the global pool (if they are there)

  for( Index k = 0 ; k < NrFi ; ++k )
   for( auto i : Rmvd[ k ] ) {
    auto h = InvItemVcblr[ k ][ i ];
    if( h < vBPar2.back() )
     Delete( h , true );
    else
     remove_from_global_pool( k , i , true );
    }
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 4th loop: handle exclusively Variable addition and removal, in forward
 // order. all removals of existing Variable are immediately performed,
 // diminishing NumVar. additions of Variable are "cached", and only performed
 // once (if ever) after the end of the cycle. removals of Variable that have
 // not been added yet just decreases the number of new Variable to be added
 // at the end. all corresponding Modification are removed from the list
 //
 // this loop may also force some linearization errors to be reset since they
 // depend not only on the initial constant (\alpha), but also from the
 // linearization itself (g) and the current stability centre (Lambda); if
 // any of those changes, the linearization error need be recomputed. in this
 // loop we consider changes of Lambda due to the removal of variables, while
 // additions never create this problem since new variables are always
 // initialized to 0, and therefore they never change the existing
 // linearization error per-se (recall that changes in \alpha and g were
 // already handled by the 2nd loop). note, however, that not all changes of
 // Lambda necessarily change the linearization error: if a Lambda_i == 0 is
 // removed, this has no impact. because this loop (unlike the 2nd one) is
 // forward and does deletions immediately, Lambda is always "in synch" with
 // the indices of the Modification and therefore the check can be easily
 // done. this is all the more important since variables are removed in
 // parallel for all components, hence if a variable with Lambda_i != 0 is
 // removed then all the linearization errors must be reset
 Index to_add = 0;
 bool addd_vars = false;  // if any Variable has ever been added
 bool rmvd_vars = false;  // if any Variable has ever been removed

 for( auto imod = v_mod_tmp.begin() ; imod != v_mod_tmp.end() ;
      // note the iterator_expression of the for() obtained by defining
      // a lambda and then immediately applying it to imod
      [ & to_delete , & v_mod_tmp ]( decltype( imod ) & it ) {
       if( to_delete )
	it =  v_mod_tmp.erase( it );
       else
	++it;
       }( imod ) ) {
  to_delete = false;
  auto mod = *imod;

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // not that we do not distinguish C05FunctionModVars* from "plain"
  // FunctionModVars*, since the only difference is whether or not the
  // operation is strongly quasi-additive, i.e., it implies or not a "hard"
  // reset, but this has already been acted upon

  {
   // a "naked" FunctionModVars
   auto tmod = std::dynamic_pointer_cast< FunctionModVars >( mod );
   if( ! tmod ) {
    // if it is not a "naked" FunctionModVars, it can still be a group of
    // identical *FunctionModVars* "dressed" into a GroupModification
    if( const auto gmod =
	std::dynamic_pointer_cast< GroupModification >( mod ) )
     // if so, pick the first one and act on it
     tmod = std::static_pointer_cast< FunctionModVars >(
					 gmod->sub_Modifications().front() );
    }

   if( tmod ) {
    // if we have a *FunctionModVars*, we have to distinguish its exact type
    // and add/delete Variable accordingly; in all cases, however, the
    // Modification is processed and can be deleted
    to_delete = true;

    if( const auto ttmod =
	std::dynamic_pointer_cast< FunctionModVarsAddd >( tmod ) ) {
     addd_vars = true;
     if( ! to_add ) {
      // the first time, check that the Modification data agrees with what
      // we expect
      if( ttmod->first() != NumVar )
       throw( std::logic_error( "wrong Variable names in FunctionModVars" ) );
      }

     to_add += ttmod->vars().size();
     continue;

     } // end( if( tmod == FunctionModVarsAddd ) )

    if( const auto ttmod =
	std::dynamic_pointer_cast< FunctionModVarsRngd >( tmod ) ) {
     rmvd_vars = true;
     Range rng = ttmod->range();

     // if any of the deleted Variable is nonzero, the linearization errors
     // will have to be recomputed for all components
     for( Index i = std::min( rng.first , NumVar ) ;
	  i < std::min( rng.second , NumVar ) ; ++i )
      if( std::abs( Lambda[ i ] ) > 1e-12 ) {
       std::fill( AlphaC.begin() , AlphaC.end() , true );
       break;
       }

     if( ( rng.first == 0 ) && ( rng.second >= NumVar ) ) {
      NumVar = 0;                      // deleting *all* Variable
      Lambda.clear();
      Lambda1.clear();
      LmbdBst.clear();
      Master->RmvVars( nullptr , 0 );  // remove from MP
      continue;                        // nothing else to do
      }

     if( rng.first >= NumVar ) {  // all the Variable are deleted already
      auto nr = rng.second - rng.first;
      if( nr > to_add )
       throw( std::logic_error( "removing non-existing Variable" ) );
      to_add -= nr;               // "virtually" remove them
      continue;                   // nothing else to do
      }
     if( rng.second >= NumVar ) {  // some of the Variable deleted already
      auto nr = rng.second - NumVar;
      if( nr > to_add )
       throw( std::logic_error( "removing non-existing Variable" ) );
      to_add -= nr;               // "virtually" remove them
      rng.second = NumVar;
      }
     if( rng.second < NumVar ) {
      // if deleting the last range of Variable nothing has to be done,
      // but deleting Variable "in the middle" rather requires moving
      // down the remaining range of values in Lambda
      std::copy( Lambda.begin() + rng.second , Lambda.end() ,
		 Lambda.begin() + rng.first );
      }
     Subset tdlt( rng.second - rng.first );
     NumVar -= tdlt.size();
     Lambda.resize( NumVar );  // adjust Lambda
     Lambda1.resize( NumVar );
     if( MaxSol > 1 )
      LmbdBst.resize( NumVar );
     std::iota( tdlt.begin() , tdlt.end() , rng.first );
     Master->RmvVars( tdlt.data() , tdlt.size() );  // remove from MP
     continue;

     }  // end( if( ttmod == FunctionModVarsRngd ) )

    if( const auto ttmod =
	std::dynamic_pointer_cast< FunctionModVarsSbst >( tmod ) ) {
     rmvd_vars = true;

     if( ttmod->subset().empty() ) {  // deleting *all* Variable
      // if any Variable is nonzero, the linearization errors
      // will have to be recomputed for all components
      for( auto & el : Lambda )
       if( std::abs( el ) > 1e-12 ) {
	std::fill( AlphaC.begin() , AlphaC.end() , true );
	break;
        }

      NumVar = 0;
      Lambda.clear();
      Lambda1.clear();
      LmbdBst.clear();
      Master->RmvVars( nullptr , 0 );  // remove from MP
      continue;                        // nothing else to do
      }

     if( ttmod->subset().front() >= NumVar ) {
      // all the Variable are deleted already
      if( ttmod->subset().size() > to_add )
       throw( std::logic_error( "removing non-existing Variable" ) );
      to_add -= ttmod->subset().size();  // "virtually" remove them
      continue;                          // nothing else to do
      }

     // if any of the deleted Variable is nonzero, the linearization errors
     // will have to be recomputed for all components
     for( auto el : ttmod->subset() ) {
      if( el >= NumVar )
       break;

      if( std::abs( Lambda[ el ] ) > 1e-12 ) {
       std::fill( AlphaC.begin() , AlphaC.end() , true );
       break;
       }
      }

     Subset tsbst;
     c_Subset * sbst = & tsbst;
     if( ttmod->subset().back() < NumVar )  // no Variable deleted already
      sbst = & ttmod->subset();             // delete them all
     else {                                 // construct the subset to delete
      auto sbstit = ttmod->subset().end();
      while( *(--sbstit) >= NumVar );
      tsbst = Subset( ttmod->subset().begin() , ++sbstit );
      auto nr = ttmod->subset().size() - tsbst.size();
      if( nr > to_add )
       throw( std::logic_error( "removing non-existing Variable" ) );
      to_add -= nr;               // "virtually" remove them
      }

     Compact( Lambda , *sbst );  // adjust Lambda
     NumVar -= sbst->size();
     Lambda.resize( NumVar );
     Lambda1.resize( NumVar );
     if( MaxSol > 1 )
      LmbdBst.resize( NumVar );
     Master->RmvVars( sbst->data() , sbst->size() );  // remove from MP
     continue;

     }  // end( if( ttmod == FunctionModVarsSbst ) )

    // if control reaches here, this is an unknown *FunctionModVars* (??)
    throw( std::logic_error( "unknown FunctionModVars" ) );

    }  // end( if( tmod == FunctionModVars ) )
   }  // end FunctionModVars
  }  // end( 4th loop, forward )

 // at this point, the set of Variable in the BundleSolver/Master Problem
 // coincides with the set of Variable in the C05Function(s), save for the
 // Variable to be added: in other words, the positions from 0 no NumVar - 1
 // in the linearizations corresponds to what BundleSolver expects

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // 5th loop: handle "horizontal" changes, i.e., changes of a given range
 // (subset) of entries in all the linearizations, i.e., C05FunctionMod* and
 // C05FunctionModLin* with which().empty(). note that due to limitations
 // of the MPSolver interface, subsets are anyway translated to a range,
 // thereby possibly requiring also entries that have not changed. if
 // Variable have been removed the original indices in the Modification need
 // be translated, and in fact if Variable in the range have been removed
 // the range can be shrank, up to disappearing altogether. the mapping is
 // more difficult if Variable have also been added. however we can exploit
 // the property that if a Variable ever has a name that is larger than its
 // index in the Modification, this can only mean that the Variable has been
 // removed and then re-added after that the Modification has been issued.
 // these Variable can therefore be ignored, since (unless they have been
 // re-removed) they will be added in the end
 //
 // this is the final loop, so the list must be empty at the end

 for( ; ! v_mod_tmp.empty() ; v_mod_tmp.pop_front() ) {
  auto mod = v_mod_tmp.front();  // pick the first Modification

  Index wFi;                     // the affected component
  Range range( NumVar , 0 );     // an empty range
  c_Subset * subset = nullptr;   // an empty subset
  c_Vec_p_Var * vars;            // the affected Variable

  // patiently sift through the possible Modification types to find what mod
  // exactly is and react accordingly

  // a C05FunctionModRngd, that at this point can only have which().empty()
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModRngd >( mod ) ) {
   if( ! tmod->which().empty() )
    throw( std::logic_error( "unexpected nonempty C05FunctionModRngd" ) );

   wFi = get_index_of_component( tmod->function() );
   vars = & tmod->vars();
   range = tmod->range();
   goto done;
   }

  // a C05FunctionModSbst, that at this point can only have which().empty()
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModSbst >( mod ) ) {
   if( ! tmod->which().empty() )
    throw( std::logic_error( "unexpected nonempty C05FunctionModSbst" ) );

   wFi = get_index_of_component( tmod->function() );
   vars = & tmod->vars();
   subset = & tmod->subset();
   goto done;
   }

  // a C05FunctionModLinRngd implies that a specific range in all the
  // linearizations must be changed (by adding something)
  if( const auto tmod =
      std::dynamic_pointer_cast< C05FunctionModLinRngd >( mod ) ) {
   wFi = get_index_of_component( tmod->function() );
   vars = & tmod->vars();
   range = tmod->range();
   goto done;
   }

  // a C05FunctionModLinSbst implies that a specific subset in all the
  // linearizations must be changed (by adding something)
  if( const auto tmod =
	   std::dynamic_pointer_cast< C05FunctionModLinSbst >( mod ) ) {
   wFi = get_index_of_component( tmod->function() );
   vars = & tmod->vars();
   subset = & tmod->subset();
   goto done;
   }

  // it is neither of the above: this should not happen
  throw( std::logic_error( "unexpected Modification slipped in" ) );

  // the range/subset (and component) have been identified: check if the
  // need to be translated due to addition/removals, and in case do it
  done:if( ! rmvd_vars ) {
   // Variable have never been removed, hence the names can be used directly
   if( subset ) {  // turn the subset into a range
    range.first = subset->front();
    range.second = subset->back() + 1;
    }
   }
  else {
   // Variable have been removed, and hence names need be actualised
   // this is done by directly checking vars() against the "active"
   // Variable of v_c05f[ 0 ], which is fairly taken as a representative
   // since all the C05Function have the same "active" Variable
   if( ! addd_vars ) {
    // ... but never added: names can have only decreased, but even more
    // importantly must have remained ordered, i.e., the first "active"
    // Variable in vars() is the first variable of the range, the last
    // "active" Variable vars() is the last variable of the range
    // note that we do not use subset and range here, as the range is
    // reconstructed from scratch using vars
    auto lit = vars->begin();
    for( ; lit != vars->end() ; ++lit ) {
     range.first = v_c05f[ 0 ]->is_active( *lit );
     if( range.first < v_c05f[ 0 ]->get_num_active_var() )
      break;
     }
    if( lit == vars->end() )  // no Variable in vars is still "active"
     continue;                // nothing else to do
    // since we know that here are some "active" Variable in vars(), this
    // second loop will necessarily end
     for( auto rit = vars->rbegin() ; ; ++rit ) {
      range.second = v_c05f[ 0 ]->is_active( *lit );
      if( range.second < v_c05f[ 0 ]->get_num_active_var() )
       break;
      }
     ++range.second;  // the range is [ first , second )
     }
   else {
    // the complicated case: Variable have both been removed and added
    // names can have changed in an almost arbitrary way, except that if
    // a name has increased then the Variable has been deleted and re-added
    // and therefore need not be included
    Subset newnames( vars->size() );
    auto lit = vars->begin();
    auto nni = newnames.begin();
    if( subset ) {
     auto sit = subset->begin();
     for( ; lit != vars->end() ; ++lit , ++sit ) {
      auto i = v_c05f[ 0 ]->is_active( *lit );
      if( ( i <= *sit ) && ( i < v_c05f[ 0 ]->get_num_active_var() ) )
       *(nni++) = i;
      }
     }
    else {
     for( ; lit != vars->end() ; ++lit , ++range.first ) {
      auto i = v_c05f[ 0 ]->is_active( *lit );
      if( ( i <= range.first ) && ( i < v_c05f[ 0 ]->get_num_active_var() ) )
       *(nni++) = i;
      }
     }
    if( nni == newnames.begin() )  // no Variable in vars is still "active"
     continue;                     // nothing else to do
    newnames.resize( std::distance( newnames.begin() , nni ) );
    range.first = newnames.front();
    range.second = newnames.back();
    ++range.second;  // the range is [ first , second )
    }
   }

  // now actually do it
  Master->ChgSubG( range.first , range.second , wFi + 1 );

  }  // end( 5th loop, forward )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are Variable to add, do it now in one blow
 // there is a trade-off here: doing this now causes Master->AddVars() to
 // (indirectly) call get_linearization_coefficients() on a smaller set of
 // linearizations, if additions are done, but on the other hand increases
 // NumVar and therefore the work done in later stages. hence, this is done
 // here only if no additions are done

 bool toadd = std::find_if( Addd.begin() , Addd.end() ,
			    []( Subset & Ak ) { return( ! Ak.empty() ); }
			    ) != Addd.end();
 if( to_add && ( ! toadd ) ) {
  NumVar += to_add;
  Lambda.resize( NumVar , 0 );
  Lambda1.resize( NumVar , 0 );
  if( MaxSol > 1 )
   LmbdBst.resize( NumVar , 0 );
  Master->AddVars( to_add );
  to_add = 0;  // done already
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are linearization to add/change, do it now in one blow
 //
 // note that due to limitations in the MPSolver interface, changing a
 // linearization is identical to adding one, even if the change was limited
 // to a range/subset of the entries
 //
 // yet, handling of additions and changes differs depending on BPar7.
 // in particular, if ( BPar7 & 4 ), then additions to the global pool
 // imply additions to the bundle, otherwise just marking that the
 // linearization exists
 //
 // note that if ( BPar7 & 3 ) >= 2, BundleSolver will happily delete from
 // the global pool any linearization it deletes from the bundle; if also
 // ! ( BPar7 & 4 ), one could therefore think it appropriate to delete
 // from the global pool any linearization that is added. however, we do not
 // do that, just refraining to add them to the bundle. these will likely be
 // overwritten during the optimization, and if memory is a problem then the
 // global pools should just be sized accordingly

 if( toadd ||
     ( std::find_if( Chgd.begin() , Chgd.end() ,
		     []( Subset & Ck ) { return( ! Ck.empty() ); }
		     ) != Chgd.end() ) ) {
  // at least a component has had linearizations added or changed

  for( Index k = 0 ; k < NrFi ; ++k ) {
   if( Addd[ k ].empty() && Chgd[ k ].empty() )  // nothing to see here
    continue;                                    // move on

   if( Addd[ k ].empty() && ( Chgd[ k ].size() >= NrItems[ k ] ) ) {
    // all existing items change and no new one is added
    if( ! reset[ k ] )
     ++cntreset;
    reset[ k ] = AlphaC[ k ] = true;             // this is a reset
    continue;
    }

   // first, cleanup Chgd[ k ] from linearizations not in the bundle
   if( ! Chgd[ k ].empty() ) {
    Subset tmp;
    for( auto i : Chgd[ k ] )
     if( InvItemVcblr[ k ][ i ] >= vBPar2.back() )
      tmp.push_back( i );

    if( ! tmp.empty() ) {
     if( tmp.size() >= Chgd[ k ].size() ) {
      Chgd[ k ].clear();
      if( Addd[ k ].empty() )  // nothing more to see here
       continue;               // move off
      }
     else
      set_difference_in_place( Chgd[ k ] , tmp );
     }
    }

   // now, if ! ( BPar7 & 4 ), also cleanup Addd[ k ] from linearizations
   // not in the bundle, but in doing so mark them into InvItemVcblr unless
   // ( BPar7 & 3 ) == 3, in which case BundleSolver does not care if there
   // are existing linearizations because it anyway freely overwrites them
   // note that for linearizations in the bundle, being in Addd[ k ] is the
   // same as being in Chgd[ k ]: the linearization has changed, and the
   // master problem must be changed to reflect this
   if( ( ! ( BPar7 & 4 ) ) && ( ! Addd[ k ].empty() ) ) {
    Subset tmp;
    for( auto i : Addd[ k ] )
     if( InvItemVcblr[ k ][ i ] >= vBPar2.back() ) {
      InvItemVcblr[ k ][ i ] = ( BPar7 & 3 ) == 3 ? InINF : vBPar2.back();
      tmp.push_back( i );
      }

    if( ! tmp.empty() ) {
     if( tmp.size() >= Addd[ k ].size() ) {
      Addd[ k ].clear();
      if( Chgd[ k ].empty() )  // nothing more to see here
       continue;               // move off
      }
     else
      set_difference_in_place( Addd[ k ] , tmp );
     }
    }

   // compute the union between Addd[ k ] and Chgd[ k ] into Addd[ k ]
   set_union_in_place( Addd[ k ] , std::move( Chgd[ k ] ) );

   // finally, add the resulting stuff to the bundle
   for( auto i : Addd[ k ] )
    add_to_bundle( k , i );
   }
  }  // end( if( additions or changes ) )

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if some component need be reset, reset the linearizations: since
 // reset[ k ] ==> AlphaC[ k ], later on also the constants will be reset

 if( cntreset == NrFi - NrEasy )  // all (non-easy) components have been reset
  // ... and if Fi0Chgd == true, then also the 0-th component has changed
  Master->ChgSubG( 0 , NumVar , Fi0Chgd ? InINF : NrFi + 1 );
 else {
  if( cntreset )                 // some (non-easy) components have been reset
   for( Index k = 0 ; k < NrFi ; ++k )
    if( reset[ k ] )             // reset[ k ] ==> ! IsEasy[ k ]
     Master->ChgSubG( 0 , NumVar , k + 1 );

  if( Fi0Chgd )                // ... and/or the 0-th component has changed
   Master->ChgSubG( 0 , NumVar , 0 );
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are constants to change entirely, do it now in one blow
 // note: in case of a full reset, get_linearization_coefficients() is
 //       called twice, once in ChgSubG() (via GetGi()) and once in the
 //       loop below. this has the potential to be horribly inefficient,
 //       but the only clean way out is to do away with MPSolver entirely

 if( std::find( AlphaC.begin() , AlphaC.end() , false ) == AlphaC.end() ) {
  for( auto & cchk : Cchg )  // all components have been reset
   cchk.clear();             // no need to change them individually

  ResetAlfa( NrFi );
  }
 else
  for( Index k = 0 ; k < NrFi ; ++k )
   if( AlphaC[ k ] ) {  // all the constants of this component are reset
    Cchg[ k ].clear();  // no need to change them individually

    ResetAlfa( k );
    }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there subsets of Alphas to change, do it now

 if( std::find_if( Cchg.begin() , Cchg.end() ,
		   []( Subset & Ck ) { return( ! Ck.empty() ); }
		   ) != Cchg.end() ) {
  std::vector< VarValue > Gi( NumVar );

  for( Index k = 0 ; k < NrFi ; ++k )
   for( auto i : Cchg[ k ] )
    if( InvItemVcblr[ k ][ i ] < vBPar2.back() ) {
     auto Ai = rs( v_c05f[ k ]->get_linearization_constant( i ) );

     #ifndef NDEBUG
      if( std::isnan( Ai ) )  // linearization no longer valid
       throw( std::logic_error( "inexistent linearization" ) );
     #endif

     // compute the linearization error in Lambda
     v_c05f[ k ]->get_linearization_coefficients( Gi.data() ,
						  Range( 0 , NumVar ) , i );
     if( ! f_convex )
      chgsign( Gi.data() , NumVar );

     Ai = UpRifFi[ k ] - Ai - std::inner_product( Lambda.begin() ,
						  Lambda.end() ,
						  Gi.begin() , double( 0 ) );
     Master->ChgAlfa( InvItemVcblr[ k ][ i ] , Ai );
     }
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // if there are (still) Variable to add, do it now in one blow

 if( to_add ) {
  NumVar += to_add;
  Lambda.resize( NumVar , 0 );
  Lambda1.resize( NumVar , 0 );
  if( MaxSol > 1 )
   LmbdBst.resize( NumVar , 0 );
  Master->AddVars( to_add );
  }

 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // and this, finally, is all! (save possibly some checks)

 //!! PrintBundle();
 #if CHECK_DS & 1
  CheckBundle();
 #endif
 #if CHECK_DS & 4
  CheckAlpha();
 #endif

 }  // end( BundleSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

void BundleSolver::CheckBundle( void )
{
 std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                           : f_log;
 // check vBPar2
 for( Index k = 0 ; k < NrFi ; ++k )
  if( Index( v_c05f[ k ]->get_int_par( C05Function::intGPMaxSz ) )
      != vBPar2[ k ] )
   *wlog << "size of global pool " << k << " does not match" << std::endl;

 // check ItemVcblr against InvItemVcblr and Master
 Subset tmp( NrFi , 0 );
 for( Index i = 0 ; i < Master->MaxName() ; ++i )
  if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
   ++tmp[ ItemVcblr[ i ].first ];
   if( Master->WComponent( i ) != ItemVcblr[ i ].first + 1 ) {
    *wlog << "position " << i << " in the bundle should be of component "
	  << ItemVcblr[ i ].first << " but Master says ";
    if( Master->WComponent( i ) == InINF )
     *wlog << "empty" << std::endl;
    else
     *wlog << Master->WComponent( i ) - 1 << std::endl;
    }

   if( InvItemVcblr[ ItemVcblr[ i ].first ][ ItemVcblr[ i ].second ] != i )
    *wlog << "position " << i << " in the bundle shoud be linearization "
	  << ItemVcblr[ i ].second << " of component "
	  << ItemVcblr[ i ].first << " but InvItemVcblr disagrees"
	  << std::endl;
   }
  else
   if( Master->WComponent( i ) < InINF )
    *wlog << "position " << i
	  << " in the bundle should be empty but Master says "
	  << Master->WComponent( i ) << std::endl;

 // check NrItems
 if( int diff = ( NrItems.back() - std::accumulate( NrItems.begin() ,
						    NrItems.begin() + NrFi ,
						    Index( 0 ) ) ) != 0 )
  *wlog << " NrItems[ NrFi ] = " << NrItems.back() << " but the sum is "
	<<  NrItems.back() + diff  << std::endl;

 for( Index k = 0 ; k < NrFi ; ++k )
  if( tmp[ k ] != NrItems[ k ] )
   *wlog << "counted " << tmp[ k ] << " items in the bundle for component "
	 << k << " but NrItems says " << NrItems[ k ] << std::endl;

 // check InvItemVcblr against ItemVcblr and C05Function
 for( Index k = 0 ; k < NrFi ; ++k )
  for( Index i = 0 ; i < vBPar2[ k ] ; ++i ) {
   if( ( InvItemVcblr[ k ][ i ] < InINF ) &&
       ( ! v_c05f[ k ]->is_linearization_there( i ) ) )
    *wlog << "linearization " << i << " in pool " << k
	  << " does not exist" << std::endl;

   if( ( InvItemVcblr[ k ][ i ] == InINF ) &&
       v_c05f[ k ]->is_linearization_there( i ) )
    *wlog << "linearization " << i << " in pool " << k
	  << " unaccounted for" << std::endl;

   if( ( InvItemVcblr[ k ][ i ] < vBPar2.back() ) &&
       ( ( ItemVcblr[ InvItemVcblr[ k ][ i ] ].first != k ) ||
	 ( ItemVcblr[ InvItemVcblr[ k ][ i ] ].second != i ) ) )
    *wlog << "linearization " << i << " in pool " << k
	  << " should be in bundle in position "
	  << InvItemVcblr[ k ][ i ] << " but ItemVcblr disagrees"
	  << std::endl;

   if( ( InvItemVcblr[ k ][ i ] < ( ( BPar7 & 3 ) ? vBPar2.back() : InINF ) )
       && ( i >= MaxItem[ k ] ) )
    *wlog << "item in position " << i << " of pool " << k
	  << " but MaxItem says " << MaxItem[ k ] << std::endl;

   if( ( InvItemVcblr[ k ][ i ] >= ( ( BPar7 & 3 ) ? vBPar2.back() : InINF ) )
       && ( i < FrFItem[ k ] ) )
     *wlog << "free item in position " << i << " of pool " << k
	   << " but FrFItem says " << FrFItem[ k ] << std::endl;

   }  // end( for( i ) )

 // check FreList (if there is anything to check)
 if( ! FreList.empty() ) {
  // copy FreList to a vector to check it (this only works since the
  // underlying container is a std::vector< Index >) and the topmost
  // element in a heap is the first element of the vector
  tmp.resize( FreList.size() );
  std::copy( &( FreList.top() ) , &( FreList.top() ) + FreList.size() ,
	     tmp.begin() );
  std::sort( tmp.begin() , tmp.end() );

  for( auto i : tmp )
   if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] )
    *wlog << "item " << i << " in FreList is not free" << std::endl;

  for( Index i = 0 ; i < Master->MaxName() ; ++i )
   if( ItemVcblr[ i ].second >= vBPar2[ ItemVcblr[ i ].first ] ) {
    auto it = std::lower_bound( tmp.begin() , tmp.end() , i );
    if( ( it == tmp.end() ) || ( *it != i ) )
     *wlog << "free item " << i << " not in FreList" << std::endl;
    }
  }
 }  // end( BundleSolver::CheckBundle )

/*--------------------------------------------------------------------------*/

void BundleSolver::CheckAlpha( void )
{
 std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                           : f_log;
 *wlog << def;
 cHpRow tA = Master->ReadLinErr();
 std::vector< VarValue > G( NumVar );
 const double eps = 1e-8;

 for( Index i = 0 ; i < Master->MaxName() ; ++i )
  if( ItemVcblr[ i ].second < vBPar2[ ItemVcblr[ i ].first ] ) {
   v_c05f[ ItemVcblr[ i ].first ]->get_linearization_coefficients( G.data() ,
						        Range( 0 , NumVar ) ,
						     ItemVcblr[ i ].second );
   if( ! f_convex )
    chgsign( G.data() , NumVar );
   HpNum tAi = rs( v_c05f[ ItemVcblr[ i ].first
			   ]->get_linearization_constant(
						   ItemVcblr[ i ].second ) );
   tAi = UpRifFi[ ItemVcblr[ i ].first ] - tAi -
                         std::inner_product( Lambda.begin() , Lambda.end() ,
					     G.begin() , double( 0 ) );

   if( std::abs( tAi - tA[ i ] ) >= eps *
       std::max( std::max( std::abs( tAi ) ,
			   std::abs( UpRifFi[ ItemVcblr[ i ].first ] ) ) ,
		 double( 1 ) ) )
    *wlog << std::endl << "Alfa[ " << i << " ]: F = " << tAi << " ~ M = "
	  << tA[ i ];
    }

 }  // end( BundleSolver::CheckAlpha )

/*--------------------------------------------------------------------------*/

void BundleSolver::CheckLBs( void )
{
 std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                           : f_log;
 *wlog << def;
 const double eps = 1e-8;
 const double one = 1;

 auto LB = Master->ReadLowerBound();
 auto GLB = f_convex ?   f_Block->get_valid_lower_bound( false )
                     : - f_Block->get_valid_upper_bound( false );

 if( TrueLB ) {  // a finite global lower bound is set
  if( LowerBound.back() == -INFshift ) {
   *wlog << std::endl << "TrueLB but no stored global bound";
   if( GLB > -INFshift )
   *wlog << std::endl << "finite global bound " << rs( GLB )
	  << " available but not set";
   }
  else
   if( GLB == -INFshift )
    *wlog << std::endl << "global bound = -INF but bound "
	  << rs( LowerBound.back() ) << " sett";
   else {
    if( std::abs( GLB - LowerBound.back() ) > eps * std::max( GLB , one ) )
     *wlog << std::endl << "global bound = " << rs( GLB )
	   << " != from stored = " << rs( LowerBound.back() );

  if( LB == -INFshift )
   *wlog << std::endl << "global bound = " << rs( GLB ) << " not set in MP";
  else {
   // translate it using the reference value of the hard components
   auto rf = UpRifFi.back();
   if( NrEasy )
    for( Index k = 0 ; k < NrFi ; ++k )
     if( IsEasy[ k ] )
      rf -= UpRifFi[ k ];

   GLB -= rf;
   }

  if( std::abs( LB - GLB ) >= eps * std::max( std::abs( LB ) , one ) )
   *wlog << std::endl << "translated global bound: F = " << rs( GLB )
	 << " ~ M = " << rs( LB );
   }
  }
 else {
  if( GLB > -INFshift )
    *wlog << std::endl << "finite global bound " << rs( GLB )
	  << " available but not set";

  if( LB > -INFshift )
   *wlog << std::endl << "unexpected global bound = " << rs( LB ) << " in MP";

  GLB =  f_convex ?   f_Block->get_valid_lower_bound( true )
                  : - f_Block->get_valid_upper_bound( true );
  if( GLB != LowerBound.back() )
   *wlog << std::endl << "conditional global bound = " << rs( GLB )
	 << " != from stored = " << rs( LowerBound.back() );
  }

 #if( ! USE_MPTESTER )
  // QPPenaltyMP does not allow individual lower bounds, and if a MPTester
  // is used then a QPPenaltyMP is involved anyway

  if( MPName & 1 )
   for( Index k = 0 ; k < NrFi ; ++k ) {
    if( NrEasy && IsEasy[ k ] )  // skip easy components
     continue;

    LB = Master->ReadLowerBound( k + 1 );
    auto C05LB = f_convex ?   v_c05f[ k ]->get_global_lower_bound()
                          : - v_c05f[ k ]->get_global_upper_bound();
    if( C05LB != LowerBound[ k ] )
     *wlog << std::endl << "bound( " << k << " ) = " << rs( C05LB )
	   << " != from stored = " << rs( LowerBound[ k ] );

    if( LB == -INFshift ) {
     if( C05LB > -INFshift )
      *wlog << std::endl << "bound( " << k << " ) = " << rs( C05LB )
	    << " not set in MP";
     }
    else
     if( C05LB == -INFshift )
      *wlog << std::endl << "unexpected bound( " << k << " ) = "
	    << rs( C05LB ) << " in MP";
     else {
      LB += UpRifFi[ k ];
      if( std::abs( LB - C05LB ) >= 1e-9 *
	  std::max( std::max( LB , UpRifFi[ k ] ) , double( 1 ) ) )
       *wlog << std::endl << "bound( " << k << " ): F = " << rs( C05LB )
	     << " ~ M = " << rs( LB );
      }
    }
 #endif

 }  // end( BundleSolver::CheckLBs )

/*--------------------------------------------------------------------------*/

void BundleSolver::PrintBundle( void )
{
 std::ostream * wlog = ( ( ! f_log ) || ( LogVerb <= 1 ) ) ? & std::cerr
                                                           : f_log;
 *wlog << def << std::endl << "Lambda = [ ";
 for( Index h = 0 ; h < NumVar - 1 ; ++h )
  *wlog << Lambda[ h ] << ", ";
 *wlog << Lambda.back() << " ]";

 auto Alfa = Master->ReadLinErr();
 std::vector< VarValue > G( NumVar );

 *wlog << std::endl;
 for( Index i = 0 ; i < Master->MaxName() ; ++i ) {
  *wlog << i << "\t";
  if( ItemVcblr[ i ].second >= vBPar2[ ItemVcblr[ i ].first ]
	  || ItemVcblr[ i ].second < 0 ) {
   *wlog << "[empty]" << std::endl;
   continue;
   }

  auto wFi = ItemVcblr[ i ].first;
  auto j = ItemVcblr[ i ].second;
  *wlog << wFi << "\t" << j << "\t[ ";

  v_c05f[ wFi ]->get_linearization_coefficients( G.data() ,
						 Range( 0 , NumVar ) , j );
  if( ! f_convex )
   chgsign( G.data() , NumVar );

  for( Index h = 0 ; h < NumVar - 1 ; ++h )
   *wlog << G[ h ] << ", ";

  *wlog << G.back() << " ]\t"
	 << rs( v_c05f[ wFi ]->get_linearization_constant( j ) )
	 << "\t" << Alfa[ i ] << std::endl;
  }
 }

#endif

/*--------------------------------------------------------------------------*/
/*--------------- METHODS OF BundleSolver::FakeFiOracle --------------------*/
/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetNumVar( void ) const
{
 return( bslv->NumVar );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetNrFi( void ) const
{
 return( bslv->NrFi );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetMaxName( void ) const
{
 return( bslv->vBPar2[ bslv->NrFi ] );
 }

/*--------------------------------------------------------------------------*/

bool BundleSolver::FakeFiOracle::GetUC( cIndex i )
{
 const auto var = bslv->LamVcblr[ i ];
 const auto lb = var->get_lb();
 if( lb > -Inf< ColVariable::VarValue >() ) {
  if( lb != ColVariable::VarValue( 0 ) )
   throw( std::logic_error( "finite lhs different from zero not allowed" ) );
  return( false );
  }

 for( Index j = 0 ; j < var->get_num_active() ; ++j ) {
  const auto cj = var->get_active( j );
  if( dynamic_cast< NNConstraint * >( cj ) )
   return( false );
  if( const auto bx = dynamic_cast< BoxConstraint * >( cj ) ) {
   const auto lhs = bx->get_lhs();
   if( lhs == -Inf< BoxConstraint::RHSValue >() )
    return( true );
   if( lhs == BoxConstraint::RHSValue( 0 ) )
    return( false );
   throw( std::logic_error( "finite lhs different from zero not allowed" ) );
   }
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/

LMNum BundleSolver::FakeFiOracle::GetUB( cIndex i )
{
 const auto var = bslv->LamVcblr[ i ];
 const auto ub = var->get_ub();
 if( ub < Inf< ColVariable::VarValue >() )
  return( LMNum( ub ) );

 for( Index j = 0 ; j < var->get_num_active() ; ++j )
  if( const auto bx = dynamic_cast< BoxConstraint * >( var->get_active( j )
						       ) )
   return( LMNum( bx->get_rhs() ) );

 return( Inf< LMNum >() );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNC( cIndex wFi )
{
 if( bslv->NrEasy && bslv->IsEasy[ wFi - 1 ] )
  return( bslv->IsEasy[ wFi - 1 ]->get_numcols() );
 else
  return( 0 );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNR( cIndex wFi )
{
 return( bslv->IsEasy[ wFi - 1 ]->get_numrows() );
 }

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetBNZ( cIndex wFi )
{
 return( bslv->IsEasy[ wFi - 1 ]->get_nzelements() );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::GetBDesc( cIndex wFi , int * Bbeg ,
					   int * Bind , double * Bval ,
					   double * lhs , double * rhs ,
					   double * cst ,
					   double * lbd , double * ubd )
{
 auto MILPSlv = bslv->IsEasy[ wFi - 1 ];

 if( Bbeg && Bind && Bval ) {
  // these three parameters can only be either all nullptr or all non-nullptr
  std::copy( MILPSlv->get_matbeg().begin() ,
	     MILPSlv->get_matbeg().end() , Bbeg );
  std::copy( MILPSlv->get_matind().begin() ,
	     MILPSlv->get_matind().end() , Bind );
  std::copy( MILPSlv->get_matval().begin() ,
	     MILPSlv->get_matval().end() , Bval );
  }

 if( cst ) {
  // note that in the concave case all the master problem objective changes
  // sign, so must do this part
  std::copy( MILPSlv->get_objective().begin() ,
	     MILPSlv->get_objective().end() , cst );
  if( ! bslv->f_convex )
   chgsign( cst , MILPSlv->get_numcols() );
  }

 if( lbd )
  std::copy( MILPSlv->get_var_lb().begin() ,
	     MILPSlv->get_var_lb().end() , lbd );
 if( ubd )
  std::copy( MILPSlv->get_var_ub().begin() ,
	     MILPSlv->get_var_ub().end() , ubd );

 if( lhs && rhs ) {
  // although the FiOracle interface allows setting all the parameters to
  // nullptr (save the first three) individually, OSIMPSolver never
  // requires lhs without rhs, so we don't handle the case
  for( int i = 0 ; i < MILPSlv->get_numrows() ; ++i )
   switch( MILPSlv->get_sense()[ i ] ) {
    case( 'L' ):  // <= constraint
     rhs[ i ] = MILPSlv->get_rhs()[ i ];
     lhs[ i ] = -INFshift;
     break;
    case( 'E' ):  // == constraint
     rhs[ i ] = MILPSlv->get_rhs()[ i ];
     lhs[ i ] = MILPSlv->get_rhs()[ i ];
     break;
    case( 'G' ):  // >= constraint
     rhs[ i ] = INFshift;
     lhs[ i ] = MILPSlv->get_rhs()[ i ];
     break;
    default: {    // that's a Ranged constraint
     double rngval = MILPSlv->get_rngval()[ i ];
     if( rngval > 0 ) {
      rhs[ i ] = MILPSlv->get_rhs()[ i ] + rngval;
      lhs[ i ] = MILPSlv->get_rhs()[ i ];
      }
     else {
      rhs[ i ] = MILPSlv->get_rhs()[ i ];
      lhs[ i ] = MILPSlv->get_rhs()[ i ] + rngval;
      }
     }
    }
  }
 }  // end( BundleSolver::FakeFiOracle::GetBDesc )

/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetANZ( cIndex wFi ,
					  cIndex strt , Index stp )
{
 return( static_cast< LagBFunction * >( bslv->v_c05f[ wFi - 1 ]
					)->get_A_nz() );
 }

/*--------------------------------------------------------------------------*/

void BundleSolver::FakeFiOracle::GetADesc( cIndex wFi , int * Abeg ,
					   int * Aind , double * Aval ,
					   cIndex strt , Index stp )
{
 auto LagB = static_cast< LagBFunction * >( bslv->v_c05f[ wFi - 1 ] );
 auto MILPSlv = bslv->IsEasy[ wFi - 1 ];
 c_Index nc = MILPSlv->get_numcols();
 Index count = 0;
 for( Index j = 0 ; j < nc ; ++j ) {
  Abeg[ j ] = count;
  if( auto cp = LagB->get_A_by_col( MILPSlv->variable_with_index( j ) ) ) {
   auto & mons = cp->second;
   auto it = mons.begin();
   if( strt )
    it = std::lower_bound( it , mons.end() , std::make_pair( strt , 0 ) ,
			   []( const auto & a , const auto & b )
			     { return( a.first < b.first ); } );

   if( stp < LagB->get_num_active_var() ) {
    for( ; it != mons.end() ; ++it ) {
     if( it->first >= stp )
      break;
     Aind[ count ] = it->first;
     Aval[ count++ ] = - it->second;
     }
    }
   else
    for( ; it != mons.end() ; ++it ) {
     Aind[ count ] = it->first;
     Aval[ count++ ] = - it->second;
     }

   }  // end( if( the variable has a Lagrangian term ) )
  }  // end( for( all columns ) )

 Abeg[ nc ] = count;  // end marker

 // like everything that goes in the objective, in the concave case A must be
 // changed sign
 if( ! bslv->f_convex )
  chgsign( Aval , count );

 }  // end( BundleSolver::FakeFiOracle::GetANZ )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/

Index BundleSolver::FakeFiOracle::GetGi( SgRow SubG , cIndex_Set & SGBse ,
					 cIndex Name , cIndex strt ,
					 Index stp )
{
 if( stp > bslv->NumVar )
  stp = bslv->NumVar;

 auto range = make_pair( strt , stp );

 if( Name == bslv->vBPar2[ bslv->NrFi ] )  // the 0-th component
  bslv->f_lf->get_linearization_coefficients( SubG , range );
 else                                      // a "normal" component
  bslv->v_c05f[ bslv->ItemVcblr[ Name ].first ]->
   get_linearization_coefficients( SubG , range ,
				   bslv->ItemVcblr[ Name ].second );

 // if one is maximising a concave function, change the sign
 if( ! bslv->f_convex )
  chgsign( SubG , stp - strt );

 SGBse = nullptr;
 return( stp - strt );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BundleSolverState -------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

void BundleSolverState::deserialize( const netCDF::NcGroup & group )
{
 auto nv = group.getDim( "BundleSolver_NumVar" );
 if( nv.isNull() )
  throw( std::logic_error(
		      "BundleSolverState::deserialize: missing NumVar" ) );

 NumVar = nv.getSize();

 auto l = group.getVar( "BundleSolver_Lambda" );
 if( l.isNull() )
  throw( std::logic_error(
		      "BundleSolverState::deserialize: missing Lambda" ) );

 Lambda.resize( nv.getSize() );
 l.getVar( Lambda.data() );

 auto tt = group.getVar( "BundleSolver_t" );
 if( tt.isNull() )
  throw( std::logic_error( "BundleSolverState::deserialize: missing t" ) );

 tt.getVar( &t );

 auto nf = group.getDim( "BundleSolver_NrFi" );
 if( nf.isNull() )
  throw( std::logic_error(
		      "BundleSolverState::deserialize: missing NrFi" ) );

 NrFi = nf.getSize();
 --NrFi;
 
 auto nup = group.getDim( "BundleSolver_UpFiLmbdef" );
 if( nup.isNull() ) {
  UpFiLmbdef = 0;
  UpFiLmb.clear();
  }
 else {
  UpFiLmbdef = nup.getSize();
  auto vup = group.getVar( "BundleSolver_UpFiLmb" );
  if( vup.isNull() )
   throw( std::logic_error(
		      "BundleSolverState::deserialize: missing UpFiLmb" ) );

  UpFiLmb.resize( NrFi + 1 );
  vup.getVar( UpFiLmb.data() );
  }

 auto nlw = group.getDim( "BundleSolver_LwFiLmbdef" );
 if( nlw.isNull() ) {
  LwFiLmbdef = 0;
  LwFiLmb.clear();
  }
 else {
  LwFiLmbdef = nlw.getSize();
  auto vup = group.getVar( "BundleSolver_LwFiLmb" );
  if( vup.isNull() )
   throw( std::logic_error(
		      "BundleSolverState::deserialize: missing LwFiLmb" ) );

  LwFiLmb.resize( NrFi + 1 );
  vup.getVar( LwFiLmb.data() );
  }

 auto f0 = group.getVar( "BundleSolver_Fi0Lmb" );
 if( f0.isNull() )
  Fi0Lmb = 0;
 else
  f0.getVar( &Fi0Lmb );

 auto lb = group.getVar( "BundleSolver_global_LB" );
 if( lb.isNull() )
  global_LB = -BundleSolver::INFshift;
 else
  lb.getVar( &global_LB );
 
 
 v_comp_State.resize( nf.getSize() , nullptr );

 for( BundleSolver::Index i = 0 ; i < v_comp_State.size() ; ++i ) {
  auto gi = group.getGroup( "Component_State_" + std::to_string( i ) );
  if( ! gi.isNull() )
    v_comp_State[ i ] = State::new_State( gi );
  }
 }  // end( BundleSolverState::deserialize )

/*--------------------------------------------------------------------------*/

void BundleSolverState::serialize( netCDF::NcGroup & group ) const
{
 // always call the method of the base class first
 State::serialize( group );

 auto nv = group.addDim( "BundleSolver_NumVar" , NumVar );

 ( group.addVar( "BundleSolver_Lambda" , netCDF::NcDouble() , nv ) ).putVar(
			               { 0 } , { NumVar } , Lambda.data() );

 ( group.addVar( "BundleSolver_t" , netCDF::NcDouble() ) ).putVar( &t );

 auto nfi = group.addDim( "BundleSolver_NrFi" , NrFi + 1 );

 if( UpFiLmbdef ) {
  group.addDim( "BundleSolver_UpFiLmbdef" , UpFiLmbdef );
  ( group.addVar( "BundleSolver_UpFiLmb" , netCDF::NcDouble() , nfi )
    ).putVar( { 0 } , { NrFi + 1 } , UpFiLmb.data() );
  }

 if( LwFiLmbdef ) {
  group.addDim( "BundleSolver_LwFiLmbdef" , LwFiLmbdef );
  ( group.addVar( "BundleSolver_LwFiLmb" , netCDF::NcDouble() , nfi )
    ).putVar( { 0 } , { NrFi + 1 } , LwFiLmb.data() );
  }

 if( Fi0Lmb != 0 )
  ( group.addVar( "BundleSolver_Fi0Lmb" , netCDF::NcDouble() )
    ).putVar( & Fi0Lmb );

 if( global_LB > - BundleSolver::INFshift )
  ( group.addVar( "BundleSolver_global_LB" , netCDF::NcDouble() )
    ).putVar( & global_LB );

 for( Index i = 0 ; i < NrFi ; ++i ) {
  if( ! v_comp_State[ i ] )
   continue;

  auto gi = group.addGroup( "Component_State_" + std::to_string( i ) );
  v_comp_State[ i ]->serialize( gi );
  }
 }  // end( BundleSolverState::serialize )


/*--------------------------------------------------------------------------*/
/*----------------------- End File BundleSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
