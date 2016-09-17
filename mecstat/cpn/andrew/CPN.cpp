/*****************************************************************************/
/* Name: CPN_vector.cpp                                                      */
/* Purpose: Monte Carlo simulation of a 2-dimensional CPN model using an     */
/*          over-heat bath algorithm.                                        */
/*****************************************************************************/
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <ratio>
#include <dSFMT.h>
#include <jump/dSFMT-jump.h>
#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <NTL/ZZ.h>
#include <omp.h>
#include <immintrin.h>
#include <x86intrin.h>

#ifdef __KNC__
  #define __AVX512__ 1
  #define MAX_VECTOR_SIZE 512
  #include "vectorclass.h"
  #define VECTORMATH 3
  #include "vectormath_lib.h"
#else
  #if defined(__SSE2__)
    #include "vectorclass.h"
    #ifdef VECTORMATH
      #include "vectormath_lib.h"
    #else
      #include "vectormath_trig.h"
      #include "vectormath_exp.h"
      #include "vectormath_hyp.h"
    #endif
  #endif
#endif

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif
#define HALF_PI (M_PI * 0.5)
#define TRUE 1
#define FALSE 0
#define CACHE_LINE 64
#define DOUBLE_CACHE_LINE 0

/*****************************************************************************/
/* Simulation parameters.                                                    */
/*****************************************************************************/
#ifndef N
  #define N 10
#endif
#ifndef ETA
  #define ETA 0.999
#endif
#define Z_LEN (2 * N)
#define POWER (2 * (N - 1))
#ifdef LEN
  #define L LEN
#else
  #define L 64
#endif
#define N_COR 1
#ifndef CORR_N_COR
  #define CORR_N_COR 5
#endif
#ifndef N_CF
  #define N_CF 10000
#endif
#define N_BINS 100
#ifndef N_THERM
  #define N_THERM 10000
#endif
#ifndef BETA
  #define BETA 0.8
#endif
#ifndef LOW_THETA_SKIP
  #define LOW_THETA_SKIP 0.0001
#endif
#define GEN_NEW_CFG 1
#define MEASURE_AUTO 0
#ifndef MEASUREMENTS
  #ifdef NO_BC
    #define MEASUREMENTS 25 // 11001
  #else
    #define MEASUREMENTS 31 // 01111
  #endif
#endif
#define LINK_ENERGY 1     // 000001
#define MAG_SUS 2         // 000010
#define G_k_1 4           // 000100
#define TOP_CHARGE 8      // 001000
#define WALL_CORR 16      // 010000
#define DIAG_CORR 32      // 100000

#define NUM_OBS 4
//#define NUM_OBS 5 /*((MEASUREMENTS & LINK_ENERGY) + \
                 (MEASUREMENTS & MAG_SUS) + \
                 (MEASUREMENTS & G_k_1) + \
                 (MEASUREMENTS & TOP_CHARGE))
/*
#ifdef NO_BC
  #define NUM_OBS 4
#else
  #define NUM_OBS 4
#endif*/
#define OBS_STRING 10
#ifndef OMP_NUM_THREADS
  #define OMP_NUM_THREADS 4
#endif
#ifndef SMT
  #define SMT 1
#endif
#define SMT_THREADS (SMT - 1)
#define RED 0
#define BLACK 1
#ifndef RESTART
  #define RESTART 0
#endif
#ifndef SAVE_MEASUREMENTS
  #define SAVE_MEASUREMENTS 1
#endif
#ifndef SAVE_THERM
  #define SAVE_THERM 1
#endif
#ifndef LOAD_THERM
  #define LOAD_THERM 0
#endif

/*****************************************************************************/
/* Allow the preprocessor to determine further useful simulation parameters. */
/*****************************************************************************/
#define V (L * L) // Lattice volume.
#define SITE_WIDTH (2 * (N + 2)) // Number of doubles for one lattice site.
#define CORR_LEN ((L / 2) + 1) // Length of folded correlator

#if OMP_NUM_THREADS >= N*N
  #define GK_THREADS N*N
  #define GK_STEP 1
#else
  #define GK_THREADS OMP_NUM_THREADS
  #if (N*N) % OMP_NUM_THREADS == 0
    #define GK_STEP (N*N/OMP_NUM_THREADS)
  #else
    #define GK_STEP ((N*N/OMP_NUM_THREADS) + 1)
  #endif
#endif

#define P_LEN (2 * N * N)
#define K_UNIT (2.0 * M_PI / L)

#define NUM_CORES (OMP_NUM_THREADS / SMT)
#if L % OMP_NUM_THREADS == 0
  #define L_STEP (L/NUM_CORES)
#else
  #define L_STEP ((L/NUM_CORES) + 1)
#endif

// Number of random numbers generated for each update per site.
#define Z_RANDS_PER_SITE 2
#define LAM_RANDS_PER_SITE 4

/*****************************************************************************/
/* Define vector type definitions and vector versions of mathematical        */
/* functions.                                                                */
/*****************************************************************************/
#ifdef __KNC__
  /***************************************************************************/
  /* If running on a KNC Xeon Phi, define KNC operations.                    */
  /***************************************************************************/
  typedef Vec8d VEC_DOUBLE;
  typedef Vec8q VEC_INT;
  typedef Vec16i VEC_INT_2;
  typedef Vec8db VEC_BOOL;

  #define N_d 8 // Number of doubles that fit into a vector.
  VEC_DOUBLE INDEX_VECTOR = {0, 2, 4, 6, 8, 10, 12, 14}; // x index template of a site vector.

  /***************************************************************************/
  /* Code to load an unaligned phi vector in the backward x direction.       */
  /***************************************************************************/
  static inline VEC_DOUBLE load_phi_xd(VEC_DOUBLE* phi_xd, VEC_DOUBLE* phi_xu)
  {
    VEC_INT_2 indices = {0, N_d * SITE_WIDTH + 1, N_d * SITE_WIDTH + 2,
                         N_d * SITE_WIDTH + 3, N_d * SITE_WIDTH + 4,
                         N_d * SITE_WIDTH + 5, N_d * SITE_WIDTH + 6,
                         N_d * SITE_WIDTH + 7, 0, 0, 0, 0, 0, 0, 0, 0};
    return _mm512_i32logather_pd(indices, ((double *)phi_xd) + N_d - 1, 1);
  }

  /***************************************************************************/
  /* Code to load an unaligned phi vector in the forward x direction.        */
  /***************************************************************************/
  static inline VEC_DOUBLE load_phi_xu(VEC_DOUBLE* phi_xd, VEC_DOUBLE* phi_xu)
  {
    VEC_INT_2 indices = {0, 1, 2, 3, 4, 5, 6, 7, 8 + N_d * SITE_WIDTH,
                         0, 0, 0, 0, 0, 0, 0};
    return _mm512_i32logather_pd(indices, ((double *)phi_xd) + 1, 1);
  }

  static inline VEC_DOUBLE selectd(VEC_BOOL const & mask,
                                   VEC_DOUBLE const & a,
                                   VEC_DOUBLE const & b)
  {
    return _mm512_mask_blend_pd(mask, b, a);
  }

  //VEC_DOUBLE sqrt(VEC_DOUBLE const & a) {return _mm512_sqrt_pd(a);}
  VEC_DOUBLE sin(VEC_DOUBLE const & a) {return _mm512_sin_pd(a);}
  VEC_DOUBLE asin(VEC_DOUBLE const & a) {return _mm512_asin_pd(a);}
  VEC_DOUBLE cos(VEC_DOUBLE const & a) {return _mm512_cos_pd(a);}
  VEC_DOUBLE acos(VEC_DOUBLE const & a) {return _mm512_acos_pd(a);}
  VEC_DOUBLE tan(VEC_DOUBLE const & a) {return _mm512_tan_pd(a);}
  VEC_DOUBLE atan(VEC_DOUBLE const & a) {return _mm512_atan_pd(a);}
  VEC_DOUBLE atan2(VEC_DOUBLE const & a, VEC_DOUBLE const & b)
  {return _mm512_atan2_pd(a, b);}
  VEC_DOUBLE sinh(VEC_DOUBLE const & a) {return _mm512_sinh_pd(a);}
  VEC_DOUBLE cosh(VEC_DOUBLE const & a) {return _mm512_cosh_pd(a);}
  VEC_DOUBLE tanh(VEC_DOUBLE const & a) {return _mm512_tanh_pd(a);}
  VEC_DOUBLE exp(VEC_DOUBLE const & a) {return _mm512_exp_pd(a);}
  VEC_DOUBLE log(VEC_DOUBLE const & a) {return _mm512_log_pd(a);}
#else
#ifdef __AVX__
  /***************************************************************************/
  /* If the processor supports AVX, define AVX operations.                   */
  /***************************************************************************/
  typedef Vec4d VEC_DOUBLE;
  typedef Vec4q VEC_INT;
  typedef Vec4db VEC_BOOL;

  #define N_d 4 // Number of doubles that fit into a vector.
  VEC_DOUBLE INDEX_VECTOR = {0, 2, 4, 6}; // x index template of a site vector.

  /***************************************************************************/
  /* Code to load an unaligned phi vector in the backward x direction.       */
  /***************************************************************************/
  static inline VEC_DOUBLE load_phi_xd(VEC_DOUBLE* phi_xd, VEC_DOUBLE* phi_xu)
  {
    __m256d temp = _mm256_loadu_pd(((double *)phi_xd) + N_d - 1);
    __m256d temp2 = _mm256_loadu_pd(((double *)phi_xu) - 1);
    return _mm256_blend_pd(temp, temp2, 0x0E); // 0x0E has bits 1110
  }

	/***************************************************************************/
	/* Code to load an unaligned phi vector in the forward x direction.        */
	/***************************************************************************/
  static inline VEC_DOUBLE load_phi_xu(VEC_DOUBLE* phi_xd, VEC_DOUBLE* phi_xu)
  {
    __m256d temp = _mm256_loadu_pd(((double *)phi_xd) + 1);
    __m256d temp2 = _mm256_loadu_pd(((double *)phi_xu) - N_d + 1);
    return _mm256_blend_pd(temp, temp2, 0x08); // 0x08 has bits 1000
  }

  /***************************************************************************/
  /* Code to blend two vectors in each way that is required to construct     */
  /* the phi vectors of adjacent sites.                                      */
  /***************************************************************************/
  static inline VEC_DOUBLE blend_phi(VEC_DOUBLE* phi_1, VEC_DOUBLE* phi_2,
                                     int type)
  {
    switch (type)
    {
      case 0: // 0000
        return *phi_1;
      case 1: // 1000
        return _mm256_blend_pd(_mm256_loadu_pd(((double *)phi_1) + type),
                               _mm256_loadu_pd(((double *)phi_2) + type - N_d),
                               0x08);
      case 2: // 1100
        return _mm256_blend_pd(_mm256_loadu_pd(((double *)phi_1) + type),
                               _mm256_loadu_pd(((double *)phi_2) + type - N_d),
                               0x0C);
      case 3: // 1110
        return _mm256_blend_pd(_mm256_loadu_pd(((double *)phi_1) + type),
                               _mm256_loadu_pd(((double *)phi_2) + type - N_d),
                               0x0E);
    }
  }

  /***************************************************************************/
  /* Blend test code.                                                        */
  /***************************************************************************/
  static inline void test_blend()
  {
    VEC_DOUBLE v1 = {0.0, 1.0, 2.0, 3.0};
    VEC_DOUBLE v2 = {4.0, 5.0, 6.0, 7.0};

    VEC_DOUBLE v3;
    for (int type = 0; type < N_d; ++type)
    {
      std::cout << "Testing case " << type << std::endl;
      v3 = blend_phi(&v1, &v2, type);
      for (int i = 0; i < N_d; ++i)
      {
        std::cout << "v[" << i << "] = " << v3[i] << std::endl;
      }
    }
    return;
  }
#else
#ifdef __SSE2__
  /***************************************************************************/
  /* If the processor supports SSE2, define SSE2 operations.                 */
  /***************************************************************************/
  typedef Vec2d VEC_DOUBLE;
  typedef Vec2q VEC_INT;
  typedef Vec2db VEC_BOOL;

  #define N_d 2 // Number of doubles that fit into a vector.
  VEC_DOUBLE INDEX_VECTOR = {0, 2}; // x index template of a site vector.

  /***************************************************************************/
  /* Code to load an unaligned phi vector in the backward x direction.       */
  /***************************************************************************/
  static inline VEC_DOUBLE load_phi_xd(VEC_DOUBLE* phi_xd, VEC_DOUBLE* phi_xu)
  {
     __m128d temp = _mm_loadu_pd(((double *)phi_xd) + N_d - 1);
     __m128d temp2 = _mm_loadu_pd(((double *)phi_xu) - 1);
    return _mm_move_sd(temp2, temp);
  }

  /***************************************************************************/
  /* Code to load an unaligned phi vector in the forward x direction.        */
  /***************************************************************************/
  static inline VEC_DOUBLE load_phi_xu(VEC_DOUBLE* phi_xd, VEC_DOUBLE* phi_xu)
  {
     __m128d temp = _mm_loadu_pd(((double *)phi_xd) + 1);
     __m128d temp2 = _mm_loadu_pd(((double *)phi_xu) - N_d + 1);
    return _mm_move_sd(temp2, temp);
  }

  /***************************************************************************/
  /* Code to blend the first entry of one vector with the second entry of    */
  /* another.                                                                */
  /***************************************************************************/
  static inline VEC_DOUBLE blend_phi(VEC_DOUBLE* phi_1, VEC_DOUBLE* phi_2,
                                     int type)
  {
    switch (type)
    {
      case 0: // No blend. Return phi_1
        return *phi_1;
        break;
      case 1: // 1:1 blend.
        return _mm_move_sd(_mm_loadu_pd(((double *)phi_2) - N_d + 1),
                           _mm_loadu_pd(((double *)phi_1) + 1)); // 1000
        break;
      default:
        return *phi_2; // No blend. Return phi_2
        break;
    }
  }
#else
  /***************************************************************************/
  /* No vectorisation - define variables for compatibility. Warning: without */
  /* vectorisation code is very slow.                                        */
  /***************************************************************************/
  #define N_d 1
  typedef double VEC_DOUBLE;
  typedef int VEC_INT;
  typedef bool VEC_BOOL;
  VEC_DOUBLE INDEX_VECTOR = {0};

  /***************************************************************************/
  /* Compatibility with load functions                                       */
  /***************************************************************************/
  static inline VEC_DOUBLE load_phi_xd(VEC_DOUBLE* phi_xd, VEC_DOUBLE* phi_xu)
  {
    return *phi_xd;
  }

  static inline VEC_DOUBLE load_phi_xu(VEC_DOUBLE* phi_xd, VEC_DOUBLE* phi_xu)
  {
    return *phi_xu;
  }

  /***************************************************************************/
  /* Compatibility with vector component manipulation functions.             */
  /***************************************************************************/
  static inline VEC_DOUBLE selectd(VEC_BOOL err, VEC_DOUBLE replace_vec,
                                   VEC_DOUBLE orig_vec)
  {
    if (err)
    {
      return replace_vec;
    }
    else
    {
      return orig_vec;
    }
  }
  static inline double horizontal_add(VEC_DOUBLE x)
  {
    return x;
  }

  /***************************************************************************/
  /* Boolean operation compatibility.                                        */
  /***************************************************************************/
  static inline VEC_BOOL horizontal_or(VEC_BOOL b)
  {
    return b;
  }
  static inline VEC_BOOL horizontal_count(VEC_BOOL b)
  {
    if(b)
    {
      return N_d;
    }
    else
    {
      return 0;
    }
  }

  /***************************************************************************/
  /* Blend compatibility function.                                           */
  /***************************************************************************/
  static inline VEC_DOUBLE blend_phi(VEC_DOUBLE* phi_1, VEC_DOUBLE* phi_2,
                                     int type)
  {
    return *phi_1;
  }

  /***************************************************************************/
  /* Min/max functions.                                                      */
  /***************************************************************************/
  static inline VEC_DOUBLE min(VEC_DOUBLE a, VEC_DOUBLE b)
  {
    return std::min(a, b);
  }
  static inline VEC_DOUBLE max(VEC_DOUBLE a, VEC_DOUBLE b)
  {
    return std::max(a, b);
  }
#endif
#endif
#endif

/*****************************************************************************/
/* Function to return a vector containing entries for the sites n sites      */
/* away from the reference vector.                                           */
/*****************************************************************************/
/*static inline VEC_DOUBLE load_phi_n_up(VEC_DOUBLE* phi_red_1,
                                       VEC_DOUBLE* phi_red_2,
                                       VEC_DOUBLE* phi_black_1,
                                       VEC_DOUBLE* phi_black_2,
                                       int shift, long dx)
{
  VEC_DOUBLE *temp1;
  VEC_DOUBLE *temp2;
  int case_no = ((dx + shift) / 2) % N_d;
  if ((dx + shift) % 2 == 0)
  {
    return blend_phi(phi_red_1, phi_red_2, case_no);
  }
}*/

/*****************************************************************************/
/* Function to return the indices for a site vector.
/*****************************************************************************/
static inline VEC_DOUBLE x_index_vector(long index, int shift)
{
  return ((1 + (index / 2)) * 2 * N_d) + shift + INDEX_VECTOR;
}


/*****************************************************************************/
/* Define variables needed for interpreting vectorized memory layout.        */
/*****************************************************************************/
// Number of red/black sites leftover after partitioning lattice into vectors.
#define EXTRA ((L / 2) % N_d)

// The value of L which fits all the complete vectors
#define MIN_L (((L / 2) / N_d) * 2 * N_d)

#if EXTRA > 0
  // The value of L to use when partioning lattice into vectors
  #define MAX_L (MIN_L + 2 * N_d)
#else
  #define MAX_L MIN_L
#endif

// Code to support lattice sizes where L/2 != N_d (still under development)
#if EXTRA > 0
 #define END_SHIFT EXTRA
#else
 #define END_SHIFT N_d
#endif

// The number of vectors in a single row of the lattice
#define NUM_VEC_BLOCKS (MAX_L / N_d)

// The number of red (or black) vectors in a single row of the lattice
#define NUM_RB_BLOCKS (NUM_VEC_BLOCKS / 2)

// The vector block width in terms of number of doubles
#define VEC_WIDTH (2 * (N + 2) * N_d)
#define VEC_Z_WIDTH (2 * N * N_d)
#define VEC_LAM_WIDTH (4 * N_d)
#define Z_STEP (2 * N_d)
#define LAM_XD_SHIFT (VEC_Z_WIDTH + (2 * N_d) + N_d - 1)

// Define how many blocks each thread updates / measures.
#define NUM_BLOCK_UPS (L * NUM_RB_BLOCKS / OMP_NUM_THREADS)
#define SPARE (L * NUM_RB_BLOCKS - (NUM_BLOCK_UPS * OMP_NUM_THREADS))

// Settings for simulations with a rectangular lattice, but square fiducial
// box where measurements are taken.
#ifdef RECTANGLE
  // Additional space to add to either end of side to make long direction.
  #ifndef L_GAP
    #define L_GAP 40
  #endif

  // Total length along long direction
  #define L_T (L + (2 * L_GAP))

  // Total volume
  #define V_TOT (L * L_T)

  // Define how many blocks each thread updates.
  #define NUM_BLOCK_UPS_RECT (L_T * NUM_RB_BLOCKS / OMP_NUM_THREADS)
  #define SPARE_RECT (L_T * NUM_RB_BLOCKS - (NUM_BLOCK_UPS_RECT * OMP_NUM_THREADS))
#else
  // Auxiliary setting for memory allocation and pointer arithmetic.
  #define L_GAP 0
  #define L_T L
  #define V_TOT V
#endif

// Define how many FLOPs each measurement performs
#define E_FLOPS (((16.0 * N + 4.5) * N_d * NUM_VEC_BLOCKS * L) + 5 * OMP_NUM_THREADS * N_d)
#define G_FLOPS ((8.0 * L * (11.25 + (N * N * NUM_VEC_BLOCKS * (N_d + 0.5)))) + (4.0 * N * N * (OMP_NUM_THREADS + 1.0)) + 4.0)
#define TOP_FLOPS (((32.0 * N) + 65.0) * N_d * L * NUM_VEC_BLOCKS + OMP_NUM_THREADS * (N_d + 1.0) + 3.0)
#define CORR_FLOPS (L * ((2 * L * ((4 * N * N) + (2 * N) + 1)) + OMP_NUM_THREADS + 2))

/*****************************************************************************/
/* Structure containing information regarding update limits.                 */
/*****************************************************************************/
struct update_limits
{
  /***************************************************************************/
  /* Start and end indices in each dimension for updates.                    */
  /***************************************************************************/
  unsigned int u_min_t;
  unsigned int u_max_t;
  unsigned int u_min_x;
  unsigned int u_max_x;
  long num_block_ups;

  /***************************************************************************/
  /* Start and end indices in each dimension for measurements.               */
  /***************************************************************************/
  unsigned int m_min_t;
  unsigned int m_max_t;
  unsigned int m_min_x;
  unsigned int m_max_x;
  long num_block_meas;
};

// Define the type of the P_t_avg operator for computing the wall-wall
// correlation function.
/*#ifdef NO_BC
  typedef VEC_DOUBLE P_type;
#else
  typedef double P_type;
#endif*/
typedef double P_type;

/*****************************************************************************/
/* Function prototypes.                                                      */
/*****************************************************************************/
static inline void F_phi_z(VEC_DOUBLE*, VEC_DOUBLE*, VEC_DOUBLE*, VEC_DOUBLE*,
                           VEC_DOUBLE*, VEC_DOUBLE*, VEC_DOUBLE*, VEC_DOUBLE*);
static inline void F_phi_lam(VEC_DOUBLE*, VEC_DOUBLE*, VEC_DOUBLE*);
static inline int update_z(VEC_DOUBLE*, VEC_DOUBLE*, dsfmt_t*);
static inline int update_lambda(VEC_DOUBLE*, VEC_DOUBLE*, dsfmt_t*);
static inline unsigned long update_loop(VEC_DOUBLE***, VEC_DOUBLE***,
                                        dsfmt_t***, unsigned int, unsigned int,
                                        unsigned int, unsigned int, int, int);
static inline unsigned long heatbath_update(VEC_DOUBLE***, VEC_DOUBLE***,
                                            dsfmt_t****, update_limits*,
                                            long*, long*);
static inline double S_1(VEC_DOUBLE***, VEC_DOUBLE***);
static inline void G_k_loop(VEC_DOUBLE***, double*, double*,
                            update_limits*);
static inline double G_k(VEC_DOUBLE***, double**, double**, double*,
                         update_limits*);
static inline void top_loop(VEC_DOUBLE***, VEC_DOUBLE*, VEC_DOUBLE*,
                            VEC_DOUBLE*, VEC_DOUBLE*, VEC_DOUBLE*,
                            update_limits*);
static inline double topological_charge(VEC_DOUBLE***, update_limits*);
static inline double comp_MC_ests(VEC_DOUBLE***, VEC_DOUBLE***, double**, double**, double*, P_type**,
                                  double **, const long, const long, update_limits*);
static inline double jackknife(double*, long, double);
void monte_carlo(std::string**, std::string, double*, double*);
static inline void generate_start_config(VEC_DOUBLE***, VEC_DOUBLE***,
                                         dsfmt_t****, VEC_DOUBLE*,
                                         std::string, long*, long*);
static inline void set_z(VEC_DOUBLE***, VEC_DOUBLE***, dsfmt_t****,
                         unsigned int, unsigned int, int, long*, long*);
static inline void initialize_lattice(VEC_DOUBLE***, VEC_DOUBLE***,
                                      dsfmt_t****, long*, long*);
static inline void initialize_counters(long **, long **, std::string);
static inline void initialize_rng(dsfmt_t*****, dsfmt_t**, unsigned long,
                                  long*, long*, std::string);
static inline void construct_lattice(VEC_DOUBLE****, VEC_DOUBLE****,
                                     VEC_DOUBLE*);
static inline void define_update_limits(update_limits*);
static inline double thermalise_lattice(VEC_DOUBLE***, VEC_DOUBLE***,
                                        dsfmt_t****, update_limits*, long*,
                                        long*, VEC_DOUBLE*, std::string);
static inline double run_simulation(VEC_DOUBLE***, VEC_DOUBLE***, dsfmt_t****,
                                    update_limits*, dsfmt_t*, long*, long*, VEC_DOUBLE*,
                                    std::string, double**, double**, double**,
                                    std::string**);

/*****************************************************************************/
/* Execute the Monte Carlo algorithm.                                        */
/*****************************************************************************/
int //(void)
{
  /***************************************************************************/
  /* Initialize variables.                                                   */
  /***************************************************************************/
  std::string *obs_list[NUM_OBS];
  int jj = 0;
#if (MEASUREMENTS & LINK_ENERGY)
  obs_list[jj++] = new std::string("E");
#endif
//#ifndef NO_BC
#if (MEASUREMENTS & MAG_SUS)
  obs_list[jj++] = new std::string("mag_sus");
#endif
#if (MEASUREMENTS & G_k_1)
  obs_list[jj++] = new std::string("G_1");
#endif
#if (MEASUREMENTS & TOP_CHARGE)
  obs_list[jj++] = new std::string("top_charge");
#endif
  std::stringstream path;
  //path << "/home/al1g13/Documents/CPN_Data/N=" << N << ",L=" << L << ",beta=" << BETA;
  //path << "/home/User/CPN_parallel/Data/
  path << "N=" << N << ",L=" << L << ",beta=" << BETA;
#ifdef NO_BC
  path << "_no_bc";
#endif
#ifdef RECTANGLE
  path << "_" << L_GAP;
#endif

#ifdef HEATBATH
 path << "_heatbath";
#endif
  std::cout << path.str() << std::endl;

  double* MC_avg;
  double* MC_err;
  MC_avg = new double[NUM_OBS];
  MC_err = new double[NUM_OBS];
  std::cout << "N = " << N << ", beta = " << BETA << ", L = " << L << "\n";
  std::cout << "N_therm = " << N_THERM << ", N_cf = " << N_CF << std::endl;
  std::cout << "Using " << OMP_NUM_THREADS << " threads." << std::endl;
  std::cout << "E = " << E_FLOPS << ", G = " << G_FLOPS << ", t = " << TOP_FLOPS;
  std::cout << ", corr = " << CORR_FLOPS << std::endl;

  /***************************************************************************/
  /* Execute the main function, time it, and print the results.              */
  /***************************************************************************/
  auto t1 = std::chrono::high_resolution_clock::now();
  monte_carlo(obs_list, path.str(), MC_avg, MC_err);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "Took ";
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()/1000.0;
  std::cout << "s to run simulation." << std::endl;

  for (int n = 0; n < NUM_OBS; ++n)
  {
    std::cout << std::setprecision(10) << *(obs_list[n]) << " ";
    std::cout << MC_avg[n] << " ";
    std::cout << MC_err[n] << std::endl;
  }
  double ms = MC_avg[1];
  double ms_err = MC_err[1];
  double G_1 = MC_avg[2];
  double G_1_err = MC_err[2];

#if (MEASUREMENTS & (MAG_SUS + G_k_1))
  double corr_len = (pow(2 * sin(M_PI / L), -2.0) * (ms/G_1 - 1.0));
  corr_len = sqrt(corr_len);
  double cl_err = (pow(2 * sin(M_PI / L), -2.0) * (ms/G_1) *
                   sqrt(pow(G_1_err/G_1, 2.0) + pow(ms_err/ms, 2.0)));
  cl_err = 0.5 * cl_err / corr_len;
  std::cout << "xi " << corr_len << " " << cl_err << "\n";
#endif

  return 1;
} /* main */


/*****************************************************************************/
/* Utility function for test purposes.                                       */
/*****************************************************************************/
static inline void print_lattice(VEC_DOUBLE*** z, int tt)
{
  long counter = 0;
  double *ptr;
  for (int t = tt; t < tt + 1; ++t)
  {
    std::cout << "_____" << t << "_____" << std::endl;
    for (int x = 0; x < NUM_VEC_BLOCKS; ++x)
    {
      std::cout << "_____BLOCK_" << x << "_____" << std::endl;
      ptr = (double *)z[t][x];
      int counter_2 = 0;
      for (int n = 0; n < VEC_WIDTH; ++n)
      {
        std::cout << ptr[n] << " ";
        ++counter_2;
        if (counter_2 >= N_d)
        {
          std::cout << "\n";
          counter_2 = 0;
        }
      }
    }
  }
}

/*****************************************************************************/
/* Force term for z variable with shift 0.                                   */
/*****************************************************************************/
static inline void F_phi_z_0(VEC_DOUBLE * phi_z_tu,
                             VEC_DOUBLE * phi_z_td,
                             VEC_DOUBLE * phi_z_xu,
                             VEC_DOUBLE * phi_z_xd,
                             VEC_DOUBLE * phi_l,
                             VEC_DOUBLE * phi_l_td,
                             VEC_DOUBLE * phi_l_xd,
                             VEC_DOUBLE * phi_l_xu,
                             VEC_DOUBLE * F_phi)
{
   VEC_DOUBLE phi_l_xd_r = load_phi_xd(&phi_l_xd[0], &phi_l_xu[0]);
   VEC_DOUBLE phi_l_xd_i = load_phi_xd(&phi_l_xd[1], &phi_l_xu[1]);
  for (int n = 0; n < (2 * N); n += 2)
  {
    /*************************************************************************/
    /* Construct the correct phi_z_xd and phi_l_xd vectors.                  */
    /*************************************************************************/
     VEC_DOUBLE phi_z_xd_r = load_phi_xd(&phi_z_xd[n], &phi_z_xu[n]);
     VEC_DOUBLE phi_z_xd_i = load_phi_xd(&phi_z_xd[n + 1],
                                                 &phi_z_xu[n + 1]);

    /*************************************************************************/
    /* Construct phi_z_xd terms.                                             */
    /*************************************************************************/
    VEC_DOUBLE accumulator_r0 = phi_l_xd_r * phi_z_xd_r;
    VEC_DOUBLE accumulator_i0 = phi_l_xd_r * phi_z_xd_i;
    VEC_DOUBLE accumulator_r1 = -(phi_l_xd_i * phi_z_xd_i);
    VEC_DOUBLE accumulator_i1 = phi_l_xd_i * phi_z_xd_r;

    /*************************************************************************/
    /* Construct phi_z_xu terms.                                             */
    /*************************************************************************/
    VEC_DOUBLE accumulator_r2 = phi_l[3] * phi_z_xu[n + 1];
    VEC_DOUBLE accumulator_i2 = -(phi_l[3] * phi_z_xu[n]);
    VEC_DOUBLE accumulator_r3 = phi_l[2] * phi_z_xu[n];
    VEC_DOUBLE accumulator_i3 = phi_l[2] * phi_z_xu[n + 1];

    /*************************************************************************/
    /* Construct phi_z_td terms.                                             */
    /*************************************************************************/
    accumulator_r0 -= phi_l_td[1] * phi_z_td[n + 1];
    accumulator_i0 += phi_l_td[1] * phi_z_td[n];
    accumulator_r1 += phi_l_td[0] * phi_z_td[n];
    accumulator_i1 += phi_l_td[0] * phi_z_td[n + 1];

    /*************************************************************************/
    /* Construct phi_z_tu terms.                                             */
    /*************************************************************************/
    accumulator_r2 += phi_l[0] * phi_z_tu[n];
    accumulator_i2 += phi_l[0] * phi_z_tu[n + 1];
    accumulator_r3 += phi_l[1] * phi_z_tu[n + 1];
    accumulator_i3 -= phi_l[1] * phi_z_tu[n];

    F_phi[n] = 2 * (accumulator_r0 + accumulator_r1 + accumulator_r2 +
                    accumulator_r3);
    F_phi[n + 1] = 2 * (accumulator_i0 + accumulator_i1 + accumulator_i2 +
                        accumulator_i3);
  }
  return;
} /* F_phi_z_0 */

/*****************************************************************************/
/* Force term for z variable with shift 1.                                   */
/*****************************************************************************/
static inline void F_phi_z_1(VEC_DOUBLE * phi_z_tu,
                             VEC_DOUBLE * phi_z_td,
                             VEC_DOUBLE * phi_z_xu,
                             VEC_DOUBLE * phi_z_xd,
                             VEC_DOUBLE * phi_l,
                             VEC_DOUBLE * phi_l_td,
                             VEC_DOUBLE * phi_l_xd,
                             VEC_DOUBLE * F_phi)
{
  for (int n = 0; n < Z_LEN; n = n + 2)
  {
    /*************************************************************************/
    /* Construct the correct phi_z_xu vector.                                */
    /*************************************************************************/
    VEC_DOUBLE phi_z_xu_r = load_phi_xu(&phi_z_xd[n], &phi_z_xu[n]);
    VEC_DOUBLE phi_z_xu_i = load_phi_xu(&phi_z_xd[n + 1], &phi_z_xu[n + 1]);

    /*************************************************************************/
    /* Construct phi_z_xu terms.                                             */
    /*************************************************************************/
    VEC_DOUBLE accumulator_r0 = phi_l[2] * phi_z_xu_r;
    VEC_DOUBLE accumulator_i0 = phi_l[2] * phi_z_xu_i;
    VEC_DOUBLE accumulator_r1 = phi_l[3] * phi_z_xu_i;
    VEC_DOUBLE accumulator_i1 = -(phi_l[3] * phi_z_xu_r);

    /*************************************************************************/
    /* Construct phi_z_xd terms.                                             */
    /*************************************************************************/
    VEC_DOUBLE accumulator_r2 = -(phi_l_xd[1] * phi_z_xd[n + 1]);
    VEC_DOUBLE accumulator_i2 = phi_l_xd[1] * phi_z_xd[n];
    VEC_DOUBLE accumulator_r3 = phi_l_xd[0] * phi_z_xd[n];
    VEC_DOUBLE accumulator_i3 = phi_l_xd[0] * phi_z_xd[n + 1];

    /*************************************************************************/
    /* Construct phi_z_td terms.                                             */
    /*************************************************************************/
    accumulator_r0 -= phi_l_td[1] * phi_z_td[n + 1];
    accumulator_i0 += phi_l_td[1] * phi_z_td[n];
    accumulator_r1 += phi_l_td[0] * phi_z_td[n];
    accumulator_i1 += phi_l_td[0] * phi_z_td[n + 1];

    /*************************************************************************/
    /* Construct phi_z_tu terms.                                             */
    /*************************************************************************/
    accumulator_r2 += phi_l[0] * phi_z_tu[n];
    accumulator_i2 += phi_l[0] * phi_z_tu[n + 1];
    accumulator_r3 += phi_l[1] * phi_z_tu[n + 1];
    accumulator_i3 -= phi_l[1] * phi_z_tu[n];

    F_phi[n] = 2 * (accumulator_r0 + accumulator_r1 + accumulator_r2 +
                    accumulator_r3);
    F_phi[n + 1] = 2 * (accumulator_i0 + accumulator_i1 + accumulator_i2 +
                        accumulator_i3);
  }
  return;
} /* F_phi_z_1 */

#ifdef NO_BC
/*****************************************************************************/
/* 3-point Force term for z variable with shift 0.                           */
/*****************************************************************************/
static inline void F_phi_z_3pt_0(VEC_DOUBLE * phi_z_t,
                                 VEC_DOUBLE * phi_z_xu,
                                 VEC_DOUBLE * phi_z_xd,
                                 VEC_DOUBLE * phi_l,
                                 VEC_DOUBLE * phi_l_t,
                                 VEC_DOUBLE * phi_l_xd,
                                 VEC_DOUBLE * phi_l_xu,
                                 VEC_DOUBLE * F_phi)
{
   VEC_DOUBLE phi_l_xd_r = load_phi_xd(&phi_l_xd[0], &phi_l_xu[0]);
   VEC_DOUBLE phi_l_xd_i = load_phi_xd(&phi_l_xd[1], &phi_l_xu[1]);
  for (int n = 0; n < (2 * N); n += 2)
  {
    /*************************************************************************/
    /* Construct the correct phi_z_xd and phi_l_xd vectors.                  */
    /*************************************************************************/
     VEC_DOUBLE phi_z_xd_r = load_phi_xd(&phi_z_xd[n], &phi_z_xu[n]);
     VEC_DOUBLE phi_z_xd_i = load_phi_xd(&phi_z_xd[n + 1],
                                                 &phi_z_xu[n + 1]);

    /*************************************************************************/
    /* Construct phi_z_xd terms.                                             */
    /*************************************************************************/
    VEC_DOUBLE accumulator_r0 = phi_l_xd_r * phi_z_xd_r;
    VEC_DOUBLE accumulator_i0 = phi_l_xd_r * phi_z_xd_i;
    VEC_DOUBLE accumulator_r1 = -(phi_l_xd_i * phi_z_xd_i);
    VEC_DOUBLE accumulator_i1 = phi_l_xd_i * phi_z_xd_r;

    /*************************************************************************/
    /* Construct phi_z_xu terms.                                             */
    /*************************************************************************/
    VEC_DOUBLE accumulator_r2 = phi_l[3] * phi_z_xu[n + 1];
    VEC_DOUBLE accumulator_i2 = -(phi_l[3] * phi_z_xu[n]);
    VEC_DOUBLE accumulator_r3 = phi_l[2] * phi_z_xu[n];
    VEC_DOUBLE accumulator_i3 = phi_l[2] * phi_z_xu[n + 1];

    /*************************************************************************/
    /* Construct phi_z_t terms.                                              */
    /*************************************************************************/
    accumulator_r0 -= phi_l_t[1] * phi_z_t[n + 1];
    accumulator_i0 += phi_l_t[1] * phi_z_t[n];
    accumulator_r1 += phi_l_t[0] * phi_z_t[n];
    accumulator_i1 += phi_l_t[0] * phi_z_t[n + 1];

    F_phi[n] = 2 * (accumulator_r0 + accumulator_r1 + accumulator_r2 +
                    accumulator_r3);
    F_phi[n + 1] = 2 * (accumulator_i0 + accumulator_i1 + accumulator_i2 +
                        accumulator_i3);
  }
  return;
} /* F_phi_z_3pt_0 */


/*****************************************************************************/
/* 3 point Force term for z variable with shift 1.                           */
/*****************************************************************************/
static inline void F_phi_z_3pt_1(VEC_DOUBLE * phi_z_t,
                                 VEC_DOUBLE * phi_z_xu,
                                 VEC_DOUBLE * phi_z_xd,
                                 VEC_DOUBLE * phi_l,
                                 VEC_DOUBLE * phi_l_t,
                                 VEC_DOUBLE * phi_l_xd,
                                 VEC_DOUBLE * F_phi)
{
  for (int n = 0; n < Z_LEN; n = n + 2)
  {
    /*************************************************************************/
    /* Construct the correct phi_z_xu vector.                                */
    /*************************************************************************/
    VEC_DOUBLE phi_z_xu_r = load_phi_xu(&phi_z_xd[n], &phi_z_xu[n]);
    VEC_DOUBLE phi_z_xu_i = load_phi_xu(&phi_z_xd[n + 1], &phi_z_xu[n + 1]);

    /*************************************************************************/
    /* Construct phi_z_xu terms.                                             */
    /*************************************************************************/
    VEC_DOUBLE accumulator_r0 = phi_l[2] * phi_z_xu_r;
    VEC_DOUBLE accumulator_i0 = phi_l[2] * phi_z_xu_i;
    VEC_DOUBLE accumulator_r1 = phi_l[3] * phi_z_xu_i;
    VEC_DOUBLE accumulator_i1 = -(phi_l[3] * phi_z_xu_r);

    /*************************************************************************/
    /* Construct phi_z_xd terms.                                             */
    /*************************************************************************/
    VEC_DOUBLE accumulator_r2 = -(phi_l_xd[1] * phi_z_xd[n + 1]);
    VEC_DOUBLE accumulator_i2 = phi_l_xd[1] * phi_z_xd[n];
    VEC_DOUBLE accumulator_r3 = phi_l_xd[0] * phi_z_xd[n];
    VEC_DOUBLE accumulator_i3 = phi_l_xd[0] * phi_z_xd[n + 1];

    /*************************************************************************/
    /* Construct phi_z_td terms.                                             */
    /*************************************************************************/
    accumulator_r0 -= phi_l_t[1] * phi_z_t[n + 1];
    accumulator_i0 += phi_l_t[1] * phi_z_t[n];
    accumulator_r1 += phi_l_t[0] * phi_z_t[n];
    accumulator_i1 += phi_l_t[0] * phi_z_t[n + 1];

    F_phi[n] = 2 * (accumulator_r0 + accumulator_r1 + accumulator_r2 +
                    accumulator_r3);
    F_phi[n + 1] = 2 * (accumulator_i0 + accumulator_i1 + accumulator_i2 +
                        accumulator_i3);
  }
  return;
} /* F_phi_z_3pt_1 */

#endif

/*****************************************************************************/
/* Force term for lambda variable with shift 0.                              */
/*****************************************************************************/
static inline void F_phi_lam(VEC_DOUBLE * phi_z,
                             VEC_DOUBLE * phi_z_up,
                             VEC_DOUBLE * F_l)
{
  F_l[0] = 0.0;
  F_l[1] = 0.0;
  for (int n = 0; n < (2 * N); n += 2)
  {
    VEC_DOUBLE accumulator_r0 = phi_z[n] * phi_z_up[n];
    VEC_DOUBLE accumulator_i0 = phi_z[n] * phi_z_up[n + 1];
    VEC_DOUBLE accumulator_r1 = phi_z[n + 1] * phi_z_up[n + 1];
    VEC_DOUBLE accumulator_i1 = phi_z[n + 1] * phi_z_up[n];
    F_l[0] += accumulator_r0 + accumulator_r1;
    F_l[1] += accumulator_i0 - accumulator_i1;
  }
  F_l[0] *= 2.0;
  F_l[1] *= 2.0;

  return;
} /* F_phi_lam */


/*****************************************************************************/
/* Force term for lambda variable with shift 1.                              */
/*****************************************************************************/
static inline void F_phi_lam_1(VEC_DOUBLE * phi_z,
                               VEC_DOUBLE * phi_z_down,
                               VEC_DOUBLE * phi_z_up,
                               VEC_DOUBLE * F_l)
{
  F_l[0] = 0;
  F_l[1] = 0;
  for (int n = 0; n < (2 * N); n += 2)
  {
    VEC_DOUBLE phi_z_up_r = load_phi_xu(&phi_z_down[n], &phi_z_up[n]);
    VEC_DOUBLE accumulator_r0 = phi_z[n] * phi_z_up_r;
    VEC_DOUBLE phi_z_up_i = load_phi_xu(&phi_z_down[n + 1], &phi_z_up[n + 1]);
    VEC_DOUBLE accumulator_i0 = phi_z[n] * phi_z_up_i;
    VEC_DOUBLE accumulator_r1 = phi_z[n + 1] * phi_z_up_i;
    VEC_DOUBLE accumulator_i1 = phi_z[n + 1] * phi_z_up_r;
    F_l[0] += accumulator_r0 + accumulator_r1;
    F_l[1] += accumulator_i0 - accumulator_i1;
  }
  F_l[0] *= 2;
  F_l[1] *= 2;

  return;
} /* F_phi_lam */


/*****************************************************************************/
/* Given a 2D vector unit perp_v that isn't perpendicular to unit vector     */
/* unit_v, adjust it so that it is indeed perpendicular to this vector.      */
/*****************************************************************************/
static inline void adjust_2D_perp_vector(VEC_DOUBLE *perp_v, VEC_DOUBLE *unit_v)
{
  VEC_DOUBLE perp_mod(0);
  VEC_DOUBLE perp_dot(0);
  perp_dot = unit_v[0] * perp_v[0];
  perp_dot += unit_v[1] * perp_v[1];
  perp_v[0] -= unit_v[0] * perp_dot;
  perp_v[1] -= unit_v[1] * perp_dot;
  perp_mod += perp_v[0] * perp_v[0];
  perp_mod += perp_v[1] * perp_v[1];
  VEC_DOUBLE inv_perp_mod = 1 / perp_mod;
  perp_v[0] *= inv_perp_mod;
  perp_v[1] *= inv_perp_mod;
}

/*****************************************************************************/
/* Given a 2N-D vector unit perp_v that isn't perpendicular to unit vector   */
/* unit_v, adjust it so that it is indeed perpendicular to this vector.      */
/*****************************************************************************/
static inline void adjust_2ND_perp_vector(VEC_DOUBLE *perp_v,
                                          VEC_DOUBLE *unit_v)
{
  VEC_DOUBLE perp_mod = 0;
  VEC_DOUBLE perp_dot = 0;
  for (int n = 0; n < Z_LEN; n += 2)
  {
    perp_dot += unit_v[n] * perp_v[n];
    perp_dot += unit_v[n + 1] * perp_v[n + 1];
    perp_v[n] -= unit_v[n] * perp_dot;
    perp_v[n + 1] -= unit_v[n + 1] * perp_dot;
    perp_mod += perp_v[n] * perp_v[n];
    perp_mod += perp_v[n + 1] * perp_v[n + 1];
  }
  VEC_DOUBLE inv_perp_mod = 1 / perp_mod;

  for (int n = 0; n < Z_LEN; n += 2)
  {
    perp_v[n] *= inv_perp_mod;
    perp_v[n + 1] *= inv_perp_mod;
  }
}


/*****************************************************************************/
/* Generate a vector of random numbers.                                      */
/*****************************************************************************/
static inline VEC_DOUBLE gen_rand(dsfmt_t** r, long* stream_length,
                                  VEC_BOOL &accepted)
{
    double p_arr[N_d] __attribute__((aligned(N_d * sizeof(double))));
    VEC_DOUBLE p;
    #if defined(__SSE2__)
        for (int i = 0; i < N_d; ++i)
        {
          if (!accepted[i])
          {
            p_arr[i] = dsfmt_genrand_open_close(r[i]);
            ++stream_length[i];
          }
        }
        p = VEC_DOUBLE().load_a(p_arr);
    #else
        p = dsfmt_genrand_open_close(*r);
        ++(*stream_length);
    #endif

    return p;
}


/*****************************************************************************/
/* Generate a 2D vector, perp_v, that is perpendicular to the unit vector    */
/* unit_v, for components specified by gen.                                  */
/*****************************************************************************/
#if defined(__SSE2__)
static inline void gen_2D_perp_vector(VEC_DOUBLE *perp_v, VEC_DOUBLE *unit_v,
                                      dsfmt_t **r, VEC_BOOL gen)
{
  /***************************************************************************/
  /* Initialise variables.                                                   */
  /***************************************************************************/
  double perp_v_arr_0[N_d];
  double perp_v_arr_1[N_d];
  perp_v[0].store(perp_v_arr_0);
  perp_v[1].store(perp_v_arr_1);

  for (int i = 0; i < N_d; ++i)
  {
    if (gen[i])
    {
      double p = dsfmt_genrand_open_close(r[i]);
      double vec_sign = 1.0;
      /***********************************************************************/
      /* Only two choices are possible - rotate the vector 90 degree either  */
      /* clockwise or anticlockwise.                                         */
      /***********************************************************************/
      if (p > 0.5)
      {
        vec_sign = -1.0;
      }
      perp_v_arr_0[i] = -vec_sign * unit_v[1][i];
      perp_v_arr_1[i] = vec_sign * unit_v[0][i];
    }
  }

  perp_v[0].load(perp_v_arr_0);
  perp_v[1].load(perp_v_arr_1);
}
#else
static inline void gen_2D_perp_vector(VEC_DOUBLE *perp_v, VEC_DOUBLE *unit_v,
                                      dsfmt_t **r, VEC_BOOL gen)
{
  /***************************************************************************/
  /* Initialise variables.                                                   */
  /***************************************************************************/
  double perp_v_arr_0;
  double perp_v_arr_1;

  if (gen)
    {
    double perp_mod = 0;
    double perp_dot = 0;

    /***********************************************************************/
    /* Construct a unit vector parallel to unit_v.                         */
    /***********************************************************************/
    do
    {
      /*********************************************************************/
      /* Generate a new random unit vector.                                */
      /*********************************************************************/
      double p = 2 * M_PI * dsfmt_genrand_open_close(*r);
      perp_v_arr_0 = cos(p);
      perp_v_arr_1 = sin(p);

      /*********************************************************************/
      /* Subtract the component parallel to unit_v and renormalise the     */
      /* vector.                                                           */
      /*********************************************************************/
      perp_dot = unit_v[0] * perp_v_arr_0;
      perp_dot += unit_v[1] * perp_v_arr_1;
      perp_v_arr_0 -= unit_v[0] * perp_dot;
      perp_v_arr_1 -= unit_v[1] * perp_dot;
      perp_mod += perp_v_arr_0 * perp_v_arr_0;
      perp_mod += perp_v_arr_1 * perp_v_arr_1;
      double inv_perp_mod = 1 / perp_mod;
      perp_v_arr_0 *= inv_perp_mod;
      perp_v_arr_1 *= inv_perp_mod;
    }
    while (abs(perp_dot) > 1e-12);
  }

  perp_v[0] = perp_v_arr_0;
  perp_v[1] = perp_v_arr_1;
}
#endif


/*****************************************************************************/
/* Generate a 2ND vector, perp_v, that is perpendicular to the unit vector   */
/* unit_v, for vectorisation components specified by gen.                    */
/*****************************************************************************/
#if defined(__SSE2__)
static inline void gen_2ND_perp_vector(VEC_DOUBLE *perp_v, VEC_DOUBLE *unit_v,
                                       dsfmt_t **r, VEC_BOOL gen)
{
  /***************************************************************************/
  /* Initialise variables.                                                   */
  /***************************************************************************/
  double perp_v_arr[Z_LEN * N_d];
  for (int n = 0; n < Z_LEN; n += 2)
  {
    perp_v[n].store(&perp_v_arr[n * N_d]);
    perp_v[n + 1].store(&perp_v_arr[(n + 1) * N_d]);
  }

  for (int i = 0; i < N_d; ++i)
  {
    if (gen[i])
    {
      /***********************************************************************/
      /* Construct a unit vector parallel to unit_v.                         */
      /***********************************************************************/
      double perp_dot;
      do
      {
        /*********************************************************************/
        /* Generate a new random unit vector.                                */
        /*********************************************************************/
        double norm = 0;

        for (int n = 0; n < Z_LEN; n += 2)
        {
          double phase = 2 * M_PI * dsfmt_genrand_open_close(r[i]);
          double mag = dsfmt_genrand_open_close(r[i]);
          perp_v_arr[N_d * n + i] = mag * cos(phase);
          perp_v_arr[N_d * (n + 1) + i] = mag * sin(phase);
          norm += mag * mag;
        }
        norm = 1/sqrt(norm);

        for (int n = 0; n < Z_LEN; n += 2)
        {
          perp_v_arr[N_d * n + i] *= norm;
          perp_v_arr[N_d * (n + 1) + i] *= norm;
        }

        /*********************************************************************/
        /* Subtract the component parallel to unit_v and renormalise the     */
        /* vector.                                                           */
        /*********************************************************************/
        perp_dot = 0;
        for (int n = 0; n < Z_LEN; n += 2)
        {
          perp_dot += unit_v[n][i] * perp_v_arr[N_d * n + i];
          perp_dot += unit_v[n + 1][i] * perp_v_arr[N_d * (n + 1) + i];
        }
        double perp_mod = 0;
        for (int n = 0; n < Z_LEN; n += 2)
        {
          perp_v_arr[N_d * n + i] -= unit_v[n][i] * perp_dot;
          perp_v_arr[N_d * (n + 1) + i] -= unit_v[n + 1][i] * perp_dot;
          perp_mod += perp_v_arr[N_d * n + i] * perp_v_arr[N_d * n + i];
          perp_mod += perp_v_arr[N_d * (n + 1) + i] * perp_v_arr[N_d * (n + 1) + i];
        }

        double inv_perp_mod = 1 / sqrt(perp_mod);
        for (int n = 0; n < Z_LEN; n += 2)
        {
          perp_v_arr[N_d * n + i] *= inv_perp_mod;
          perp_v_arr[N_d * (n + 1) + i] *= inv_perp_mod;
        }
      }
      while (abs(perp_dot) > 1e-12);
    }
  }

  for (int n = 0; n < Z_LEN; n += 2)
  {
    perp_v[n].load_a(&perp_v_arr[N_d * n]);
    perp_v[n + 1].load_a(&perp_v_arr[N_d * (n + 1)]);
  }
}
#else

#endif


static inline VEC_DOUBLE get_sign(VEC_DOUBLE *phi, VEC_DOUBLE *F_phi,
                                  VEC_DOUBLE cos_theta, VEC_DOUBLE theta_old)
{
  /***************************************************************************/
  /* Determine the sign to use for the perpendicular component.              */
  /***************************************************************************/
  VEC_DOUBLE phi_pos = cos_theta * F_phi[0] - sin(theta_old) * F_phi[1];
  VEC_DOUBLE phi_neg = cos_theta * F_phi[0] + sin(theta_old) * F_phi[1];

  VEC_BOOL sign_1 = (abs(phi_pos - phi[0]) < 3e-8);
  VEC_BOOL sign_2 = (abs(phi_neg - phi[0]) < 3e-8);

  double sign_arr[N_d];
  for (int i = 0; i < N_d; ++i)
  {
    if (sign_1[i])
    {
      if (sign_2[i])
      {
        int rn = rand() % 2;
        if (rn == 0)
        {
          sign_arr[i] = -1.0;
        }
      }
      else
      {
        sign_arr[i] = 1.0;
      }
    }
    else
    {
      sign_arr[i] = -1.0;
    }
  }
  VEC_DOUBLE sin_sign = VEC_DOUBLE().load_a(sign_arr);
  return sin_sign;
}


/*****************************************************************************/
/* Variable Update Function - z fields.                                      */
/*****************************************************************************/
static inline int update_z(VEC_DOUBLE * __restrict__ phi,
                           VEC_DOUBLE * __restrict__ F_phi,
                           dsfmt_t** r, long * stream_length)
{
  /***************************************************************************/
  /* Calculate the old angle between phi and F and the modulus of F.         */
  /***************************************************************************/
  VEC_DOUBLE F_phi_dot_prod(0);
  VEC_DOUBLE F_mod(0);

  for (int n = 0; n < Z_LEN; n += 2)
  {
    F_phi_dot_prod += phi[n] * F_phi[n];
    F_mod += F_phi[n] * F_phi[n];
    F_phi_dot_prod += phi[n + 1] * F_phi[n + 1];
    F_mod += F_phi[n + 1] * F_phi[n + 1];
  }

  F_mod = sqrt(F_mod);
  VEC_DOUBLE inv_F_mod = 1.0 / F_mod;
  VEC_DOUBLE cos_theta = F_phi_dot_prod * inv_F_mod;

  for (int n = 0; n < Z_LEN; n += 2)
  {
    F_phi[n] *= inv_F_mod;
    F_phi[n + 1] *= inv_F_mod;
  }

  /***************************************************************************/
  /* Catch any numerical errors that will return an invalid value for acos.  */
  /***************************************************************************/
  VEC_BOOL err = (cos_theta < -1.0);
  VEC_DOUBLE unit = 1.0;
  cos_theta = selectd(err, -unit, cos_theta);
  err = (cos_theta > 1.0);
  cos_theta = selectd(err, unit, cos_theta);
  VEC_DOUBLE theta_old = acos(cos_theta);

  /***************************************************************************/
  /* If theta is close to the extremes of 0 or pi, perform a heatbath update */
  /* to avoid introducing numerical errors.                                  */
  /***************************************************************************/
  VEC_BOOL smooth = ((theta_old < LOW_THETA_SKIP) ||
                                       ((M_PI - theta_old) <= LOW_THETA_SKIP));
#ifdef HEATBATH
  VEC_BOOL skip = true;
#else
  VEC_BOOL skip = (theta_old < 1e-16);
#endif

  /***************************************************************************/
  /* Make a vector perpendicular to F_phi using the sin of the angle between */
  /* F_phi and z. By default create the vector necessary to maximise the     */
  /* scalar product between the old value for phi and new one.               */
  /***************************************************************************/
  VEC_DOUBLE perp_vector[2 * N];
  for (int n = 0; n < Z_LEN; n += 2)
  {
    perp_vector[n] = (phi[n] - cos_theta * F_phi[n]) / sin(theta_old);
    perp_vector[n + 1] = (phi[n + 1] - cos_theta * F_phi[n + 1]) / sin(theta_old);
  }

  /***************************************************************************/
  /* If theta_old is zero within machine precision, generate a new random    */
  /* perpendicular vector.                                                   */
  /***************************************************************************/
  if (horizontal_or(skip))
  {
    gen_2ND_perp_vector(perp_vector, F_phi, r, skip);
  }

  /***************************************************************************/
  /* Smooth out numerical errors in the generate of this vector.             */
  /***************************************************************************/
  if (horizontal_or(smooth))
  {
    adjust_2ND_perp_vector(perp_vector, F_phi);
  }

  /***************************************************************************/
  /* Generate PDF parameters.                                                */
  /***************************************************************************/
  VEC_DOUBLE a = F_mod * N * BETA;
  VEC_DOUBLE z = (N - 1.0) / a;
  VEC_DOUBLE b = sqrt(1 + z*z);
  VEC_DOUBLE cos_x_0 = b - z;
  VEC_DOUBLE x_0 = acos(cos_x_0);
  VEC_DOUBLE c = sqrt(a * b);

  VEC_BOOL accepted(false);
  VEC_DOUBLE sin_theta_new;
  VEC_DOUBLE cos_theta_new;

  /*************************************************************************/
  /* Keep generating a new trial angle until one is accepted.              */
  /*************************************************************************/
  VEC_DOUBLE theta_new;
  int its = 0;
  double p_arr[N_d] __attribute__((aligned(N_d * sizeof(double))));
  double q_arr[N_d] __attribute__((aligned(N_d * sizeof(double))));

  do
  {
#if defined(__SSE2__)
    for (int i = 0; i < N_d; ++i)
    {
      if (!accepted[i])
      {
        p_arr[i] = dsfmt_genrand_open_close(r[i]);
        ++stream_length[i];
      }
    }
    VEC_DOUBLE p = VEC_DOUBLE().load_a(p_arr);
#else
    VEC_DOUBLE p = dsfmt_genrand_open_close(*r);
    ++(*stream_length);
#endif
    VEC_DOUBLE a1 = c * ((double)M_PI - x_0);
    VEC_DOUBLE a2 = p * atan(a1);
    VEC_DOUBLE a3 = (p - 1.0) * atan(c * x_0);
    VEC_DOUBLE tan_expr = tan(a2 + a3);
    theta_new = x_0 + tan_expr/c;
    sin_theta_new = sin(theta_new);
    cos_theta_new = cos(theta_new);
    VEC_DOUBLE sin_term = sin_theta_new / sin(x_0);
    VEC_DOUBLE exponent = a * (cos_theta_new - cos_x_0);
    VEC_DOUBLE sin_pow_term = pow(sin_term, POWER);
    VEC_DOUBLE exp_term = exp(exponent);
    VEC_DOUBLE density = 1.0 + pow(c * (theta_new - x_0), 2);
    VEC_DOUBLE ratio = sin_pow_term * exp_term;
    VEC_DOUBLE acc = ratio * density * ETA;
    
#if defined(__SSE2__)
    for (int i = 0; i < N_d; ++i)
    {
      if (!accepted[i])
      {
        q_arr[i] = dsfmt_genrand_open_close(r[i]);
        ++stream_length[i];
      }
    }
    VEC_DOUBLE q = VEC_DOUBLE().load_a(q_arr);
#else
    VEC_DOUBLE q = dsfmt_genrand_open_close(*r);
    ++(*stream_length);
#endif
    accepted = (acc >= q);
    ++its;
  }
  while (horizontal_count(accepted) != N_d);

  /***************************************************************************/
  /* Once the trial angle has been accepted, construct the new phi vector.   */
  /***************************************************************************/
  VEC_DOUBLE phi_mod(0.0);
  VEC_DOUBLE sin_term = sin(theta_new);

  for (int n = 0; n < Z_LEN; n += 2)
  {
    phi[n] = F_phi[n] * cos_theta_new - perp_vector[n] * sin_term;
    phi_mod += phi[n] * phi[n];
    phi[n + 1] = F_phi[n + 1] * cos_theta_new - perp_vector[n + 1] * sin_term;
    phi_mod += phi[n + 1] * phi[n + 1];
  }
  phi_mod = 1 / sqrt(phi_mod);

  /***************************************************************************/
  /* Normalise the new vector.                                               */
  /***************************************************************************/
  for (int n = 0; n < Z_LEN; n += 2)
  {
    phi[n] *= phi_mod;
    phi[n + 1] *= phi_mod;
  }

  return its;
} /* update_variable */


/*****************************************************************************/
/* Variable Update Function - U(1) lambda variable.                          */
/*****************************************************************************/
static inline int update_lambda(VEC_DOUBLE * __restrict__ phi,
                                VEC_DOUBLE * __restrict__ F_phi,
                                dsfmt_t** r, long* stream_length)
{
  /***************************************************************************/
  /* Calculate the old angle between phi and F and the modulus of F.         */
  /***************************************************************************/
  VEC_DOUBLE F_phi_dot_prod = phi[0] * F_phi[0];
  VEC_DOUBLE F_mod = F_phi[0] * F_phi[0];
  F_phi_dot_prod += phi[1] * F_phi[1];
  F_mod += F_phi[1] * F_phi[1];

  F_mod = sqrt(F_mod);
  VEC_DOUBLE inv_F_mod = 1.0 / F_mod;
  VEC_DOUBLE cos_theta = F_phi_dot_prod * inv_F_mod;
  F_phi[0] *= inv_F_mod;
  F_phi[1] *= inv_F_mod;

  /***************************************************************************/
  /* Catch any numerical errors that will return an invalid value for acos.  */
  /***************************************************************************/
  VEC_BOOL err = (cos_theta < -1.0);
  VEC_DOUBLE unit = 1.0;
  cos_theta = selectd(err, -unit, cos_theta);
  err = (cos_theta > 1.0);
  cos_theta = selectd(err, unit, cos_theta);
  VEC_DOUBLE theta_old = acos(cos_theta);

  /***************************************************************************/
  /* If theta is close to the extremes of 0 or pi, skip the update to avoid  */
  /* introducing numerical errors. Set theta_old to a safe value so that the */
  /* update procedure does not return a nan.                                 */
  /***************************************************************************/
  VEC_BOOL smooth = ((theta_old < LOW_THETA_SKIP) ||
                                       ((M_PI - theta_old) <= LOW_THETA_SKIP));
#ifdef HEATBATH
  VEC_BOOL skip = true;
#else
  VEC_BOOL skip = (theta_old < 1e-16);
#endif

  /***************************************************************************/
  /* Generate PDF parameters.                                                */
  /***************************************************************************/
  VEC_DOUBLE a = F_mod * N * BETA;
  VEC_DOUBLE sin_theta_new;
  VEC_DOUBLE cos_theta_new;

  /***************************************************************************/
  /* Use the standard algorithm to generate a new angle.                     */
  /***************************************************************************/
#ifdef SLOW_METHOD
  VEC_DOUBLE c = sqrt(a);

  VEC_BOOL accepted(false);
  double p_arr[N_d] __attribute__((aligned(N_d * sizeof(double))));
  double q_arr[N_d] __attribute__((aligned(N_d * sizeof(double))));
  VEC_DOUBLE eta_vec = 0.78125;
  VEC_DOUBLE eta_small(0);
  VEC_BOOL small = (a < 1.122);

  if (horizontal_or(small))
  {
    eta_small = exp(-2*a) * (1 + a * M_PI);
  }
  selectd(small, eta_vec, eta_small);

  /***************************************************************************/
  /* Keep generating a new trial angle until one is accepted.                */
  /***************************************************************************/
  VEC_DOUBLE theta_new;
  int its = 0;

  do
  {
#if defined(__SSE2__)
    for (int i = 0; i < N_d; ++i)
    {
      if (!accepted[i])
      {
        p_arr[i] = dsfmt_genrand_open_close(r[i]);
        ++stream_length[i];
      }
    }
    VEC_DOUBLE p = VEC_DOUBLE().load_a(p_arr);
#else
    VEC_DOUBLE p = dsfmt_genrand_open_close(*r);
    ++(*stream_length);
#endif
    VEC_DOUBLE a1 = c * ((double)M_PI);
    VEC_DOUBLE a2 = p * atan(a1);
    VEC_DOUBLE tan_expr = tan(a2);
    theta_new = tan_expr/c;
    sin_theta_new = sin(theta_new);
    cos_theta_new = cos(theta_new);
    VEC_DOUBLE exponent = a * (cos_theta_new - 1.0);
    VEC_DOUBLE exp_term = exp(exponent);
    VEC_DOUBLE density = 1.0 + a * pow(theta_new, 2);
    VEC_DOUBLE acc = exp_term * density * eta_vec;
#if defined(__SSE2__)
    for (int i = 0; i < N_d; ++i)
    {
      if (!accepted[i])
      {
        q_arr[i] = dsfmt_genrand_open_close(r[i]);
        ++stream_length[i];
      }
    }
    VEC_DOUBLE q = VEC_DOUBLE().load_a(q_arr);
#else
    VEC_DOUBLE q = dsfmt_genrand_open_close(*r);
    ++(*stream_length);
#endif
    accepted = (acc >= q);
    ++its;
  }
  while (horizontal_count(accepted) != N_d);

  /***************************************************************************/
  /* Use the optimised cosh method to generate a new angle.                  */
  /***************************************************************************/
#else
//static inline VEC_DOUBLE optimised_cosh(VEC_DOUBLE a, VEC_DOUBLE)

  /***************************************************************************/
  /* Initialise variables.                                                   */
  /***************************************************************************/
  double eps = 0.001;
  double as = 0.798953686083986;
  VEC_DOUBLE dap = max(0.0, a - as);
  VEC_DOUBLE del = 0.35 * dap + 1.03 * sqrt(dap);
  VEC_DOUBLE alp = min(sqrt(a * (2 - eps)), max(sqrt(eps * a), del));
  VEC_DOUBLE cosh_term = cosh(M_PI * alp);
  VEC_DOUBLE exp_term = exp(2 * a);
  VEC_DOUBLE bet = max(alp * alp / a, (cosh_term - 1)/(exp_term - 1)) - 1;
  VEC_DOUBLE bt1 = sqrt((1 + bet) / (1 - bet));

  VEC_BOOL acc(false); // accepted or not
  double p_arr[N_d] __attribute__((aligned(N_d * sizeof(double))));
  double q_arr[N_d] __attribute__((aligned(N_d * sizeof(double))));
  VEC_DOUBLE p, h1, theta_new, g, q;
  int its(0);

  /***************************************************************************/
  /* Keep generating a new trial angle until one is accepted.                */
  /***************************************************************************/
  do
  {

#if defined(__SSE2__)
    for (int i = 0; i < N_d; ++i)
    {
      if (!acc[i])
      {
        p_arr[i] = dsfmt_genrand_open_close(r[i]);
        ++stream_length[i];
      }
    }
    p = VEC_DOUBLE().load_a(p_arr);
#else
    p = dsfmt_genrand_open_close(*r);
    ++(*stream_length);
#endif
    h1 = bt1 * tan((2.0 * p - 1.0) * atan(tanh(HALF_PI * alp) / bt1));
    theta_new = log((1.0+h1)/(1.0-h1)) / alp;
    cos_theta_new = cos(theta_new);
    // decide if accept or reject
    g = (exp(-a * (1.0 - cos_theta_new)) * (cosh(alp * theta_new) + bet) /
                                                                (1.0 + bet));
#if defined(__SSE2__)
    for (int i = 0; i < N_d; ++i)
    {
      if (!acc[i])
      {
        q_arr[i] = dsfmt_genrand_open_close(r[i]);
        ++stream_length[i];
      }
    }
    q = VEC_DOUBLE().load_a(q_arr);
#else
    q = dsfmt_genrand_open_close(*r);
    ++(*stream_length);
#endif
    acc = (q < g);
    ++its;
  }
  while(horizontal_count(acc) != N_d);
#endif

  /***************************************************************************/
  /* Once the trial angle has been accepted, construct the new phi vector by */
  /* minimising the scalar product between the new and old phi.              */
  /***************************************************************************/
  VEC_DOUBLE phi_mod(0.0);
  VEC_DOUBLE zero = 0.0;

  VEC_DOUBLE perp_vector[2];
  perp_vector[0] = (phi[0] - cos_theta * F_phi[0]) / sin(theta_old);
  perp_vector[1] = (phi[1] - cos_theta * F_phi[1]) / sin(theta_old);

  /***************************************************************************/
  /* If theta is very small, then the perpendicular vector is not well       */
  /* determined, so we must make sure it is indeed perpendicular to F_phi    */
  /* to an acceptable precision.                                             */
  /***************************************************************************/
  if (horizontal_or(smooth))
  {
    /*************************************************************************/
    /* If this vector contains a value for theta that is zero to machine     */
    /* precision, generate a new one randomly.                               */
    /*************************************************************************/
    if (horizontal_or(skip))
    {
      /*std::cout << "Before:" << std::endl;

      for (int i = 0; i < N_d; ++i)
      {
        if (skip[i])
        {
          std::cout << "F[0] = " << F_phi[0][i] << ", F[1] = " << F_phi[1][i] << std::endl;
          std::cout << "v[0] = " << -F_phi[1][i] << ", v[1] = " << F_phi[0][i] << std::endl;
        }
      }*/
      gen_2D_perp_vector(perp_vector, F_phi, r, skip);
      /*for (int i = 0; i < N_d; ++i)
      {
        if (skip[i])
        {
          std::cout << "p[0] = " << perp_vector[0][i] << ", p[1] = " << perp_vector[1][i] << std::endl;
        }
      }*/
    }
    adjust_2D_perp_vector(perp_vector, F_phi);

    /*************************************************************************/
    /* If theta is not too small, the the existing vector will be slightly   */
    /* out alignment, and so we readjust it so that it is indeed             */
    /* perpendicular to F_phi.                                               */
    /*************************************************************************/
    //adjust_2D_perp_vector(perp_vector, F_phi);
  }

  /***************************************************************************/
  /* We may now safely construct the new phi vector.                         */
  /***************************************************************************/
  VEC_DOUBLE sin_sign = get_sign(phi, F_phi, cos_theta, theta_old);
  VEC_DOUBLE sin_term = sin(theta_new);
  //cos_theta_new = selectd(skip, cos_theta, cos_theta_new);
  //phi[0] = cos_theta_new * F_phi[0] - sin_sign * sin_term * F_phi[1];
  //phi[1] = cos_theta_new * F_phi[1] + sin_sign * sin_term * F_phi[0];

  /*VEC_DOUBLE phi_r1 = F_phi[0] * cos_theta_new - F_phi[1] * sin_term;
  VEC_DOUBLE phi_r2 = F_phi[0] * cos_theta_new + F_phi[1] * sin_term;
  VEC_DOUBLE phi_i1 = F_phi[1] * cos_theta_new + F_phi[0] * sin_term;
  VEC_DOUBLE phi_i2 = F_phi[1] * cos_theta_new - F_phi[0] * sin_term;*/
  phi[0] = (F_phi[0] * cos_theta_new - perp_vector[0] * sin_term);
  phi[1] = (F_phi[1] * cos_theta_new - perp_vector[1] * sin_term);
  phi_mod += phi[0] * phi[0];
  phi_mod += phi[1] * phi[1];
  phi_mod = 1/sqrt(phi_mod);

  /***************************************************************************/
  /* Normalise the new vector.                                               */
  /***************************************************************************/
  phi[0] *= phi_mod;
  phi[1] *= phi_mod;

  //for (int i = 0; i < N_d; ++i)
  //{
  //  std::cout << "old: " << std::setprecision(16) << sin_sign[i] * theta_old[i] << "  new : " << std::setprecision(16) << theta_new[i] << std::endl;
  //}
  return its;
} /* update_lambda */


static inline unsigned long update_func(VEC_DOUBLE ***z, VEC_DOUBLE *F_phi,
                                        dsfmt_t**** r,
                                        unsigned int t, unsigned int x,
                                        unsigned int td, unsigned int xd,
                                        unsigned int tu, unsigned int xu,
                                        int shift,
                                        long* z_stream_counter,
                                        long* lam_stream_counter)
{
  unsigned long flops = 0;
  int its;
  VEC_DOUBLE *phi_z;
  VEC_DOUBLE *phi_z_td;
  VEC_DOUBLE *phi_z_xd;
  VEC_DOUBLE *phi_z_tu;
  VEC_DOUBLE *phi_z_xu;
  VEC_DOUBLE *phi_l_td;
  VEC_DOUBLE *phi_l_xd;
  VEC_DOUBLE *phi_l;

  /***************************************************************************/
  /* Construct the phi vectors for this and neighbouring sites. If we do     */
  /* not want to impose periodic boundary conditions, replace links on the   */
  /* boundary with zero vectors.                                             */
  /***************************************************************************/
  phi_z = z[t][x];
  phi_l = &z[t][x][Z_LEN];
  phi_z_td = z[td][x];
  phi_l_td = &z[td][x][Z_LEN];
  phi_l_xd = &z[t][xd][Z_LEN + 2];
  phi_z_tu = z[tu][x];
  phi_z_xd = z[t][xd];
  phi_l_xd = &z[t][xd][Z_LEN + 2];
  phi_z_xu = z[t][xu];

#ifdef NO_BC
  if (shift == 0)
  {
    if (t == 0)
    {
      F_phi_z_3pt_0(phi_z_tu, phi_z_xu, phi_z_xd, phi_l, phi_l,
                    phi_l_xd, phi_z_xu + Z_LEN + 2, F_phi);
    }
    else if (tu == 0)
    {
      F_phi_z_3pt_0(phi_z_td, phi_z_xu, phi_z_xd, phi_l, phi_l_td,
                    phi_l_xd, phi_z_xu + Z_LEN + 2, F_phi);
    }
    else
    {
      F_phi_z_0(phi_z_tu, phi_z_td, phi_z_xu, phi_z_xd, phi_l, phi_l_td,
                phi_l_xd, phi_z_xu + Z_LEN + 2, F_phi);
    }
  }
  else
  {
    if (t == 0)
    {
      F_phi_z_3pt_1(phi_z_tu, phi_z_xu, phi_z_xd, phi_l, phi_l,
                    phi_l_xd, F_phi);
    }
    else if (tu == 0)
    {
      F_phi_z_3pt_1(phi_z_td, phi_z_xu, phi_z_xd, phi_l, phi_l_td,
                    phi_l_xd, F_phi);
    }
    else
    {
      F_phi_z_1(phi_z_tu, phi_z_td, phi_z_xu, phi_z_xd, phi_l, phi_l_td,
                phi_l_xd, F_phi);
    }
  }
#else
  if (shift == 0)
  {
    F_phi_z_0(phi_z_tu, phi_z_td, phi_z_xu, phi_z_xd, phi_l, phi_l_td,
              phi_l_xd, phi_z_xu + Z_LEN + 2, F_phi);
  }
  else
  {
    F_phi_z_1(phi_z_tu, phi_z_td, phi_z_xu, phi_z_xd, phi_l, phi_l_td,
              phi_l_xd, F_phi);
  }
#endif

  /***************************************************************************/
  /* Construct the force term to determine the local action, then            */
  /* update the variable corresponding to this site.                         */
  /***************************************************************************/
  its = update_z(phi_z, F_phi, r[t][x], z_stream_counter);
  flops += (58 * N) + (307 * its) + 138;

  /***************************************************************************/
  /* Update the gauge link in the t direction if we are on the edge of the   */
  /* lattice and have no boundary conditions, do not update the link as it   */
  /* is fixed at zero.                                                       */
  /***************************************************************************/
#ifdef NO_BC
  if (t != 0)
  {
    F_phi_lam(phi_z, phi_z_tu, F_phi);
    its = update_lambda(phi_l, F_phi, r[t][x], lam_stream_counter);
    flops += (8 * N) + (306 * its) + 193;
  }
#else
  F_phi_lam(phi_z, phi_z_tu, F_phi);
  its = update_lambda(phi_l, F_phi, r[t][x], lam_stream_counter);
  flops += (8 * N) + (306 * its) + 193;
#endif

  /***************************************************************************/
  /* Repeat for the gauge link in the x direction.                           */
  /***************************************************************************/
  if (shift == 0)
  {
    F_phi_lam(phi_z, phi_z_xu, F_phi);
  }
  else
  {
    F_phi_lam_1(phi_z, phi_z_xd, phi_z_xu, F_phi);
  }
  its = update_lambda(phi_l + 2, F_phi, r[t][x], lam_stream_counter);
  flops += (8 * N) + (306 * its) + 193;

  return flops * N_d;
}


static inline unsigned long update_loop(VEC_DOUBLE ***z, VEC_DOUBLE ***lam,
                            dsfmt_t**** r,  unsigned int min_t,
                            unsigned int max_t, unsigned int min_x,
                            unsigned int max_x, int step, int shift,
                            long* z_stream_counter, long* lam_stream_counter)
{
  /***************************************************************************/
  /* Initialize variables.                                                   */
  /***************************************************************************/
  unsigned long flops = 0;
  unsigned int ind;
  unsigned int x;
  int start = min_x;
  int end = NUM_VEC_BLOCKS - 1;
  bool skip_last(false);
  VEC_DOUBLE F_phi[2 * N];

  if (shift == 1)
  {
    ++start;
  }

  unsigned int td = min_t - 1;
  unsigned int tu = min_t + 1;
  long index = (min_t * NUM_VEC_BLOCKS + start) * N_d;

  for (unsigned int t = min_t; t < max_t; ++t)
  {
    /*************************************************************************/
    /* If this is the last t iteration, set the last x iteration.            */
    /*************************************************************************/
    if (t == (max_t - 1))
    {
      end = max_x;
      if (max_x < NUM_VEC_BLOCKS - 1)
      {
        skip_last = true;
      }
      else if (max_x > NUM_VEC_BLOCKS - 1)
      {
        end = NUM_VEC_BLOCKS - 1;
      }
    }
    if (t == 0)
    {
      td = L_T - 1;
    }
    else if (t == 1)
    {
      td = 0;
    }
    else if (t == (L_T - 1))
    {
      tu = 0;
    }

    /*************************************************************************/
    /* Update the first entry on this row.                                   */
    /*************************************************************************/
    if (start == 0)
    {
      flops += update_func(z, F_phi, r, t, 0, td, NUM_VEC_BLOCKS - 1,
                           tu, 1, shift,
                           &z_stream_counter[index],
                           &lam_stream_counter[index]);
      start += step;
      index += step * N_d;
    }

    /*************************************************************************/
    /* Update all remaining entries on this row except the last.             */
    /*************************************************************************/
    for (x = start; x < end; x += step)
    {
      flops += update_func(z, F_phi, r, t, x, td, x - 1, tu, x + 1, shift,
                           &z_stream_counter[index],
                           &lam_stream_counter[index]);
      index += step * N_d;
    }

    /*************************************************************************/
    /* Update the last entry on this row if specified.                       */
    /*************************************************************************/
    if (!skip_last)
    {
      if (shift == 1)
      {
        flops += update_func(z, F_phi, r, t, x, td, x - 1, tu, 0, shift,
                             &z_stream_counter[index],
                             &lam_stream_counter[index]);
      }

      /***********************************************************************/
      /* Set variables for the next loop iteration.                          */
      /***********************************************************************/
      index += N_d;
      tu++;
      td++;
      shift = (shift + 1) % 2;
      start = shift;
    }
  }

  return flops;
} /* update_loop */

static inline unsigned long heatbath_update(VEC_DOUBLE ***z, VEC_DOUBLE ***lam,
                              dsfmt_t**** r, update_limits* limits,
                              long * z_stream_counter,
                              long * lam_stream_counter)
{
  unsigned long flops[OMP_NUM_THREADS] = {0};
  unsigned long total_flops = 0;

  /***************************************************************************/
  /* Loop over the whole lattice and update the lattice variables at each    */
  /* point sequentially. Update in a checkerboard pattern: update all the    */
  /* "black" points in parallel, followed by all of the "red" points.        */
  /***************************************************************************/
  #pragma omp parallel num_threads(OMP_NUM_THREADS)
  {
    /*************************************************************************/
    /* Update black.                                                         */
    /*************************************************************************/
    #pragma omp for schedule(static)
    for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
    {
      int shift = limits[ID].u_min_t % 2;
      int step = 2;
      flops[ID] += update_loop(z, lam, r, limits[ID].u_min_t,
                               limits[ID].u_max_t, limits[ID].u_min_x,
                               limits[ID].u_max_x, step, shift,
                               z_stream_counter, lam_stream_counter);
    }

    /*************************************************************************/
    /* Update red.                                                           */
    /*************************************************************************/
    #pragma omp for schedule(static)
    for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
    {
      int shift = (limits[ID].u_min_t + 1) % 2;
      int step = 2;
      flops[ID] += update_loop(z, lam, r, limits[ID].u_min_t,
                               limits[ID].u_max_t, limits[ID].u_min_x,
                               limits[ID].u_max_x, step, shift,
                               z_stream_counter, lam_stream_counter);
    }
  }

  for (int i = 0; i < OMP_NUM_THREADS; ++i)
  {
    total_flops += flops[i];
  }
  return total_flops;
} /* heatbath_update */


/*****************************************************************************/
/* Action.                                                                   */
/*****************************************************************************/
static inline double S_1(VEC_DOUBLE ***z, VEC_DOUBLE ***lam)
{
  //VEC_DOUBLE S[OMP_NUM_THREADS];
  double sum(0.0);

  #pragma omp parallel num_threads(OMP_NUM_THREADS / SMT) reduction(+:sum)
  {
    #pragma omp for
    for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
    {
      unsigned int min_t = ID * L_STEP;
      unsigned int max_t = (ID + 1) * L_STEP;
      // HACKY TEST, L -> L_T (3 instances)
      int shift = 0;
      if (max_t > L)
      {
        max_t = L;
      }

      VEC_DOUBLE acc1 = 0;
      VEC_DOUBLE acc2 = 0;
      VEC_DOUBLE acc3 = 0;
      VEC_DOUBLE acc4 = 0;

      for (unsigned int t = min_t; t < max_t; ++t)
      {
        // tu will wrap around only if feducial box is whole volume.
        unsigned int tu = (t + 1) % L_T;
        for (unsigned int x = 0; x < NUM_VEC_BLOCKS; x += 2)
        {
          /*******************************************************************/
          /* Construct array indices and phi vectors.                        */
          /*******************************************************************/
          unsigned int xu = x + 1;
          unsigned int xu_2 = (xu + 1) % NUM_VEC_BLOCKS;
          VEC_DOUBLE * phi_z = z[t][x];
          VEC_DOUBLE * phi_z_tu = z[tu][x];
          VEC_DOUBLE * phi_z_xu = z[t][xu];
          VEC_DOUBLE * lambda = &z[t][x][Z_LEN];
          VEC_DOUBLE * phi_z_tu_2 = z[tu][xu];
          VEC_DOUBLE * phi_z_xu_2 = z[t][xu_2];
          VEC_DOUBLE * lambda_2 = &z[t][xu][Z_LEN];
          VEC_DOUBLE prod_1 = 0;
          VEC_DOUBLE prod_2 = 0;
          VEC_DOUBLE prod_3 = 0;
          VEC_DOUBLE prod_4 = 0;
          VEC_DOUBLE prod_5 = 0;
          VEC_DOUBLE prod_6 = 0;
          VEC_DOUBLE prod_7 = 0;
          VEC_DOUBLE prod_8 = 0;

          /*********************************************************************/
          /* For each given link, compute the action at the specified lattice  */
          /* point.                                                            */
          /*********************************************************************/
          for (int n = 0; n < Z_LEN; n += 2)
          {
            prod_1 += (phi_z_tu[n] * phi_z[n] + phi_z_tu[n + 1] * phi_z[n + 1]);
            prod_2 += (phi_z_tu[n + 1] * phi_z[n] - phi_z_tu[n] * phi_z[n + 1]);
            prod_3 += (phi_z_xu[n] * phi_z[n] + phi_z_xu[n + 1] * phi_z[n + 1]);
            prod_4 += (phi_z_xu[n + 1] * phi_z[n] - phi_z_xu[n] * phi_z[n + 1]);
            VEC_DOUBLE phi_z_xu_r = load_phi_xu(&phi_z[n], &phi_z_xu_2[n]);
            VEC_DOUBLE phi_z_xu_i = load_phi_xu(&phi_z[n + 1], &phi_z_xu_2[n + 1]);
            prod_5 += (phi_z_tu_2[n] * phi_z_xu[n] + phi_z_tu_2[n + 1] * phi_z_xu[n + 1]);
            prod_6 += (phi_z_tu_2[n + 1] * phi_z_xu[n] - phi_z_tu_2[n] * phi_z_xu[n + 1]);
            prod_7 += (phi_z_xu_r * phi_z_xu[n] + phi_z_xu_i * phi_z_xu[n + 1]);
            prod_8 += (phi_z_xu_i * phi_z_xu[n] - phi_z_xu_r * phi_z_xu[n + 1]);
          }
          acc1 += prod_1 * lambda[0];
          acc2 += prod_2 * lambda[1];
          acc3 += prod_3 * lambda[2];
          acc4 += prod_4 * lambda[3];
          acc1 += prod_5 * lambda_2[0];
          acc2 += prod_6 * lambda_2[1];
          acc3 += prod_7 * lambda_2[2];
          acc4 += prod_8 * lambda_2[3];
          acc1 -= 4;
        }
      }
      sum = horizontal_add(acc1 + acc2 + acc3 + acc4);
    }
  }
  // HACKY TEST - V -> L * L_T
  return -(sum / V);
} /*   Action   */


/* Average link - force angle */
static inline void link_angle(VEC_DOUBLE ***z, VEC_DOUBLE ***lam,
                              VEC_DOUBLE *angles)
{
  //VEC_DOUBLE S[OMP_NUM_THREADS];
  //double sum(0.0);

  #pragma omp parallel num_threads(OMP_NUM_THREADS / SMT) //reduction(+:sum)
  {
    #pragma omp for
    for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
    {
      VEC_DOUBLE F_phi_1[2];
      VEC_DOUBLE F_phi_2[2];
      unsigned int min_t = ID * L_STEP;
      unsigned int max_t = (ID + 1) * L_STEP;
      // HACKY TEST, L -> L_T (3 instances)
      int shift = 0;
      if (max_t > L)
      {
        max_t = L;
      }

      VEC_DOUBLE acc1 = 0;
      VEC_DOUBLE acc2 = 0;

      for (unsigned int t = min_t; t < max_t; ++t)
      {
        // tu will wrap around only if the fiducial volume is the whole one.
        unsigned int tu = (t + 1) % L_T;
        for (unsigned int x = 0; x < NUM_VEC_BLOCKS; x += 2)
        {
          /*******************************************************************/
          /* Construct array indices and phi vectors.                        */
          /*******************************************************************/
          unsigned int xu = x + 1;
          unsigned int xu_2 = (xu + 1) % NUM_VEC_BLOCKS;
          VEC_DOUBLE * phi_z = z[t][x];
          VEC_DOUBLE * phi_z_tu = z[tu][x];
          VEC_DOUBLE * phi_z_xu = z[t][xu];
          VEC_DOUBLE * lambda = &z[t][x][Z_LEN];
          VEC_DOUBLE * phi_z_tu_2 = z[tu][xu];
          VEC_DOUBLE * phi_z_xu_2 = z[t][xu_2];
          VEC_DOUBLE * lambda_2 = &z[t][xu][Z_LEN];

          long V_ind = 2 * (NUM_VEC_BLOCKS * t + x);

          /*******************************************************************/
          /* For each given link, compute the angle between lambda and its   */
          /* associated force vector.                                        */
          /*******************************************************************/
          // t direction, red & black
          F_phi_lam(phi_z, phi_z_tu, F_phi_1);
          F_phi_lam(phi_z_xu, phi_z_tu_2, F_phi_2);

          VEC_DOUBLE F_phi_dot_prod = lambda[0] * F_phi_1[0];
          VEC_DOUBLE F_mod = F_phi_1[0] * F_phi_1[0];
          F_phi_dot_prod += lambda[1] * F_phi_1[1];
          F_mod += F_phi_1[1] * F_phi_1[1];
          F_mod = sqrt(F_mod);
          F_phi_1[0] /= F_mod;
          F_phi_1[1] /= F_mod;

          VEC_DOUBLE cos_theta = F_phi_dot_prod / F_mod;
          VEC_BOOL err = (cos_theta < -1.0);
          VEC_DOUBLE unit = 1.0;
          cos_theta = selectd(err, -unit, cos_theta);
          err = (cos_theta > 1.0);
          cos_theta = selectd(err, unit, cos_theta);
          VEC_DOUBLE theta_old = acos(cos_theta);
          //angles[V_ind] = theta_old * get_sign(lambda, F_phi_1, cos_theta, theta_old);
          angles[V_ind] = F_mod * F_phi_1[0];
          //acc1 += get_sign(lambda, F_phi_1, cos_theta, theta_old);

          F_phi_dot_prod = lambda_2[0] * F_phi_2[0];
          F_mod = F_phi_2[0] * F_phi_2[0];
          F_phi_dot_prod += lambda_2[1] * F_phi_2[1];
          F_mod += F_phi_2[1] * F_phi_2[1];
          F_mod = sqrt(F_mod);
          F_phi_2[0] /= F_mod;
          F_phi_2[1] /= F_mod;
          cos_theta = F_phi_dot_prod / F_mod;
          err = (cos_theta < -1.0);
          cos_theta = selectd(err, -unit, cos_theta);
          err = (cos_theta > 1.0);
          cos_theta = selectd(err, unit, cos_theta);
          theta_old = acos(cos_theta);
          //angles[V_ind + 2] = theta_old * get_sign(lambda_2, F_phi_2, cos_theta, theta_old);
          angles[V_ind + 2] = F_mod * F_phi_2[0];
          //acc2 += get_sign(lambda_2, F_phi_2, cos_theta, theta_old);

          // x direction, red & black
          F_phi_lam(phi_z, phi_z_xu, F_phi_1);
          F_phi_lam_1(phi_z_xu, phi_z, phi_z_xu_2, F_phi_2);

          F_phi_dot_prod = lambda[2] * F_phi_1[0];
          F_mod = F_phi_1[0] * F_phi_1[0];
          F_phi_dot_prod += lambda[3] * F_phi_1[1];
          F_mod += F_phi_1[1] * F_phi_1[1];
          F_mod = sqrt(F_mod);
          F_phi_1[0] /= F_mod;
          F_phi_1[1] /= F_mod;
          cos_theta = F_phi_dot_prod / F_mod;
          err = (cos_theta < -1.0);
          cos_theta = selectd(err, -unit, cos_theta);
          err = (cos_theta > 1.0);
          cos_theta = selectd(err, unit, cos_theta);
          theta_old = acos(cos_theta);
          //angles[V_ind + 1] = theta_old * get_sign(&lambda[2], F_phi_1, cos_theta, theta_old);
          angles[V_ind + 1] = F_mod * F_phi_1[0];
          //acc2 += get_sign(&lambda[2], F_phi_1, cos_theta, theta_old);

          F_phi_dot_prod = lambda_2[2] * F_phi_2[0];
          F_mod = F_phi_2[0] * F_phi_2[0];
          F_phi_dot_prod += lambda_2[3] * F_phi_2[1];
          F_mod += F_phi_2[1] * F_phi_2[1];
          F_mod = sqrt(F_mod);
          F_phi_2[0] /= F_mod;
          F_phi_2[1] /= F_mod;
          cos_theta = F_phi_dot_prod / F_mod;
          err = (cos_theta < -1.0);
          cos_theta = selectd(err, -unit, cos_theta);
          err = (cos_theta > 1.0);
          cos_theta = selectd(err, unit, cos_theta);
          theta_old = acos(cos_theta);
          //angles[V_ind + 3] = theta_old * get_sign(&lambda_2[2], F_phi_2, cos_theta, theta_old);
          angles[V_ind + 3] = F_mod * F_phi_2[0];
          //acc2 += get_sign(&lambda_2[2], F_phi_2, cos_theta, theta_old);

          /*if (V_ind != 0)
          {
            angles[V_ind] = angles[0];
          }
          angles[V_ind + 1] = angles[0];
          angles[V_ind + 2] = angles[0];
          angles[V_ind + 3] = angles[0];*/
        }
      }
      //sum = horizontal_add(acc1 + acc2);
    }
  }
  // HACKY TEST - V -> L * L_T
//  return sum / (2 * V);
  return;
} /*   Link angle   */


/*****************************************************************************/
/* Compuation loop for G_k.                                                  */
/*****************************************************************************/
static inline void G_k_loop(VEC_DOUBLE ***z, double * __restrict__ P_0,
                            double * __restrict__ P_1,
                            update_limits* __restrict__ limits)
{
  /***************************************************************************/
  /* Initialise variables.                                                   */
  /***************************************************************************/
  unsigned int start = limits->m_min_x;
  unsigned int end = NUM_VEC_BLOCKS;

  for (unsigned int t = limits->m_min_t; t < limits->m_max_t; ++t)
  {
    double cos_term = cos(t * K_UNIT);
    double sin_term = sin(t * K_UNIT);

    if ((t + 1) == limits->m_max_t)
    {
      end = limits->m_max_x;
    }

    for (unsigned int x = start; x < end; x += 2)
    {
      unsigned int ii = 0;
      unsigned int jj = 0;
      VEC_DOUBLE real_part_0;
      VEC_DOUBLE imag_part_0;
      VEC_DOUBLE real_part_1;
      VEC_DOUBLE imag_part_1;
      VEC_DOUBLE real_part_2;
      VEC_DOUBLE imag_part_2;
      VEC_DOUBLE real_part_3;
      VEC_DOUBLE imag_part_3;

      for (long kk = 0; kk < P_LEN; kk += 2)
      {
        real_part_0 = z[t][x][ii] * z[t][x][jj];
        real_part_1 = z[t][x][ii + 1] * z[t][x][jj + 1];
        imag_part_0 = z[t][x][ii] * z[t][x][jj + 1];
        imag_part_1 = z[t][x][ii + 1] * z[t][x][jj];
        real_part_2 = z[t][x + 1][ii] * z[t][x + 1][jj];
        real_part_3 = z[t][x + 1][ii + 1] * z[t][x + 1][jj + 1];
        imag_part_2 = z[t][x + 1][ii] * z[t][x + 1][jj + 1];
        imag_part_3 = z[t][x + 1][ii + 1] * z[t][x + 1][jj];
        double sum_r_1 = horizontal_add(real_part_0 + real_part_1);
        double sum_i_1 = horizontal_add(imag_part_0 - imag_part_1);
        double sum_r_2 = horizontal_add(real_part_2 + real_part_3);
        double sum_i_2 = horizontal_add(imag_part_2 - imag_part_3);
        double temp1 = sum_r_1 + sum_r_2;
        double temp2 = sum_i_1 + sum_i_2;

        P_0[kk] += temp1;
        P_0[kk + 1] += temp2;
        P_1[kk] += temp1 * cos_term - temp2 * sin_term;
        P_1[kk + 1] += temp1 * sin_term + temp2 * cos_term;

        jj += 2;
        if (jj == (Z_LEN))
        {
          ii += 2;
          jj = 0;
        }
      }
    }
    start = 0;
  }

  return;
}

static inline void G_k_loop_alt(VEC_DOUBLE ***z, double * __restrict__ P_0,
                                double * __restrict__ P_1,
                                int ID, update_limits* limits)
{
  /***************************************************************************/
  /* Initialise variables.                                                   */
  /***************************************************************************/
  unsigned int start = limits->m_min_x;
  unsigned int end = NUM_VEC_BLOCKS;

  for (unsigned int t = limits->m_min_t; t < limits->m_max_t; ++t)
  {
    if ((t + 1) == limits->m_max_t)
    {
      end = limits->m_max_x;
    }

    for (unsigned int x = start; x < end; x += 2)
    {
      unsigned int ii = 0;
      unsigned int jj = 0;
      VEC_DOUBLE real_part_0;
      VEC_DOUBLE imag_part_0;
      VEC_DOUBLE real_part_1;
      VEC_DOUBLE imag_part_1;
      VEC_DOUBLE real_part_2;
      VEC_DOUBLE imag_part_2;
      VEC_DOUBLE real_part_3;
      VEC_DOUBLE imag_part_3;

      VEC_DOUBLE ind_vec = x_index_vector(x, 0);
      VEC_DOUBLE cos_term_1 = cos(K_UNIT * ind_vec);
      VEC_DOUBLE sin_term_1 = sin(K_UNIT * ind_vec);
      ind_vec += 1;
      VEC_DOUBLE cos_term_2 = cos(K_UNIT * ind_vec);
      VEC_DOUBLE sin_term_2 = sin(K_UNIT * ind_vec);

      for (long kk = 0; kk < P_LEN; kk += 2)
      {
        real_part_0 = z[t][x][ii] * z[t][x][jj];
        real_part_1 = z[t][x][ii + 1] * z[t][x][jj + 1];
        imag_part_0 = z[t][x][ii] * z[t][x][jj + 1];
        imag_part_1 = z[t][x][ii + 1] * z[t][x][jj];
        real_part_2 = z[t][x + 1][ii] * z[t][x + 1][jj];
        real_part_3 = z[t][x + 1][ii + 1] * z[t][x + 1][jj + 1];
        imag_part_2 = z[t][x + 1][ii] * z[t][x + 1][jj + 1];
        imag_part_3 = z[t][x + 1][ii + 1] * z[t][x + 1][jj];
        VEC_DOUBLE sum_r_1 = real_part_0 + real_part_1;
        VEC_DOUBLE sum_i_1 = imag_part_0 - imag_part_1;
        VEC_DOUBLE sum_r_2 = real_part_2 + real_part_3;
        VEC_DOUBLE sum_i_2 = imag_part_2 - imag_part_3;
        P_0[kk] += horizontal_add(sum_r_1 + sum_r_2);
        P_0[kk + 1] += horizontal_add(sum_i_1 + sum_i_2);

        VEC_DOUBLE real = sum_r_1 * cos_term_1 + sum_r_2 * cos_term_2;
        real -= sum_i_1 * sin_term_1 + sum_i_2 * sin_term_2;
        VEC_DOUBLE imag = sum_r_1 * sin_term_1 + sum_r_2 * sin_term_2;
        imag += sum_i_1 * cos_term_1 + sum_i_2 * cos_term_2;
        P_1[kk] += horizontal_add(real);
        P_1[kk + 1] += horizontal_add(imag);

        jj += 2;
        if (jj == (Z_LEN))
        {
          ii += 2;
          jj = 0;
        }
      }
    }
    start = 0;
  }

  return;
}


#ifdef NO_BC
/*****************************************************************************/
/* Compuation loop for G_k.                                                  */
/*****************************************************************************/
static inline void G_k_loop_no_bc(VEC_DOUBLE ***z, double * __restrict__ P_0,
                                double * __restrict__ P_1,
                                update_limits* __restrict__ limits)
{
  unsigned int start = limits->m_min_x;
  unsigned int end = NUM_VEC_BLOCKS;
  for (unsigned int t = limits->m_min_t; t < limits->m_max_t; ++t)
  {
    if ((t + 1) == limits->m_max_t)
    {
      end = limits->m_max_x;
    }
    for (unsigned int x = start; x < end; x += 2)
    {
      unsigned int ii = 0;
      unsigned int jj = 0;
      VEC_DOUBLE real_part_0;
      VEC_DOUBLE imag_part_0;
      VEC_DOUBLE real_part_1;
      VEC_DOUBLE imag_part_1;
      VEC_DOUBLE real_part_2;
      VEC_DOUBLE imag_part_2;
      VEC_DOUBLE real_part_3;
      VEC_DOUBLE imag_part_3;

      VEC_DOUBLE ind_vec = x_index_vector(x, 0);
      VEC_DOUBLE cos_term_1 = cos(K_UNIT * ind_vec);
      VEC_DOUBLE sin_term_1 = sin(K_UNIT * ind_vec);
      VEC_DOUBLE cos_term_1a = cos(K_UNIT * (ind_vec + t));
      VEC_DOUBLE sin_term_1a = sin(K_UNIT * (ind_vec + t));
      ind_vec += 1;
      VEC_DOUBLE cos_term_2 = cos(K_UNIT * ind_vec);
      VEC_DOUBLE sin_term_2 = sin(K_UNIT * ind_vec);
      VEC_DOUBLE cos_term_2a = cos(K_UNIT * (ind_vec + t));
      VEC_DOUBLE sin_term_2a = sin(K_UNIT * (ind_vec + t));

      for (long kk = 0; kk < P_LEN; kk += 2)
      {
        real_part_0 = z[t][x][ii] * z[t][x][jj];
        real_part_1 = z[t][x][ii + 1] * z[t][x][jj + 1];
        imag_part_0 = z[t][x][ii] * z[t][x][jj + 1];
        imag_part_1 = z[t][x][ii + 1] * z[t][x][jj];
        real_part_2 = z[t][x + 1][ii] * z[t][x + 1][jj];
        real_part_3 = z[t][x + 1][ii + 1] * z[t][x + 1][jj + 1];
        imag_part_2 = z[t][x + 1][ii] * z[t][x + 1][jj + 1];
        imag_part_3 = z[t][x + 1][ii + 1] * z[t][x + 1][jj];
        VEC_DOUBLE sum_r_1 = real_part_0 + real_part_1;
        VEC_DOUBLE sum_i_1 = imag_part_0 - imag_part_1;
        VEC_DOUBLE sum_r_2 = real_part_2 + real_part_3;
        VEC_DOUBLE sum_i_2 = imag_part_2 - imag_part_3;

        VEC_DOUBLE real = sum_r_1 * cos_term_1 + sum_r_2 * cos_term_2;
        real -= sum_i_1 * sin_term_1 + sum_i_2 * sin_term_2;
        VEC_DOUBLE imag = sum_r_1 * sin_term_1 + sum_r_2 * sin_term_2;
        imag += sum_i_1 * cos_term_1 + sum_i_2 * cos_term_2;
        P_0[kk] += horizontal_add(real);
        P_0[kk + 1] += horizontal_add(imag);

        real = sum_r_1 * cos_term_1a + sum_r_2 * cos_term_2a;
        real -= sum_i_1 * sin_term_1a + sum_i_2 * sin_term_2a;
        imag = sum_r_1 * sin_term_1a + sum_r_2 * sin_term_2a;
        imag += sum_i_1 * cos_term_1a + sum_i_2 * cos_term_2a;
        P_1[kk] += horizontal_add(real);
        P_1[kk + 1] += horizontal_add(imag);

        jj += 2;
        if (jj == (Z_LEN))
        {
          ii += 2;
          jj = 0;
        }
      }
    }
    start = 0;
  }

  return;
}
#endif


/*****************************************************************************/
/* Determine the simple wall-wall correlation function.                      */
/* Computed by making the object P(t) to perform convolution over t.         */
/*****************************************************************************/
static inline void wall_corr(VEC_DOUBLE ***z,
                             double** __restrict__ P_t_avg,
                             double* __restrict__ correlator)
{
  /***************************************************************************/
  /* Construct the object P_t_avg; that is sum the composite operator P(t,x) */
  /* over the t direction for each value of x.                               */
  /***************************************************************************/
  #pragma omp parallel for num_threads(OMP_NUM_THREADS / SMT) schedule(static)
  for (unsigned int t = 0; t < L; ++t)
  {
    /*************************************************************************/
    /* Assign temporary vector memory for this average P operator.           */
    /*************************************************************************/
    VEC_DOUBLE P_op[P_LEN];
    for (int kk = 0; kk < P_LEN; kk += 2)
    {
      P_op[kk] = 0;
      P_op[kk + 1] = 0;
    }

    /*************************************************************************/
    /* Average the composite operator over this dimension.                   */
    /*************************************************************************/
    for (unsigned int x = 0; x < NUM_VEC_BLOCKS; ++x)
    {
      /***********************************************************************/
      /* Add the tensor product of the z vector at this site with its        */
      /* complex conjugate to the corresponding P_t_avg entry.               */
      /***********************************************************************/
      int ii = 0;
      int jj = 0;

      for (long kk = 0; kk < P_LEN; kk += 2)
      {
        /*********************************************************************/
        /* The matrices in P_t_avg are symmetric, so compute only the upper  */
        /* triangular entries, and the rest can be recovered from the        */
        /* symmetry.                                                         */
        /*********************************************************************/
        if (jj >= ii)
        {
          P_op[kk] += (z[t][x][ii] * z[t][x][jj] +
                                            z[t][x][ii + 1] * z[t][x][jj + 1]);
          P_op[kk + 1] += (z[t][x][ii] * z[t][x][jj + 1] -
                                                z[t][x][ii + 1] * z[t][x][jj]);
        }
        else
        {
          long kk_sym = jj * N + ii;
          P_op[kk] = P_op[kk_sym];
          P_op[kk + 1] = - P_op[kk_sym + 1];
        }

        /*********************************************************************/
        /* Construct indices for the next pass.                              */
        /*********************************************************************/
        jj += 2;
        if (jj == Z_LEN)
        {
          ii += 2;
          jj = 0;
        }
      }
    }

    /*************************************************************************/
    /* Finish the computation of the P_t_avg operator by summing across      */
    /* the vector.                                                           */
    /*************************************************************************/
    for (long kk = 0; kk < P_LEN; kk += 2)
    {
      P_t_avg[t][kk] = horizontal_add(P_op[kk]);
      P_t_avg[t][kk + 1] = horizontal_add(P_op[kk + 1]);
    }
  }

  /***************************************************************************/
  /* Allocate memory resources for computation.                              */
  /***************************************************************************/
  double correlator_pool[NUM_CORES * (CORR_LEN + DOUBLE_CACHE_LINE)] = {0};

  /***************************************************************************/
  /* Compute the propagator from the object P_t_avg.                         */
  /***************************************************************************/
  #pragma omp parallel num_threads(OMP_NUM_THREADS / SMT)
  {
    int ID = omp_get_thread_num();
    double *corr;
    corr = &correlator_pool[ID * (CORR_LEN + DOUBLE_CACHE_LINE)];
    VEC_DOUBLE P_op[P_LEN];

    #pragma omp for schedule(static)
    for (unsigned int y = 0; y < L; y += N_d)
    {
      /***********************************************************************/
      /* Construct a vectorised P operator to speed up the calculation.      */
      /***********************************************************************/
    #if defined(__SSE2__)
      for (int i = 0; i < N_d; ++i)
      {
        for (long kk = 0; kk < P_LEN; kk += 2)
        {
          P_op[kk].insert(i, P_t_avg[y + i][kk]);
          P_op[kk + 1].insert(i, P_t_avg[y + i][kk + 1]);
        }
      }
    #else
      for (long kk = 0; kk < P_LEN; kk += 2)
      {
        P_op[kk] = P_t_avg[y][kk];
        P_op[kk + 1] = P_t_avg[y][kk + 1];
      }
    #endif

      /***********************************************************************/
      /* Add to each relevant propagator entry with respect to site y.       */
      /***********************************************************************/
      for (unsigned int x = 0; x < L; ++x)
      {
        /*********************************************************************/
        /* Add to the relevant propagator entries by computing the trace of  */
        /* the product of the two time averaged composite operators.         */
        /*********************************************************************/
        VEC_DOUBLE trace_1 = 0;
        VEC_DOUBLE trace_2 = 0;
        for (long kk = 0; kk < P_LEN; kk += 2)
        {
          trace_1 += P_t_avg[x][kk] * P_op[kk];
          trace_2 -= P_t_avg[x][kk + 1] * P_op[kk + 1];
        }
        trace_1 += trace_2;

        /*********************************************************************/
        /* Add this term to the corresponding correlator entries.            */
        /*********************************************************************/
        #if defined(__SSE2__)
        for (int i = 0; i < N_d; ++i)
        {
          unsigned int step = abs(y + i - x);
          if (step > (L / 2))
          {
            step = L - step;
          }
          corr[step] += trace_1[i];
        }
        #else
        unsigned int step = y;
        if (step > (L / 2))
        {
          step = L - step;
        }
        corr[step] += trace_1;
        #endif
      }
    }
  }

  /***************************************************************************/
  /* Combine the correlators of each thread and subtract the disconnected    */
  /* contribution.                                                           */
  /***************************************************************************/
  #pragma omp parallel for num_threads(NUM_CORES) schedule(static)
  for (unsigned int x = 0; x < CORR_LEN; ++x)
  {
    correlator[x] = 0;
    for (int ID = 0; ID < NUM_CORES; ++ID)
    {
      correlator[x] += correlator_pool[ID * (CORR_LEN + DOUBLE_CACHE_LINE) + x];
    }
    correlator[x] /= (V * L);

    // Average the folded portion of the correlator.
    if ((x > 0) && (x < (L / 2)))
    {
      correlator[x] *= 0.5;
    }
    correlator[x] -= 1.0 / N;
  }

  return;
}

/*****************************************************************************/
/* Wall-wall correlator, performed without convulation over time.            */
/*****************************************************************************/
static inline void wall_corr_alt(VEC_DOUBLE ***z,
                                 P_type** __restrict__ P_t_avg,
                                 double* __restrict__ correlator)
{
  /***************************************************************************/
  /* Allocate memory resources for computation.                              */
  /***************************************************************************/
  double correlator_pool[NUM_CORES * CORR_LEN] = {0};

  /***************************************************************************/
  /* For each timeslice in the feducial box, computer the correlator in the  */
  /* x direction.                                                            */
  /***************************************************************************/
  #pragma omp parallel for schedule(static)
  for (unsigned int t = 0; t < L; ++t)
  {
    /*************************************************************************/
    /* Get a pointer to the correlator memory for this thread.               */
    /*************************************************************************/
    int ID = omp_get_thread_num();
    double *corr;
    corr = &correlator_pool[ID * CORR_LEN];

    /*************************************************************************/
    /* For each x index in a vector, compute the correlation function.       */
    /*************************************************************************/
    for (unsigned int x = 0; x < NUM_VEC_BLOCKS; x += 2)
    {
      /***********************************************************************/
      /* Initialise variables to facilitate the construction of vectors with */
      /* points n steps away from this one.                                  */
      /***********************************************************************/
      int vec_no = x / 2;
      long x1 = vec_no * (2 * N_d);

      VEC_DOUBLE* phi_ref_r = z[0][x];
      VEC_DOUBLE* phi_ref_b = z[0][x + 1];
      VEC_DOUBLE* phi_1;
      VEC_DOUBLE* phi_2;

      /***********************************************************************/
      /* For each step away from this reference vector, compute the          */
      /* correlation.                                                        */
      /***********************************************************************/
      for (long dx = 0; dx < L; ++dx)
      {
        /*********************************************************************/
        /* Calculate the necessary variables to construct a z vector that is */
        /* dx sites away from the reference red site vector, and dx - 1 from */
        /* the black site vector.                                            */
        /*********************************************************************/
        long dx2 = dx - 1;
        if (dx2 < 0)
        {
          dx2 = L - 1;
        }

        long ind = 2 * (((x1 + dx) / (2 * N_d)) % NUM_RB_BLOCKS) + (dx % 2);
        long ind_2 = (ind + 2) % NUM_VEC_BLOCKS;
        int vec_type = (dx / 2) % N_d;
        phi_1 = z[t][ind];
        phi_2 = z[t][ind_2];

        /*********************************************************************/
        /* The correlator entry is the scalar product of the phi vector at   */
        /* the site x and the phi vector dx sites away.                      */
        /*********************************************************************/
        VEC_DOUBLE real_part_r1 = 0;
        VEC_DOUBLE real_part_r2 = 0;
        VEC_DOUBLE imag_part_r1 = 0;
        VEC_DOUBLE imag_part_r2 = 0;

        VEC_DOUBLE real_part_b1 = 0;
        VEC_DOUBLE real_part_b2 = 0;
        VEC_DOUBLE imag_part_b1 = 0;
        VEC_DOUBLE imag_part_b2 = 0;

        for (int ii = 0; ii < Z_LEN; ii += 2)
        {
          VEC_DOUBLE phi_dx_r = blend_phi(&phi_1[ii], &phi_2[ii], vec_type);
          VEC_DOUBLE phi_dx_i = blend_phi(&phi_1[ii + 1],
                                          &phi_2[ii + 1], vec_type);
          real_part_r1 += phi_dx_r * phi_ref_r[ii];
          real_part_r2 += phi_dx_i * phi_ref_r[ii + 1];
          imag_part_r1 += phi_dx_r * phi_ref_r[ii + 1];
          imag_part_r2 += phi_dx_i * phi_ref_r[ii];
          real_part_b1 += phi_dx_r * phi_ref_b[ii];
          real_part_b2 += phi_dx_i * phi_ref_b[ii + 1];
          imag_part_b1 += phi_dx_r * phi_ref_b[ii + 1];
          imag_part_b2 += phi_dx_i * phi_ref_b[ii];
        }

        VEC_DOUBLE real_part_r = real_part_r1 + real_part_r2;
        VEC_DOUBLE imag_part_r = imag_part_r1 - imag_part_r2;
        VEC_DOUBLE corr_entry_r = real_part_r * real_part_r +
                                  imag_part_r * imag_part_r;

        VEC_DOUBLE real_part_b = real_part_b1 + real_part_b2;
        VEC_DOUBLE imag_part_b = imag_part_b1 - imag_part_b2;
        VEC_DOUBLE corr_entry_b = real_part_b * real_part_b +
                                  imag_part_b * imag_part_b;

        /*********************************************************************/
        /* Assemble the index of the correlator for this value of dx and add */
        /* to the master array.                                              */
        /*********************************************************************/
        long corr_ind_r, corr_ind_b;
        if (dx <= (L / 2))
        {
          corr_ind_r = dx;
        }
        else
        {
          corr_ind_r = L - dx;
        }
        if (dx2 <= (L / 2))
        {
          corr_ind_b = dx2;
        }
        else
        {
          corr_ind_b = L - dx2;
        }
        corr[corr_ind_r] += horizontal_add(corr_entry_r);
        corr[corr_ind_b] += horizontal_add(corr_entry_b);
      }
    }
  }
  /***************************************************************************/
  /* Combine the correlators of each thread and subtract the disconnected    */
  /* contribution.                                                           */
  /***************************************************************************/
  #pragma omp parallel for schedule(static)
  for (unsigned int x = 0; x < CORR_LEN; ++x)
  {
    correlator[x] = 0;
    for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
    {
      correlator[x] += correlator_pool[ID * CORR_LEN + x];
    }
    correlator[x] /= V;

    // Average the folded part of the correlator.
    if ((x > 0) && (x < (L / 2)))
    {
      correlator[x] *= 0.5;
    }
    correlator[x] -= 1.0 / N;
  }

  return;
} /* wall_corr_alt */

#ifdef NO_BCS // This #ifdef is mostly depricated.
/*****************************************************************************/
/* Determine the simple wall-wall correlation function.                      */
/*****************************************************************************/
static inline void wall_corr_no_bc(VEC_DOUBLE ***z,
                                   P_type** __restrict__ P_t_avg,
                                   double* __restrict__ correlator)
{
  /***************************************************************************/
  /* Zero the P_t_avg object.                                                */
  /***************************************************************************/
  #pragma omp parallel for schedule(static)
  for (unsigned int x = 0; x < NUM_VEC_BLOCKS; ++x)
  {
    for (long kk = 0; kk < P_LEN; kk += 2)
    {
      P_t_avg[x][kk] = 0;
      P_t_avg[x][kk + 1] = 0;
    }
  }

  /***************************************************************************/
  /* Construct the object P_t_avg; that is sum the composite operator P(t,x) */
  /* over the t direction for each value of x.                               */
  /***************************************************************************/
  #pragma omp parallel for schedule(static)
  for (unsigned int x = 0; x < NUM_VEC_BLOCKS; ++x)
  {
    /*************************************************************************/
    /* Average the composite operator over each timeslice.                   */
    /*************************************************************************/
    for (unsigned int t = L_GAP; t < L - L_GAP; ++t)
    {
      /***********************************************************************/
      /* Add the tensor product of the z vector at this site with its        */
      /* complex conjugate to the corresponding P_t_avg entry.               */
      /***********************************************************************/
      int ii = 0;
      int jj = 0;

      for (long kk = 0; kk < P_LEN; kk += 2)
      {
        /*********************************************************************/
        /* The matrices in P_t_avg are symmetric, so compute only the upper  */
        /* triangular entries, and the rest can be recovered from the        */
        /* symmetry.                                                         */
        /*********************************************************************/
        if (jj >= ii)
        {
          P_t_avg[x][kk] += (z[t][x][ii] * z[t][x][jj] +
                                            z[t][x][ii + 1] * z[t][x][jj + 1]);
          P_t_avg[x][kk + 1] += (z[t][x][ii] * z[t][x][jj + 1] -
                                                z[t][x][ii + 1] * z[t][x][jj]);
        }
        else
        {
          long kk_sym = jj * N + ii;
          P_t_avg[x][kk] = P_t_avg[x][kk_sym];
          P_t_avg[x][kk + 1] = - P_t_avg[x][kk_sym + 1];
        }

        /*********************************************************************/
        /* Construct indices for the next pass.                              */
        /*********************************************************************/
        jj += 2;
        if (jj == Z_LEN)
        {
          ii += 2;
          jj = 0;
        }
      }
    }
  }

  /***************************************************************************/
  /* Allocate memory resources for computation.                              */
  /***************************************************************************/
  double correlator_pool[OMP_NUM_THREADS * CORR_LEN] = {0};

  /***************************************************************************/
  /* Compute the propagator from the object P_t_avg.                         */
  /***************************************************************************/
  #pragma omp for schedule(static)
  for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
  {
    unsigned int min_y = ID * L_STEP;
    unsigned int max_y = (ID + 1) * L_STEP;
    if (max_y > L)
    {
      max_y = L;
    }
    double *corr;
    corr = &correlator_pool[ID * CORR_LEN];

    for (unsigned int y = min_y; y < max_y; ++y)
    {
      /***********************************************************************/
      /* Save vector indices corresponding to the matrix we will use         */
      /* repeatedly.                                                         */
      /***********************************************************************/
      unsigned int y_v = (y / (2 * N_d)) * 2 + (y % 2);
      unsigned int y_i = y % (2 * N_d);
      y_i /= 2;

      /***********************************************************************/
      /* Add to each relevant propagator entry with respect to site y.       */
      /***********************************************************************/
      for (unsigned int x = 0; x < NUM_VEC_BLOCKS; ++x)
      {
        /*********************************************************************/
        /* Add to the relevant propagator entries by computing the trace of  */
        /* the product of the two time averaged composite operators.         */
        /*********************************************************************/
        VEC_DOUBLE trace_1 = 0;
        VEC_DOUBLE trace_2 = 0;
        for (long kk = 0; kk < P_LEN; kk += 2)
        {
          trace_1 += P_t_avg[x][kk] * P_t_avg[y_v][kk][y_i];
          trace_2 -= P_t_avg[x][kk + 1] * P_t_avg[y_v][kk + 1][y_i];
        }
        trace_1 += trace_2;

        /*********************************************************************/
        /* Add this term to the corresponding correlator entries.            */
        /*********************************************************************/
        int x_start = x - (x % 2);
        x_start *= N_d;
        x_start += x % 2;
        x_start -= y;

        for (int i = 0; i < N_d; ++i)
        {
          unsigned int step = abs(x_start + 2 * i);
          if (step > (L / 2))
          {
            step = L - step;
          }
          corr[step] += trace_1[i];
        }
      }
    }
  }

  /***************************************************************************/
  /* Combine the correlators of each thread and subtract the disconnected    */
  /* contribution.                                                           */
  /***************************************************************************/
  #pragma omp parallel for schedule(static)
  for (unsigned int x = 0; x < CORR_LEN; ++x)
  {
    correlator[x] = 0;
    for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
    {
      correlator[x] += correlator_pool[ID * CORR_LEN + x];
    }
    correlator[x] /= (V_TOT * L_T);

    // Average the folded part of the correlator.
    if ((x > 0) && (x < (L / 2)))
    {
      correlator[x] *= 0.5;
    }
    correlator[x] -= 1.0 / N;
  }

  return;
} /* wall_corr_no_bc */


static inline void wall_corr_no_bc_alt(VEC_DOUBLE ***z,
                                       P_type** __restrict__ P_t_avg,
                                       double* __restrict__ correlator)
{
  /***************************************************************************/
  /* Allocate memory resources for computation.                              */
  /***************************************************************************/
  double correlator_pool[NUM_CORES * CORR_LEN] = {0};

  /***************************************************************************/
  /* For each timeslice in the feducial box, computer the correlator in the  */
  /* x direction.                                                            */
  /***************************************************************************/
  #pragma omp parallel for schedule(static)
  for (unsigned int t = L_GAP; t < L - L_GAP; ++t)
  {
    /*************************************************************************/
    /* Get a pointer to the correlator memory for this thread.               */
    /*************************************************************************/
    int ID = omp_get_thread_num();
    double *corr;
    corr = &correlator_pool[ID * CORR_LEN];

    /*************************************************************************/
    /* For each x index in a vector, compute the correlation function.       */
    /*************************************************************************/
    for (unsigned int x = 0; x < NUM_VEC_BLOCKS; x += 2)
    {
      /***********************************************************************/
      /* Initialise variables to facilitate the construction of vectors with */
      /* points n steps away from this one. Keep the reference vector fixed  */
      /* in the centre of the time direction of the lattice to minimise      */
      /* influence from the open boundary condition.                         */
      /***********************************************************************/
      int vec_no = x / 2;
      long x1 = vec_no * (2 * N_d);

      VEC_DOUBLE* phi_ref_r = z[L / 2][x];
      VEC_DOUBLE* phi_ref_b = z[L / 2][x + 1];
      VEC_DOUBLE* phi_1;
      VEC_DOUBLE* phi_2;

      /***********************************************************************/
      /* For each step away from this reference vector, compute the          */
      /* correlation.                                                        */
      /***********************************************************************/
      for (long dx = 0; dx < L; ++dx)
      {
        /*********************************************************************/
        /* Calculate the necessary variables to construct a z vector that is */
        /* dx sites away from the reference red site vector, and dx - 1 from */
        /* the black site vector.                                            */
        /*********************************************************************/
        long dx2 = dx - 1;
        if (dx2 < 0)
        {
          dx2 = L - 1;
        }

        long ind = 2 * (((x1 + dx) / (2 * N_d)) % NUM_RB_BLOCKS) + (dx % 2);
        long ind_2 = (ind + 2) % NUM_VEC_BLOCKS;
        int vec_type = (dx / 2) % N_d;
        phi_1 = z[t][ind];
        phi_2 = z[t][ind_2];

        /*********************************************************************/
        /* The correlator entry is the scalar product of the phi vector at   */
        /* the site x and the phi vector dx sites away.                      */
        /*********************************************************************/
        VEC_DOUBLE real_part_r1 = 0;
        VEC_DOUBLE real_part_r2 = 0;
        VEC_DOUBLE imag_part_r1 = 0;
        VEC_DOUBLE imag_part_r2 = 0;

        VEC_DOUBLE real_part_b1 = 0;
        VEC_DOUBLE real_part_b2 = 0;
        VEC_DOUBLE imag_part_b1 = 0;
        VEC_DOUBLE imag_part_b2 = 0;

        for (int ii = 0; ii < Z_LEN; ii += 2)
        {
          VEC_DOUBLE phi_dx_r = blend_phi(&phi_1[ii], &phi_2[ii], vec_type);
          VEC_DOUBLE phi_dx_i = blend_phi(&phi_1[ii + 1],
                                          &phi_2[ii + 1], vec_type);
          real_part_r1 += phi_dx_r * phi_ref_r[ii];
          real_part_r2 += phi_dx_i * phi_ref_r[ii + 1];
          imag_part_r1 += phi_dx_r * phi_ref_r[ii + 1];
          imag_part_r2 += phi_dx_i * phi_ref_r[ii];
          real_part_b1 += phi_dx_r * phi_ref_b[ii];
          real_part_b2 += phi_dx_i * phi_ref_b[ii + 1];
          imag_part_b1 += phi_dx_r * phi_ref_b[ii + 1];
          imag_part_b2 += phi_dx_i * phi_ref_b[ii];
        }

        VEC_DOUBLE real_part_r = real_part_r1 + real_part_r2;
        VEC_DOUBLE imag_part_r = imag_part_r1 - imag_part_r2;
        VEC_DOUBLE corr_entry_r = real_part_r * real_part_r +
                                  imag_part_r * imag_part_r;

        VEC_DOUBLE real_part_b = real_part_b1 + real_part_b2;
        VEC_DOUBLE imag_part_b = imag_part_b1 - imag_part_b2;
        VEC_DOUBLE corr_entry_b = real_part_b * real_part_b +
                                  imag_part_b * imag_part_b;

        /*********************************************************************/
        /* Assemble the index of the correlator for this value of dx and add */
        /* to the master array.                                              */
        /*********************************************************************/
        long corr_ind_r, corr_ind_b;
        if (dx <= (L / 2))
        {
          corr_ind_r = dx;
        }
        else
        {
          corr_ind_r = L - dx;
        }
        if (dx2 <= (L / 2))
        {
          corr_ind_b = dx2;
        }
        else
        {
          corr_ind_b = L - dx2;
        }
        corr[corr_ind_r] += horizontal_add(corr_entry_r);
        corr[corr_ind_b] += horizontal_add(corr_entry_b);
      }
    }
  }

  /***************************************************************************/
  /* Combine the correlators of each thread and subtract the disconnected    */
  /* contribution.                                                           */
  /***************************************************************************/
  #pragma omp parallel for schedule(static)
  for (unsigned int x = 0; x < CORR_LEN; ++x)
  {
    correlator[x] = 0;
    for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
    {
      correlator[x] += correlator_pool[ID * CORR_LEN + x];
    }
    correlator[x] /= (L * L_T);

    // Average the folded part of the correlator.
    if ((x > 0) && (x < (L / 2)))
    {
      correlator[x] *= 0.5;
    }
    correlator[x] -= 1.0 / N;
  }

  return;
} /* wall_corr_no_bc_alt */
#endif


/*****************************************************************************/
/* First two components of the Fourier Transform of the correlation          */
/* function.                                                                 */
/*****************************************************************************/
static inline double G_k(VEC_DOUBLE ***z, double **P_0, double **P_1,
                         double *G_components, update_limits* limits)
{
  /***************************************************************************/
  /* First find the first two components of the Fourier Transform of the     */
  /* projection operator.                                                    */
  /***************************************************************************/
  G_components[0] = 0.0;
  G_components[1] = 0.0;

  /***************************************************************************/
  /* Compute the entries of the projection operator in parallel.             */
  /***************************************************************************/
  #pragma omp parallel for num_threads(OMP_NUM_THREADS) schedule(static)
  for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
  {
    /*************************************************************************/
    /* Zero all the entries of the projection operator before computing      */
    /* them.                                                                 */
    /*************************************************************************/
    for (long kk = 0; kk < P_LEN; kk += 2)
    {
      P_0[ID][kk] = 0.0;
      P_1[ID][kk] = 0.0;
      P_0[ID][kk + 1] = 0.0;
      P_1[ID][kk + 1] = 0.0;
    }
//#ifdef NO_BC
//    G_k_loop_alt(z, P_0[ID], P_1[ID], &limits[ID]);
//#else
    G_k_loop(z, P_0[ID], P_1[ID], &limits[ID]);
//#endif
  }

  /***************************************************************************/
  /* Gather the results produced by the different threads.                   */
  /***************************************************************************/
  #pragma omp parallel for num_threads(OMP_NUM_THREADS) schedule(static)
  for (long kk = 0; kk < P_LEN; kk += 2)
  {
    for (int ID = 1; ID < OMP_NUM_THREADS; ++ID)
    {
      P_0[0][kk] += P_0[ID][kk];
      P_0[0][kk + 1] += P_0[ID][kk + 1];
      P_1[0][kk] += P_1[ID][kk];
      P_1[0][kk + 1] += P_1[ID][kk + 1];
    }
  }

  /***************************************************************************/
  /* Next find the trace of the matrix product of the projection operator    */
  /* with its complex conjugate.                                             */
  /***************************************************************************/
  for (int kk = 0; kk < P_LEN; kk += 2)
  {
    G_components[0] += P_0[0][kk] * P_0[0][kk];
    G_components[0] += P_0[0][kk + 1] * P_0[0][kk + 1];
    G_components[1] += P_1[0][kk] * P_1[0][kk];
    G_components[1] += P_1[0][kk + 1] * P_1[0][kk + 1];
  }

  G_components[0] /= V;
  G_components[0] -= (double)V / (double)N;
  G_components[1] /= V;

  return G_FLOPS;
} /* G_k */


#ifdef NO_BC
/*****************************************************************************/
/* First two components of the Fourier Transform of the correlation          */
/* function.                                                                 */
/*****************************************************************************/
static inline double G_k_no_bc(VEC_DOUBLE ***z, double **P_0, double **P_1,
                               double *G_components, update_limits* limits)
{
  /***************************************************************************/
  /* First find the first two components of the Fourier Transform of the     */
  /* projection operator.                                                    */
  /***************************************************************************/
  G_components[0] = 0.0;
  G_components[1] = 0.0;

  /***************************************************************************/
  /* Compute the entries of the projection operator in parallel.             */
  /***************************************************************************/
  #pragma omp parallel for num_threads(OMP_NUM_THREADS) schedule(static)
  for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
  {
    /*************************************************************************/
    /* Zero all the entries of the projection operator before computing      */
    /* them.                                                                 */
    /*************************************************************************/
    for (long kk = 0; kk < P_LEN; kk += 2)
    {
      P_0[ID][kk] = 0.0;
      P_1[ID][kk] = 0.0;
      P_0[ID][kk + 1] = 0.0;
      P_1[ID][kk + 1] = 0.0;
    }
    G_k_loop_no_bc(z, P_0[ID], P_1[ID], &limits[ID]);
  }

  /***************************************************************************/
  /* Gather the results produced by the different threads.                   */
  /***************************************************************************/
  #pragma omp parallel for num_threads(OMP_NUM_THREADS) schedule(static)
  for (long kk = 0; kk < P_LEN; kk += 2)
  {
    for (int ID = 1; ID < OMP_NUM_THREADS; ++ID)
    {
      P_0[0][kk] += P_0[ID][kk];
      P_0[0][kk + 1] += P_0[ID][kk + 1];
      P_1[0][kk] += P_1[ID][kk];
      P_1[0][kk + 1] += P_1[ID][kk + 1];
    }
  }

  /***************************************************************************/
  /* Next find the trace of the matrix product of the projection operator    */
  /* with its complex conjugate.                                             */
  /***************************************************************************/
  for (int kk = 0; kk < P_LEN; kk += 2)
  {
    G_components[0] += P_0[0][kk] * P_0[0][kk];
    G_components[0] += P_0[0][kk + 1] * P_0[0][kk + 1];
    G_components[1] += P_1[0][kk] * P_1[0][kk];
    G_components[1] += P_1[0][kk + 1] * P_1[0][kk + 1];
  }

  G_components[0] /= V;
  //G_components[0] -= V / N;
  G_components[1] /= V;

  return G_FLOPS;
} /* G_k_no_bc */
#endif

/*****************************************************************************/
/* Compuation loop for topological susceptibility.                           */
/*****************************************************************************/
void top_loop(VEC_DOUBLE ***z, VEC_DOUBLE __restrict__ *top_charge,
              update_limits * __restrict__ limits)
{
  unsigned int t = limits->m_min_t;
  unsigned int x = limits->m_min_x;
  long max_k = limits->num_block_meas;
  // tu only wraps around if fiducial volume is equivalent to full volume.
  unsigned int tu = (t + 1) % L_T;
  unsigned int xu = (x + 1);
  unsigned int xu2 = (xu + 1) % NUM_VEC_BLOCKS;
  for (long kk = 0; kk < max_k; kk += 2)
  {
    /*************************************************************************/
    /* Compute the trace of the product of three projection operators; to    */
    /* save time use the definition in terms of scalar products of z         */
    /* vectors. Start by zeroing the 8 necessary components of each          */
    /* product.                                                              */
    /*************************************************************************/
    //VEC_DOUBLE x0[2];
    VEC_DOUBLE x1[2];
    VEC_DOUBLE x2[2];
    VEC_DOUBLE x3[2];
    VEC_DOUBLE x4[2];
    VEC_DOUBLE x5[2];
    VEC_DOUBLE x6[2];
    VEC_DOUBLE x7[2];
    VEC_DOUBLE x8[2];
    //VEC_DOUBLE x9[2];

    for (int k = 0; k < 2; ++k)
    {
      //x0[k] = 0;
      x1[k] = 0;
      x2[k] = 0;
      x3[k] = 0;
      x4[k] = 0;
      x5[k] = 0;
      x6[k] = 0;
      x7[k] = 0;
      x8[k] = 0;
      //x9[k] = 0;
    }

    VEC_DOUBLE * phi_z = z[t][x];
    VEC_DOUBLE * phi_z_tu = z[tu][x];
    VEC_DOUBLE * phi_z_xu = z[t][xu];
    VEC_DOUBLE * phi_z_txu = z[tu][xu];

    /*************************************************************************/
    /* The first phi_up vector is aligned correctly; we can proceed safely.  */
    /* Compute the 4 scalar products between adjacent lattice sites in a     */
    /* square.                                                               */
    /*************************************************************************/
    for (int n = 0; n < Z_LEN; n += 2)
    {
      /*x0[0] += phi_z_txu[n] * phi_z[n];
      x0[0] += phi_z_txu[n + 1] * phi_z[n + 1];
      x0[1] += phi_z_txu[n] * phi_z[n + 1];
      x0[1] -= phi_z_txu[n + 1] * phi_z[n];*/
      x1[0] += phi_z_tu[n] * phi_z_txu[n];
      x1[0] += phi_z_tu[n + 1] * phi_z_txu[n + 1];
      x1[1] += phi_z_tu[n] * phi_z_txu[n + 1];
      x1[1] -= phi_z_tu[n + 1] * phi_z_txu[n];
      x2[0] += phi_z[n] * phi_z_tu[n];
      x2[0] += phi_z[n + 1] * phi_z_tu[n + 1];
      x2[1] += phi_z[n] * phi_z_tu[n + 1];
      x2[1] -= phi_z[n + 1] * phi_z_tu[n];
      x3[0] += phi_z_xu[n] * phi_z[n];
      x3[0] += phi_z_xu[n + 1] * phi_z[n + 1];
      x3[1] += phi_z_xu[n] * phi_z[n + 1];
      x3[1] -= phi_z_xu[n + 1] * phi_z[n];
      x4[0] += phi_z_txu[n] * phi_z_xu[n];
      x4[0] += phi_z_txu[n + 1] * phi_z_xu[n + 1];
      x4[1] += phi_z_txu[n] * phi_z_xu[n + 1];
      x4[1] -= phi_z_txu[n + 1] * phi_z_xu[n];
    }

    /*************************************************************************/
    /* For the next site we must load the aligned vectors in each iteration. */
    /*************************************************************************/
    VEC_DOUBLE * phi_z_xd = phi_z;
    VEC_DOUBLE * phi_z_tuxd = phi_z_tu;
    phi_z = phi_z_xu;
    phi_z_tu = phi_z_txu;
    phi_z_xu = z[t][xu2];
    phi_z_txu = z[tu][xu2];
    for (int n = 0; n < Z_LEN; n += 2)
    {
      VEC_DOUBLE phi_z_txu_r = load_phi_xu(&phi_z_tuxd[n],

                                                    &phi_z_txu[n]);
      VEC_DOUBLE phi_z_txu_i = load_phi_xu(&phi_z_tuxd[n + 1],
                                                    &phi_z_txu[n + 1]);
      x5[0] += phi_z_tu[n] * phi_z_txu_r;
      x5[1] += phi_z_tu[n] * phi_z_txu_i;
      x6[0] += phi_z[n] * phi_z_tu[n];
      x6[1] += phi_z[n] * phi_z_tu[n + 1];
      x5[0] += phi_z_tu[n + 1] * phi_z_txu_i;
      x5[1] -= phi_z_tu[n + 1] * phi_z_txu_r;
      x6[0] += phi_z[n + 1] * phi_z_tu[n + 1];
      x6[1] -= phi_z[n + 1] * phi_z_tu[n];
      VEC_DOUBLE phi_z_xu_r = load_phi_xu(&phi_z_xd[n],
                                                   &phi_z_xu[n]);
      VEC_DOUBLE phi_z_xu_i = load_phi_xu(&phi_z_xd[n + 1],
                                                   &phi_z_xu[n + 1]);
      x7[0] += phi_z_xu_r * phi_z[n];
      x7[1] += phi_z_xu_r * phi_z[n + 1];
      x8[0] += phi_z_txu_r * phi_z_xu_r;
      x8[1] += phi_z_txu_r * phi_z_xu_i;
      x7[0] += phi_z_xu_i * phi_z[n + 1];
      x7[1] -= phi_z_xu_i * phi_z[n];
      x8[0] += phi_z_txu_i * phi_z_xu_i;
      x8[1] -= phi_z_txu_i * phi_z_xu_r;
      /*x9[0] += phi_z_txu_r * phi_z[n];
      x9[0] += phi_z_txu_i * phi_z[n + 1];
      x9[1] += phi_z_txu_r * phi_z[n + 1];
      x9[1] -= phi_z_txu_i * phi_z[n];*/
    }

    /***********************************************************************/
    /* Compute the topological charge from the scalar products, skipping   */
    /* the factor of 1/2pi for efficiency.                                 */
    /***********************************************************************/
    VEC_DOUBLE real_part_1 = x1[0] * x2[0] - x1[1] * x2[1];
    VEC_DOUBLE imag_part_1 = x1[0] * x2[1] + x1[1] * x2[0];
    VEC_DOUBLE real_part_2 = x3[0] * x4[0] - x3[1] * x4[1];
    VEC_DOUBLE imag_part_2 = x3[0] * x4[1] + x3[1] * x4[0];
    VEC_DOUBLE real1 = real_part_1 * real_part_2 - imag_part_1 * imag_part_2;
    VEC_DOUBLE imag1 = real_part_1 * imag_part_2 + imag_part_1 * real_part_2;
    VEC_DOUBLE real_part_3 = x5[0] * x6[0] - x5[1] * x6[1];
    VEC_DOUBLE imag_part_3 = x5[0] * x6[1] + x5[1] * x6[0];
    VEC_DOUBLE real_part_4 = x7[0] * x8[0] - x7[1] * x8[1];
    VEC_DOUBLE imag_part_4 = x7[0] * x8[1] + x7[1] * x8[0];
    VEC_DOUBLE real2 = real_part_3 * real_part_4 - imag_part_3 * imag_part_4;
    VEC_DOUBLE imag2 = real_part_3 * imag_part_4 + imag_part_3 * real_part_4;
    VEC_DOUBLE charges_1 = atan2(imag1, real1);
    VEC_DOUBLE charges_2 = atan2(imag2, real2);

    /*#pragma omp critical
    {
      for (int i = 0; i < N_d; ++i)
      {
        std::cout << charges_1[i] / (2 * M_PI) << std::endl;
        std::cout << charges_2[i] / (2 * M_PI) << std::endl;
      }
    }*/

    *top_charge += atan2(imag1, real1) + atan2(imag2, real2);

    /***********************************************************************/
    /* Calculate the indices for the next lattice point.                   */
    /***********************************************************************/
    x += 2;
    xu += 2;
    xu2 += 2;
    if (xu2 == NUM_VEC_BLOCKS)
    {
      xu2 = 0;
    }
    else if (x == NUM_VEC_BLOCKS)
    {
      ++t;
      ++tu;
      if (tu == L)
      {
        tu = 0;
      }
      x = 0;
      xu = 1;
      if (t == L)
      {
        kk = max_k;
      }
    }
  }
}

/*****************************************************************************/
/* Topological charge.                                                       */
/*****************************************************************************/
double topological_charge(VEC_DOUBLE ***z, update_limits *limits)
{
  /***************************************************************************/
  /* Compute the topological charge at each lattice point.                   */
  /***************************************************************************/
  double top_charge(0);
  #pragma omp parallel for num_threads(OMP_NUM_THREADS) reduction(+:top_charge)
  for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
  {
    VEC_DOUBLE top_charges = 0.0;
    top_loop(z, &top_charges, &limits[ID]);
    top_charge += horizontal_add(top_charges);
  }

  top_charge /= 2 * M_PI;
  //top_charge /= (2 * M_PI * L);
  //top_charge *= top_charge;
  // HACKY TEST -> V FACTOR -> L * L_T
  //top_charge /= (L * L_T);

  if (top_charge != top_charge)
  {
    std::cout << "Error in topological charge measurement\n";
    //exit(0);
  }
  return top_charge;
} /* topological_charge */


/*****************************************************************************/
/* Compute Monte Carlo estimates.                                            */
/*****************************************************************************/
static inline double comp_MC_ests(VEC_DOUBLE ***z, VEC_DOUBLE ***lam,
       double **P_0, double **P_1, double* correlator,
       P_type** P_t_avg, double **MC_bin, const long bin_size,
       const long ii, update_limits* limits, VEC_DOUBLE *angles)
{
  int n = 0;
  double G_ests[2];
  double FLOPs = 0;
  double FLOPS(0);
  long jj = ii % bin_size;

  /* Action */
#if (MEASUREMENTS & LINK_ENERGY)
  MC_bin[n++][jj] = S_1(z, lam);
  FLOPs += E_FLOPS;
#endif

  /* Wall-wall correlation */
#if (MEASUREMENTS & WALL_CORR)
  //link_angle(z, lam, &angles[2 * L * NUM_VEC_BLOCKS * jj]);
  long corr_jj = jj % CORR_N_COR;
  if (corr_jj == 0)
  {
    corr_jj = jj / CORR_N_COR;
    wall_corr_alt(z, P_t_avg, &correlator[corr_jj * CORR_LEN]);
    FLOPS += CORR_FLOPS;
  }
#endif

  /* G_k */
#if (MEASUREMENTS & (MAG_SUS + G_k_1))
  FLOPs += G_k(z, P_0, P_1, G_ests, limits);
  #if (MEASUREMENTS & MAG_SUS)
    MC_bin[n++][jj] = G_ests[0];
  #endif
  #if (MEASUREMENTS & G_k_1)
    MC_bin[n++][jj] = G_ests[1];
  #endif
#endif
  /*
#ifdef NO_BC
  #if (MEASUREMENTS & (G_1_0 + G_1_1))
    FLOPs += G_k_no_bc(z, P_0, P_1, G_ests, limits);
    #if (MEASUREMENTS & G_1_0)
    MC_bin[n++][jj] = G_ests[0];
    #endif
    #if (MEASUREMENTS & G_1_1)
    MC_bin[n++][jj] = G_ests[1];
    #endif
  #endif
#else
  #if (MEASUREMENTS & (MAG_SUS + G_k_1))
    FLOPs += G_k_no_bc(z, P_0, P_1, G_ests, limits);
    #if (MEASUREMENTS & MAG_SUS)
    MC_bin[n++][jj] = G_ests[0];
    #endif
    #if (MEASUREMENTS & G_k_1)
    MC_bin[n++][jj] = G_ests[1];
    #endif
  #endif
#endif
*/
  /* Topological charge */
#if (MEASUREMENTS & TOP_CHARGE)
  MC_bin[n][jj] = topological_charge(z, limits);
  FLOPs += TOP_FLOPS;
#endif

  return FLOPs;
} /* comp_MC_ests */


/*****************************************************************************/
/* Jackknife error analysis.                                                 */
/*****************************************************************************/
static inline double jackknife(double *MC_ests, long len, double MC_avg)
{
  double err = 0;
  double *arr;
  arr = new double[len];

  /***************************************************************************/
  /* Compute averages for copies of the array, dropping one element at a     */
  /* time.                                                                   */
  /***************************************************************************/
  double sum = 0;
  for(long ii = 0; ii < len; ii++)
  {
    sum += MC_ests[ii];
  }
  
  for (long ii = 0; ii < len; ii++)
  {
    arr[ii] = (sum - MC_ests[ii]) / (len - 1);
  }

  /***************************************************************************/
  /* Compute the error from this array of averages.                          */
  /***************************************************************************/
  for(int ii = 0; ii < len; ii++)
  {
    double square = (arr[ii] - MC_avg) * (arr[ii] - MC_avg);
    err += square;
  }

  return sqrt(err * (len - 1) / len);
} /* jackknife */

/*****************************************************************************/
/* Thermalise the lattice.                                                   */
/*****************************************************************************/
static inline double thermalise_lattice(VEC_DOUBLE ***z, VEC_DOUBLE ***lam,
                                        dsfmt_t**** r, update_limits* limits,
                                        long * z_stream_counter,
                                        long * lam_stream_counter,
                                        VEC_DOUBLE *lattice_memory,
                                        std::string path)
{
  long FLOPs;
  double total_FLOPs(0);
  if ((!RESTART) && (!LOAD_THERM))
  {
    for (long ii = 0; ii < N_THERM; ii++)
    {
      FLOPs = heatbath_update(z, lam, r, limits, z_stream_counter,
                              lam_stream_counter);
      total_FLOPs += (double)FLOPs;
    }

    /*************************************************************************/
    /* If the flag is set, save off the thermalised config.                  */
    /*************************************************************************/
    if (SAVE_THERM != 0)
    {
      std::stringstream lattice_file;
      std::ofstream outfile;
      lattice_file << path << "_thermalised_z_config";
      outfile.open(lattice_file.str(), std::ios::out | std::ios::binary);
      if (outfile.is_open())
      {
        outfile.write(reinterpret_cast<const char*>(lattice_memory),
                      std::streamsize(N_d * L * NUM_VEC_BLOCKS *
                                      SITE_WIDTH * sizeof(double)));
        outfile.close();
      }
    }
  }
  return total_FLOPs;
}

void print_update_limits(update_limits* limits)
{
  for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
  {
    std::cout << "num block updates = " << limits[ID].num_block_ups << std::endl;
    std::cout << "min t = " << limits[ID].u_min_t;
    std::cout << ", max t = " << limits[ID].u_max_t;
    std::cout << ", min x = " << limits[ID].u_min_x;
    std::cout << ", max x = " << limits[ID].u_max_x << std::endl;
    std::cout << "num block meas = " << limits[ID].num_block_meas << std::endl;
    std::cout << "min t = " << limits[ID].m_min_t;
    std::cout << ", max t = " << limits[ID].m_max_t;
    std::cout << ", min x = " << limits[ID].m_min_x;
    std::cout << ", max x = " << limits[ID].m_max_x << std::endl;
  }
}


/*****************************************************************************/
/* Calculate the acceptance rate of the heatbath update procedure in the     */
/* simulation by comparing the minimum required number of generated random   */
/* numbers to the actual numbers provided by the stream counters.            */
/*****************************************************************************/
double calculate_acc_rate(long* stream_counter, long len, int n_per_site)
{
  double acc_rate = 0;
  for (long ii = 0; ii < len; ++ii)
  {
    acc_rate += (double)(n_per_site * N_CF) / (double)stream_counter[ii];
  }
  acc_rate /= len;

  double acc_rate_err = 0;
  for (long ii = 0; ii < len; ++ii)
  {
    acc_rate_err += pow(((double)(n_per_site * N_CF) /
                        (double)stream_counter[ii]) - acc_rate, 2.0);
  }
  acc_rate_err /= len;
  acc_rate_err = sqrt(acc_rate_err);

  std::cout << acc_rate << " (" << acc_rate_err << ")" << std::endl;
  return acc_rate;
}


/*****************************************************************************/
/* Monte Carlo algorithm.                                                    */
/*****************************************************************************/
void monte_carlo(std::string **obs_list, std::string path, double *MC_avg,
                 double *MC_err)
{
  /***************************************************************************/
  /* Initialize variables.                                                   */
  /***************************************************************************/
  long index, ind;
  unsigned long flops = 0;
  double flops_per_sec = 0;

  /***************************************************************************/
  /* Initialize the random number generator and allocate memory for an array */
  /* of random number generator instances: one for each lattice site.        */
  /***************************************************************************/
  unsigned long seed = 1000000;
  dsfmt_t* gen_mem;
  dsfmt_t**** r;
  long *z_stream_counter;
  long *lam_stream_counter;
  initialize_counters(&z_stream_counter, &lam_stream_counter, path);
  initialize_rng(&r, &gen_mem, seed, z_stream_counter, lam_stream_counter, path);

  /***************************************************************************/
  /* Allocate memory for the complex scalar field z at each lattice point    */
  /* and the gauge links between each point in both the time and space       */
  /* directions. Construct the array such that the forward gauge links       */
  /* between site n and n+mu are stored at site n.                           */
  /***************************************************************************/
  #define ALLOCATE_DYNAMICALLY 1
#ifdef ALLOCATE_DYNAMICALLY
  VEC_DOUBLE *lattice_memory;
  void *aligned_memory;
  int rc = posix_memalign(&aligned_memory, N_d * sizeof(double),
                          L_T * NUM_VEC_BLOCKS * SITE_WIDTH *
                          sizeof(VEC_DOUBLE));
  if (rc != 0)
  {
    std::cout << "Could not allocate memory\n" << std::endl;
    exit(0);
  }
  lattice_memory = (VEC_DOUBLE *)aligned_memory;
#else
  VEC_DOUBLE lattice_memory[L_T * NUM_VEC_BLOCKS * SITE_WIDTH];
#endif
  VEC_DOUBLE ***z, ***lam;
  construct_lattice(&z, &lam, lattice_memory);

  /***************************************************************************/
  /* Define the limits for the blocks that each thread will update.          */
  /***************************************************************************/
  omp_set_num_threads(OMP_NUM_THREADS);
  update_limits limits[OMP_NUM_THREADS];
  define_update_limits(limits);

  /***************************************************************************/
  /* Allocate memory for the 2D Monte Carlo measurement array.               */
  /***************************************************************************/
  double **MC_bin, **MC_est;
  MC_bin = new double*[NUM_OBS];
  MC_est = new double*[NUM_OBS];
  long bin_size = N_CF / N_BINS;
  for (int n = 0; n < NUM_OBS; ++n)
  {
    MC_bin[n] = new double[bin_size];
    MC_est[n] = new double[N_BINS];
  }

  /***************************************************************************/
  /* Allocate memory for projection operator matrices. Pad the array in      */
  /* order to eliminate false sharing between threads.                       */
  /***************************************************************************/
  double *P_0[OMP_NUM_THREADS];
  double *P_1[OMP_NUM_THREADS];

  for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
  {
    P_0[ID] = new double[P_LEN];
    P_1[ID] = new double[P_LEN + CACHE_LINE];
  }

  /***************************************************************************/
  /* Generate a starting configuration for the simulation; either by loading */
  /* a saved configuration or by randomly generating a new one.              */
  /***************************************************************************/
  generate_start_config(z, lam, r, lattice_memory, path, z_stream_counter,
                        lam_stream_counter);

// HACKY TEST - DISABLE NO_BC
 /* #undef NO_BC
  for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
  {
    limits[ID].num_block_meas = limits[ID].num_block_ups;
    limits[ID].m_min_t = limits[ID].u_min_t;
    limits[ID].m_max_t = limits[ID].u_max_t;
    limits[ID].m_min_x = limits[ID].u_min_x;
    limits[ID].m_max_x = limits[ID].u_max_x;
  }
  #if L_T % OMP_NUM_THREADS == 0
    #define L_STEP (L_T/NUM_CORES)
  #else
    #define L_STEP ((L_T/NUM_CORES) + 1)
  #endif
  //#define L_GAP 0
  std::cout << "L_T = " << L_T << std::endl;
// END OF HACKY TEST

  /***************************************************************************/
  /* Update the lattice N_therm times to thermalise the lattice.             */
  /***************************************************************************/
  auto t1 = std::chrono::high_resolution_clock::now();
  flops_per_sec = thermalise_lattice(z, lam, r, limits, z_stream_counter,
                                     lam_stream_counter, lattice_memory, path);
  auto t2 = std::chrono::high_resolution_clock::now();
  double time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()*1e-6;
  flops_per_sec = flops_per_sec / time;
  std::cout << "Took ";
  std::cout << time << "s to thermalise." << std::endl;
  std::cout << "Thermalisation GFLOPS = " << flops_per_sec*1e-9 << std::endl;

  /***************************************************************************/
  /* Update the lattice N_cor times then use it to make an estimate of the   */
  /* lattice observable specified by obs. Repeat N_cf times.                 */
  /***************************************************************************/
  double total_FLOPs;
  auto start = std::chrono::high_resolution_clock::now();
  total_FLOPs = run_simulation(z, lam, r, limits, gen_mem, z_stream_counter,
                               lam_stream_counter, lattice_memory, path,
                               P_0, P_1, MC_est, obs_list);
  auto end = std::chrono::high_resolution_clock::now();
  time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()*1e-3;
  total_FLOPs /= time;
  std::cout << "Simulation GFLOPS = " << total_FLOPs*1e-9 << std::endl;
  std::cout << "z acceptance rate = ";
  double z_acc_rate = calculate_acc_rate(z_stream_counter,
                                         L_T * NUM_VEC_BLOCKS * N_d,
                                         Z_RANDS_PER_SITE);
  std::cout << "lam acceptance rate = ";
  double lam_acc_rate = calculate_acc_rate(lam_stream_counter,
                                           L_T * NUM_VEC_BLOCKS * N_d,
                                           LAM_RANDS_PER_SITE);

  /***************************************************************************/
  /* Compute Monte Carlo averages and estimates for each observable.         */
  /***************************************************************************/
  for (int n = 0; n < NUM_OBS; n++)
  {
    MC_avg[n] = 0;
    for (long ii = 0; ii < N_BINS; ii++)
    {
      MC_avg[n] += MC_est[n][ii];
    }
    MC_avg[n] /= N_BINS;
    MC_err[n] = jackknife(MC_est[n], N_BINS, MC_avg[n]);
  }

  /***************************************************************************/
  /* Free all the memory we are finished with.						                   */
  /***************************************************************************/
  for (int n = 0; n < NUM_OBS; ++n)
  {
    delete[] MC_est[n];
    delete[] MC_bin[n];
  }
  delete[] MC_est;
  delete[] MC_bin;

  for (int ID = OMP_NUM_THREADS - 1; ID >= 0; --ID)
  {
    delete[] P_1[ID];
    delete[] P_0[ID];
  }

  for (int t = L - 1; t >= 0; --t)
  {
    for (int x = NUM_VEC_BLOCKS - 1; x >= 0; --x)
    {
      delete[] r[t][x];
    }
    delete[] r[t];
    delete[] lam[t];
    delete[] z[t];
  }
   delete[] lam;
   delete[] z;
   delete[] r;
   delete[] gen_mem;
} /* monte_carlo */

/*****************************************************************************/
/* Generate a configuration from which to start the simulation.              */
/*****************************************************************************/
void generate_start_config(VEC_DOUBLE ***z, VEC_DOUBLE ***lam, dsfmt_t ****r,
                           VEC_DOUBLE *lattice_memory, std::string path,
                           long *z_stream_counter, long *lam_stream_counter)
{
  if (!RESTART)
  {
    if (LOAD_THERM)
    {
      std::ifstream infile;
      std::stringstream lattice_file;
      lattice_file << path << "_thermalised_z_config";
      infile.open(lattice_file.str(), std::ios::in | std::ios::binary);
      if (infile.is_open())
      {
        infile.read((char *) lattice_memory, L_T * NUM_VEC_BLOCKS * SITE_WIDTH
                                                     * N_d * (sizeof(double)));
        std::cout << infile.gcount() << " bytes read." << std::endl;
      }
      else
      {
        std::cout << "Could not open file!" << std::endl;
        exit(0);
      }
    }
    /*************************************************************************/
    /* Set the initial z vector to be a random complex vector of modulus 1 & */
    /* each gauge link to be a random complex number of modulus 1. Group the */
    /* values together in multiples of the number of doubles that will fit   */
    /* in a vector register.                                                 */
    /***************************************************************************/
    else
    {
      initialize_lattice(z, lam, r, z_stream_counter, lam_stream_counter);
    }
  }
  else
  {
    std::ifstream infile;
    std::stringstream lattice_file;
    lattice_file << path << "_z_config";
    infile.open(lattice_file.str(), std::ios::in | std::ios::binary);
    if (infile.is_open())
    {
      infile.read((char *) lattice_memory, L_T * NUM_VEC_BLOCKS * SITE_WIDTH
                                                     * N_d * (sizeof(double)));
      std::cout << infile.gcount() << " bytes read." << std::endl;
    }
    else
    {
      std::cout << "Could not open file!" << std::endl;
      exit(0);
    }
  }
} /* generate_start_config */

static inline void set_z(VEC_DOUBLE ***z, VEC_DOUBLE ***lam, dsfmt_t**** r,

                         unsigned int t, unsigned int x, int max_n,
                         long *z_stream_counter, long *lam_stream_counter)
{
  double z_vec[VEC_Z_WIDTH] __attribute__((aligned(N_d * \
                                                       sizeof(double)))) = {0};
  double lam_vec[VEC_LAM_WIDTH] __attribute__((aligned(N_d * \
                                                       sizeof(double)))) = {0};
  double mag, norm, phase;

  for (int n = 0; n < max_n; ++n)
  {
    /*************************************************************************/
    /* Initialize z vector.                                                  */
    /*************************************************************************/
    norm = 0;
    for (int i = 0; i < VEC_Z_WIDTH; i += Z_STEP)
    {
      phase = 2 * M_PI * dsfmt_genrand_open_close(r[t][x][n]);
      mag = dsfmt_genrand_open_close(r[t][x][n]);
      z_vec[i + n] = mag * cos(phase);
      z_vec[i + N_d + n] = mag * sin(phase);
      norm += mag * mag;
    }
    norm = 1/sqrt(norm);
    for (int i = 0; i < VEC_Z_WIDTH; i += Z_STEP)
    {
      z_vec[i + n] *= norm;
      z_vec[i + N_d + n] *= norm;
    }

    /*************************************************************************/
    /* Initialize lambda vector.                                             */
    /*************************************************************************/
#ifdef NO_BC
    // HACKY TEST - DISABLING NO_BC
    if (t == 0)
    {
      phase = 0;
    }
    else
    {
      phase = 2 * M_PI * dsfmt_genrand_open_close(r[t][x][n]);
    }
#else
    phase = 2 * M_PI * dsfmt_genrand_open_close(r[t][x][n]);
#endif
    lam_vec[n] = cos(phase);
    lam_vec[n + N_d] = sin(phase);
#ifdef NO_BC
    // HACKY TEST - DISABLING NO_BC
    if (t == 0)
    {
      phase = 0;
    }
    else
    {
      phase = 2 * M_PI * dsfmt_genrand_open_close(r[t][x][n]);
    }
#else
    phase = 2 * M_PI * dsfmt_genrand_open_close(r[t][x][n]);
#endif
    lam_vec[n + 2 * N_d] = cos(phase);
    lam_vec[n + 3 * N_d] = sin(phase);
  }

  /***************************************************************************/
  /* Copy the vectors across into the lattice memory space.                  */
  /***************************************************************************/
  for (int i = 0, j = 0; j < VEC_Z_WIDTH; ++i, j += N_d)
  {
#if defined (__SSE2__)
    z[t][x][i] = VEC_DOUBLE().load_a(&z_vec[j]);
#else
    z[t][x][i] = z_vec[j];
#endif
  }
  for (int i = 0, j = 0; j < VEC_LAM_WIDTH; ++i, j += N_d)
  {
#if defined (__SSE2__)
    lam[t][x][i] = VEC_DOUBLE().load_a(&lam_vec[j]);
#else
    lam[t][x][i] = lam_vec[j];
#endif
  }
} /* set_z */

/*****************************************************************************/
/* Set the initial z vector to be a random complex vector of modulus 1 and   */
/* each gauge link to be a random complex number of modulus 1. Group the     */
/* values together in multiples of the number of doubles that will fit in    */
/* a vector register.                                                        */
/*****************************************************************************/
static inline void initialize_lattice(VEC_DOUBLE ***z, VEC_DOUBLE ***lam,
                                      dsfmt_t**** r, long *z_stream_counter,
                                      long *lam_stream_counter)
{
  for (unsigned int t = 0; t < L_T; ++t)
  {
    /*************************************************************************/
    /* Set all the entries on this row.                                      */
    /*************************************************************************/
    for (unsigned int x = 0; x < NUM_VEC_BLOCKS; ++x)
    {
      set_z(z, lam, r, t, x, N_d, z_stream_counter, lam_stream_counter);
    }
  }
} /* initialize_lattice */


static inline void initialize_counters(long **z_stream_counter,
                                       long **lam_stream_counter,
                                       std::string path)
{
  /***************************************************************************/
  /* Allocate memory for counters.                                           */
  /***************************************************************************/
  *z_stream_counter = new long[2 * L_T * NUM_VEC_BLOCKS * N_d];

  /***************************************************************************/
  /* Initialize counters to track how many random numbers we have generated. */
  /* If we are restarting from a saved state, load these counters from a     */
  /* file.                                                                   */
  /***************************************************************************/
  if (RESTART != 0)
  {
    std::ifstream infile;
    std::stringstream counter_file;
    counter_file << path << "_stream_counters";
    infile.open(counter_file.str(), std::ios::in | std::ios::binary);
    if (infile.is_open())
    {
      infile.read((char *) (*z_stream_counter),
                  2 * L_T * NUM_VEC_BLOCKS * N_d * (sizeof(long)));
    }
  }

  /***************************************************************************/
  /* If this is a new run, zero all the counters.                            */
  /***************************************************************************/
  else
  {
    for (long index = 0; index < (2 * V_TOT); ++index)
    {
      (*z_stream_counter)[index] = 0;
    }
  }

  *lam_stream_counter = &((*z_stream_counter)[V_TOT]);
} /* initialize_counters */


/*****************************************************************************/
/* Utility function to jump ahead in RNG stream.                             */
/*****************************************************************************/
// static inline void jump_ahead_by_long(dsfmt_t *r, long stream_length)
// {
  // /***************************************************************************/
  // /* Initialise variables.                                                   */
  // /***************************************************************************/
  // NTL::GF2X lcmpoly;
  // read_file(lcmpoly, 0, "poly.1279.txt");
  // NTL::ZZ step;
  // std::stringstream ss;
  // ss << stream_length;
  // ss >> step;

  // /***************************************************************************/
  // /* Calculate the input long as a polynomial jump string.                   */
  // /***************************************************************************/
  // std::string jump_str;
  // calc_jump(jump_str, step, lcmpoly);

  // /***************************************************************************/
  // /* Jump ahead by the specified length.                                     */
  // /***************************************************************************/
  // dSFMT_jump(r, jump_str);

  // return;
// }

/*****************************************************************************/
/* Initialize the random number generator.                                   */
/*****************************************************************************/
static inline void initialize_rng(dsfmt_t***** r, dsfmt_t** gen_mem,
                                  unsigned long seed,
                                  long *z_stream_counter,
                                  long *lam_stream_counter,
                                  std::string path)
{
  static const char * jumppoly = "2c6b058dca1fbfb57ebf41e67fec066c8828f2bf"
      "9414331d2767fa740aba89987685a79114b3543edbc83476e35fd1e52b8b2436932"
      "b9d18a728d8c7d7009a31aabed9bf4646909b8138f3e2a05e611c48dd1ce58a4618"
      "3adabf3314da38599af92720efca1535872e7f85ef916b2c1e41dfe8ea764730f6b"
      "a2654ab287a55214bcf08cfae416906e4979108d606819b5d9e5b2f11ce028577f6"
      "c3788a9c688f1d64f5ae341eb95169824954";
  *gen_mem = new dsfmt_t[(L_T / 2) * NUM_VEC_BLOCKS * N_d];
  *r = new dsfmt_t***[L_T];
  long gen_index = 0;
  long stream_index = 0;

  /***************************************************************************/
  /* If restarting a run, load the RNG state from a file.                    */
  /***************************************************************************/
  if (RESTART != 0)
  {
    std::ifstream infile;
    std::stringstream rng_file;
    rng_file << path << "_rng_state";
    infile.open(rng_file.str(), std::ios::in | std::ios::binary);
    if (infile.is_open())
    {
      infile.read((char *) (*gen_mem), (L_T / 2) * NUM_VEC_BLOCKS * N_d *
                                                            (sizeof(dsfmt_t)));
    }
  }

  /***************************************************************************/
  /* Else if we are starting a new run, initialize the RNG.                  */
  /***************************************************************************/
  else
  {
    for (gen_index = 0; gen_index < (V_TOT / 2); ++gen_index)
    {
      if (gen_index == 0)
      {
        dsfmt_init_gen_rand(&(*gen_mem)[gen_index], seed);
      }
      else
      {
        (*gen_mem)[gen_index] = (*gen_mem)[gen_index - 1];
        dSFMT_jump(&(*gen_mem)[gen_index], jumppoly);
      }
    }
  }

  /***************************************************************************/
  /* Construct the pointer array that points to the RNG at each site.        */
  /***************************************************************************/
  gen_index = 0;
  for (unsigned int t = 0; t < L_T; ++t)
  {
    (*r)[t] = new dsfmt_t**[NUM_VEC_BLOCKS];
    for(unsigned int x = 0; x < NUM_VEC_BLOCKS; x += 2)
    {
      (*r)[t][x] = new dsfmt_t*[N_d];
      (*r)[t][x + 1] = new dsfmt_t*[N_d];

      /***********************************************************************/
      /* As pairs of red and black sites are always updated by the same      */
      /* thread, they may share a random number generator.                   */
      /***********************************************************************/
      for (int n = 0; n < N_d; ++n)
      {
        (*r)[t][x][n] = &(*gen_mem)[gen_index];
        (*r)[t][x + 1][n] = &(*gen_mem)[gen_index];
        ++gen_index;
      }
    }
  }
}


/*****************************************************************************/
/* Construct an array of pointers to access the lattice memory.              */
/*****************************************************************************/
static inline void construct_lattice(VEC_DOUBLE**** z, VEC_DOUBLE**** lam,
                                     VEC_DOUBLE* lattice_memory)
{
  *z = new VEC_DOUBLE**[L_T];
  *lam = new VEC_DOUBLE**[L_T];

  for (unsigned int t = 0; t < L_T; ++t)
  {
    (*z)[t] = new VEC_DOUBLE*[L_T];
    (*lam)[t] = new VEC_DOUBLE*[L_T];
    for (unsigned int x = 0; x < NUM_VEC_BLOCKS; ++x)
    {
      long z_index = ((t * NUM_VEC_BLOCKS) + x) * SITE_WIDTH;
      long lam_index = z_index + 2 * N;
      (*z)[t][x] = &lattice_memory[z_index];
      (*lam)[t][x] = &lattice_memory[lam_index];
    }
  }
}


/*****************************************************************************/
/* Define the limits for the blocks that each thread will update.            */
/*****************************************************************************/
static inline void define_update_limits(update_limits *limits)
{
#ifdef RECTANGLE
  unsigned int spare_counter_rect(SPARE_RECT);
  unsigned int rect_min_x(0);
  unsigned int rect_max_x(0);
  unsigned int rect_min_t(0);
  unsigned int rect_max_t(0);
  long tot_updates_rect = 0;
#endif
  unsigned int spare_counter(SPARE);
  unsigned int min_x(0);
  unsigned int max_x(0);
  unsigned int min_t(0);
  unsigned int max_t(0);
  long tot_meas_sites = 0;

  for (int ID = 0; ID < OMP_NUM_THREADS; ++ID)
  {
    /*************************************************************************/
    /* Define how many sites this thread will measure.                       */
    /*************************************************************************/
    long num_block_meas = NUM_BLOCK_UPS;
    if (spare_counter > 0)
    {
      num_block_meas++;
      spare_counter--;
    }

    /*************************************************************************/
    /* Define how many sites will be updated overall after this thread's     */
    /* contribution. Use this information to calculate the start and end     */
    /* coordinates of our 2D lattice.                                        */
    /*************************************************************************/
    tot_meas_sites += num_block_meas;
    max_t = tot_meas_sites / NUM_RB_BLOCKS;
    max_x = 2 * (tot_meas_sites % NUM_RB_BLOCKS);

    if (max_x == 0)
    {
      max_x = NUM_VEC_BLOCKS;
    }
    else
    {
      max_t += 1;
    }
    limits[ID].m_min_t = min_t;
    limits[ID].m_max_t = max_t;
    limits[ID].m_min_x = min_x;
    limits[ID].m_max_x = max_x;
    limits[ID].num_block_meas = 2 * num_block_meas;

    if (max_x == NUM_VEC_BLOCKS)
    {
      min_x = 0;
    }
    else
    {
      min_x = max_x;
      max_t -= 1;
    }
    min_t = max_t;

    /*************************************************************************/
    /* If we are using a longer extent on one edge, we must define update    */
    /* limits within an extended box.                                        */
    /*************************************************************************/
    #ifdef RECTANGLE
      /***********************************************************************/
      /* Define how many sites this thread will measure.                    */
      /***********************************************************************/
      long num_block_ups = NUM_BLOCK_UPS_RECT;
      if (spare_counter_rect > 0)
      {
        num_block_ups++;
        spare_counter_rect--;
      }

      /***********************************************************************/
      /* Define how many sites will be used for measurement by this thread's */
      /* contribution. Use this information to calculate the start and end   */
      /* coordinates of our 2D lattice.                                      */
      /***********************************************************************/
      tot_updates_rect += num_block_ups;
      rect_max_t = (tot_updates_rect / NUM_RB_BLOCKS);
      rect_max_x = 2 * (tot_updates_rect % NUM_RB_BLOCKS);

      if (rect_max_x == 0)
      {
        rect_max_x = NUM_VEC_BLOCKS;
      }
      else
      {
        rect_max_t += 1;
      }
      limits[ID].u_min_t = rect_min_t;
      limits[ID].u_max_t = rect_max_t;
      limits[ID].u_min_x = rect_min_x;
      limits[ID].u_max_x = rect_max_x;
      limits[ID].num_block_ups = 2 * num_block_ups;
      if (rect_max_x == NUM_VEC_BLOCKS)
      {
        rect_min_x = 0;
      }
      else
      {
        rect_min_x = rect_max_x;
        rect_max_t -= 1;
      }
      rect_min_t = rect_max_t;

    /*************************************************************************/
    /* If we are not using open boundary conditions, the update and          */
    /* measurement limits should be the same.                                */
    /*************************************************************************/
    #else
      limits[ID].u_min_t = limits[ID].m_min_t;
      limits[ID].u_max_t = limits[ID].m_max_t;
      limits[ID].u_min_x = limits[ID].m_min_x;
      limits[ID].u_max_x = limits[ID].m_max_x;
      limits[ID].num_block_ups = limits[ID].num_block_meas;
    #endif
    /*std::cout << "Thread " << ID << " update limits." << std::endl;
    std::cout << limits[ID].u_min_t << std::endl;
    std::cout << limits[ID].u_max_t << std::endl;
    std::cout << limits[ID].u_min_x << std::endl;
    std::cout << limits[ID].u_max_x << std::endl;
    std::cout << "Thread " << ID << " measurement limits." << std::endl;
    std::cout << limits[ID].m_min_t << std::endl;
    std::cout << limits[ID].m_max_t << std::endl;
    std::cout << limits[ID].m_min_x << std::endl;
    std::cout << limits[ID].m_max_x << std::endl;
    std::cout << std::endl;*/
  }
} /* define_update_limits */

static inline double run_simulation(VEC_DOUBLE ***z, VEC_DOUBLE ***lam,
                                    dsfmt_t**** r, update_limits* limits,
                                    dsfmt_t* gen_mem,
                                    long * z_stream_counter,
                                    long * lam_stream_counter,
                                    VEC_DOUBLE *lattice_memory,
                                    std::string path,
                                    double **P_0, double **P_1,
                                    double **MC_est, std::string **obs_list)
{
  /***************************************************************************/
  /* Initialize variables.                                                   */
  /***************************************************************************/
  long bin_num = 0;
  long FLOPs;
  double total_FLOPs;
  double meas_FLOPs;
  double meas_FLOPS(0);

  /***************************************************************************/
  /* Allocate memory for the bin that we will put measurements into.         */
  /***************************************************************************/
  double **MC_bin;
  MC_bin = new double*[NUM_OBS];
  long bin_size = N_CF / N_BINS;
  for (int n = 0; n < NUM_OBS; ++n)
  {
    MC_bin[n] = new double[bin_size];
  }

  /***************************************************************************/
  /* Allocate memory for the correlator and for the time averaged composite  */
  /* operator structure used to compute it.                                  */
  /***************************************************************************/
  double *correlator;
  long corr_bin = bin_size / CORR_N_COR;
  correlator = new double[CORR_LEN * corr_bin];

  void *angle_mem;
  int rc = posix_memalign(&angle_mem, N_d * sizeof(double),
                          2 * V * bin_size * sizeof(double));
  if (rc != 0)
  {
    std::cout << "Could not allocate memory\n" << std::endl;
    exit(0);
  }
  VEC_DOUBLE *angles = (VEC_DOUBLE *)angle_mem;

  void *P_op_mem;
  rc = posix_memalign(&P_op_mem, N_d * sizeof(double),
                                                   L * P_LEN * sizeof(double));
  if (rc != 0)
  {
    std::cout << "Could not allocate memory\n" << std::endl;
    exit(0);
  }

  double *P_op_mem_ptr = (double *)P_op_mem;
/*#ifdef NO_BC
  VEC_DOUBLE **P_t_avg;
  P_t_avg = new VEC_DOUBLE*[NUM_VEC_BLOCKS];
  for (unsigned int x = 0; x < NUM_VEC_BLOCKS; ++x)
  {
    P_t_avg[x] = (VEC_DOUBLE *)(P_op_mem_ptr + x * N_d * P_LEN);
  }
#else*/
  double **P_t_avg;
  P_t_avg = new double*[L];
  for (unsigned int x = 0; x < L; ++x)
  {
    P_t_avg[x] = P_op_mem_ptr + x * P_LEN;
  }
//#endif
  bool save_configs = false;
  long first_save_config = 0;

  /***************************************************************************/
  /* Update the lattice configuration and extract measurements. Repeat       */
  /* N_cf times.                                                             */
  /***************************************************************************/
  for(long ii = 0; ii < N_CF; ++ii)
  {
    for(long jj = 0; jj < N_COR; ++jj)
    {
      FLOPs = heatbath_update(z, lam, r, limits, z_stream_counter,
                              lam_stream_counter);
      total_FLOPs += (double)FLOPs;
    }

    /*************************************************************************/
    /* Take measurements within the feducial box of side length L. For       */
    /* periodic BCs this is the full lattice.                                */
    /*************************************************************************/
    auto t1 = std::chrono::high_resolution_clock::now();
    // HACKY TEST, [L_GAP] -> [0]
    meas_FLOPs = comp_MC_ests(&z[L_GAP], &lam[L_GAP], P_0, P_1, correlator,
                              P_t_avg, MC_bin, bin_size, ii, limits,
                              angles);
    auto t2 = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()*1e-9;
    meas_FLOPS += meas_FLOPs / time;
    total_FLOPs += meas_FLOPs;

    /*************************************************************************/
    /* If we've reached the end of a bin, average the estimates and save     */
    /* the configurations.                                                   */
    /*************************************************************************/
    if ((ii + 1) % bin_size == 0)
    {
      for (int n = 0; n < NUM_OBS; ++n)
      {
        /*********************************************************************/
        /* Save off data for this bin.                                       */
        /*********************************************************************/
        if (SAVE_MEASUREMENTS != 0)
        {
          /*******************************************************************/
          /* Save off the measurements.                                      */
          /*******************************************************************/
          std::string filename = path + "_" + *(obs_list[n]);
          std::ofstream outfile;
          outfile.open(filename, std::ios::out | std::ios::app | std::ios::binary);
          if (outfile.is_open())
          {
            outfile.write(reinterpret_cast<const char*>(MC_bin[n]),
                          std::streamsize(bin_size * sizeof(double)));
            outfile.close();
          }
        }

        MC_est[n][bin_num] = 0;
        for (long k = 0; k < bin_size; k++)
        {
          MC_est[n][bin_num] += MC_bin[n][k];
        }
        MC_est[n][bin_num] /= bin_size;
      }

      /***********************************************************************/
      /* Save off the correlator if specified.                               */
      /***********************************************************************/
      #if (MEASUREMENTS & WALL_CORR)
        if (SAVE_MEASUREMENTS != 0)
        {
          /*******************************************************************/
          /* Save off the correlator.                                        */
          /*******************************************************************/
          std::string filename = path + "_correlator_alt";
          std::ofstream outfile;
          outfile.open(filename, std::ios::out | std::ios::app | std::ios::binary);
          if (outfile.is_open())
          {
            outfile.write(reinterpret_cast<const char*>(correlator),
                          std::streamsize(CORR_LEN * corr_bin *
                                                              sizeof(double)));
            outfile.close();
          }
        }
        /*
        if (SAVE_MEASUREMENTS != 0)
        {
          /*******************************************************************/
          /* Save off the correlator.                                        */
          /*******************************************************************/
          /*std::string filename = path + "_angles";
          std::ofstream outfile;
          outfile.open(filename, std::ios::out | std::ios::app | std::ios::binary);
          if (outfile.is_open())
          {
            outfile.write(reinterpret_cast<const char*>(angles),
                          std::streamsize(bin_size * 2 * V * sizeof(double)));
            outfile.close();
          }
        }*/
      #endif

      /***********************************************************************/
      /* Save off the latest lattice configuration and RNG state.            */
      /***********************************************************************/
      if (SAVE_MEASUREMENTS != 0)
      {
        std::ofstream outfile;
        std::stringstream lattice_file;
        lattice_file << path << "_z_config";
        outfile.open(lattice_file.str(), std::ios::out | std::ios::binary);
        if (outfile.is_open())
        {
          outfile.write(reinterpret_cast<const char*>(lattice_memory),
                        std::streamsize(N_d * L_T * NUM_VEC_BLOCKS *
                                        SITE_WIDTH * sizeof(double)));
          outfile.close();
        }

        std::stringstream stream_file;
        stream_file << path << "_rng_state";
        outfile.open(stream_file.str(), std::ios::out | std::ios::binary);
        if (outfile.is_open())
        {
          outfile.write(reinterpret_cast<const char*>(gen_mem),
                        std::streamsize(N_d * (L_T / 2) * NUM_VEC_BLOCKS *
                                                             sizeof(dsfmt_t)));
          outfile.close();
        }
      }
      bin_num++;
    }
  }

  /***************************************************************************/
  /* Free the memory we are finished with.                                   */
  /***************************************************************************/
  delete[] P_t_avg;
  free(P_op_mem);
  delete[] correlator;
  for (int n = 0; n < NUM_OBS; ++n)
  {
    delete[] MC_bin[n];
  }
  delete[] MC_bin;

  std::cout << "Measurement GFLOPS = " << meas_FLOPS*1e-9/N_CF << std::endl;
  return total_FLOPs;
}
