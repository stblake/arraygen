/***************************************
 *
 *         Array Generator
 *
 ***************************************/

/* Written by Sam Blake */

/* Started on 25 Feb 2014 */

/* Version 1.0  */

/* Version 1.0 frozen on 2 April 2015 */

/* 

Notes from SAMB in 0722: 

This program was originally designed to be a one-stop-shop for 
constructions of sequences and arrays with good autocorrelation 
and pairwise cross correlation. However, we (Andrew and I) quickly 
realised that constructions alone are not a great commercial offering 
and this program morphed into a research platform. 

This program has been, more or less, unused since around the end of 2015. 

Subsequent to the creation of this program, some of these constructions 
were patented by Optimark in the US and should not be used without their
permission. However, subsequent constructions outperform these patented 
constructions making their use obsolete. 
*/

#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h> 
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <pwd.h>
#include <unistd.h> 
#include <termios.h>
/* #include <crypt.h> */ /* Not needed on mac OS */


typedef int bool;
#define TRUE 1
#define FALSE 0

#define DEBUG 1

#define CHECK_EXPIRY 1
#define CHECK_PASSWORD 0

#define VERSION "1.0"

#define EXPIRY_YEAR  2030
#define EXPIRY_MONTH 12
#define EXPIRY_DAY   31

#define NEGINF -9999

#ifndef MAX_LU_PRIME
#define MAX_LU_PRIME 0
#endif

#define MAX_NUMBER_OF_PRIME_FACTORS 64

#define MEMCHECK(p) \
  if (p == NULL)    \
    {\
      printf("\n\nERROR: Out of memory.\n\n");\
      exit(1);\
    }

#define INTERNAL_ERROR(k) printf("\nINTERNAL ERROR %i: Please contact Sam Blake or Andrew Tirkel.\n\n", k)

typedef struct polynomial_struct polynomial;

/* We use a dense polynomial representation. */
struct polynomial_struct 
  {
    int degree;
    int *coefficients; // Coefficient list with constant term first.
  };

#define poly_deg(e) e.degree
#define poly_list(e) e.coefficients
#define poly_free(e) free(e.coefficients)

typedef struct bivariate_polynomial_struct bivariate_polynomial;

/* 
   Notes on bivariate_polynomial_struct: Stores a bivariate polynomial using a 
   dense representation such that a dot product with:
     
   1       y     y^2 ...     y^degY
   x     x*y   x*y^2 ...   x*y^degY
   x^2 x^2*y x^2*y^2 ... x^2*y^degY
   ...                        .
   .                          .
   .                          .
   x^degX x^degX*y x^degX*y^2 ... x^degX*y^degY

   gives the bivariate polynomial of total degree degX + degY. */

struct bivariate_polynomial_struct
  {
    int degX;
    int degY;
    int **coeff_array;
};

#define poly_degX(e) e.degX
#define poly_degY(e) e.degY
#define poly_array(e) e.coeff_array

typedef struct galois_field gf;

struct galois_field 
  {
    int p; // p is the prime characteristic
    int d; // Galois field with p^d elements (Depreciated. Should extract from poly_deg(prim_poly))
    polynomial prim_poly; // primitive polynomial
    polynomial poly; // polynomial in the Galois field
  };

#define gf_prime(e) e.p
#define gf_deg(e) e.prim_poly.degree

#define gf_prim_poly(e) e.prim_poly
#define gf_prim_poly_list(e) e.prim_poly.coefficients

#define gf_poly(e) e.poly
#define gf_poly_list(e) e.poly.coefficients
#define gf_poly_deg(e) e.poly.degree

void run_tests();
void display_help();
void readpass(char password[128]);
bool checkexpirydate();
bool gf_equal(gf a, gf b);
bool gf_field_equal(gf a, gf b);
gf gf_create(int p, int degree, int *polycoeffs, int *elems);
gf gf_copy(gf a);
void gf_print(gf a);
void gf_delete(gf a);
gf gf_plus(gf a, gf b);
gf gf_times(gf a, gf b);
bool primep(int n);
polynomial poly_mult(polynomial a, polynomial b);
polynomial poly_rem(polynomial u, polynomial v, int mod);
int mult_inv(int a, int b);
void gf_power_table(int prime, polynomial primpoly);
gf gf_power(gf a, int n);
void print_array_1D(int *seq, int d0);
int* sidelnikov1D(int p, int d, polynomial primpoly, polynomial perm);
bool zeroq(polynomial l);
int power_representation(polynomial l);
void free_power_table();
polynomial gen_rand_poly(int p, int d);
int poly_print(polynomial poly);
int random_in_range(unsigned int min, unsigned int max);
void print_array_2D(int **array, int d0, int d1);
polynomial gen_rand_primitive_poly(int p, int d);
int** legendre2Dgeneralised(int p, polynomial primpoly, polynomial numpoly, 
			    polynomial denpoly, bool generalised);
void print_array_3D(int ***array, int d0, int d1, int d2);
int*** legendre3Dgeneralised(int p, polynomial primpoly, polynomial numpoly, 
			     polynomial denpoly, bool generalised);
bool gf_equal_to_1(gf a);
void write_array_1D(char *filename, bool csvp, int *array, int d0);
void write_array_2D(char *filename, bool csvp, int **array, int d0, int d1);
void write_array_3D(char *filename, bool csvp, int ***array, int d0, int d1, int d2);
void binarize1D(int *array, int d0);
void binarize2D(int **array, int d0, int d1);
void binarize3D(int ***array, int d0, int d1, int d2);
void construction_Sidelnikov1D(char *filename, int p, polynomial primpoly, bool binarize, bool acv, 
			       bool pipe, bool binaryfile, bool verbose, bool csvp);
void construction_Legendre2D(char *filename, int p, polynomial numpoly, polynomial denpoly, polynomial primpoly, 
			     bool binarize, bool acv, bool pipe, bool binaryfile, bool verbose, bool csvp);
void construction_Legendre3D(char *filename, int p, polynomial numpoly, polynomial denpoly, polynomial primpoly, 
			     bool binarize, bool acv, bool pipe, bool binaryfile, bool zip, bool verbose, bool csvp);
int*** D2_array_3D(int p, polynomial prim_poly, polynomial numpoly);
void construction_d2(char *filename, int p, polynomial primpoly, polynomial numpoly, 
		     bool binarize, bool acv, bool pipe, bool binaryfile, bool zip, bool verbose, bool csvp);
void construction_d2_family(char *filename, int p, polynomial primpoly, int shift_poly_degree, 
			    bool binarize, bool acv, bool pipe, bool binaryfile, bool zip, bool verbose, bool csvp);
polynomial gf_poly_eval(int p, int k, polynomial poly, polynomial primpoly);
polynomial poly_copy(polynomial poly);
void construction_c1(char *filename, int p, polynomial primpoly, polynomial numpoly, 
		     bool binarize, bool acv, bool pipe, bool binaryfile, bool zip, bool verbose, bool csvp, 
		     bool lincomp);
int* acv_1d_slow(int *array, int d0);
int** acv_2d_slow(int **array, int d0, int d1);
int*** acv_3d_slow(int ***array, int d0, int d1, int d2);
void pretty_print_array_2D(int **array, int d0, int d1);
int intlength(int k);
void print_power_table(int p, polynomial primpoly);
void pretty_print_array_3D(int ***array, int d0, int d1, int d2);
void increment(int *counter, int len, int b);
void reverse(int arr[], int start, int end);
void NextPermutation(int *array, int length);
void swap(int *i, int *j);
void write_binary_array_1D(char *filename, int *array, int d0);
void write_binary_array_2D(char *filename, int **array, int d0, int d1);
void write_binary_array_3D(char *filename, int ***array, int d0, int d1, int d2);
int* legendre1D(int p);
int* quadratic_residues(int p);
void construction_Legendre1D(char *filename, int p, bool acv, bool pipe, bool binary, bool verbose, bool csvp);
int power_mod(int base, int exponent, int mod);
int poly_eval(polynomial poly, int k, int mod);
int** Tirkel_Hall_2D(int p, polynomial shiftpoly, int map);
void construction_Tirkel_Hall2D(char *filename, int p, int zeromap, polynomial shiftpoly, 
				bool acv, bool pipe, bool binaryfile, bool verbose, bool csvp);
void construction_Tirkel_Hall2D_family(char *filename, int p, int shift_poly_degree, bool acv,
				       int zeromap, bool pipe, bool binaryfile, bool verbose, bool csvp);
void generate_power_table(int prime, polynomial poly);
void io_1d(char *filename, int *array, int d0, bool csvp, bool acv, bool pipe, bool binaryfile, bool verbose);
void io_2d(char *filename, int **array, int d0, int d1, bool csvp, bool acv, bool pipe, bool binaryfile, bool verbose);
void io_3d(char *filename, int ***array, int d0, int d1, int d2, bool csvp, bool binarize, bool acv, bool pipe, 
	   bool binaryfile, bool zip, bool verbose, bool lincomp);
bool system_solver_2D(int p, int c3, int c5, int c6);
int poly_eval_2D(bivariate_polynomial poly, int x, int y, int mod);
int*** Blake_Tirkel_3D(int p, bivariate_polynomial shiftpoly, int map);
int power(int base, int exp);
void construction_Blake_Tirkel3D(char *filename, int p, int zeromap, bivariate_polynomial shiftpoly, 
				 bool acv, bool pipe, bool binaryfile, bool zip, bool verbose, bool csvp);
void poly_print_2D(bivariate_polynomial poly);
bool already_near_binary_3D(int ***array, int d0, int d1, int d2);
bivariate_polynomial  gen_special_rand_poly_2D(int p);
int gcd(int u, int v);
void lower(char *s);
void construction_Blake_Tirkel3D_special_family(char *filename, int p, int zeromap, bool pipe, bool acv, 
						bool binaryfile, bool zip, bool verbose, bool csvp);
bivariate_polynomial gen_rand_poly_2D(int p, int polydegX, int polydegY);
void construction_d3(char *filename, int p, polynomial primpoly, polynomial numpoly, 
		     polynomial denpoly, bool binarize, bool acv, bool pipe, bool binaryfile, 
		     bool zip, bool verbose, bool csvp);
polynomial gf_rational_eval(int p, int K, polynomial numpoly, polynomial denpoly, polynomial primpoly);
int** legendre2D(int p, polynomial primpoly);
int*** legendre3D(int p, polynomial primpoly);
int*** D3_array_3D(int p, polynomial primpoly, polynomial numpoly, polynomial denpoly);
void construction_c2(char *filename, int p, polynomial primpoly, polynomial numpoly, 
		     polynomial denpoly, bool binarize, bool acv, bool pipe, bool binaryfile, 
		     bool zip, bool verbose, bool csvp);
int*** C1_array_3D(int p, polynomial prim_poly, polynomial shift_poly);
int*** C2_array_3D(int p, polynomial prim_poly, polynomial numpoly, polynomial denpoly);
bool neginfq(polynomial p);
int* berlenkamp_massey(int *seq, int N);
void copy_vec(int *src, int srcpos, int *dest, int destpos, int len);
int lc(int *seq, int n);
bool coeffs_equal(int *c1, int *c2, int len);
void construction_legendre_projection_2D(char *filename, int p, bool acv, bool pipe, 
					 bool binaryfile, bool zip, bool verbose, bool csvp);
int** legendre_projection_2D(int p, int shift);
void euler_phi_disp(int n);
int euler_phi(int n);
void primitive_root_disp(int n);
int primitive_root(int n);
void factor_integer_disp(int n);
int* factor_integer(int n);
int** A_array_2D(int p, polynomial shiftpoly, int map);
void construction_a(char *filename, int p, int zeromap, bool enum_diag, polynomial shiftpoly, 
		    bool acv, bool pipe, bool binary, bool verbose, bool csv);
void construction_a_family(char *filename, int p, int shift_poly_degree, int zeromap,  
			   bool enum_diag, bool acv, bool pipe, bool binaryfile, 
			   bool verbose, bool csv);
int* CRTenumerate2D(int **array, int d0, int d1);
void primitive_polynomial_test(int p, int d, int n);
int* m_sequence(polynomial primpoly, int p);
void construction_m_sequence(char *filename, int p, polynomial primpoly,
			    bool acv, bool pipe, bool binaryfile, 
			     bool zip, bool verbose, bool csv);
int ipatov_psi(int a, int p);


int malloc_count = 0; // Used for tracking memory usage.
int gf_alloc_count = 0; // Total number of GF expression allocations. 

int **power_table;
bool power_table_init = FALSE;
int power_table_length = 0;
int power_table_degree = 0;
bool primitive_polynomial_check = TRUE;

char *default_filename = "arraygen";
char *default_txt_extension = ".txt";
char *default_bin_extension = ".bin";

#define SIDELNIKOV_1D_FLAG 5
#define LEGENDRE_1D_FLAG 11
#define LEGENDRE_2D_FLAG 6
#define LEGENDRE_3D_FLAG 7
#define A_FLAG 24
#define C1_FLAG 3
#define C2_FLAG 4
#define D2_FLAG 1
#define D3_FLAG 2
#define TIRKEL_HALL_FLAG 13
#define BLAKE_TIRKEL_FLAG 17
#define LEGENDREPROJECTION_2D_FLAG 19
#define MSEQUENCE_FLAG 27

int main(int argc, char **argv) {
  char pwdin[128], filename[128], *construction_name;
#if CHECK_PASSWORD
  char *pwdcheck;
  const char *const pass = "$1iz8QbXkEakU"; // white@ssnow
  int pwdokay;
#endif
  const int Sidelnikov1D_flag =  SIDELNIKOV_1D_FLAG, 
    Legendre1D_flag = LEGENDRE_1D_FLAG,
    Legendre2D_flag = LEGENDRE_2D_FLAG, 
    Legendre3D_flag = LEGENDRE_3D_FLAG, 
    LegendreProjection2D_flag = LEGENDREPROJECTION_2D_FLAG,
    A_flag = A_FLAG,
    D2_flag = D2_FLAG, 
    D3_flag = D3_FLAG, 
    C1_flag = C1_FLAG, 
    C2_flag = C2_FLAG,
    Tirkel_Hall_flag = TIRKEL_HALL_FLAG,
    Blake_Tirkel_flag = BLAKE_TIRKEL_FLAG,
    Msequence_flag = MSEQUENCE_FLAG;
  int k, n, p = -1, array_flag = -1, polydegree = -1,
    numpolydegree = -1, denpolydegree = -1, polydegX = -1, 
    polydegY = -1, **shift_poly_coeffs_2D, zeromap = 0,
    x, y, den_coeffs[1] = {1};
  double start, stop, elapsed;
  bool expired, autofilename = TRUE, genfamily = FALSE, 
    verbose = FALSE, binary = FALSE, polyfromuser = FALSE, binarize = FALSE,
    pipe = FALSE, randomise = FALSE, numpolyfromuser = FALSE, 
    denpolyfromuser = FALSE, acv = FALSE, xcv = FALSE, disp_pow_table = FALSE,
    zip = FALSE, special = FALSE, random_prim_poly_generated = FALSE,
    poly_2D_specified = FALSE, csv = FALSE, linear_complexity = FALSE,
    enum_diagonals = FALSE;
  polynomial poly, numpoly, denpoly;
  bivariate_polynomial shiftpoly2D;

  
  /* Set random seed. */
  srand(time(NULL));

  /* Incase these are not initialised. */
  poly_deg(poly) = 0;
  poly_deg(numpoly) = 0;
  poly_deg(denpoly) = 0;
  poly_list(denpoly) = den_coeffs; // By default denominator poly is 1. 

  memset(pwdin, '\0', sizeof(pwdin));
  memset(filename, '\0', sizeof(filename));

  /* Convert command line arguments to lower case. */
  for (n = 1; n < argc; n++) lower(argv[n]);

  /* Options requiring a higher priority. */
  for (n = 1; n < argc; n++) {
    if (strcmp(argv[n], "-version") == 0) {
      printf("Version %s. Copyright Sam Blake and Andrew Tirkel, 2015.\n", VERSION);
      return 0;
    }
    else if (strcmp(argv[n], "-help") == 0 || strcmp(argv[n], "-h") == 0) {
      display_help(); 
      return 0;
      }
    else if (strcmp(argv[n], "-pipe") == 0) 
      pipe = TRUE;
  }

  if (!pipe) {
    printf("\n               * * * * * * * * * * * * * * * * * * * * * *   \n");
    printf(  "             * *                                       * *   \n");
    printf(  "           *   *                                      *  *   \n");
    printf(  "         *     * * * * * * * * * * * * * * * * * * * * * *   \n");
    printf(  "       *     *                                      *   *    \n");
    printf(  "     * * * * * * * * * * * * * * * * * * * * * * * *   *     \n");
    printf(  "     *   *            -- ARRAYGEN --               *  *      \n");
    printf(  "     * *             array generator               * *       \n");
    printf(  "     * * * * * * * * * * * * * * * * * * * * * * * *        \n\n");
    printf(  "Version %s. Copyright Sam Blake and Andrew Tirkel, 2015.\n", VERSION);
  }

#if CHECK_EXPIRY
  /* Check expiry date. */
  expired = checkexpirydate();
  if (expired) return 0;
#endif

#if CHECK_PASSWORD
  /* Check users password. */
  readpass(pwdin);
  /* printf("%s\n\n", pwdin); */
  pwdcheck = crypt(pwdin, pass);
  /* printf("crypt of pwd is %s\n", pwdcheck); */
  pwdokay = strcmp(pwdcheck, pass) == 0;

  if (! pwdokay) 
  {
    printf("The password is incorrect. Please contact Sam Blake or Andrew Tirkel.\n\n");
    return 0;
  }
#else
  if (!pipe) printf("\n    *** UNPROTECTED VERSION  ***\n");
#endif
  
  /* Set default file name. */
  strcpy(filename, default_filename); 
  strcat(filename, default_txt_extension); // Default is to output text file.

  /* Process command line inputs. */
  for (n = 1; n < argc; n++)  /* Skip argv[0] (program name). */
  {
    if ( strcmp(argv[n], "-help") == 0 )
      continue;
    // Pipe.
    else if (strcmp(argv[n], "-pipe") == 0) 
      continue;
    // Verbose displays information as the program executes.
    else if ( strcmp(argv[n], "-verbose") == 0 || strcmp(argv[n], "-v") == 0 ) 
      {
	// -pipe overrides -verbose
	if (pipe) 
	  verbose = FALSE;
	else
	  verbose = TRUE;
      }
    // Randomize (or Randomise) generated a primitive polynomial at random. 
    else if ( strcmp(argv[n], "-randomise") == 0 || strcmp(argv[n], "-randomize") == 0 )
      {
	randomise = TRUE;
      }
    // Compute the linear complexity in the x, y, and z directions. 
    else if ( strcmp(argv[n], "-linearcomplexity") == 0 || strcmp(argv[n], "-lc") == 0 )
      {
	linear_complexity = TRUE;
	if (verbose) printf("\nCOMPUTING THE LINEAR COMPLEXITY.");
      }
    // Factor an integer.
    else if ( strcmp(argv[n], "-factor") == 0 )
      {
	k = atoi(argv[++n]);
	if (verbose) printf("\nCOMPUTING THE PRIME FACTORISATION OF %i:", k);
	factor_integer_disp(k);
	return 0;
      }
    // Compute Euler's Phi function.
    else if ( strcmp(argv[n], "-eulerphi") == 0 )
      {
	k = atoi(argv[++n]);
	if (verbose) printf("\nCOMPUTING EULER'S PHI OF %i:", k);
	euler_phi_disp(k);
	return 0;
      }
    // Enumerate the array along the diagonals. (Chinese remainder theorem.)
    else if ( strcmp(argv[n], "-crt") == 0 )
      {
	if (verbose) printf("\nENUMERATING ARRAY ALONG ALL DIAGONALS USING CRT.");
	enum_diagonals = TRUE;
      }
    // Compute the primitive root.
    else if ( strcmp(argv[n], "-primitiveroot") == 0 )
      {
	k = atoi(argv[++n]);
	if (verbose) printf("\nCOMPUTING THE PRIMITIVE ROOT OF %i:", k);
	primitive_root_disp(k);
	return 0;
      }
    // Generate primitive polynomials. 
    else if ( strcmp(argv[n], "-primitivepolynomials") == 0 )
      {
	k = atoi(argv[++n]);
	if (p == -1) {
	  printf("\nERROR: -prime must be specified before -primitivepolynomials.\n");
	  exit(1);
	}
	if (polydegree == -1) {
	  printf("\nERROR: -degree must be specified before -primitivepolynomials.\n");
	  exit(1);
	}
	if (verbose) printf("\nGENERATING %i PRIMITIVE POLYNOMIALS.", k);
	primitive_polynomial_test(p, polydegree, k);
	return 0;
      }
    // Construction name. 
    else if ( strcmp(argv[n], "-sidelnikov") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    SIDELNIKOV1D");
	array_flag = Sidelnikov1D_flag;
	construction_name = "SIDELNIKOV1D";	
      }
    else if ( strcmp(argv[n], "-legendre1d") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    LEGENDRE SEQUENCE");
	array_flag = Legendre1D_flag;
	construction_name = "LEGENDRE1D";
      }
    else if ( strcmp(argv[n], "-legendre2d") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    LEGENDRE2D");
	array_flag = Legendre2D_flag;
	construction_name = "LEGENDRE2D";
      }
    else if ( strcmp(argv[n], "-legendreprojection2d") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    LEGENDREPROJECTION2D");
	array_flag = LegendreProjection2D_flag;
	construction_name = "LEGENDREPROJECTION2D";
      }
    else if ( strcmp(argv[n], "-legendre3d") == 0 ) 
      {
	if (verbose) printf("\nCONSTRUCTION:    LEGENDRE3D");
	array_flag = Legendre3D_flag;
	construction_name = "LEGENDRE3D";
      }
    else if ( strcmp(argv[n], "-a") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    A");
	array_flag = A_flag;
	construction_name = "A";
      }
    else if ( strcmp(argv[n], "-d2") == 0)
      {
	if (verbose) printf("\nCONSTRUCTION:    D2");
	array_flag = D2_flag;
	construction_name = "D2";
      }
    else if ( strcmp(argv[n], "-d3") == 0 ) 
      {
	if (verbose) printf("\nCONSTRUCTION:    D3");
        array_flag = D3_flag;
	construction_name = "D3";
      }
    else if ( strcmp(argv[n], "-c1") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    C1");
        array_flag = C1_flag;
	construction_name = "C1";
      }
    else if ( strcmp(argv[n], "-c2") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    C2");
	array_flag = C2_flag;
	construction_name = "C2";
      }
    else if ( strcmp(argv[n], "-tirkel2d") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    TIRKEL-HALL");
	array_flag = Tirkel_Hall_flag;
	construction_name = "TIRKEL-HALL2D";
      }
    else if ( strcmp(argv[n], "-blake3d") == 0 )
      {
	if (verbose) printf("\nCONSTRUCTION:    BLAKE-TIRKEL");
	array_flag = Blake_Tirkel_flag;
	construction_name = "BLAKE-TIRKEL3D";
      }
    else if ( strcmp(argv[n], "-special") == 0 )
      {
	if (verbose) printf("\nGENERATING SPECIAL SHIFT POLYNOMIAL.");
	special = TRUE;
      }
    // Get prime. 
    else if ( strcmp(argv[n], "-prime") == 0 || strcmp(argv[n], "-p") == 0)
      {
	p = atoi(argv[++n]);
	if (! primep(p))
	  {
	    printf("\nERROR: %i is not prime.\n", p);
	    exit(1);
	  }
	if (verbose) printf("\nPRIME:    %i", p);
      }
    // Get primitive polynomial degree.
    else if ( strcmp(argv[n], "-degree") == 0 || strcmp(argv[n], "-deg") == 0 || 
	      strcmp(argv[n], "-d") == 0 )
      {
	polydegree = atoi(argv[++n]);
	if (verbose) printf("\nPOLYNOMIAL DEGREE:    %i", polydegree);
      }
    // Get primitive polynomial.
    else if ( strcmp(argv[n], "-polynomial") == 0 ||  strcmp(argv[n], "-poly") == 0 )
      {
	polyfromuser = TRUE;
	if (polydegree == -1) 
	  {
	    printf("\nERROR: polynomial degree must be specified before the polynomial coefficients.\n");
	    exit(1);
	  }
	// Read in primitive polynomial.
	poly_deg(poly) = polydegree;
	poly_list(poly) = (int*) malloc((1 + polydegree)*sizeof(int));
	for (k = 0; k <= polydegree; k++)
	  if ( strcmp(argv[++n], "-oo") == 0 )
	    poly_list(poly)[k] = NEGINF;
	  else
	    poly_list(poly)[k] = atoi(argv[n]);
	if (verbose) 
	  {
	     printf("\nPOLYNOMIAL:    ");
	     poly_print(poly);
	  }
      }
    // Get the degree of the numerator (transformation) polynomial. 
    else if ( strcmp(argv[n], "-numeratordegree") == 0 || strcmp(argv[n], "-ndeg") == 0 )
      {
	numpolydegree = atoi(argv[++n]);
	if (verbose) printf("\nNUMERATOR POLYNOMIAL DEGREE:    %i", numpolydegree);
      }
    // Get the numerator (transformation) polynomial.
    else if ( strcmp(argv[n], "-numeratorpolynomial") == 0 || strcmp(argv[n], "-npoly") == 0 )
      {
	numpolyfromuser = TRUE;
	if (numpolydegree == -1) 
	  {
	    printf("\nERROR: numerator polynomial degree must be specified \
before the numerator polynomial coefficients.\n");
	    exit(1);
	  }
	// Read in shift sequence polynomial.
	poly_deg(numpoly) = numpolydegree;
	poly_list(numpoly) = (int*) malloc((1 + numpolydegree)*sizeof(int));
	for (k = 0; k <= numpolydegree; k++)
	  if ( strcmp(argv[++n], "-oo") == 0 )
	    poly_list(numpoly)[k] = NEGINF;
	  else
	    poly_list(numpoly)[k] = atoi(argv[n]);
	if (verbose) 
	  {
	     printf("\nNUMERATOR TRANSFORMATION POLYNOMIAL:    ");
	     poly_print(numpoly);
	  }
      }
    // Get the degree of the denominator (transformation) polynomial. 
    else if ( strcmp(argv[n], "-denominatordegree") == 0 || strcmp(argv[n], "-ddeg") == 0 )
      {
	denpolydegree = atoi(argv[++n]);
	if (verbose) printf("\nDENOMINATOR POLYNOMIAL DEGREE:    %i", denpolydegree);
      }
    // Get the denominator (transformation) polynomial.
    else if ( strcmp(argv[n], "-denominatorpolynomial") == 0 || strcmp(argv[n], "-dpoly") == 0 )
      {
	denpolyfromuser = TRUE;
	if (denpolydegree == -1) 
	  {
	    printf("\nERROR: denominator polynomial degree must be specified \
before the denominator polynomial coefficients.\n");
	    exit(1);
	  }
	// Read in shift sequence polynomial.
	poly_deg(denpoly) = denpolydegree;
	poly_list(denpoly) = (int*) malloc((1 + denpolydegree)*sizeof(int));
	for (k = 0; k <= denpolydegree; k++)
	  if ( strcmp(argv[++n], "-oo") == 0 )
	    poly_list(denpoly)[k] = NEGINF;
	  else
	    poly_list(denpoly)[k] = atoi(argv[n]);
	if (verbose) 
	  {
	     printf("\nDENOMINATOR TRANSFORMATION POLYNOMIAL:    ");
	     poly_print(denpoly);
	  }
      } 
    // Get the X degree of the bivariate shift polynomial. 
    else if ( strcmp(argv[n], "-degx") == 0 || strcmp(argv[n], "-polydegx") == 0 ) 
      {
	polydegX = atoi(argv[++n]);
	if (verbose) printf("\nDEGREE IN X OF BIVARIATE SHIFT POLYNOMIAL:    %i", polydegX);
      }
    // Get the Y degree of the bivariate shift polynomial. 
    else if ( strcmp(argv[n], "-degy") == 0 || strcmp(argv[n], "-polydegy") == 0 ) 
      {
	polydegY = atoi(argv[++n]);
	if (verbose) printf("\nDEGREE IN Y OF BIVARIATE SHIFT POLYNOMIAL:    %i", polydegY);
      }
    // Get the bivariate shift polynomial. 
    else if ( strcmp(argv[n], "-shiftpoly2d") == 0 ) 
      {
	if (polydegX == -1 || polydegY == -1) 
	  {
	     printf("\nERROR: polydegX and polydegY must be specified before the \
bivariate shift polynomial.\n\n");
	     exit(1);
	  }
	
	poly_2D_specified = TRUE;

	// Define polynomial degrees.
	poly_degX(shiftpoly2D) = polydegX;
	poly_degY(shiftpoly2D) = polydegY;

	// Allocate memory for the coefficients. 
	shift_poly_coeffs_2D = (int**) malloc((1 + polydegX)*sizeof(int*));
	for (x = 0; x <= polydegX; x++) 
	  shift_poly_coeffs_2D[x] = (int*) malloc((1 + polydegY)*sizeof(int));

	// Read in bivariate shift sequence polynomial. 
	for (x = 0; x <= polydegX; x++) {
	  for (y = 0; y <= polydegY; y++)
	    shift_poly_coeffs_2D[x][y] = atoi(argv[++n]);
	}
	poly_array(shiftpoly2D) = shift_poly_coeffs_2D;
 
	if (verbose) {
	  printf("\nBIVARIATE TRANSFORMATION POLYNOMIAL IS: \n");
	  poly_print_2D(shiftpoly2D);
	}
      }
    // Get output file name. 
    else if ( strcmp(argv[n], "-o") == 0 || strcmp(argv[n], "-output") == 0 )
      {
	autofilename = FALSE;
	strcpy(filename,argv[++n]);
	if (verbose) printf("\nFILENAME:    %s", filename);
      }
    // Generate binary output file. 
    else if ( strcmp(argv[n], "-b") == 0 || strcmp(argv[n], "-bin") == 0 || 
	      strcmp(argv[n], "-binary") == 0 )
      {
	binary = TRUE;
	strcpy(filename, ""); // Clear file name.
	strcpy(filename, default_filename);
	strcat(filename, default_bin_extension);
	if (verbose) printf("\nCREATING BINARY ARRAY FILE.");
      }
    // Generate family. 
    else if (  strcmp(argv[n], "-f") == 0 ||  strcmp(argv[n], "-family") == 0 )
      {
	genfamily = TRUE;
	if (verbose) printf("\nGENERATING FULL FAMILY.");
      }
    else if ( strcmp(argv[n], "-binarize") == 0 ||  strcmp(argv[n], "-binarise") == 0 )
      {
	binarize = TRUE;
	if (verbose) printf("\nCONVERTING ARRAY TO NEAR-BINARY.");
      }
    else if ( strcmp(argv[n], "-autocorrelate") == 0 ||  strcmp(argv[n], "-acv") == 0 )
      {
	acv = TRUE;
	if (verbose) printf("\nAUTOCORRELATIONS WILL BE COMPUTED.");
      }
    else if ( strcmp(argv[n], "-cross-correlate") == 0 || strcmp(argv[n], "-ccv") == 0 || 
	      strcmp(argv[n], "-xcv") == 0 )
      {
	xcv = TRUE;
	if (verbose) printf("\nCROSS-CORRELATIONS WILL BE COMPUTED.");
      }
    else if (strcmp(argv[n], "-zeromap") == 0) 
      {
	zeromap = atoi(argv[++n]);
	if (verbose) printf("\nMAPPING ZERO TO %i.", zeromap);	
      }
    else if ( strcmp(argv[n], "-powertable") == 0 ) 
      {
	disp_pow_table = TRUE;
	if (verbose) printf("\nDISPLAYING TABLE OF POWERS OF ALPHA.");
      }
    else if ( strcmp(argv[n], "-gz") == 0 || strcmp(argv[n], "-gzip") == 0 || 
	      strcmp(argv[n], "-zip") == 0) {
      zip = TRUE;
      if (verbose) printf("\nGZIPPING OUTPUT FILE.");
    }
    else if ( strcmp(argv[n], "-csv") == 0 )
      {
	csv = TRUE;
	if (verbose) printf("\nWRITING OUTPUT IN CSV FORMAT.");
      }
    else if ( strcmp(argv[n], "-msequence") == 0 )
      {
	array_flag = Msequence_flag;
	if (verbose) printf("\nCONSTRUCTION:    M-SEQUENCE");
	construction_name = "M-SEQUENCE";
      }
    else
      {
	printf("\nERROR: unknown option:  %s\n\n", argv[n]);
	exit(1);
      }
  }
  
  /* User must specify a polynomial degree if the construction is a m-sequence. */
  if (array_flag == Msequence_flag && polydegree == -1) {
    printf("\nERROR: primitive polynomial degree not specified. \n\n");
    exit(1);
  }

  /* Check the user has specified the polynomial degree. */
  if (polydegree == -1 && array_flag != Legendre1D_flag && array_flag != Tirkel_Hall_flag && 
      array_flag != Blake_Tirkel_flag && array_flag != LegendreProjection2D_flag && 
      array_flag != A_flag) {
    polydegree = (array_flag == Legendre3D_flag) ? 3 : 2;
    if (verbose) printf("\nDEFAULT PRIMITIVE POLYNOMIAL DEGREE:    %i", polydegree);
  }

  start = clock(); 

  /* If lookup table is not loaded, then randomly generate a primitive polynomial. */
  if ( MAX_LU_PRIME == 0 ) randomise = TRUE;

  /* If the user has not specified a primitive polynomial, then create one at random. */
  if (array_flag != Legendre1D_flag && array_flag != Tirkel_Hall_flag && 
      array_flag != Blake_Tirkel_flag && array_flag != LegendreProjection2D_flag && 
      array_flag != A_flag) {
    if (polyfromuser == FALSE && randomise == TRUE) {
      if (verbose) printf("\nGENERATING RANDOM PRIMITIVE POLYNOMIAL...");
      poly = gen_rand_primitive_poly(p, polydegree);
      if (verbose) {
	printf("\nPRIMITIVE POLYNOMIAL IS:    ");
	poly_print(poly);
      }
      random_prim_poly_generated = TRUE;
    }
    else if (polyfromuser == FALSE && randomise == FALSE) {
      if (verbose) printf("\nRANDOMLY SELECTING PRIMITIVE POLYNOMIAL FROM TABLE LOOKUP...");  
      // UNDER CONSTRUCTION...
      if (verbose) printf("\nPRIMITIVE POLYNOMIAL IS:    ");
      poly_print(poly);
    }
  }

  /* Generate random univariate (numerator) shift polynomial. */
  if (numpolyfromuser == FALSE && genfamily == FALSE && array_flag != Sidelnikov1D_flag && 
      array_flag != Legendre2D_flag && array_flag != Legendre3D_flag && array_flag != -1 && 
      array_flag != Legendre1D_flag && array_flag != Blake_Tirkel_flag && 
      array_flag != LegendreProjection2D_flag && array_flag != Msequence_flag) {
    if (numpolydegree == -1) {
      if (array_flag != D3_flag && array_flag != C2_flag) {
	numpolydegree = 2;
	numpoly = gen_rand_poly(p, numpolydegree); /* Quadratic by default. */
      } else {
	numpolydegree = 1;
	numpoly = gen_rand_poly(p, numpolydegree); /* Linear by default. */
      }
    } else {
      numpoly = gen_rand_poly(p, numpolydegree);
    }
    if (verbose) {
      if (array_flag != D3_flag && array_flag != C2_flag)
	printf("\nRANDOM SHIFT SEQUENCE POLYNOMIAL IS:    ");
      else
	printf("\nRANDOM NUMERATOR SHIFT SEQUENCE POLYNOMIAL IS:    ");
      poly_print(numpoly);
    }
  }

  /* Generate random univariate (denominator) shift polynomial. */
  if (denpolyfromuser == FALSE && genfamily == FALSE && 
        (array_flag == D3_flag || array_flag == C2_flag)) {
    if (denpolydegree == -1) {
      denpolydegree = 1;
      denpoly = gen_rand_poly(p, denpolydegree); /* Linear by default. */
    } else {
      denpoly = gen_rand_poly(p, denpolydegree);
    }
    if (verbose) {
	printf("\nRANDOM DENOMINATOR SHIFT SEQUENCE POLYNOMIAL IS:    ");
      poly_print(denpoly);
    }
  }

  /* Generate random bivariate shift polynomial. */
  if ( array_flag == Blake_Tirkel_flag && ! poly_2D_specified && ! genfamily ) {
    if (special) {
      shiftpoly2D = gen_special_rand_poly_2D(p);
    } else if (! (polydegX == -1 || polydegY == -1) ) {
      shiftpoly2D = gen_rand_poly_2D(p, polydegX, polydegY);
    } else {
      printf("\nError: bivariate polynomial degree not specified.\n\n");
      exit(1);
    }

    if (verbose) printf("\nRANDOM BIVARIATE SHIFT SEQUENCE POLYNOMIAL IS:    ");
    if (verbose) poly_print_2D(shiftpoly2D);    
  }

  stop = clock();
  elapsed = ((double) (stop - start)) / CLOCKS_PER_SEC;
  if (verbose && random_prim_poly_generated) printf("\nElapsed time: %g seconds\n", elapsed);

  if (verbose) printf("\n"); /* Don't remove me! */

  if (disp_pow_table)
    print_power_table(p, poly);

  /* Construct array(s). */
  switch (array_flag)
    {
        case SIDELNIKOV_1D_FLAG:
	  construction_Sidelnikov1D(filename, p, poly, binarize, acv, pipe, binary, verbose, csv);
          break;
        case LEGENDRE_1D_FLAG:
	  construction_Legendre1D(filename, p, acv, pipe, binary, verbose, csv);
	  break;
        case LEGENDRE_2D_FLAG:
	  if (genfamily) {
	    construction_legendre_projection_2D(filename, p, acv, pipe, binary, zip, verbose, csv);
	  } else {
	    construction_Legendre2D(filename, p, numpoly, denpoly, poly, binarize, 
				    acv, pipe, binary, verbose, csv);
	  }
	  break;
        case LEGENDREPROJECTION_2D_FLAG:
	  construction_legendre_projection_2D(filename, p, acv, pipe, binary, zip, verbose, csv);
	  break;
        case LEGENDRE_3D_FLAG:
	  construction_Legendre3D(filename, p, numpoly, denpoly, poly, binarize, 
				  acv, pipe, binary, zip, verbose, csv);
	  break;
        case A_FLAG:
	  if (genfamily) {
	    if (numpolydegree == -1) numpolydegree = 2;
	    construction_a_family(filename, p, numpolydegree, zeromap, enum_diagonals, 
				  acv, pipe, binary, verbose, csv);
	  } else {
	    construction_a(filename, p, zeromap, enum_diagonals, numpoly, 
			   acv, pipe, binary, verbose, csv);
	  }
	  break;
        case  TIRKEL_HALL_FLAG:
	  if (genfamily) {
	    if (numpolydegree == -1) numpolydegree = 2;
	    construction_Tirkel_Hall2D_family(filename, p, numpolydegree, acv,
					      zeromap, pipe, binary, verbose, csv);
	  } else {
	    construction_Tirkel_Hall2D(filename, p, zeromap, numpoly, acv, 
				       pipe, binary, verbose, csv);
	  }
	  break;
        case BLAKE_TIRKEL_FLAG:
	  if (genfamily) {
	    construction_Blake_Tirkel3D_special_family(filename, p, zeromap, pipe, acv,
						       binary, zip, verbose, csv);
	  } else {
	    construction_Blake_Tirkel3D(filename, p, zeromap, shiftpoly2D, 
					acv, pipe, binary, zip, verbose, csv);
	  }
	  break;
        case C1_FLAG:
	  if (genfamily) {
	    printf("\nERROR: The array generator cannot construct the given family.\n\n");
	    exit(1);
          } else {
	    construction_c1(filename, p, poly, numpoly, binarize, acv, pipe, 
			  binary, zip, verbose, csv, linear_complexity);
	  }
	  break;
        case D2_FLAG:
          if (genfamily) {
	    if (numpolydegree == -1) numpolydegree = 2;
	    construction_d2_family(filename, p, poly, numpolydegree, binarize, 
				   acv, pipe, binary, zip, verbose, csv);
          } else {
	    construction_d2(filename, p, poly, numpoly, binarize, acv, pipe, 
			    binary, zip, verbose, csv);
	  }
	  break;
        case D3_FLAG:
          if (genfamily) {
	    printf("\nERROR: The array generator cannot construct the given family.\n\n");
	    exit(1);
	  } else {
	    construction_d3(filename, p, poly, numpoly, denpoly, binarize, 
			    acv, pipe, binary, zip, verbose, csv);
	  }
	  break;
        case C2_FLAG:
          if (genfamily) 
	    printf("\nERROR: The array generator cannot construct the given family.\n\n");
          else
	    construction_c2(filename, p, poly, numpoly, denpoly, binarize, 
			    acv, pipe, binary, zip, verbose, csv);
	  break;
        case MSEQUENCE_FLAG:
          construction_m_sequence(filename, p, poly,
				 acv, pipe, binary, zip, verbose, csv);
	  break;
    }

  /* Free polynomials. */
  // TODO

  stop = clock();
  elapsed = ((double) (stop - start)) / CLOCKS_PER_SEC;
  if (verbose) printf("\nElapsed time: %g seconds\n", elapsed); 
#if DEBUG
  if (verbose) printf("\nmalloc_count = %d", malloc_count);
  if (verbose) printf("\ngf_alloc_count = %d\n", gf_alloc_count);
#endif
  return 0;
}

/*
 *
 *    PROGRAM INFORMATION
 *    -------------------
 */
void display_help() {
  printf("\nNAME\n\
    arraygen - This is our experimental array generator platform.\n\n");

  printf("SYNTAX\n\
    arraygen -help -pipe -verbose -prime -binary -output -family -binarise\n\
              -randomise -deg -poly -ndeg -npoly -ddeg -dpoly -sidelnikov\n\
              -legendre2d -legendre3d -A -C1 -C2 -D2 -D3 -tirkel2d -blake3d\n\
              -msequence -autocorrelate -csv -crt \n\n");

  printf("DESCRIPTION\n\
    This autocorrelation/cross-correlation family array generator constructs 1, 2, and 3 dimensional\n\
    arrays using both known, published and patented array constructions.\n\
        Output formats:\n\
            At present the output formats supported are text files containing space \n\
            separated arrays and csv files. Some experimental support for binary files exists.\n\
        Arithmetic: \n\
            arraygen implements arithmetic in GF(p^d). Initially, the program\n\
    will create a table of powers of alpha, the primitive element of the Galois field.\n\
    This table is used for subsequent calculations in GF(p^d). The program can\n\
    pseudo-randomly generate a primitive polynomial in GF(p^d). Otherwise, if the user\n\
    specifies a primitive polynomial, the program will check that the input polynomial \n\
    is indeed primitive before proceeding with any array constructions.\n\n\
    The univariate shift polynomials are entered in \"powers of \\alpha\" notation. Thus,\n\
    the polynomial specified by 4 3 7 is stored internally as \n\
    \\alpha^4 + \\alpha^3*x + \\alpha^7*x^2.\n\n\
    The bivariate shift polynomials are not stored as powers of \\alpha. These polynomials \n\
    are entered into the program as a flattened array of coefficients, such that a \n\
    dot product with\n\n\
          1       y     y^2 ...     y^degY\n\
          x     x*y   x*y^2 ...   x*y^degY\n\
          x^2 x^2*y x^2*y^2 ... x^2*y^degY\n\
          ...                        .\n\
          .                          .\n\
          .                          .\n\
          x^degX x^degX*y x^degX*y^2 ... x^degX*y^degY\n\n\
    will generate a bivariate polynomial of total degree degX + degY.\n\n\
    The following options are available (options are not case sensitive):\n\n\
        -help\n\
            Displays help messages and returns. (This message!)\n\n\
        -version\n\
            Returns the version of the program.\n\n\
        -p or prime [integer]\n\
            Specifies the prime to use in the Galois field.\n\n\
        -v or -verbose\n\
            Prints messages to stdout about the programs execution.\n\n\
        -pipe\n\
            Tell the program the array generated will be piped to another program.\n\n\
        -o or -output [file name]\n\
            Specify the output file name. By default the filename is arraygen.out for\n\
            ASCII files and arraygen.bin for binary files.\n\n\
        -csv\n\
            Writes a csv file.\n\n\
        -binary\n\
            Specify that the output file should be a binary file.\n\n\
        -family\n\
            Specifies that the full family of p^d arrays is to be generated.\n\n\
        -binarise\n\
            Converts the array to -1, +1, with -oo -> 0.\n\n\
        -randomise\n\
            Specifies that a primitive polynomial should be generated randomly.\n\n\
        -crt\n\
            Instructs the program to enumerate the array along the diagonals via the \n\
            Chinese remainder theorem.\n\n\
        -autocorrelate or -acv\n\
            Specifies that the periodic autocorrelation of the array(s) should be generated. \n\
            The array of autocorrelations are written to arraygen.acv. Note that computing \n\
            correlations of large arrays significantly slows the execution of arraygen.\n\n\
        -powertable\n\
            Prints the power table to stdout when -verbose is specified.\n\n\
        -deg or -degree [integer]\n\
            Degree of the primitive polynomial.\n\n\
        -poly or -polynomial [list of coefficients]\n\
            Specifies the coefficient list of the primitive polynomial with the\n\
            constant term first.\n\n\
        -ndeg or -numeratordegree [integer]\n\
            Degree of the numerator shift polynomial.\n\n\
        -npoly or -numeratorpolynomial [list of coefficients]\n\
            Specifies the coefficient list of the numerator shift polynomial with\n\
            the constant term first.\n\n\
        -ddeg or -denominatordegree [integer]\n\
            Degree of the denominator shift polynomial.\n\n\
        -dpoly or -denominatorpolynomial [list of coefficients]\n\
            Specifies the coefficient list of the denominator shift polynomial with\n\
            the constant term first.\n\n\
        -factor [integer]\n\
            Returns the prime factorisation of the integer input.\n\n\
        -eulerphi [integer]\n\
            Returns Euler's phi of an integer input. \n\n\
        -primitiveroot [integer]\n\
            Returns the smallest primitive root of the integer input.\n\n\
        -sidelnikov\n\
            Instructs the program to generate a Sidelnikov sequence.\n\n\
        -legendre1d\n\
            Instructs the program to generate a Legendre sequence.\n\n\
        -legendre2d\n\
            Instructs the program to generate a 2D Legendre array.\n\n\
        -legendre3d\n\
            Instructs the program to generate a 3D Legendre array.\n\n\
        -tirkel2d\n\
            Instructs the program to generate a 2D Tirkel array of size p x p over 0, -1, 1.\n\n\
        -blake3d\n\
            Instructs the program to generate a 3D Blake-Tirkel array of size p x p x p \n\
            over 0, -1, 1.\n\n\
        -a\n\
            Instructs the program to generate a A-type array.\n\n\
        -c1\n\
            Instructs the program to generate a C1-type array.\n\n\
        -c2\n\
            Instructs the program to generate a C2-type array.\n\n\
        -d1\n\
            Instructs the program to generate a D1-type array.\n\n\
        -d2\n\
            Instructs the program to generate a D2-type array.\n\n\
        -d3\n\
            Instructs the program to generate a D3-type array.\n\n");

  printf("EXAMPLES\n\
    Compute Euler's totient function of 9999191:\n\n\
        $ ./arraygen -eulerphi 9999191\n\n\
    Factor 2^29 - 1: \n\n\
        $ ./arraygen -factor 536870911\n\n\
    Compute the smallest primitive root of 17449:\n\n\
        $ ./arraygen -primitiveroot 17449\n\n\
    To generate a length 11^2 - 1 Sidelnikov sequence (p = 11) with a randomly\n\
    generated primitive polynomial:\n\n\
        $ ./arraygen -v -p 11 -sidelnikov\n\n\
    To generate a length 29^2 - 1 Sidelnikov sequence (p = 29) with the primitive\n\
    polynomial 21 + 6 x + x^2:\n\n\
        $ ./arraygen -v -p 29 -deg 2 -poly 21 6 1 -sidelnikov\n\n\
    To generate a 2D Legendre array of size 17 by 17 with a randomly generated\n\
    primitive polynomial:\n\n\
        $ ./arraygen -v -p 17 -legendre2d\n\n\
    To generate a 2D legendre array with the primitive polynomial 3 + 3 x + x^2: \n\n\
        $ ./arraygen -v -p 5 -legendre2d -deg 2 -poly 3 3 1\n\n\
    To generate a family of 17, 17 x 17 Legendre arrays via projections and \n\
    compute the autocorrelation of each array:\n\n\
        $ ./arraygen -v -p 17 -legendreprojection2d -autocorrelate\n\n\
    To generate a 3D Legendre array of size 5 by 5 by 5 with a randomly generated\n\
    primitive polynomial:\n\n\
        $ ./arraygen -v -p 5 -legendre3d\n\n\
    To generate a 3D Legendre array of size 19 x 19 x 19 with the\n\
    primitive polynomial 17 + 4 x + 3 x^2 + x^3:\n\n\
        $ ./arraygen -v -p 19 -legendre3d -deg 3 -poly 17 4 3 1\n\n\
    To generate a 17 x 16 array from construction A with shift polynomial\n\
    11 + 6 x + 13 x^2:\n\n\
        $ ./arraygen -v -p 17 -ndeg 2 -npoly 11 6 13 -a\n\n\
    To generate a 11 x 10 array from construction A with random shift polynomial\n\
    and enumerate the array along the diagonal via the Chinese remainder theorem:\n\n\
        $ ./arraygen -v -p 11 -a -crt -acv\n\n\
    To generate a 3D array of size 11 x 11 x 120 using construction C1 with a randomly\n\
    generated primitive polynomial and shift polynomial:\n\n\
        $ ./arraygen -v -p 11 -c1\n\n\
    To generate a 3D array of size 37 x 37 x 1368 using construction C1 with the primitive\n\
    polynomial 17 + 32 x + x^2 and a random shift polynomial:\n\n\
        $ ./arraygen -v -p 37 -c1 -deg 2 -poly 17 32 1\n\n\
    Similarly, one can specify the shift polynomial and not the primitive polynomial:\n\n\
        $ ./arraygen -v -p 37 -c1 -ndeg 2 -npoly 1 2 3\n\n\
    And, of course, both the primitive polynomial and the shift polynomial can be specified:\n\n\
        $ ./arraygen -v -p 37 -c1 -deg 2 -poly 17 32 1 -ndeg 2 -npoly 1 2 3\n\n\
    To generate a 3D array of size 7 x 7 x 48 from construction C2 with randomly generated\n\
    primitive and shift polynomials:\n\n\
        $ ./arraygen -v -p 7 -c2\n\n\
    Similarly, to generate an array of size 7 x 7 x 48 from construction D2 with a randomly\n\
    generated primitive polynomial and shift polynomial:\n\n\
        $ ./arraygen -v -p 7 -d2\n\n\
    The program can also compute the autocorrelation of an array. For example:\n\n\
        $ ./arraygen -v -p 23 -legendre2d -autocorrelate\n\n\
    computes the autocorrelation of a 23 x 23 Legendre array.\n\n\
    Two-dimensional Tirkel-Hall arrays can be generated. For example, the following generates\n\
    a family of 48 arrays of size 7 x 7 using a cubic shift polynomial:\n\n\
        $ ./arraygen -v -p 7 -ndeg 3 -tirkel2d -family -autocorrelate\n\n\
    Three-dimensional Blake-Tirkel arrays can be generated. For example: \n\n\
        $ ./arraygen -v -p 7 -blake3D -polydegx 2 -polydegy 2 -shiftpoly2D 0 0 1 0 1 0 1 0 0\n\n\
    computes a 3D Blake-Tirkel array of size 7 x 7 x 7 with shift polynomial y^2 + x y + x^2.\n\n\
    Here we generate a Blake-Tirkel array of size 11 x 11 x 11 with a randomly selected special \n\
    shift polynomial:\n\n\
        $ ./arraygen -v -p 11 -blake3D -special\n\n\
    To generate a 3D D3-type array of size 7 x 7 x 48, with numerator polynomial randomly \n\
    constructed of degree 5 and denominator polynomial randomly constructed of degree 2, \n\
    binarize the array and compute the autocorrelation:\n\n\
        $ ./arraygen -v -p 7 -d3 -binarize -ndeg 5 -ddeg 2 -acv\n\n\
    To compute the table of powers of alpha for a random primitive polynomial of degree 5 \n\
    in GF(3^5):\n\n\
        $ ./arraygen -v -p 3 -deg 5 -powertable\n\n\
    note that no sequence or array construction has been specified.\n\n\
    The following command generates 32 primitive polynomials of degree 12 mod 3:\n\n\
        $ ./arraygen -v -p 3 -d 12 -primitivepolynomials 32\n\n");

  printf("DEBUGGING\n\
    Some timing and debugging information is printed to screen when -verbose (or -v) is given \n\
    to arraygen. The variables malloc_count and gf_alloc_count are displayed at the end of the \n\
    programs execution. malloc_count displays the total number of calls to malloc which have not \n\
    been freed (ie. memory leaks). gf_alloc_count counts the total number of Galois field symbolic \n\
    expressions generated during the entire execution of the program. The internal arithmetic of \n\
    the program can be inspected by specifying -powertable, which displays the full table of powers \n\
    of \\alpha in polynomial form.\n\n\
    If something completely unexpected has occured during the programs execution, the message\n\
    INTERNAL ERROR /error code/ may appear. The program inputs and the error code should be reported \n\
    to Sam Blake (samuel.thomas.blake [at] gmail.com).\n\n");

}


/*
 *
 *    PROGRAM INTERFACE TO THE CONSTRUCTION OF M SEQUENCES
 *    ----------------------------------------------------
 */

void construction_m_sequence(char *filename, int p, polynomial primpoly,
			    bool acv, bool pipe, bool binaryfile, 
			    bool zip, bool verbose, bool csv) {
  int *mseq, deg, len;
  
  deg = poly_deg(primpoly); // Degree of the primitive polynomial.
  len = power(p, deg) - 1;  // Length of the m-sequence. 

  mseq = m_sequence(primpoly, p); 

  io_1d(filename, mseq, len, csv, acv, pipe, binaryfile, verbose);
}

/*
 *
 *    PROGRAM INTERFACE TO CONSTRUCTION A 2D ARRAYS FAMILY GENERATOR
 *    --------------------------------------------------------------
 */

void construction_a_family(char *filename, int p, int shift_poly_degree, int zeromap,  
			   bool enum_diag, bool acv, bool pipe, bool binaryfile, 
			   bool verbose, bool csv) {
  int k, j, *cl, npad,
    family_size = power(p, shift_poly_degree - 1);
  char temp[128], textfile[128], binfile[128], index[16], sprime[20], pattern[20];
  polynomial shift_poly;

  strcpy(textfile, default_filename);
  strcat(textfile, default_txt_extension);

  strcpy(binfile, default_filename);
  strcat(binfile, default_bin_extension);

  if (verbose) printf("FAMILY SIZE = %i\n", family_size);

  memset(index, '\0', sizeof(index));
  memset(temp, '\0', sizeof(temp));

  /* Initialise array of coefficients of the shift polynomial to 1's. */
  cl = (int*) malloc((shift_poly_degree + 1)*sizeof(int));
  MEMCHECK(cl);
  for (k = 0; k <= shift_poly_degree; k++) 
    cl[k] = 1;

  poly_list(shift_poly) = cl;

  /* Set the degree of the shift polynomial. */
  poly_deg(shift_poly) = shift_poly_degree;

  sprintf(sprime, "_A_p%d_deg%d", p, shift_poly_degree);

  for (k = 1; k <= family_size; k++) {
    /* Generate file name. */  
    npad = intlength(family_size) - intlength(k);
    strcpy(pattern, "_");
    for (j = 0; j < npad; j++) 
      strcat(pattern, "0");
    strcat(pattern, "%d");
    sprintf(index, pattern, k);

    if ( strcmp(filename, textfile) == 0 ) {
      /* Write file name as, for example: arraygen_p7_deg2_392.txt */
      strcpy(temp, default_filename);
      strcat(temp, sprime);
      strcat(temp, index);
      strcat(temp, default_txt_extension);
    } else if ( strcmp(filename, binfile) == 0 ) {
      /* Write file name as, for example: arraygen_p7_deg2_392.bin */
      strcpy(temp, default_filename);
      strcat(temp, sprime);
      strcat(temp, index);
      strcat(temp, default_bin_extension);
    } else {
      strcpy(temp, filename);
      strcat(temp, index);
    }

    if (verbose) printf("\nFILENAME = %s\n", temp);
    if (verbose) {
      printf("SHIFT POLYNOMIAL = ");
      poly_print(shift_poly);
      printf("\n");
    }

    /* Create a A-type array. */
    construction_a(temp, p, zeromap, enum_diag, shift_poly, acv, pipe, binaryfile, verbose, csv);
    
    /* Increment the coefficient list of the polynomial. NOTE: We do not touch the 
       last two coefficients, these are fixed. That is the coefficients of x^1 and x^0
       are held constant. */
    reverse(cl, 0, shift_poly_degree + 1);
    increment(cl, shift_poly_degree, p); 
    reverse(cl, 0, shift_poly_degree + 1);
    poly_list(shift_poly) = cl;
  }
}

/*
 *
 *    PROGRAM INTERFACE TO CONSTRUCTION A 2D ARRAYS
 *    ---------------------------------------------
 */

void construction_a(char *filename, int p, int zeromap, bool enum_diag, polynomial shiftpoly, 
		    bool acv, bool pipe, bool binary, bool verbose, bool csv) {
  int **array2d, *seq;
  char temp[128];

  memset(temp, '\0', 128);

  /* Construct the array. */
  array2d = A_array_2D(p, shiftpoly, zeromap);

  /* If requested, enumerate the array along all diagonals. */
  if (enum_diag) seq = CRTenumerate2D(array2d, p, p - 1);

  io_2d(filename, array2d, p, p - 1, csv, acv, pipe, binary, verbose);

  if (enum_diag) {
    strncpy(temp, filename, strlen(filename) - 4);
    /* This is a hack! */
    if (binary)
      strcat(temp, "_crt.bin");
    else if (csv)
      strcat(temp, "_crt.csv");
    else 
      strcat(temp, "_crt.txt");
    if (verbose) printf("\n\n");
    io_1d(temp, seq, p*(p - 1), csv, acv, pipe, binary, verbose);
    if (verbose) printf("\n");
  }
}

/*
 *
 *    PROGRAM INTERFACE TO THE LEGENDRE 2D FAMILY GENERATOR (VIA PROJECTIONS)
 *    -----------------------------------------------------------------------
 */

void construction_legendre_projection_2D(char *filename, int p, bool acv, bool pipe, 
				      bool binaryfile, bool zip, bool verbose, bool csvp) {
  int *legendre, **lp, i, j, k, s, l, c, npad;
  char temp[128], textfile[128], pattern[32], sprime[16], index[16];

  memset(index, '\0', sizeof(index));
  memset(temp, '\0', sizeof(temp));

  sprintf(sprime, "_p%d", p);

  strcpy(textfile, default_filename);
  strcat(textfile, default_txt_extension);

  if (! primep(p))
    {
      printf("\nERROR: %i is not prime.\n", p);
      exit(1);
    }

  /* Construct the Legendre sequence. */
  legendre = legendre1D(p);

  for (s = 0; s < p; s++) {
    /* Allocate memory for the p x p array. */
    lp = (int**) malloc(p*sizeof(int*));
    MEMCHECK(lp);
    
    for (k = 0; k < p; k++) {
      lp[k] = (int*) malloc(p*sizeof(int));
      MEMCHECK(lp[k]);
      for (l = 0; l < p; l++) lp[k][l] = 0;
    }

    /* Generate array. */ 
    for (l = 0; l < p; l++) {
      c = legendre[(l + s) % p];
      for (i = 0; i < p; i++) {
	lp[i][(l*i)%p] = c;
	lp[0][i] = -1;
      }
    }

    for (l = 0; l < p; l++) {
      for (i = 0; i < p; i++) {
    	if (lp[l][i] == 0) lp[l][i] = 1;
      }
    }
    
    lp[0][0] = 0;

    /* Create file name. */
    npad = intlength(p) - intlength(s);
    strcpy(pattern, "_");
    for (j = 0; j < npad; j++) 
      strcat(pattern, "0");
    strcat(pattern, "%d");
    sprintf(index, pattern, s);
    if ( strcmp(filename, textfile) == 0 ) {
      strcpy(temp, default_filename);
      strcat(temp, sprime);
      strcat(temp, index);      
    } else {
      strcpy(temp, filename);
      strcat(temp, sprime);
      strcat(temp, index);
    }
    if (csvp)
      strcat(temp, ".csv");
    else if (binaryfile)
      strcat(temp, ".bin");
    else
      strcat(temp, ".txt");

    if (verbose) printf("\n%s\n", temp);

    /* Export array to file. */
    io_2d(temp, lp, p, p, csvp, acv, pipe, binaryfile, verbose);
  }

  free(legendre);
}

int** legendre_projection_2D(int p, int s) {
  int *legendre, **lp, n, l, i, c;

  if (! primep(p))
    {
      printf("\nERROR: %i is not prime.\n", p);
      exit(1);
    }

  legendre = legendre1D(p);

  /* Allocate memory for p x p array. */
  lp = (int**) malloc(p*sizeof(int*));
  MEMCHECK(lp);

  for (n = 0; n < p; n++) {
    lp[n] = (int*) malloc(p*sizeof(int));
    MEMCHECK(lp[n]);
  }  

  for (l = 0; l < p; l++) {
    c = legendre[(l + s) % p];
    for (i = 0; i < p; i++) {
      lp[i][(l*i)%p] = c;
      lp[0][i] = -1;
    }
  }

  for (l = 0; l < p; l++) {
    for (i = 0; i < p; i++) {
      if (lp[l][i] == 0) lp[l][i] = 1;
    }
  }
    
  lp[0][0] = 0;

  free(legendre);

  return lp;
}

/*
 *
 *    PROGRAM INTERFACE TO CONSTRUCTION C2 CODE
 *    -----------------------------------------
 */
void construction_c2(char *filename, int p, polynomial primpoly, polynomial numpoly, 
		     polynomial denpoly, bool binarize, bool acv, bool pipe, bool binaryfile, 
		     bool zip, bool verbose, bool csvp) {
  int ***array, len = p*p - 1;

  if (poly_deg(numpoly) == -1 || poly_deg(denpoly) == -1) 
    {
      printf("\nERROR: shift polynomials not specified or randomly generated.\n");
      exit(1);
    }

  /* Create a C2 array. */
  array = C2_array_3D(p, primpoly, numpoly, denpoly);

  /* Binarize array. */
  if (binarize == TRUE)
    binarize3D(array, p, p, len);

  /* Export results. */
  io_3d(filename, array, p, p, len, csvp, binarize, acv, pipe, binaryfile, zip, verbose, FALSE);
}

/*
 *
 *    PROGRAM INTERFACE TO CONSTRUCTION D3 CODE
 *    -----------------------------------------
 */

void construction_d3(char *filename, int p, polynomial primpoly, polynomial numpoly, 
		     polynomial denpoly, bool binarize, bool acv, bool pipe, bool binaryfile, 
		     bool zip, bool verbose, bool csvp) {
  int ***array, len = p*p - 1;

  if (poly_deg(numpoly) == -1 || poly_deg(denpoly) == -1) 
    {
      printf("\nERROR: shift polynomials not specified or randomly generated.\n");
      exit(1);
    }

  /* Create a D3 array. */
  array = D3_array_3D(p, primpoly, numpoly, denpoly);

  /* Binarize array. */
  if (binarize == TRUE)
    binarize3D(array, p, p, len);

  /* Export results. */
  io_3d(filename, array, p, p, len, csvp, binarize, acv, pipe, binaryfile, zip, verbose, FALSE);
}

/*
 *
 *    PROGRAM INTERFACE TO CONSTRUCTION D2 FAMILY GENERATOR
 *    -----------------------------------------------------
 */
/* Construct a family of p^d arrays from construction D2. */
void construction_d2_family(char *filename, int p, polynomial primpoly, 
			    int shift_poly_degree, bool binarize, bool acv, bool pipe, 
			    bool binaryfile, bool zip, bool verbose, bool csvp) {
  int k, j, *cl, npad,
    d = poly_deg(primpoly), 
    pd = power(p,d),
    family_size = power(pd, shift_poly_degree - 1);
  char temp[128], textfile[128], binfile[128], index[16], sprime[20], pattern[20];
  polynomial shift_poly;

  strcpy(textfile, default_filename);
  strcat(textfile, default_txt_extension);

  strcpy(binfile, default_filename);
  strcat(binfile, default_bin_extension);

  if (verbose) printf("FAMILY SIZE = %i\n", family_size);

  memset(index, '\0', sizeof(index));
  memset(temp, '\0', sizeof(temp));

  /* Initialise array of coefficients of the shift polynomial to 1's. */
  cl = (int*) malloc((shift_poly_degree + 1)*sizeof(int));
  MEMCHECK(cl);
  for (k = 0; k <= shift_poly_degree; k++) 
    cl[k] = 0;

  poly_list(shift_poly) = cl;

  /* Set the degree of the shift polynomial. */
  poly_deg(shift_poly) = shift_poly_degree;

  sprintf(sprime, "_D2_p%d_deg%d", p, shift_poly_degree);

  for (k = 1; k <= family_size; k++) {
    /* Generate file name. */  
    npad = intlength(family_size) - intlength(k);
    strcpy(pattern, "_");
    for (j = 0; j < npad; j++) 
      strcat(pattern, "0");
    strcat(pattern, "%d");
    sprintf(index, pattern, k);

    if ( strcmp(filename, textfile) == 0 ) {
      /* Write file name as, for example: arraygen_p7_deg2_392.txt */
      strcpy(temp, default_filename);
      strcat(temp, sprime);
      strcat(temp, index);
      strcat(temp, default_txt_extension);
    } else if ( strcmp(filename, binfile) == 0 ) {
      /* Write file name as, for example: arraygen_p7_deg2_392.bin */
      strcpy(temp, default_filename);
      strcat(temp, sprime);
      strcat(temp, index);
      strcat(temp, default_bin_extension);
    } else {
      strcpy(temp, filename);
      strcat(temp, index);
    }

    if (verbose) printf("\nFILENAME = %s\n", temp);
    if (verbose) {
      printf("SHIFT POLYNOMIAL = ");
      poly_print(shift_poly);
      printf("\n");
    }

    /* Create a D2 array. */
    construction_d2(temp, p, primpoly, shift_poly, binarize, acv, pipe, binaryfile, zip, verbose, csvp);

    /* Increment the coefficient list of the polynomial. NOTE: We do not touch the 
       last two coefficients, these are fixed. That is the coefficients of x^1 and x^0
       are held constant. */
    reverse(cl, 0, shift_poly_degree + 1);
    increment(cl, shift_poly_degree, pd); 
    reverse(cl, 0, shift_poly_degree + 1);
    poly_list(shift_poly) = cl;
  }
}

/*
 *
 *    PROGRAM INTERFACE TO CONSTRUCTION C1 CODE
 *    -----------------------------------------
 */
void construction_c1(char *filename, int p, polynomial primpoly, polynomial numpoly, 
		     bool binarize, bool acv, bool pipe, bool binaryfile, bool zip, bool verbose, 
		     bool csvp, bool lincomp) {
  int ***array, len = p*p - 1;

  if (poly_deg(numpoly) == -1) 
    {
      printf("\nERROR: shift polynomial not specified or randomly generated.\n");
      exit(1);
    }

  /* Create a C1 array. */
     array = C1_array_3D(p, primpoly, numpoly);

  /* Binarize array. */
  if (binarize == TRUE)
    binarize3D(array, p, p, len);

  /* Export results. */
  io_3d(filename, array, p, p, len, csvp, binarize, acv, pipe, binaryfile, zip, verbose, lincomp);
}

/*
 *
 *    PROGRAM INTERFACE TO CONSTRUCTION D2 CODE
 *    -----------------------------------------
 */
void construction_d2(char *filename, int p, polynomial primpoly, polynomial numpoly, 
		     bool binarize, bool acv, bool pipe, bool binaryfile, bool zip, bool verbose, bool csvp) {
  int ***array, len = p*p - 1;

  if (poly_deg(numpoly) == -1) 
    {
      printf("\nERROR: shift polynomial not specified or randomly generated.\n");
      exit(1);
    }

  /* Create a D2 array. */
     array = D2_array_3D(p, primpoly, numpoly);

  /* Binarize array. */
  if (binarize == TRUE)
    binarize3D(array, p, p, len);

  /* Export results. */
  io_3d(filename, array, p, p, len, csvp, binarize, acv, pipe, binaryfile, zip, verbose, FALSE);
}

/*
 *
 *    PROGRAM INTERFACE TO LEGENGRE 3D CODE
 *    -------------------------------------
 */
void construction_Legendre3D(char *filename, int p, polynomial numpoly, 
			     polynomial denpoly, polynomial primpoly, 
			     bool binarize, bool acv, bool pipe, bool binaryfile, 
			     bool zip, bool verbose, bool csvp) {
  int ***array;

  if (poly_deg(numpoly) == 0) {
    /* Create traditional 3D Legendre array. */
    array = legendre3D(p, primpoly);
  } else {
    /* Create generalised 3D Legendre array. */
    array = legendre3Dgeneralised(p, primpoly, numpoly, denpoly, TRUE);
  }

  /* Binarize array. */
  if (binarize == TRUE)
    binarize3D(array, p, p, p);

  /* Export results. */
  io_3d(filename, array, p, p, p, csvp, binarize, acv, pipe, binaryfile, zip, verbose, FALSE);
}

/*
 *
 *    PROGRAM INTERFACE TO LEGENGRE 2D CODE
 *    -------------------------------------
 */
void construction_Legendre2D(char *filename, int p, polynomial numpoly, polynomial denpoly, 
			     polynomial primpoly, bool binarize, bool acv, bool pipe, 
			     bool binaryfile, bool verbose, bool csvp) {
  int **array; 

  if (poly_deg(numpoly) == 0) {
    /* Create traditional 2D Legendre array. */
    array = legendre2D(p, primpoly);
  } else {
    /* Create generalised 2D Legendre array. */
    array = legendre2Dgeneralised(p, primpoly, numpoly, denpoly, TRUE);
  }

  /* Binarize array. */
  if (binarize == TRUE)
    binarize2D(array, p, p);

  /* Export results. */
  io_2d(filename, array, p, p, csvp, acv, pipe, binaryfile, verbose);
}

/*
 *
 *    PROGRAM INTERFACE TO SIDELNIKOV1D CODE
 *    --------------------------------------
 */
void construction_Sidelnikov1D(char *filename, int p, polynomial primpoly, bool binarize, bool acv, 
			       bool pipe, bool binaryfile, bool verbose, bool csvp) {
  int n, len, *seq;
  polynomial perm;

  /* Create permutation polynomial: poly = 1 + 0 x + ... 0 x^poly.degree */
  poly_deg(perm) = poly_deg(primpoly) - 1;
  poly_list(perm) = (int*) malloc(poly_deg(primpoly)*sizeof(int));
  MEMCHECK(poly_list(perm));
  poly_list(perm)[0] = 1;
  for (n = 1; n <= poly_deg(perm); n++) poly_list(perm)[n] = 0;

  /* Generate sequence. */
  seq = sidelnikov1D(p, poly_deg(primpoly), primpoly, perm);
  len = power(p, poly_deg(primpoly)) - 1;

  /* Binarize result. */
  if (binarize == TRUE) 
    binarize1D(seq, len);

  io_1d(filename, seq, len, csvp, acv, pipe, binaryfile, verbose);
}

/*
 *
 *    PROGRAM INTERFACE TO LEGENDRE SEQUENCE CONSTRUCTION
 *    ---------------------------------------------------
 */
void construction_Legendre1D(char *filename, int p, bool acv, bool pipe, bool binaryfile, 
			     bool verbose, bool csvp) {
  int *seq;

  /* Generate sequence. */
  seq = legendre1D(p);

  /* Export results. */
  io_1d(filename, seq, p, csvp, acv, pipe, binaryfile, verbose);
}

/*
 *
 *    PROGRAM INTERFACE TO TIRKEL-HALL 2D ARRAY FAMILY GENERATOR
 *    ----------------------------------------------------------
 */
void construction_Tirkel_Hall2D_family(char *filename, int p, int shift_poly_degree, bool acv,
				       int zeromap, bool pipe, bool binaryfile, 
				       bool verbose, bool csvp) {
  int k, j, *cl, *clp, npad, family_size;
  char temp[128], textfile[128], binfile[128], index[16], sprime[20], pattern[20];
  polynomial shift_poly;

  if (shift_poly_degree <= 0) 
    shift_poly_degree = 2;

  family_size = power(p, shift_poly_degree - 1) - 1;

  strcpy(textfile, default_filename);
  strcat(textfile, default_txt_extension);

  strcpy(binfile, default_filename);
  strcat(binfile, default_bin_extension);

  if (verbose) printf("FAMILY SIZE = %i\n", family_size);

  memset(index, '\0', sizeof(index));
  memset(temp, '\0', sizeof(temp));

  /* Initialise the n-tuples counter to zeros. */
  clp = (int*) malloc((shift_poly_degree - 1)*sizeof(int));
  MEMCHECK(clp);
  for (k = 0; k < shift_poly_degree - 1; k++) 
    clp[k] = 0;

  clp[0] = 1;

  /* Initialise array of coefficients of the shift polynomial to 1. */
  cl = (int*) malloc((shift_poly_degree + 1)*sizeof(int));
  MEMCHECK(cl);
  for (k = 0; k < shift_poly_degree; k++) 
    cl[k] = 0;

  cl[0] = 1; 
  cl[shift_poly_degree - 1] = 1;
  cl[shift_poly_degree] = 1;

  poly_list(shift_poly) = cl;

  /* Set the degree of the shift polynomial. */
  poly_deg(shift_poly) = shift_poly_degree;

  sprintf(sprime, "_Tirkel2D_p%d_deg%d", p, shift_poly_degree);

  for (k = 1; k <= family_size; k++) {
    /* Generate file name. */  
    npad = intlength(family_size) - intlength(k);
    strcpy(pattern, "_");
    for (j = 0; j < npad; j++) 
      strcat(pattern, "0");
    strcat(pattern, "%d");
    sprintf(index, pattern, k);

    if ( strcmp(filename, textfile) == 0 ) {
      /* Write file name as, for example: arraygen_p7_deg2_392.txt */
      strcpy(temp, default_filename);
      strcat(temp, sprime);
      strcat(temp, index);
      strcat(temp, default_txt_extension);
    } else if ( strcmp(filename, binfile) == 0 ) {
      /* Write file name as, for example: arraygen_p7_deg2_392.bin */
      strcpy(temp, default_filename);
      strcat(temp, sprime);
      strcat(temp, index);
      strcat(temp, default_bin_extension);
    } else {
      strcpy(temp, filename);
      strcat(temp, index);
    }

    if (verbose) printf("\nFILENAME = %s\n", temp);
    if (verbose) {
      printf("SHIFT POLYNOMIAL = ");
      poly_print(shift_poly);
      printf("\n");
    }

    /* Create the Tirkel Hall 2D array. */
    construction_Tirkel_Hall2D(temp, p, zeromap, shift_poly, acv, pipe, binaryfile, verbose, csvp);

    /* Increment the coefficient list of the polynomial. NOTE: We do not touch the 
       coefficients of x^0 and x^(n-1) (where the degree of the polynomial is n). */
    reverse(clp, 0, shift_poly_degree - 1);
    increment(clp, shift_poly_degree, p); 
    reverse(clp, 0, shift_poly_degree - 1);

#if 0
    printf("\nclp = ");
    for (j = 0; j < shift_poly_degree - 1; j++) printf("%i ", clp[j]);
    printf("\n");
#endif

    /* Coefficients of x^0 and x^(n-1) are fixed. */
    for (j = 1; j < shift_poly_degree - 1; j++)
      cl[j] = clp[j - 1];
    cl[shift_poly_degree] = clp[shift_poly_degree - 2];
    poly_list(shift_poly) = cl;
  }
}

/*
 *
 *    PROGRAM INTERFACE TO TIRKEL-HALL 2D CODE
 *    ----------------------------------------
 */
void construction_Tirkel_Hall2D(char *filename, int p, int zeromap, polynomial shiftpoly, 
				bool acv, bool pipe, bool binaryfile, bool verbose, bool csvp) {
  int **array; 

  /* Compute array. */
  array = Tirkel_Hall_2D(p, shiftpoly, zeromap);

  /* Export results. */
  io_2d(filename, array, p, p, csvp, acv, pipe, binaryfile, verbose);
}

/*
 *
 *    PROGRAM INTERFACE TO BLAKE-TIRKEL 3D ARRAY FAMILY GENERATOR
 *    -----------------------------------------------------------
 */

/* The coefficients c_3, c_5, and c_6 in 

        q_2(i,j) = c_1 + c_2 j + c_3 j^2 + c_4 i + c_5 i j + c_6 j^2

   take all combinations in Z_p. The coefficients c_1, c_2 and c_4 are fixed. */

void construction_Blake_Tirkel3D_special_family(char *filename, int p, int zeromap, bool pipe, bool acv,
						bool binaryfile, bool zip, bool verbose, bool csvp) {
  int **coeffs, k = 0, j, npad, c3, c5, c6, x, polydegX = 2, polydegY = 2, 
    conj_fam_size = power(p,3) - power(p,2);
  bivariate_polynomial shift_poly;
  char temp[128], textfile[128], binfile[128], index[100], sprime[100], pattern[100];
  bool solnp;

  /* Set random seed. */
  srand(time(NULL));

  strcpy(textfile, default_filename);
  strcat(textfile, default_txt_extension);

  strcpy(binfile, default_filename);
  strcat(binfile, default_bin_extension);

  /* Initialise shift polynomial. */
  poly_degX(shift_poly) = polydegX;
  poly_degY(shift_poly) = polydegY;

  /* Allocate memory for array of coefficients. */
  coeffs = (int**) malloc((1 + polydegX)*sizeof(int*));
  MEMCHECK(coeffs);
  for (x = 0; x <= polydegX; x++) {
    coeffs[x] = (int*) malloc((1 + polydegY)*sizeof(int));
    MEMCHECK(coeffs[x]);
  }

  poly_array(shift_poly) = coeffs;

  /* The coefficient array is of the form
   
      c_1  c_2  c_3
      c_4  c_5  0
      c_6  0    0
*/  
  
  /* Assign the fixed coefficients. */
  coeffs[1][2] = 0;
  coeffs[2][1] = 0;
  coeffs[2][2] = 0;

  sprintf(sprime, "_Blake_Tirkel_Special_3D_p%d", p);

  /* Iterate through all 3-tuples of c3, c5, c6. */
  for (c3 = 0; c3 < p; c3++) {
    for (c5 = 0; c5 < p; c5++) {
      for (c6 = 0; c6 < p; c6++) {
	// printf("\n c3 = %i, c5 = %i, c6 = %i", c3, c5, c6);
	solnp = system_solver_2D(p, c3, c5, c6);
	if ( ! solnp ) {

          /* Lexicographic array number. */
	  k++; 

	  /* Randomly assign arbitrary coefficints. */
	  coeffs[0][0] = random_in_range(0, p - 1); /* c_1 */
	  coeffs[0][1] = random_in_range(0, p - 1); /* c_2 */
	  coeffs[1][0] = random_in_range(0, p - 1); /* c_4 */

	  /* Generate file name. */
	  npad = intlength(conj_fam_size) - intlength(k);
	  strcpy(pattern, "_");
	  for (j = 0; j < npad; j++) 
	    strcat(pattern, "0");
	  strcat(pattern, "%d");
	  sprintf(index, pattern, k);

	  if ( strcmp(filename, textfile) == 0 ) {
	    /* .txt */
	    strcpy(temp, default_filename);
	    strcat(temp, sprime);
	    strcat(temp, index);
	    strcat(temp, default_txt_extension);
	  } else if ( strcmp(filename, binfile) == 0 ) {
	    /* .bin */
	    strcpy(temp, default_filename);
	    strcat(temp, sprime);
	    strcat(temp, index);
	    strcat(temp, default_bin_extension);
	  } else {
	    strcpy(temp, filename);
	    strcat(temp, index);
	  }

	  if (verbose) printf("\nFILENAME = %s\n", temp);
	  
	  /* Update coefficients of shift polynomial */
	  coeffs[0][2] = c3;
	  coeffs[1][1] = c5;
	  coeffs[2][0] = c6;

	  if (verbose) {
	    printf("ARRAY NUMBER: %i", k);
	    printf("\nSHIFT POLYNOMIAL = ");
	    poly_print_2D(shift_poly);
	    printf("\n\n");
	  }

	  /* Create the Blake-Tirkel 3D array. */
	  construction_Blake_Tirkel3D(temp, p, zeromap, shift_poly, 
				      acv, pipe, binaryfile, zip, verbose, csvp);
	}
      }
    }
  }

  if (verbose) printf("FAMILY SIZE FOR P = %i IS: %i\n", p, k);
}

/*
 *
 *    PROGRAM INTERFACE TO THE BLAKE 3D CODE
 *    --------------------------------------
 */

void construction_Blake_Tirkel3D(char *filename, int p, int zeromap, bivariate_polynomial shiftpoly, 
				 bool acv, bool pipe, bool binaryfile, bool zip, bool verbose, bool csvp) {
  int ***array; 

  /* Compute array. */
  array = Blake_Tirkel_3D(p, shiftpoly, zeromap);

  /* Export results. */
  io_3d(filename, array, p, p, p, csvp, FALSE, acv, pipe, binaryfile, zip, verbose, FALSE);
}

/*
 *
 *    IO FOR 1D ARRAYS
 *    ----------------
 */
void io_1d(char *filename, int *array, int d0, bool csvp, bool acv, bool pipe, 
	   bool binaryfile, bool verbose) {
  int *corr;
  double start, stop, elapsed;

  /* Compute autocorrelations. */
  if (acv == TRUE) {
    if (verbose) 
      printf("\nCOMPUTING THE AUTOCORRELATION...\n");
    start = clock();
    corr = acv_1d_slow(array, d0);
    stop = clock();
    elapsed = ((double) (stop - start)) / CLOCKS_PER_SEC;
    if (verbose) 
      printf("FINISHED COMPUTING THE AUTOCORRELATION AFTER %6.3f (seconds)\n\n", elapsed);
  }
  /* Pipe results? */
  if (pipe == TRUE) {
      print_array_1D(array, d0); /* Write arrays to stdout */
      if (acv == TRUE) {
	printf("\n\n");
	print_array_1D(corr, d0);
	printf("\n");
	free(corr);
      }
  } else {
    /* Write sequence to file. */
    if (binaryfile == TRUE)
      write_binary_array_1D(filename, array, d0); 
    else
      write_array_1D(filename, csvp, array, d0); 
    if (verbose == TRUE) 
      print_array_1D(array, d0);
    /* Write autocorrelations to file. */
    if (acv == TRUE) {
      write_array_1D("arraygen.acv", csvp, corr, d0);
      if (verbose == TRUE) {
	printf("\n\n");
	print_array_1D(corr, d0);
	printf("\n");
      }
      free(corr);
    }
  }

  free(array);
}

/*
 *
 *    IO FOR 2D ARRAYS
 *    ----------------
 */
void io_2d(char *filename, int **array, int d0, int d1, bool csvp, 
	   bool acv, bool pipe, bool binaryfile, bool verbose) {
  int n, m, **corr, maxcorr;
  double start, stop, elapsed;

  /* Compute autocorrelations */
  if (acv == TRUE) {
    if (verbose) 
      printf("\nCOMPUTING THE AUTOCORRELATION...\n");
    start = clock();
    corr = acv_2d_slow(array, d0, d1);
    stop = clock();
    elapsed = ((double) (stop - start)) / CLOCKS_PER_SEC;
    if (verbose) 
      printf("FINISHED COMPUTING THE AUTOCORRELATION AFTER %6.3f (seconds)\n\n", elapsed);
  }

  /* Output results. */
  if (pipe == TRUE) {
    print_array_2D(array, d0, d1); /* Write arrays to stdout */
    if (acv == TRUE) {
      printf("\n\n");
      print_array_2D(corr, d0, d1);
      printf("\n");
      /* Free array of correlations. */
      for(n = 0; n < d0; n++) 
	free(corr[n]);
      free(corr);
    }
  } else {
    /* Write array to file. */
    if (binaryfile == TRUE)
      write_binary_array_2D(filename, array, d0, d1); 
    else
      write_array_2D(filename, csvp, array, d0, d1); 

    if (verbose == TRUE) pretty_print_array_2D(array, d0, d1);
    /* Write autocorrelations to file. */
    if (acv == TRUE) {
      write_array_2D("arraygen.acv", csvp, corr, d0, d1);
      if (verbose == TRUE) {
	printf("\n\n");
	pretty_print_array_2D(corr, d0, d1);
	printf("\n");

	/* Compute the maximum off-peak non-zero autocorrelation. */
	maxcorr = INT_MIN;

	for (n = 0; n < d0; n++) {
	  for (m = 0; m < d1; m++) {
	    if (corr[n][m] > maxcorr && !(n == 0 && m == 0))
		maxcorr = corr[n][m];
	  }
	}

	printf("\n");
	printf("MAX OFF-PEAK CORRELATION = %i, %.2f%% OF PEAK    (PEAK = %i)", 
	       maxcorr, 100.0*((float) maxcorr)/((float) corr[0][0]), corr[0][0]);
	printf("\n");
      }
      /* Free array of correlations. */
      for(n = 0; n < d0; n++) 
	free(corr[n]);
      free(corr);
    }
  }

  /* Free memory. */
  for(n = 0; n < d0; n++) 
    free(array[n]);
  free(array);
}

/*
 *
 *    IO FOR 3D ARRAYS
 *    ----------------
 */
void io_3d(char *filename, int ***array, int d0, int d1, int d2, 
	   bool csvp, bool binarize, bool acv, bool pipe, 
	   bool binaryfile, bool zip, bool verbose, bool lincomp) {
  int n, m, k, ***corr, maxcorr;
  char gzip_str[128];
  double start, stop, elapsed;

  /* Compute autocorrelations */
  if (acv == TRUE) {
    if (binarize == FALSE) 
      binarize3D(array, d0, d1, d2);
    if (verbose) 
      printf("\nCOMPUTING THE AUTOCORRELATION...\n");
    start = clock();
    corr = acv_3d_slow(array, d0, d1, d2);
    stop = clock();
    elapsed = ((double) (stop - start)) / CLOCKS_PER_SEC;
    if (verbose) 
      printf("FINISHED COMPUTING THE AUTOCORRELATION AFTER %6.3f (seconds)\n\n", elapsed);
  }

  /* Output results. */
  if (pipe == TRUE) {
    print_array_3D(array, d0, d1, d2);
    if (acv == TRUE) {
      printf("\n\n");
      pretty_print_array_3D(corr, d0, d1, d2);
      printf("\n");
      /* Free array of correlations. */
      for(n = 0; n < d0; n++) {
	for (m = 0; m < d1; m++) 
	  free(corr[n][m]);
	free(corr[n]);
      }
      free(corr);
    }
  } else {
    /* Write array to file. */
    if (binaryfile == TRUE)
      write_binary_array_3D(filename, array, d0, d1, d2); 
    else
      write_array_3D(filename, csvp, array, d0, d1, d2); 

    /* Compress output file? */
    if (zip == TRUE) {
      strcpy(gzip_str, "gzip -f -9 ");
      strcat(gzip_str, filename);
      if (verbose) printf("\nRUNNING THE EXTERNAL COMMAND: %s\n", gzip_str);
      system(gzip_str);
    }
    /* Write array to stdout? */
    if (verbose) pretty_print_array_3D(array, d0, d1, d2);
    /* Write autocorrelations to file. */
    if (acv == TRUE) {
      write_array_3D("arraygen.acv", csvp, corr, d0, d1, d2);
      /* Compress output file? */
      if (zip == TRUE)
	system("gzip -f -9 arraygen.acv");
      /* Print array of correlations to stdout? */
      if (verbose == TRUE) {
	printf("\n\n");
	pretty_print_array_3D(corr, d0, d1, d2);
	
	/* Compute the maximum off-peak non-zero autocorrelation. */
	maxcorr = INT_MIN;

	for (n = 0; n < d0; n++) {
	  for (m = 0; m < d1; m++) {
	    for (k = 0; k < d2; k++) {
	      if (corr[n][m][k] > maxcorr && !(n == 0 && m == 0 && k == 0))
		maxcorr = corr[n][m][k];
	    }
	  }
	}

	printf("\n");
	printf("MAX OFF-PEAK CORRELATION = %i, %.2f%% OF PEAK    (PEAK = %i)", 
	       maxcorr, 100.0*((float) maxcorr)/((float) corr[0][0][0]), corr[0][0][0]);
	printf("\n");
      }
    }
    
    /* Compute the 1D linear complexity of the array. */
    if (lincomp == TRUE) {
      if (binarize == FALSE) 
	binarize3D(array, d0, d1, d2);
      /* Replace any zeros with 1's, then map the array to a 0,1 array. */
      for (n = 0; n < d0; n++) {
	for (m = 0; m < d1; m++) {
	  for (k = 0; k < d2; k++) {
	    if (array[n][m][k] == 0) array[n][m][k] = 1;
	    if (array[n][m][k] == 1) 
	      array[n][m][k] = 0;
	    else
	      array[n][m][k] = 1;
	  }
	}
      }
      /* linear complexity with y, z constant. */
      
    }
  }

  /* Free memory. */
  for(n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) 
      free(array[n][m]);
    free(array[n]);
  }
  free(array);
}

/*
 *
 *    CRT-TYPE ARRAY ENUMERATION
 *    --------------------------
 */
int* CRTenumerate2D(int **array, int d0, int d1) {
  int n, *seq, len = d0*d1; 

  if ( gcd(d0, d1) != 1 ) {
    printf("ERROR: CRT-type array enumeration is only valid when\n\
the array dimensions are coprime.");
    exit(1);
  }

  seq = (int*) malloc(len*sizeof(int));
  MEMCHECK(seq);

  for (n = 0; n < len; n++) 
    seq[n] = array[n%d0][n%d1];

  return seq;
}

/*
 *
 *    CONSTRUCTION OF M SEQUENCES
 *    ---------------------------
 */

int* m_sequence(polynomial primpoly, int p) {
  int *state, *recurr, *mseq, deg, len, n, m, b;

  /* Initial state of the LFSR. */
  deg = poly_deg(primpoly);
  state = (int*) malloc(deg*sizeof(int));

  state[0] = 1;
  for(n = 1; n < deg; n++)
    state[n] = 0;

  /* construct recurrence coefficients.  */
  recurr = (int*) malloc(deg*sizeof(int));

  for (n = 0; n < deg; n++) {
    recurr[n] = (-1*poly_list(primpoly)[deg - n - 1])%p;
    if (recurr[n] < 0) recurr[n] += p;
  }

  /* Length of the M sequence. */
  len = power(p, deg) - 1;

  /* Allocate memory for the sequence. */
  mseq = (int*) malloc(len*sizeof(int));

  /* Simulate the LFSR. */
  for (n = 0; n < len; n++) {

    b = 0;
    for (m = 0; m < deg; m++) {
      b += recurr[m]*state[m];
    }
    b = b%p;
    if (b < 0) b += p;

    /* shift */
    for (m = deg - 1; m > 0; m--) {
      state[m] = state[m - 1];
    }

    state[0] = b;
    mseq[n] = state[deg - 1];
  }

  free(state);
  free(recurr);
  return mseq;
}

/*
 *
 *    IPATOV'S PHI - EXTENDED BINARY CHARACTER
 *    ----------------------------------------
 */

/* Reference: eqn. 6.18, pp. 171, Spread Spectrum and CDMA - Principles
              and Applications, V. Ipatov, Wiley, 2005.  */

int ipatov_psi(int a, int p) {
  int *qr, n, psi;

  if (a%p == 0) 
    return -1; // Extension from page 241. 

  /* Compute quadratic residues. */
  qr = quadratic_residues(p);

  psi = -1;

  for (n = 0; n < p; n++) {
    if (qr[n] == a) {
      psi = 1;
      break ;
    }
  }

  free(qr);

  return psi;
}


/*
 *
 *    CONSTRUCTION A (2D) ARRAYS
 *    --------------------------
 */
int** A_array_2D(int p, polynomial shiftpoly, int map) {
  int n, m, shift, pr, *legendre, **array2d;
  
  if (! primep(p))
    {
      printf("\nERROR: %i is not prime.\n", p);
      exit(1);
    }

  if (map != 0 && (p + 1)%4 != 0) {
    printf("\nERROR: %i is not of the form 4k - 1.", p);
    exit(1);
  }
  
  if (! ( map == 1 || map == -1 || map == 0 )) {
    printf("\nERROR: mapping should be 0, -1, or 1.\n");
    exit(1);
  }

  /* Construct the Legendre sequence. */
  legendre = legendre1D(p);

  /* Compute the primitive root of p. */
  pr = primitive_root(p);

  /* Allocate memory for p x p - 1 array. */
  array2d = (int**) malloc(p*sizeof(int*));
  MEMCHECK(array2d);
  for (n = 0; n < p; n++) {
    array2d[n] = (int*) malloc((p - 1)*sizeof(int));
    MEMCHECK(array2d[n]);
  }

  for (n = 0; n < p - 1; n++) {
    shift = poly_eval(shiftpoly, power_mod(pr, n, p), p);
    // printf("\npower = %i, shift = %i\n", power_mod(pr, n, p), shift);
    for (m = 0; m < p; m++) {
      array2d[m][n] = legendre[(m + shift)%p];
      // printf("%i ", array2d[m][n]);
      if (map != 0 && array2d[m][n] == 0)
        array2d[m][n] = map;
    }
  }

  return array2d;
}


/*
 *
 *    BLAKE TIRKEL 3D ARRAYS
 *    ----------------------
 */

int*** Blake_Tirkel_3D(int p, bivariate_polynomial shiftpoly, int map) {
  int n, m, k, s, *legendre, ***array3d;

  if (! primep(p))
    {
      printf("\nERROR: %i is not prime.\n", p);
      exit(1);
    }

  if (map != 0 && (p + 1)%4 != 0) {
    printf("\nERROR: %i is not of the form 4k - 1.", p);
    exit(1);
  }

  if (! ( map == 1 || map == -1 || map == 0 )) {
    printf("\nERROR: mapping should be 0, -1, or 1.\n");
    exit(1);
  }

  /* Construct the Legendre sequence. */
  legendre = legendre1D(p);

  /* Allocate memory for p x p x p array. */
  array3d = (int***) malloc(p*sizeof(int**));
  MEMCHECK(array3d);

  for (n = 0; n < p; n++) {
    array3d[n] = (int**) malloc(p*sizeof(int*));
    MEMCHECK(array3d[n]);

    for (m = 0; m < p; m++) {
      array3d[n][m] = (int*) malloc(p*sizeof(int));
      MEMCHECK(array3d[n][m]);
    }
  }

  /* Construct the array. */
  for (n = 0; n < p; n++) {
    for (m = 0; m < p; m++) {
      /* Evaluate shift polynomial at x == n, y == m. */
      s = poly_eval_2D(shiftpoly, n, m, p);
      for (k = 0; k < p; k++) {
	array3d[k][m][n] = legendre[(k + s) % p]; // periodic shift.
	/* Apply the mapping. */
	if (map != 0 && array3d[k][m][n] == 0) 
	  array3d[k][m][n] = map;
      }
    }
  }

  /* Free memory. */
  free(legendre);
  
  return array3d;
}


/*
 *
 *    TIRKEL HALL 2D ARRAYS
 *    ---------------------
 */
int** Tirkel_Hall_2D(int p, polynomial shiftpoly, int map) {
  int *legendre, **array2d, n, m, s;

  if (! primep(p))
    {
      printf("\nERROR: %i is not prime.\n", p);
      exit(1);
    }

  if (map != 0 && (p + 1)%4 != 0) {
    printf("\nERROR: %i is not of the form 4k - 1.", p);
    exit(1);
  }

  if (! ( map == 1 || map == -1 || map == 0 )) {
    printf("\nERROR: mapping should be 0, -1, or 1.\n");
    exit(1);
  }

  /* Construct the Legendre sequence. */
  legendre = legendre1D(p);

  /* Allocate memory for p x p array. */
  array2d = (int**) malloc(p*sizeof(int*));
  MEMCHECK(array2d);

  for (n = 0; n < p; n++) {
    array2d[n] = (int*) malloc(p*sizeof(int));
    MEMCHECK(array2d[n]);
  }

  /* Construct the array. */
  for (n = 0; n < p; n++) {
    /* Evaluate shift polynomial at x == n. */
    s = poly_eval(shiftpoly, n, p);
    for (m = 0; m < p; m++) {
      array2d[m][n] = legendre[(m + s) % p]; // periodic shift.
      /* Apply the mapping. */
      if (map != 0 && array2d[m][n] == 0) 
	array2d[m][n] = map;
    }
  }

  free(legendre);

  return array2d;
}

/*
 *
 *    LEGENDRE SEQUENCES
 *    ------------------
 */
int* legendre1D(int p) {
  int *qr, *ls, n, m;

  qr = quadratic_residues(p);

#if 0
  /* Print quadratic residues. */
  printf("\n qr = ");
  for (n = 0; n < p; n++) 
    printf("%i ", qr[n]);
  printf("\n\n");
#endif

  ls = (int*) malloc(p*sizeof(int));

  MEMCHECK(ls);

  ls[0] = 0;

  for (n = 1; n < p; n++) {
    for (m = 0; m < p; m++) {
      if (n == qr[m]) {
	ls[n] = 1;
	goto cludge;
      }
    }
    ls[n] = -1;
  cludge: ;
  }

  free(qr);
  
  return ls;
}

/*
 *
 *    QUADRATIC RESIDUES
 *    ------------------
 */

int* quadratic_residues(int p) {
  /* Creates a (non-distinct) list of quadratic residues of p. */
  int n, *qr;

  /* p is assumed to be prime. (It should have been checked further up the tree.) */

  qr = (int*) malloc(p*sizeof(int));

  MEMCHECK(qr);

  qr[0] = 1; // hack
  for (n = 1; n < p; n++) 
    qr[n] = power_mod(n, 2, p); // Redundant for quadratics. 

  return qr;
}

/*
 *
 *    MODULAR EXPONENTIATION
 *    ----------------------
 */
int power_mod(int base, int exponent, int mod) {
  int rem = 1;

  while (exponent > 0) {
    if (exponent%2 == 1)
      rem = (rem*base) % mod;
    exponent /= 2;
    base = (base*base) % mod;
  }
  return rem;
}

/*
 *
 *    SMALLEST PRIMITIVE ROOT
 *    -----------------------
 */
void primitive_root_disp(int n) {
  int pr; 
  pr = primitive_root(n);
  printf("\ng(%i) = %i\n\n", n, pr);
}

int primitive_root(int n) {
  int *factors, s, a, k, nfacs = 0, prev; 

  if (n == 2) 
    return 1;

  if (n == 4) 
    return 3;

  s = euler_phi(n);
  factors = factor_integer(s);
  
  for (k = 0; k < MAX_NUMBER_OF_PRIME_FACTORS; k++) {
    if (factors[k] != 0) 
      nfacs++;
  }

  for (a = 2; a < n; a++) {
    prev = 0;
    for (k = 0; k < nfacs; k++) {
      if (factors[k] == prev) 
        continue;
      if (power_mod(a, s/factors[k], n) == 1)
        goto next;
      prev = factors[k];
    }
    return a; 
    next: ;
  }

  return 0;
}

/*
 *
 *    PRIME FACTORISATION
 *    -------------------
 */
void factor_integer_disp(int n) {
  int i = 0, *factors;
  factors = factor_integer(n);
  printf("\n%i = ", n);
  while (factors[i] != 0) {
    printf("%i", factors[i++]);
    if (factors[i] != 0)
      printf("*");
  }
  printf("\n\n");
}

int* factor_integer(int n) {
  int i, *factor_list, k = 0;

  factor_list = (int*) malloc(MAX_NUMBER_OF_PRIME_FACTORS*sizeof(int));

  for (i = 0; i < MAX_NUMBER_OF_PRIME_FACTORS; i++)
    factor_list[i] = 0;

  /* Factor out 2. */
  while (n%2 == 0) {
    n = n/2;
    factor_list[k++] = 2;
  }

  for (i = 3; i*i <= n; i += 2) {
    while (n%i == 0) {
      factor_list[k++] = i;
      n /= i;
    }
  }

  if (n > 2) 
    factor_list[k++] = n;

  return factor_list;
}

/*
 *
 *    EULER PHI
 *    ---------
 */

void euler_phi_disp(int n) {
  int phi;  
  phi = euler_phi(n);
  printf("\nphi(%i) = %i\n\n", n, phi);
}

int euler_phi(int n) {
  int k, phi = 1;

  /* This is slow. It should work from the prime factorisation of n. */

  if (n == 0) 
    return 0;

  for(k = 2; k < n; k++) {
    if (gcd(k, n) == 1) 
      phi++;
  }
  return phi;
}

/* 
 *
 *    GCD
 *    ---
 */

int gcd(int u, int v) {
  int t;
  while (v) {
    t = u; 
    u = v; 
    v = t % v;
  }
  return u < 0 ? -u : u;
}

/*
 *
 *    3D CONSTRUCTION C1
 *    ------------------
 */

int*** C1_array_3D(int p, polynomial prim_poly, polynomial shift_poly) {
  int one[1] = {1};
  polynomial one_poly;

  /* Set the denominator to 1. */
  poly_deg(one_poly) = 0;
  poly_list(one_poly) = one;

  return C2_array_3D(p, prim_poly, shift_poly, one_poly); 
}

/*
 *
 *    3D CONSTRUCTION C2
 *    ------------------
 */
int*** C2_array_3D(int p, polynomial prim_poly, polynomial numpoly, polynomial denpoly) {
  int ***array3d, **legendre, *sidelnikov, d = poly_deg(prim_poly), 
    n, i, j, k, tau, len;
  polynomial perm_poly;

  /* Sanity checks. */
  if (poly_deg(numpoly) < 1 || poly_deg(prim_poly) < 1) {
    printf("\nERROR: malformed polynomials in Construction C1.\n");
    exit(1);
  }

  /* Construct a 2D Legendre array of size p x p. */
  legendre = legendre2Dgeneralised(p, prim_poly, numpoly, denpoly, TRUE);
  
  /* Construct a 1D Sidenlikov sequence of length p^d - 1. */
  len = power(p,d) - 1;
  /* Create permutation polynomial. */
  poly_deg(perm_poly) = poly_deg(prim_poly);
  poly_list(perm_poly) = (int*) malloc(len*sizeof(int));
  poly_list(perm_poly)[0] = 1;
  for (n = 1; n < len; n++) 
    poly_list(perm_poly)[n] = 0;
  sidelnikov = sidelnikov1D(p, d, prim_poly, perm_poly);

  /* Allocate memory for 3D array. */
  /* p x p x p^d - 1 */
  array3d = (int***) malloc(p*sizeof(int**));
  MEMCHECK(array3d);

  for (i = 0; i < p; i++) {
    array3d[i] = (int**) malloc(p*sizeof(int*));
    MEMCHECK(array3d[i]);

    for (j = 0; j < p; j++) {
      array3d[i][j] = (int*) malloc(len*sizeof(int));
      MEMCHECK(array3d[i][j]);

      for (k = 0; k < len; k++) 
	array3d[i][j][k] = 0;
    }
  }

  /* Construct 3D array. */
  for (i = 0; i < p; i++) {
    for (j = 0; j < p; j++) {
      tau = legendre[i][j];
      if (tau == NEGINF)
	for (k = 0; k < len; k++)
	  array3d[i][j][k] = 0;
      else
	for (k = 0; k < len; k++)
	  array3d[i][j][k] = sidelnikov[(k + tau)%len];
    }
  }
  
  return array3d;
}

/*
 *
 *    3D CONSTRUCTION D2
 *    ------------------
 */
int*** D2_array_3D(int p, polynomial prim_poly, polynomial shift_poly) {
  int pp = p*p, n, i, j, k, h, v, **legendre, ***array3d, zero_coeffs[2];
  gf alpha, *alpha_coeffs, temp1, temp2, temp3, shift;

  /* Sanity checks. */
  if (poly_deg(shift_poly) < 1 || poly_deg(prim_poly) < 1) {
    printf("\nERROR: malformed polynomials in Construction D2.\n");
    exit(1);
  }

  /* Allocate memory for 3D array. */
  /* p x p x p^2 - 1 */
  array3d = (int***) malloc(p*sizeof(int**));
  MEMCHECK(array3d);

  for (i = 0; i < p; i++) {
    array3d[i] = (int**) malloc(p*sizeof(int*));
    MEMCHECK(array3d[i]);

    for (j = 0; j < p; j++) {
      array3d[i][j] = (int*) malloc(pp*sizeof(int));
      MEMCHECK(array3d[i][j]);

      for (k = 0; k < pp; k++) 
	array3d[i][j][k] = 0;
    }
  }  

  /* Construct a 2D Legendre array of size p x p. */
  legendre = legendre2D(p, prim_poly);

  /* Construct the 2D array shifts from the polynomial form of powers of 
     alpha of the form a + b \alpha mapping to a shift of 'a' in the 
     vertical and 'b' in the horizontal. */

  alpha = gf_create(p, 2, poly_list(prim_poly), power_table[1]);

  /* Compute the coefficient list of the shift sequence. */
  alpha_coeffs = (gf*) malloc((1 + poly_deg(shift_poly))*sizeof(gf));

  for (n = 0; n <= poly_deg(shift_poly); n++)
    alpha_coeffs[n] = gf_power(alpha, poly_list(shift_poly)[n]);

  for (k = 0; k < pp - 1; k++) {
    /* Set shift to zero. */
    zero_coeffs[0] = 0; zero_coeffs[1] = 0;
    shift = gf_create(p, poly_deg(prim_poly), poly_list(prim_poly), zero_coeffs);

    /* Compute the shift polynomial. */
    for (n = 0; n <= poly_deg(shift_poly); n++)
      {
	temp1 = gf_power(alpha, k*n);
	temp2 = gf_times(temp1, alpha_coeffs[n]);
	gf_delete(temp1);
	temp3 = gf_plus(shift, temp2);
	gf_delete(temp2);
	gf_delete(shift);
	shift = gf_copy(temp3);
	gf_delete(temp3);
      }

    /* Extract the horizontal and vertical shifts. */
    h = gf_poly_list(shift)[0];
    v = gf_poly_list(shift)[1];
    gf_delete(shift);

    // printf("\nh,v = %i, %i\n", h, v);

    /* 2D periodic shift */
    for (i = 0; i < p; i++)
      for (j = 0; j < p; j++) 
	array3d[i][j][k] = legendre[(i+v)%p][(j+h)%p];
  }

  /* Free powers of alpha. */
  gf_delete(alpha);
  for (n = 0; n <= poly_deg(shift_poly); n++) 
    gf_delete(alpha_coeffs[n]);

  /* Free Legendre array. */
  for (i = 0; i < p; i++) 
    free(legendre[i]);
  free(legendre);

  /* Free coefficient list of the shift sequence. */
  free(alpha_coeffs);

  return array3d;
}

/*
 *
 *    3D CONSTRUCTION D3
 *    ------------------
 */

int*** D3_array_3D(int p, polynomial prim_poly, polynomial numpoly, polynomial denpoly) {
  int pp = p*p, i, j, k, h, v, **legendre, ***array3d;
  polynomial shift;

  /* Sanity checks. */
  if (poly_deg(numpoly) < 1 || poly_deg(prim_poly) < 1) {
    printf("\nERROR: malformed polynomials in Construction D2.\n");
    exit(1);
  }

  /* Allocate memory for 3D array. */
  /* p x p x p^2 - 1 */
  array3d = (int***) malloc(p*sizeof(int**));
  MEMCHECK(array3d);

  for (i = 0; i < p; i++) {
    array3d[i] = (int**) malloc(p*sizeof(int*));
    MEMCHECK(array3d[i]);

    for (j = 0; j < p; j++) {
      array3d[i][j] = (int*) malloc(pp*sizeof(int));
      MEMCHECK(array3d[i][j]);

      for (k = 0; k < pp; k++) 
	array3d[i][j][k] = 0;
    }
  }

  /* Construct a 2D Legendre array of size p x p. */
  legendre = legendre2D(p, prim_poly);

  /* Construct the 2D array shifts from the polynomial form of powers of 
     alpha of the form a + b \alpha mapping to a shift of 'a' in the 
     vertical and 'b' in the horizontal. */

  for (k = 0; k < pp - 1; k++) {

    shift = gf_rational_eval(p, k, numpoly, denpoly, prim_poly);

    /* Extract the horizontal and vertical shifts. */
    h = poly_list(shift)[0];
    v = poly_list(shift)[1];
    poly_free(shift);

    // printf("\nh,v = %i, %i\n", h, v);

    if (h == NEGINF || v == NEGINF) {
      /* plane of -oo */
      for (i = 0; i < p; i++)
	for (j = 0; j < p; j++) 
	  array3d[i][j][k] = NEGINF;
    } else {
      /* 2D periodic shift */
      for (i = 0; i < p; i++)
	for (j = 0; j < p; j++) 
	  array3d[i][j][k] = legendre[(i+v)%p][(j+h)%p];
    }
  }

  /* Free Legendre array. */
  for (i = 0; i < p; i++) 
    free(legendre[i]);
  free(legendre);

  return array3d;
}

/*
 *
 *    3D LEGENDRE ARRAYS
 *    ------------------
 */
int*** legendre3D(int p, polynomial primpoly) {
  int one[1] = {1}, x[2] = {0,1};
  polynomial x_poly, one_poly;

  poly_deg(x_poly) = 1;
  poly_list(x_poly) = x;

  poly_deg(one_poly) = 0;
  poly_list(one_poly) = one;

  return legendre3Dgeneralised(p, primpoly, x_poly, one_poly, FALSE);
}

int*** legendre3Dgeneralised(int p, polynomial primpoly, polynomial numpoly, 
			     polynomial denpoly, bool generalised) {
  int ***array, poly_coeffs[3], i, j, k, K;
  polynomial poly, genpoly;

  if (poly_deg(primpoly) != 3)
    {
      printf("\nERROR: legendre3D requires a primitive polynomial of degree 3.\n");
      exit(1);
    }

  /* Initialise power table lookup. */
  generate_power_table(p, primpoly);

  /* Allocate memory. */
  array = (int***) malloc(p*sizeof(int**));
  MEMCHECK(array);

  for (i = 0; i < p; i++) {
    array[i] = (int**) malloc(p*sizeof(int*));
    MEMCHECK(array[i]);

    for (j = 0; j < p; j++) {
      array[i][j] = (int*) malloc(p*sizeof(int));
      MEMCHECK(array[i][j]);
    }
  }

  poly_deg(poly) = 2;

  /* Generate entries of the Legendre array. */
  for (i = 0; i < p; i++) {
    poly_coeffs[0] = i;
    for (j = 0; j < p; j++) {
      poly_coeffs[1] = j; 
      for (k = 0; k < p; k++) {
	poly_coeffs[2] = k;
	poly_list(poly) = poly_coeffs;
        /* Find position of i + j x + k x^2 in the power table. */
        K = power_representation(poly);
	
	if (generalised) {
	  /* Evaluate shift rational function (numpoly/denpoly) at x = \alpha^K. */
	  genpoly = gf_rational_eval(p, K, numpoly, denpoly, primpoly);
	  /* Convert to power of \alpha representation. */
	  K = power_representation(genpoly);
	  if ( neginfq(genpoly) )
	    K = NEGINF;
	  else
	    K = power_representation(genpoly);
	  /* Delete polynomial form of power of \alpha.  */
	  poly_free(genpoly);
	}

	if (K != NEGINF) 
	  K++;
	array[i][j][k] = K;
      }
    }
  }
  return array;
}

/*
 *
 *    2D GENERALISED LEGENDRE ARRAYS
 *    ------------------------------
 */

int** legendre2D(int p, polynomial primpoly) {
  int one[1] = {1}, x[2] = {0,1};
  polynomial x_poly, one_poly;

  poly_deg(x_poly) = 1;
  poly_list(x_poly) = x;

  poly_deg(one_poly) = 0;
  poly_list(one_poly) = one;

  return legendre2Dgeneralised(p, primpoly, x_poly, one_poly, FALSE);
}

int** legendre2Dgeneralised(int p, polynomial primpoly, polynomial numpoly,
			    polynomial denpoly, bool generalised) {
  int **array, poly_coeffs[2], i, j, K;
  polynomial poly, genpoly;


  if (poly_deg(primpoly) != 2)
    {
      printf("\nERROR: legendre2D requires a primitive polynomial of degree 2.\n");
      exit(1);
    }

  /* Initialise power table lookup. */
  generate_power_table(p, primpoly);

  /* Allocate memory. */
  array = (int**) malloc(p*sizeof(int*));
  MEMCHECK(array);

  for (i = 0; i < p; i++) {
    array[i] = (int*) malloc(p*sizeof(int));
    MEMCHECK(array[i]);
  }

  poly_deg(poly) = 1;

  /* Generate entries of the Legendre array. */
  for (i = 0; i < p; i++) {
    poly_coeffs[0] = i;
    for (j = 0; j < p; j++) {
      poly_coeffs[1] = j;
      poly_list(poly) = poly_coeffs;
      /* Find position of i + j x in the power table. */
      K = power_representation(poly);
      
      if (generalised) {
	/* Evaluate shift rational function (numpoly/denpoly) at x = \alpha^K. */
	genpoly = gf_rational_eval(p, K, numpoly, denpoly, primpoly);
	/* Convert to power of \alpha representation. */
	if ( neginfq(genpoly) )
	  K = NEGINF;
	else
	  K = power_representation(genpoly);
	/* Delete polynomial form of power of \alpha.  */
	poly_free(genpoly);
      }
      
      if (K != NEGINF) 
	K++;
      array[i][j] = K;
    }
  }
  return array;
}

/* 
 *
 *    SIDELNIKOV SEQUENCES
 *    --------------------
 */

/* sidelnikov1D generates a length p^d - 1 sequence.  */
int* sidelnikov1D(int p, int d, polynomial primpoly, polynomial perm) {
  int n, *sidelseq;
  gf gf_perm, gfpl, gfp;

  /* Initialise power table lookup. */
  generate_power_table(p, primpoly);

  /* Convert the perm to a Galois field expression. */
  gf_perm = gf_create(p, d, poly_list(primpoly), poly_list(perm));

  /* Generate Sidelnikov1D sequence term-by-term. */
  sidelseq = (int*) malloc(power_table_length*sizeof(int));
  for (n = 0; n < power_table_length; n++) {
    /* Add the permutation to the power of alpha. */
    gfpl = gf_create(p, d, poly_list(primpoly), power_table[n]);
    gfp = gf_plus(gf_perm, gfpl);
    /* Find power of \alpha representation of gfp. */
    sidelseq[n] = power_representation(gf_poly(gfp));
    gf_delete(gfpl);
    gf_delete(gfp);
  }
  gf_delete(gf_perm);
  return sidelseq;
}

/*
 *
 *    EVALUATE A POLYNOMIAL
 *    ---------------------
 */

/* UNIVARIATE */

int poly_eval(polynomial poly, int k, int mod) {
  /* Computes poly(x) at x == k. */
  int n, sum = 0;

  for (n = 0; n <= poly_deg(poly); n++)
    sum += (poly_list(poly)[n]*power(k,n)) % mod;

  return sum;
}

/* BIVARIATE */

int poly_eval_2D(bivariate_polynomial poly, int x, int y, int mod) {
  /* Computes poly(x,y). */
  int n, m, C, xd, yd, sum = 0;
  
  xd = poly_degX(poly);
  yd = poly_degY(poly);

  for (n = 0; n <= xd; n++) {
    for (m = 0; m <= yd; m++) {
      C = poly_array(poly)[n][m];
      sum += C*power(x,n)*power(y,m);
      sum = sum % mod;
    }
  }

  return sum;
}

/*
 *
 *    EVALUATE A RATIONAL FUNCTION IN A GF
 *    ------------------------------------
 */

polynomial gf_rational_eval(int p, int K, polynomial numpoly, polynomial denpoly, 
			    polynomial primpoly) {
  int d = poly_deg(primpoly), *ninf;
  polynomial numeval, deneval, result; 
  gf numgf, dengf, deninv, prod;

  /* Evaluate numerator polynomial at x = \alpha^K. */
  numeval = gf_poly_eval(p, K, numpoly, primpoly);

  /* Optimisation: check if denpoly == 1. */
  if (poly_deg(denpoly) == 0 && poly_list(denpoly)[0] == 1) 
    return numeval;

  /* Evaluate denominator polynomial at x = \alpha^K.  */
  deneval = gf_poly_eval(p, K, denpoly, primpoly);

  // printf("\nden = ");
  // poly_print(deneval);

  if (zeroq(deneval)) 
    {
      poly_deg(result) = 0;
      ninf = (int*) malloc(sizeof(int));
      ninf[0] = NEGINF;
      poly_list(result) = ninf;
      return result;
    }

  /* Express the numerator and denominator in GF representation. */
  numgf = gf_create(p, d, poly_list(primpoly), poly_list(numeval));
  dengf = gf_create(p, d, poly_list(primpoly), poly_list(deneval));

  /* denpoly^-1 == denpoly^(-1 mod p^d - 1) == denpoly^(p^d - 2) */
  deninv = gf_power(dengf, power(p, d) - 2);
  prod = gf_times(numgf, deninv);

  // printf("\nprod = ");
  // gf_print(prod);

  gf_delete(numgf);
  gf_delete(dengf);
  gf_delete(deninv);

  result = poly_copy(gf_poly(prod));
  
  gf_delete(prod);

  return result;
}

/*
 *
 *    EVALUATE A POLYNOMIAL IN A GF
 *    -----------------------------
 */
polynomial gf_poly_eval(int p, int k, polynomial poly, polynomial primpoly) {
  int n = poly_deg(poly), i, *zeros;
  polynomial result;
  gf alpha, temp0, temp1, temp2, temp3, sum;

  /* Initialise power table lookup. */
  generate_power_table(p, primpoly);

  /* Create primitive power of alpha. */
  alpha = gf_create(p, poly_deg(primpoly), poly_list(primpoly), power_table[1]);

  /* Set sum to zero. */
  zeros = (int*) malloc((1 + poly_deg(poly))*sizeof(int));
  for (i = 0; i <= poly_deg(poly); i++) 
    zeros[i] = 0;
  sum = gf_create(p, poly_deg(primpoly), poly_list(primpoly), zeros);

  /* Let poly = [c_0, c_1, ..., c_n] be a polynomial in the variable x, 
     then we compute poly evaluated at x = \alpha^k, that is: 
     \alpha^c_0 + \alpha^c_1 \alpha^k + \alpha^c_2 \alpha^(2 k) + ... + 
     \alpha^c_n \alpha^(n k) */

    for (i = 0; i <= n; i++) {
      /* \alpha^c_i */
      if (poly_list(poly)[i] == NEGINF)
	temp0 = gf_power(alpha, NEGINF); // \alpha^-oo
      else
	temp0 = gf_power(alpha, poly_list(poly)[i]); // \alpha^c_i
      /* \alpha^{k*i} */
      if (k == NEGINF)
	temp1 = gf_power(alpha, NEGINF); // \alpha^-oo
      else
	temp1 = gf_power(alpha, k*i); // \alpha^{k*i}
      /* \alpha^c_i \alpha^{k i} */
      temp2 = gf_times(temp0, temp1);
      gf_delete(temp0);
      gf_delete(temp1);
      /* sum += \alpha^c_i \alpha^{k i} */
      temp3 = gf_plus(sum, temp2);
      gf_delete(temp2);
      gf_delete(sum);
      sum = gf_copy(temp3);
      gf_delete(temp3);
    }

    gf_delete(alpha);
    result = poly_copy(gf_poly(sum)); // Copy the sum polynomial data structure to result.
    gf_delete(sum); // Delete the sum. 

    return result;
}

/*
 *
 *    DISPLAYING MULTIPLE RANDOM PRIMITIVE POLYNOMIALS
 *    ------------------------------------------------
 */

void primitive_polynomial_test(int p, int d, int n) {
  int k;
  polynomial poly;
  double start, stop, elapsed;

  printf("\n");
  for (k = 0; k < n; k++) {
    start = clock();
    poly = gen_rand_primitive_poly(p, d);
    stop = clock();
    elapsed = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("%8.2g(secs):    ", elapsed);
    poly_print(poly);
    printf("\n");
    poly_free(poly);
  }
}

/*
 *
 *    GENERATE RANDOM PRIMITIVE POLYNOMIAL
 *    ------------------------------------
 */

polynomial gen_rand_primitive_poly(int p, int d) {
  /* Search for a primitive polynomial of a user specified
   degree, d, mod p. */
  polynomial poly;

  while (TRUE) {
    /* Generate random polynomial mod p of degree d. */
    poly = gen_rand_poly(p,d);
    /* Check polynomial is primitive. */
    gf_power_table(p,poly);
    if (primitive_polynomial_check == TRUE)
      return poly;
     else
      poly_free(poly);
  }
}

/*
 *
 *    GENERATE RANDOM POLYNOMIAL
 *    --------------------------
 */
polynomial gen_rand_poly(int p, int d) {
  polynomial poly;
  int n;

  if (p < 0 || d < 0) {
      printf("\nERROR: malformed polynomial specification.\n");
      exit(1);
    }

  poly_deg(poly) = d;
  poly_list(poly) = (int*) malloc((1 + d)*sizeof(int));
  poly_list(poly)[d] = 1;

  for (n = 0; n < d; n++) 
    poly_list(poly)[n] = random_in_range((unsigned int) 0, (unsigned int) p);

  while(poly_list(poly)[0] == 0)
    poly_list(poly)[0] = random_in_range((unsigned int) 0, (unsigned int) p);
  
  return poly;
}

/*
 *
 *    GENERATE RANDOM BIVARIATE POLYNOMIAL
 *    ------------------------------------
 */
bivariate_polynomial gen_rand_poly_2D(int p, int polydegX, int polydegY) {
  bivariate_polynomial poly;
  int x, y, **coeffs;

  /* Set degree in x and y. */
  poly_degX(poly) = polydegX;
  poly_degY(poly) = polydegY;

  /* Allocate memory for array of coefficients. */
  coeffs = (int**) malloc((1 + polydegX)*sizeof(int*));
  for (x = 0; x <= polydegX; x++) 
    coeffs[x] = (int*) malloc((1 + polydegY)*sizeof(int));

  for (x = 0; x <= polydegX; x++) {
    for (y = 0; y <= polydegY; y++) {
      coeffs[x][y] = random_in_range(0, p - 1);
    }
  }
  
  poly_array(poly) = coeffs;
  return poly;
}



/*
 *
 *    GENERATE RANDOM SPECIAL BIVARIATE POLYNOMIAL
 *    --------------------------------------------
 */
/* Generates a bivariate polynomial of the form 
       c_1 + c_2 y + c_3 y^2 + c_4 x + c_5 x y + c_6 x^2, 
   where the coefficients are random integers mod p. */

bivariate_polynomial  gen_special_rand_poly_2D(int p) { 
  bivariate_polynomial poly;
  int x, c3, c5, c6, **coeffs, polydegX = 2, polydegY = 2;

  /* Set degree in x and y. */
  poly_degX(poly) = polydegX;
  poly_degY(poly) = polydegY;

  /* Allocate memory for array of coefficients. */
  coeffs = (int**) malloc((1 + polydegX)*sizeof(int*));
  MEMCHECK(coeffs);
  for (x = 0; x <= polydegX; x++) { 
    coeffs[x] = (int*) malloc((1 + polydegY)*sizeof(int));
    MEMCHECK(coeffs[x]);
  }

  /* The coefficient array is of the form
   
      c_1  c_2  c_3
      c_4  c_5  0
      c_6  0    0
*/

  coeffs[1][2] = 0;
  coeffs[2][1] = 0;
  coeffs[2][2] = 0;
  coeffs[0][0] = random_in_range(0, p - 1);
  coeffs[0][1] = random_in_range(0, p - 1);
  coeffs[1][0] = random_in_range(0, p - 1);

  do {
    coeffs[0][2] = c3 = random_in_range(0, p - 1);
    coeffs[1][1] = c5 = random_in_range(0, p - 1);
    coeffs[2][0] = c6 = random_in_range(0, p - 1);
#if 0
    printf("\nc3 = %i, c_5 = %i, c6 = %i", c3, c5, c6);
#endif
  } while ( system_solver_2D(p, c3, c5, c6) );

  poly_array(poly) = coeffs;
  return poly;
}

/* 
 The routine system_solver_2D( ... ) solves the system

    2 c_6 s_1 + c_5 s_2 = 0 \mod p

    2 c_3 s_2 + c_5 s_1 = 0 \mod p

    for s_1, s_2 \in Z_p, where c_3, c_5, c_6 \in Z_p are known. If a 
 solution exists, then system_solver_2D() returns TRUE, otherwise FALSE. 

 The algorithm used is not, in any way, clever. It should be 
 sufficient providing p is not too big. */

bool system_solver_2D(int p, int c3, int c5, int c6) {
  int s1, s2;

  for (s1 = 0; s1 < p; s1++) {
    for (s2 = 0; s2 < p; s2++) {
      if (s1 == 0 && s2 == 0) 
	continue;
      if ((2*c6*s1 + c5*s2) % p == 0 && (2*c3*s2 + c5*s1) % p == 0) {
#if 0
	printf("\nSystems solution: s1 = %i, s2 = %i\n", s1, s2);
#endif
	return TRUE;
      }
    }
  }

#if 0
  printf("\nSystem has no solution.");
#endif
  return FALSE;
}

/* 
 *
 *    RANDOM INTEGER WITHIN A RANGE
 *    -----------------------------
 */

/* Source: 
 http://stackoverflow.com/questions/2509679/how-to-generate-a-random-number-from-within-a-range */

/* Would like a semi-open interval [min, max) */
int random_in_range (unsigned int min, unsigned int max)
{
  int n = 0;
  int base_random, range, remainder, bucket;

  while(++n < 16) rand(); /* hack! STB */
  base_random = rand(); /* in [0, RAND_MAX] */
  if (RAND_MAX == base_random) return random_in_range(min, max);
  /* now guaranteed to be in [0, RAND_MAX) */
  range     = max - min;
  remainder = RAND_MAX % range;
  bucket    = RAND_MAX / range;
  /* There are range buckets, plus one smaller interval
     within remainder of RAND_MAX */
  if (base_random < RAND_MAX - remainder)
    return min + base_random/bucket;
  else
    return random_in_range(min, max);
}

/* 
 *
 *    ARITHMETIC IN GALOIS FIELD
 *    --------------------------
 */

void generate_power_table(int prime, polynomial poly) {
  double start, stop, elapsed;
  /* Initialise power table lookup. */
  if (! power_table_init) {
#if DEBUG
    printf("\nGENERATING POWER TABLE LOOKUP...\n");
    start = clock();
#endif
    gf_power_table(prime, poly);
    /* Check user input polynomial is primitive. */
    if (primitive_polynomial_check == FALSE)
      {
	printf("\nERROR: input polynomial is not primitive.\n");
	exit(1);
      }
#if DEBUG
    stop = clock();
    elapsed = ((double) (stop - start)) / CLOCKS_PER_SEC;
    printf("\nFINISHED GENERATING POWER TABLE: %6.3f (seconds)\n\n", elapsed);
#endif
    }
} 

/* Display the powers of \alpha */
void print_power_table(int p, polynomial primpoly) {
  int n, k, maxdigs, ndigs, nspaces, plen;
  polynomial poly;
  
  /* Initialise power table lookup. */
  generate_power_table(p, primpoly);

  printf(  "|----------|-----------------------------------------|");
  printf("\n| alpha^n  |       polynomial representation         |\n");
  printf(  "|----------|-----------------------------------------|\n");

  maxdigs = intlength(power_table_length);
  poly_deg(poly) = poly_deg(primpoly) - 1;

  for (n = 0; n < power_table_length; n++) {
    ndigs = intlength(n);
    nspaces = 9 - ndigs; 
    printf("|");
    for (k = 0; k < nspaces; k++) 
      printf(" ");
    printf("%i | ", n);
    poly_list(poly) = power_table[n];
    plen = poly_print(poly);
    nspaces = 40 - plen;
    for (k = 0; k < nspaces; k++)
      printf(" ");
    printf("|\n");
  }

  printf(  "|----------|-----------------------------------------|\n");
}

/* Table of powers of \alpha */

void gf_power_table(int prime, polynomial primpoly) {
  int fieldsize, deg, n, m, *alpha_poly, *p_poly;
  gf p0, p1, alpha;

  /* Free old power table. */
  if (power_table_init)
    free_power_table();

  power_table_init = TRUE; /* Global. */
  primitive_polynomial_check = TRUE; /* Global. */

  deg = poly_deg(primpoly);

  if (deg == 1) 
    {
      printf("\nLinear primitive polynomials are not supported.\n");
      exit(1);
    }

  fieldsize = power(prime,deg) - 1;
  power_table_length = fieldsize; /* Global. */
  power_table_degree = deg; /* Global. */
  power_table = (int**) malloc(fieldsize*sizeof(int*));
  for (n = 0; n < fieldsize; n++) {
    power_table[n] = (int*) malloc(deg*sizeof(int));
    /* initialise to zero */
    for (m = 0; m < deg; m++) power_table[n][m] = 0;
  }

  /* \alpha^0 */
  power_table[0][0] = 1;

  /* \alpha^1 */
  /* Create polynomial */
  alpha_poly = (int*) malloc(deg*sizeof(int));
  for(m = 0; m < deg; m++) 
    alpha_poly[m] = 0;
  alpha_poly[1] = 1;
  
  /* Create gf expr for \alpha */
  alpha = gf_create(prime, deg, poly_list(primpoly), alpha_poly);

  /* Set p to 1  */
  p_poly = (int*) malloc(deg*sizeof(int));
  for(m = 0; m < deg; m++) 
    p_poly[m] = power_table[0][m];
  p0 = gf_create(prime, deg, poly_list(primpoly), p_poly);

  for (n = 1; n < fieldsize; n++) {
    /* Multiply p by alpha. */
    p1 = gf_times(alpha, p0);

#if 0
    printf("\nalpha^%i = ", n);
    gf_print(p1);
#endif

    /* Is the power of alpha 1? If so, the polynomial is not primitive. */
    if (gf_equal_to_1(p1)) 
      {
	primitive_polynomial_check = FALSE;
	gf_delete(p0);
	gf_delete(p1);
	gf_delete(alpha);
	free_power_table();
	power_table_init = FALSE;
	return ;
      }

    /* Copy polynomial form of power of \alpha into power_table. */
    for (m = 0; m < deg; m++)
      power_table[n][m] = gf_poly_list(p1)[m];

    /* Clean up memory. */
    gf_delete(p0);
    p0 = gf_create(gf_prime(p1), gf_deg(p1), gf_prim_poly_list(p1), gf_poly_list(p1));
    gf_delete(p1);
  }

  gf_delete(p0);
  gf_delete(alpha);
}

void free_power_table() {
  int n;
  if (! power_table_init) return ;
  power_table_init = FALSE;
  for(n = 0; n < power_table_length; n++) free(power_table[n]);
  free(power_table);
}

/* Compute power of \alpha */
gf gf_power(gf a, int n) {
  int i, m, exponent, *zeros;
  gf a_pow_n;

  /* Initialise power table lookup. */
  generate_power_table(gf_prime(a), gf_prim_poly(a));

  /* a^-oo == 0 */
  if (n == NEGINF) {
    zeros = (int*) malloc(gf_deg(a)*sizeof(int));
    for (i = 0; i < gf_deg(a); i++)
      zeros[i] = 0;
    a_pow_n = gf_create(gf_prime(a), gf_deg(a), gf_prim_poly_list(a), zeros);
    return a_pow_n;
  }

  /* Lookup 'a's representation as a power of alpha, say m. Then
   a^n = (alpha^m)^n = alpha^(m * n mod p^d - 1). */
  m = power_representation(gf_poly(a));

  /* m == -oo  =>  a == 0 => a^n == 0  */
  if (m == NEGINF)
    return gf_copy(a);

  if (m < 0 || n < 0) 
    {
      (n < 0) ? INTERNAL_ERROR(737) : INTERNAL_ERROR(693);
      exit(1);
    }

  exponent = (m*n) % power_table_length;

  a_pow_n = gf_create(gf_prime(a), gf_deg(a), gf_prim_poly_list(a), power_table[exponent]);

  return a_pow_n;
}

int power_representation(polynomial l) {
  /* Find the power of alpha representation of the polynomial. */
  int i, len = poly_deg(l);

  if (! power_table_init) 
    {
      INTERNAL_ERROR(335);
      exit(1);
    }

  /* zero? */
  if (zeroq(l))
    return NEGINF;

  for (i = 0; i < power_table_length; i++) {
    if (coeffs_equal(poly_list(l), power_table[i], len + 1)) {
      return i;
    }
  }

  /* poly_print(l); */
  INTERNAL_ERROR(375);
  exit(1);
}

bool coeffs_equal(int *c1, int *c2, int len) {
  int n;

  for (n = 0; n < len; n++) {
    if (c1[n] != c2[n]) 
      return FALSE;
  }

  return TRUE;
}

bool neginfq(polynomial p) {
  if (poly_deg(p) == 0 && poly_list(p)[0] == NEGINF)
    return TRUE;
  else
    return FALSE;
}

bool zeroq(polynomial l) {
  int n, len = poly_deg(l);
  for (n = 0; n <= len; n++) {
    if (poly_list(l)[n] != 0)
      return FALSE;
  }
  return TRUE;
}

/* multiplication */

gf gf_times(gf a, gf b) {
  gf ab;
  polynomial abp, abmodpx;

  if (! gf_field_equal(a, b) )
    {
      INTERNAL_ERROR(397);
      exit(1);
    }

  /* 
  poly_deg(ap) = gf_poly_deg(a);
  poly_list(ap) = gf_poly_list(a);

  poly_deg(bp) = gf_poly_deg(b);
  poly_list(bp) = gf_poly_list(b); */

  /* r = a*b mod p(x) */
  abp = poly_mult(gf_poly(a), gf_poly(b));
  abmodpx = poly_rem(abp, gf_prim_poly(a), gf_prime(a));
  ab = gf_create(gf_prime(a), gf_deg(a), gf_prim_poly_list(a), poly_list(abmodpx));

  /* Free memory. */
  poly_free(abp);
  poly_free(abmodpx);
  return ab;
}


/* Integer polynomial multiplication.  */
polynomial poly_mult(polynomial a, polynomial b) {
  polynomial ab;
  int n, m;
  
  /* Multiply polynomials, a and b. */
  
  poly_deg(ab) = poly_deg(a) + poly_deg(b); 
  poly_list(ab) = (int*) malloc((1 + poly_deg(ab))*sizeof(int));

  /* Initialise coefficients to zero. */
  for (n = 0; n<= poly_deg(ab); n++) 
    poly_list(ab)[n] = 0;

  /* For each coefficient we compute: a x^n * b x^m = a*b*x^(n + m) */
  for (n = 0; n <= poly_deg(a); n++) {
    for (m = 0; m <= poly_deg(b); m++)
      poly_list(ab)[n + m] += poly_list(a)[n]*poly_list(b)[m];
  }
  return ab;
}

/* Integer polynomial (long) division. */
polynomial poly_rem(polynomial u, polynomial v, int p) {
  /* Ref: TAOCP, vol. 2, pp. 402, 2nd Ed, D. Knuth, 1980.  */
  /* Computes the remainder of u(x)/v(x) mod p. */
  polynomial r; // remainder polynomial
  int q, m, n, j, k, num, den;

  m = poly_deg(u);
  n = poly_deg(v);

  /* Degree of u(x) is less than the degree of v(x). */
  if (m < n) 
    return poly_copy(u);

  /* u(x) and v(x) are integers */
  if (m == 0 && n == 0) 
    {
      poly_deg(r) = 0;
      poly_list(r) = (int*) malloc(sizeof(int));
      poly_list(r)[0] = poly_list(u)[0]%poly_list(v)[0];
      return r;
    }

  /* Division by zero. */
  if (n == 0 && poly_list(v)[0] == 0)
    {
      INTERNAL_ERROR(666);
      exit(1);
    }

  /* Create remainder polynomial. */
  poly_deg(r) = n - 1;
  poly_list(r) = (int*) malloc((poly_deg(u) + 1)*sizeof(int));
  for (k = 0; k <= poly_deg(u); k++) 
    poly_list(r)[k] = poly_list(u)[k]; 

  /* D1 */
  for (k = m - n; k >= 0; k--) {
    /* D2 */
    for (j = n + k - 1; j >= k; j--) {
          num = poly_list(r)[n + k]*poly_list(v)[j - k];
	  den = poly_list(v)[n];
	  if (num%den == 0)
	    q = (num/den)%p;
	  else 
	    {
	      q = mult_inv(den,p);
	      q = num*q%p;
	    }
	  poly_list(r)[j] -= q;
	}
    }

  /* Reduce coefficients mod p. */
  for (k = 0; k <= poly_deg(r); k++) {
    poly_list(r)[k] = poly_list(r)[k]%p;
    if (poly_list(r)[k] < 0) 
      poly_list(r)[k] = poly_list(r)[k] + p;
    }

  /* Remove (ignore) zero coefficients of higher-degree terms. */
  for (k = n; k <= poly_deg(r); k++) 
    poly_list(r)[k] = 0;

  /* Compute the minimal degree. */
  for (k = n; k >= 0; k--)
    if (poly_list(r)[k] != 0) break;
  poly_deg(r) = k;

  return r;
}

/* Multiplicative inverse mod p. */
/* Source: http://rosettacode.org/wiki/Modular_inverse#C */

int mult_inv(int a, int b)
{
	int b0 = b, t, q;
	int x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

/* addition */

gf gf_plus(gf a, gf b) {
  /* a + b */
  gf result;
  int n, *apb;

  if (! gf_field_equal(a, b) ) 
    {
      INTERNAL_ERROR(519);
      exit(1);
    }

  apb = (int *) malloc((1 + gf_poly_deg(a))*sizeof(int));

  /* add */
  for (n = 0; n <= gf_poly_deg(a); n++)
      apb[n] = (gf_poly_list(a)[n] + gf_poly_list(b)[n]) % a.p;

  result = gf_create(gf_prime(a), gf_deg(a), gf_prim_poly_list(a), apb);

  free(apb);
  return result;
}

/* 
 *
 *    NAIVE PRIMALITY CHECKING
 *
 */

/* Just for small numbers. */

bool primep(int n) {
    int i;

    if (n == 2 || n == 3) return TRUE; 

    if (n <= 1 || (n % 2 == 0)) return FALSE; 

    for (i = 3; i*i <= n; i += 2)
      if (n % i == 0) return FALSE;

    return TRUE;
}


/* 
 *
 *    CREATE GALOIS FIELD DATA STRUCTURE
 *    ----------------------------------
 */

gf gf_create(int p, int degree, int *prim_poly_coeffs, int *poly_coeffs) {

  polynomial prim_poly, poly;
  gf result;
  int n, *prim_poly_coeffs_copy, *poly_coeffs_copy;

  /* Create primitive polynomial as polynomial struct. */
  prim_poly_coeffs_copy = (int *) malloc((degree + 1)*sizeof(int));
  for (n = 0; n <= degree; n++) 
    prim_poly_coeffs_copy[n] = prim_poly_coeffs[n];
  poly_deg(prim_poly) = degree;
  poly_list(prim_poly) = prim_poly_coeffs_copy;

  /* Create GF data structure. */
  gf_prime(result) = p;
  gf_deg(result) = degree; // depreciated. 
  gf_prim_poly(result) = prim_poly;

  poly_coeffs_copy = (int *) malloc(degree*sizeof(int));
  for (n = 0; n < degree; n++) 
    poly_coeffs_copy[n] = poly_coeffs[n]; 

  poly_deg(poly) = degree - 1;
  poly_list(poly) = poly_coeffs_copy;

  gf_poly(result) = poly;
  
  malloc_count += 2; /* Global variable. */
  gf_alloc_count++;  /* Global variable. */

  return result;
}

/* 
 *
 *    COPY GALOIS FIELD DATA STRUCTURE
 *    --------------------------------
 */
gf gf_copy(gf a) {
  gf copy;
  copy = gf_create(gf_prime(a), gf_deg(a), gf_prim_poly_list(a), gf_poly_list(a));
  return copy;
}

/* 
 *
 *    DELETE GALOIS FIELD DATA STRUCTURE
 *    ----------------------------------
 */

void gf_delete(gf a) {
  malloc_count -= 2; /* Global variable. */
  if (malloc_count < 0) 
    {
      INTERNAL_ERROR(721);
      exit(1);
    }
  free(gf_poly_list(a));
  free(gf_prim_poly_list(a));
}

/* 
 *
 *    DISPLAY GALOIS FIELD (debugging)
 *    --------------------------------
 */

void gf_print(gf a) {
  int n, deg = gf_deg(a);
  printf("\nGF[%d, {", gf_prime(a));
  for (n = 0; n < deg; n++)
    printf("%d, ", gf_prim_poly_list(a)[n]);
  printf("%d}][{", gf_prim_poly_list(a)[deg]);
  for (n = 0; n < deg - 1; n++)
    printf("%d, ", gf_poly_list(a)[n]);
  printf("%d}]\n", gf_poly_list(a)[deg - 1]);
}

/* 
 *
 *    EQUALITY IN GALOIS FIELD
 *    ------------------------
 */

bool gf_field_equal(gf a, gf b) {
  int n;

  /* Check primes are equal.  */
  if (gf_prime(a) != gf_prime(b)) return FALSE;
 
  /* Check equal number of elements.  */
  if (gf_poly_deg(a) != gf_poly_deg(b)) return FALSE;

  /* Check polynomial degrees are equal. */
  if (gf_deg(a) != gf_deg(b)) return FALSE;

  /* Check primitive polynomials are equal. */
  for (n = 0; n < gf_deg(a); n++) 
      if (gf_prim_poly_list(a)[n] != gf_prim_poly_list(b)[n] ) return FALSE;
  return TRUE;
}

bool gf_equal(gf a, gf b) {
  int n, nelems;

  /* Check fields are equal. */
  if ( ! gf_field_equal(a, b) ) return FALSE;

  /* Check elements are equal.  */
  nelems = gf_poly_deg(a);
  for (n = 0; n < nelems; n++) 
      if (gf_poly_list(a)[n] != gf_poly_list(b)[n]) return FALSE;
  return TRUE;
}

bool gf_equal_to_1(gf a) {
  int n;
  if (gf_poly_list(a)[0] != 1) return FALSE;
  for (n = 1; n < gf_deg(a); n++)
    if (gf_poly_list(a)[n] != 0) return FALSE;
  return TRUE;
}

/* 
 *
 *    COPY A POLYNOMIAL DATA STRUCTURE
 *    --------------------------------
 */
polynomial poly_copy(polynomial poly) {
  int n;
  polynomial copy;

  /* Set the polynomial degree. */
  poly_deg(copy) = poly_deg(poly); 
  /* Assign memory for the polynomial coefficient list. */
  poly_list(copy) = (int*) malloc((1 + poly_deg(copy))*sizeof(int)); 
  MEMCHECK(poly_list(copy));
  /* Copy the coefficient list term-by-term. */
  for (n = 0; n <= poly_deg(poly); n++) 
    poly_list(copy)[n] = poly_list(poly)[n];
  return copy;
}

/*
 *
 *    PRETTY DISPLAY OF A BIVARIATE POLYNOMIAL
 *    ----------------------------------------
 */

void poly_print_2D(bivariate_polynomial poly) {
  int n, m, degX = poly_degX(poly), degY = poly_degY(poly), coeff;
  bool first = TRUE;

  for (n = 0; n <= degX; n++) {
    for (m = 0; m <= degY; m++) {

      /* Print coefficient. */
      coeff = poly_array(poly)[n][m];

      if (coeff != 0) {
	if (coeff < 0) {
	  if (first) 
	    printf("- %i", coeff);
	  else
	    printf(" - %i", coeff);
	} else {
	  if (!first) 
	    printf(" +");
	  if (coeff != 1 || (n == 0 && m == 0))
	    printf(" %i", coeff);
	}

	/* Print power of x. */
	if (n != 0) {
	  if (n == 1) {
	    printf(" x"); 
	  } else {
	    printf(" x^%i", n);
	  }
	}

	/* Print power of y. */
	if (m != 0) {
	  if (m == 1) {
	    printf(" y");
	  } else {
	  printf(" y^%i", m);
	  }
	}
	
	first = FALSE;
      }
    }
  }
}

/*
 *
 *    PRETTY DISPLAY OF A POLYNOMIAL
 *    ------------------------------
 */
int poly_print(polynomial poly) {
  int d = poly_deg(poly), n, c, nchars = 0;
  bool allzero = TRUE;

  if (zeroq(poly))
    {
      printf("0");
      return 1;
    }

  for(n = 0; n <= d; n++) {
    c = poly_list(poly)[n];
    if (n > 0 && poly_list(poly)[n-1] != 0)
      allzero = FALSE;
    if (c == 0) continue;
    if (n == 0)
      {
	if (c == NEGINF) {
	  printf("-oo");
	  nchars += 3;
	} else {
	  printf("%i", c);
	  nchars += intlength(c);
	}
      }
    else if (n == 1)
      {
	if (c < 0) {
	    if (c == NEGINF) {
	      printf("-oo x");
	      nchars += 5;
	    } else {
	      printf(" %i x", c);
	      nchars += 3 + intlength(c);
	    }
	}
	else if (c == 1)
	  {
	      if (allzero == TRUE) {
		printf("x");
		nchars++;
	      } else {
		printf(" + x");
		nchars += 4;
	      }
	  }
	else
	  {
	    if (allzero == TRUE) {
	      if (c == NEGINF) {
		printf(" - oo x");
		nchars += 7;
	      } else {
		printf("%i x", c);
		nchars += 2 + intlength(c);
	      }
	    } else { 
	      if (c == NEGINF) {
	        printf(" - oo x");
		nchars += 7;
	      } else {
		printf(" + %i x", c);
		nchars += 5 + intlength(c);
	      }
	    }
	  }
      }
    else
      {
	if (c < 0) {
	    if (c == NEGINF) {
	      printf(" - oo x^%i", n);
	      nchars += 8 + intlength(n);
	    } else {
	      printf(" %i x^%i", c, n);
	      nchars += 4 + intlength(c) + intlength(n);
	    }
	}
	else if (c == 1)
	  {
	    if (allzero == TRUE) {
	      printf("x^%i", n);
	      nchars += 2 + intlength(n);
	    } else {
	      printf(" + x^%i", n);
	      nchars += 5 + intlength(n);
	    }
	  }
	else
	  {
	    if (allzero == TRUE) {
	      if (c == NEGINF) {
		printf("- oo x^%i", n);
		nchars += 7 + intlength(n);
	      } else {
		printf("%i x^%i", c, n);
		nchars += 3 + intlength(c) + intlength(n);
	      }
	    } else {
	      if (c == NEGINF) {
		printf(" - oo x^%i", n);
		nchars += 8 + intlength(n);
	      } else {
		printf(" + %i x^%i", c, n);
		nchars += 6 + intlength(c) + intlength(n);
	      }
	    }
	  }
      }
  }
  return nchars;
}

/*
 *
 *	PERMUTATION RELATED CODE
 *      ------------------------       
 */

void NextPermutation(int *array, int length) {
  /* Knuth, combinatorial algorithms, Algorithm L, page 2, section 7.2.1.2 */	
  int j, l, k;

  /* L2 -- find j */
  j = length - 2;
  while(array[j] >= array[j + 1]){
    if(j == 0) return;
    j--;
  }

  /* L3 -- increase a[j] */
  l = length - 1;
  while(array[j] >= array[l])
    l--;

  /* swap a[j] and a[l] */
  swap(&array[l], &array[j]);

  /* L4 -- reverse a[j] to a[length - 1] */
  k = j + 1;
  l = length - 1;
  while(k < l)
    swap(&array[k++], &array[l--]);
}

void swap(int *i, int *j) {
    int t;
    t = *i;
    *i = *j;
    *j = t;
}

void reverse(int arr[], int start, int end) {
  int temp;
  while(start < end)
  {
    temp = arr[start];   
    arr[start] = arr[end];
    arr[end] = temp;
    start++;
    end--;
  }   
}

void increment(int *counter, int len, int b) {
  int n, i;
  
  /* Increment a counter by 1 in base b. */
  for (n = len - 1; n >= 0; n--){
    if (counter[n] == b - 1){
      continue ;
    } else {
      counter[n]++;
      /* Fill trailing positions with zeros. */
      for (i = n + 1; i < len; i++)
	counter[i] = 0;
      break ;
    }
  }
}


/* 
 *
 *    BINARIZE 1D ARRAY
 *    -----------------
 */
/* Basically just does: int -> int mod 2 followed by: 1 -> -1, 0 -> 1, -oo -> 0 */
void binarize1D(int *array, int d0) {
  int n;
  for (n = 0; n < d0; n++) {
    if (array[n] == NEGINF)
      array[n] = 0;
    else
      {
	array[n] = array[n]%2;
	if (array[n] == 0)
	  array[n] = 1;
	else
	  array[n] = -1;
      }
  }
}

/* 
 *
 *    BINARIZE 2D ARRAY
 *    -----------------
 */
void binarize2D(int **array, int d0, int d1) {
  int n, m;
  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      if (array[n][m] == NEGINF)
	array[n][m] = 0;
      else
	{
	  array[n][m] = array[n][m]%2;
	  if (array[n][m] == 0)
	    array[n][m] = 1;
	  else
	    array[n][m] = -1;
	}
    }
  }
}


/* 
 *
 *    BINARIZE 3D ARRAY
 *    -----------------
 */
void binarize3D(int ***array, int d0, int d1, int d2) {
  int n, m, k;

  if (already_near_binary_3D(array, d0, d1, d2)) 
    return;

  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      for (k = 0; k < d2; k++) {
	if (array[n][m][k] == NEGINF)
	  array[n][m][k] = 0;
	else
	  {
	    array[n][m][k] = array[n][m][k]%2;
	    if (array[n][m][k] == 0)
	      array[n][m][k] = 1;
	    else
	      array[n][m][k] = -1;
	}
      }
    }
  }
}

bool already_near_binary_3D(int ***array, int d0, int d1, int d2) {
  int n, m, k;

  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      for (k = 0; k < d2; k++) {
	if (!(array[n][m][k] == 1 || array[n][m][k] == -1 || array[n][m][k] == 0))
	  return FALSE;
      }
    }
  }

  return TRUE;
}

/* 
 *
 *    DISPLAY 1D ARRAY
 *    ----------------
 */

#define print(k) if (k == NEGINF)\
      printf("-oo");\
    else\
      printf("%i", k)

void print_array_1D(int *seq, int d0) {
  int n;
  for (n = 0; n < d0; n++) {
    print(seq[n]);
    if (n != d0 - 1) printf(" ");
  }
}

/* 
 *
 *    DISPLAY 2D ARRAY
 *    ----------------
 */

void print_array_2D(int **array, int d0, int d1) {
  int n, m;
  
  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      print(array[n][m]);
      if (m != d1 - 1) printf(" ");
    }
    printf("\n");
  }
}

/* Computes the number of decimal digits in an int. */
int intlength(int k) {
  int m;
  m = (k > 0) ? k : -k;
  
  /* Inelegant, but faster than 1 + ((int) log10(((float) k))) */
  if (m < 10) return 1;
  if (m < 100) return 2;
  if (m < 1000) return 3;
  if (m < 10000) return 4;
  if (m < 100000) return 5;
  if (m < 1000000) return 6;
  if (m < 10000000) return 7;
  if (m < 100000000) return 8;
  if (m < 1000000000) return 9;
  /* The next number would be greater than 2^32 - 1 == 4294967295. */
  return ((int) 1.0 + floor(log10(m)));
}

/* This routine is significantly slower than print_array_2D. It is useful for 
displaying arrays of autocorrelations as everything is aligned. */
void pretty_print_array_2D(int **array, int d0, int d1) {
  int v, n, m, k, max = 0, ndigs, maxdigs;

  /* Compute maximum value (in magnitude) in the array. */
  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      v = array[n][m];
      if (v == NEGINF && max < 100) {
	max = 100; // three characters in -oo
	continue;
      }
      if (v > 0) {
	if (v > max) {
	  max = v;
	}
      } else {
	if (-10*v > max) {
	  max = -10*v; // add one digit for the minus sign.
	}
      }
    }
  }  
  
  maxdigs = intlength(max);

  // printf("\nmaxdigs = %i\n", maxdigs);

  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      v = array[n][m];
      if (v == NEGINF) { 
	for (k = 0; k < maxdigs - 3; k++) printf(" ");
	printf("-oo");
      } else {
	ndigs = intlength(v);
	if (v < 0) ndigs++;
	for (k = 0; k < maxdigs - ndigs; k++) printf(" ");
	printf("%i", v);
      }
      /* One space between numbers. */
      if (m != d1 - 1) printf(" ");
    }
    printf("\n");
  }
}

/* 
 *
 *    DISPLAY 3D ARRAY
 *    ----------------
 */

void pretty_print_array_3D(int ***array, int d0, int d1, int d2) {
  int v, n, m, k, l, max = 0, ndigs, maxdigs;

  /* Compute maximum value (in magnitude) in the array. */
  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      for (k = 0; k < d2; k++) {
	v = array[n][m][k];
	if (v == NEGINF && max < 100) {
	  max = 100; // three characters in -oo
	} else if (v > 0) {
	  if (v > max) {
	    max = v;
	  }
	} else {
	  if (-10*v > max) {
	    max = -10*v; // add one digit for the minus sign.
	  }
	}
      }
    }
  }
  
  maxdigs = intlength(max);
  
  // printf("\nmaxdigs = %i\n", maxdigs);

  for (k = 0; k < d2; k++) {
    for (n = 0; n < d0; n++) {
      for (m = 0; m < d1; m++) {
	v = array[n][m][k];
	if (v == NEGINF) {
	  for (l = 0; l < maxdigs - 3; l++) printf(" ");
	  printf("-oo");
	} else {
	  ndigs = intlength(v);
	  if (v < 0) ndigs++;
	  for (l = 0; l < maxdigs - ndigs; l++) printf(" ");
	  printf("%i", v);
	}
	/* One space between numbers. */
	if (m != d1 - 1) printf(" ");
      }
      printf("\n");
    }
    printf("\n");
  }
}

void print_array_3D(int ***array, int d0, int d1, int d2) {
  int n, m, k;
  
  for (k = 0; k < d2; k++) {
    for (n = 0; n < d0; n++) {
      for (m = 0; m < d1; m++) {
        print(array[n][m][k]);
        if (m != d1 - 1) printf(" ");
      }
      printf("\n");
    }
    printf("\n");
  }
}

/* 
 *
 *    WRITE 1D ARRAY TO BINARY FILE
 *    -----------------------------
 */

void write_binary_array_1D(char *filename, int *array, int d0) {
  int n, val, intinf = INT_MAX;
  
  FILE *f = fopen(filename, "wb");
  if (f == NULL) 
    {
      printf("\nERROR: cannot create output file.\n");
      exit(1);
    }

  for (n = 0; n < d0; n++) {
    val = array[n];
    if (val == NEGINF)
      fwrite(&intinf, sizeof(int), 1, f);
    else
      fwrite(&val, sizeof(int), 1, f);
  }
  fclose(f);
}

/* 
 *
 *    WRITE 2D ARRAY TO BINARY FILE
 *    -----------------------------
 */

void write_binary_array_2D(char *filename, int **array, int d0, int d1) {
  int n, m, val, intinf = INT_MAX;
  
  FILE *f = fopen(filename, "wb");
  if (f == NULL) 
    {
      printf("\nERROR: cannot create output file.\n");
      exit(1);
    }

  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
      val = array[n][m];
      if (val == NEGINF)
        fwrite(&intinf, sizeof(int), 1, f);
      else
	fwrite(&val, sizeof(int), 1, f);
    }
  }
  fclose(f);
}

/* 
 *
 *    WRITE 3D ARRAY TO BINARY FILE
 *    -----------------------------
 */

void write_binary_array_3D(char *filename, int ***array, int d0, int d1, int d2) {
  int n, m, k, val, intinf = INT_MAX;
  
  FILE *f = fopen(filename, "wb");
  if (f == NULL) 
    {
      printf("\nERROR: cannot create output file.\n");
      exit(1);
    }

  for (k = 0; k < d2; k++) {
    for (n = 0; n < d0; n++) {
      for (m = 0; m < d1; m++) {
	val = array[n][m][k];
	if (val == NEGINF)
	  fwrite(&intinf, sizeof(int), 1, f);
	else
	  fwrite(&val, sizeof(int), 1, f);
      }
    }
  }
  fclose(f);
}

/* 
 *
 *    WRITE 1D ARRAY
 *    --------------
 */

void write_array_1D(char *filename, bool csvp, int *array, int d0) {
  int n;
  
  FILE *f = fopen(filename, "w");
  if (f == NULL) 
    {
      printf("\nERROR: cannot create output file.\n");
      exit(1);
    }

  for (n = 0; n < d0; n++) {
    if (array[n] == NEGINF)
      fprintf(f,"-oo");
    else
      fprintf(f,"%i", array[n]);
    if (n != d0 - 1) {
      if (csvp)
	fprintf(f,",");
      else
	fprintf(f," ");
    }
  }
  fprintf(f,"\n");
  fclose(f);
}

/* 
 *
 *    WRITE 2D ARRAY
 *    --------------
 */

void write_array_2D(char *filename, bool csvp, int **array, int d0, int d1) {
  int n, m;
  
  FILE *f = fopen(filename, "w");
  if (f == NULL) 
    {
      printf("\nERROR: cannot create output file.\n");
      exit(1);
    }

  for (n = 0; n < d0; n++) {
    for (m = 0; m < d1; m++) {
	if (array[n][m] == NEGINF)
	  fprintf(f,"-oo");
	else
	  fprintf(f,"%i", array[n][m]);
	if (m != d1 - 1) {
	  if (csvp)
	    fprintf(f,",");
	  else
	    fprintf(f," ");
	}
    }
    fprintf(f, "\n");
  }
  fprintf(f,"\n");
  fclose(f);
}

/* 
 *
 *    WRITE 3D ARRAY
 *    --------------
 */

void write_array_3D(char *filename, bool csvp, int ***array, int d0, int d1, int d2) {
  int n, m, k;
  
  FILE *f = fopen(filename, "w");

  if (f == NULL) 
    {
      printf("\nERROR: cannot create output file.\n");
      exit(1);
    }

  for (k = 0; k < d2; k++) {
    for (n = 0; n < d0; n++) {
      for (m = 0; m < d1; m++) {
	if (array[n][m][k] == NEGINF)
	  fprintf(f,"-oo");
	else
	  fprintf(f,"%i", array[n][m][k]);
	if (m != d1 - 1) {
	  if (csvp)
	    fprintf(f,",");
	  else
	    fprintf(f," ");
	}
      }
      fprintf(f, "\n");
    }
  }
  fclose(f);
}

/* 
 *
 *    CORRELATIONS
 *    ------------
 */

/* BEWARE: All the correlations code in this program assumes the arrays are 
integers, NOT roots of unity in exponent form. Thus, there is no conjugation 
and no complex arithmetic. */

/* -- 1D PERIODIC AUTOCORRELATION  -- */

int* acv_1d_slow(int *array, int d0) {

  int *corr, n, tau;

  /* Create a 1D array to store the autocorrelations. */
  corr = (int*) malloc(d0*sizeof(int));

  /* Compute the 1D autocorrelations. */
  for(tau = 0; tau < d0; tau++) {
    corr[tau] = 0;
    for (n = 0; n < d0; n++) {
      corr[tau] += array[n]*array[(n + tau) % d0];
    }
  }

  return corr;
}

/* -- 2D PERIODIC AUTOCORRELATION  -- */

int** acv_2d_slow(int **array, int d0, int d1) {

  int **corr, n, m, s0, s1;

  /* Create a 2D array to store the autocorrelations. */
  corr = (int**) malloc(d0*sizeof(int*));
  for (n = 0; n < d0; n++) 
    corr[n] = (int*) malloc(d1*sizeof(int));

  /* Compute the 2D autocorrelations. */
  for (s0 = 0; s0 < d0; s0++) {
    for (s1 = 0; s1 < d1; s1++) {
      corr[s0][s1] = 0;
      for (n = 0; n < d0; n++) {
	for (m = 0; m < d1; m++) {
	  corr[s0][s1] += array[n][m]*array[(n + s0) % d0][(m + s1) % d1];
	}
      }
    }
  }

  return corr;
}

/* -- 3D PERIODIC AUTOCORRELATION  -- */

int*** acv_3d_slow(int ***array, int d0, int d1, int d2) {

  int ***corr, n, m, k, s0, s1, s2;

  /* Create a 3D array to store the autocorrelations. */
  corr = (int***) malloc(d0*sizeof(int**));
  for (n = 0; n < d0; n++) {
    corr[n] = (int**) malloc(d1*sizeof(int*));
    for (m = 0; m < d1; m++)
      corr[n][m] = (int*) malloc(d2*sizeof(int*));
  }

  /* Compute the 3D autocorrelations. */
  for (s0 = 0; s0 < d0; s0++) {
    for (s1 = 0; s1 < d1; s1++) {
      for (s2 = 0; s2 < d2; s2++) {
	corr[s0][s1][s2] = 0;
	for (n = 0; n < d0; n++) {
	  for (m = 0; m < d1; m++) {
	    for (k = 0; k < d2; k++) {
	      corr[s0][s1][s2] += array[n][m][k]*array[(n + s0) % d0][(m + s1) % d1][(k + s2) % d2];
	    }
	  }
	}
      }
    }
  }

  return corr;
}

/* 
 *
 *    RECURSIVE INTEGER EXPONENTIATION
 *    --------------------------------
 */

int power(int base, int exp) {

  if (base < 0 || exp < 0)
    INTERNAL_ERROR(409);

  if (exp == 0)
    return 1;
  else if (exp % 2)
    return base * power(base, exp - 1);
  else {
    int temp = power(base, exp / 2);
    return temp * temp;
  }
}

/* 
 *
 *    1D LINEAR COMPLEXITY - BERLENKAMP MASSEY ALGORITHM
 *    --------------------------------------------------
 */
int* berlenkamp_massey(int *seq, int N) {
  int *b, *c, *t, d, l, m, n, nm, i, j;

  b = (int*) malloc(N*sizeof(int));
  c = (int*) malloc(N*sizeof(int));
  t = (int*) malloc(N*sizeof(int));

  b[0] = 1;
  for (n = 1; n < N; n++) b[n] = 0;
  c[0] = 1;
  for (n = 1; n < N; n++) c[n] = 0;
  l = 0; 
  m = -1;

  for (n = 0; n < N; n++) {
    d = 0;
    for (i = 0; i <= l; i++)
      d ^= c[i]*seq[n - i];

    if (d == 1) {
      copy_vec(c, 0, t, 0, N);
      nm = n - m;
      for (j = 0; j < N - nm; j++)
	c[nm + j] ^= b[j];
      if (l <= n/2) {
	l = n + 1 - l;
	m = n;
	copy_vec(t, 0, b, 0, N);
      }
    }
    //    print_vector(c, N);
  }

  return c;
}

/* lc - Returns the degree of the connection polynomial from the
   Berlenkamp Massey algorithm. */
int lc(int *seq, int n) {
  int nzeros = 0, k; 

  for (k = n - 1; k >= 0; k--) {
    if (seq[k] == 0)
      nzeros++;
    else
      break ;
  }

  return n - nzeros - 1;
}


/* 
 *
 *    1D ARRAY MEMORY COPY
 *    --------------------
 */

void copy_vec(int *src, int srcpos, int *dest, int destpos, int len) {
  int n;
  for (n = 0; n < len; n++)
    dest[destpos + n] = src[srcpos + n];
}

/* 
 *
 *    PASSWORD PROTECTION AND EXPIRY DATES
 *    ------------------------------------
 */


#define EXPIRED  printf("\nThis version of the array generator has expired. \
Please contact Sam Blake or Andrew Tirkel. \n\n")

bool checkexpirydate() {
  time_t rawtime;
  struct tm *info;
  int year, mon, day;

  time( &rawtime );
  info = localtime( &rawtime );
  year = info->tm_year + 1900;
  mon = info->tm_mon + 1;
  day = info->tm_mday;
  if (year < EXPIRY_YEAR) 
    {
      return FALSE;
    }
  else if (year > EXPIRY_YEAR)
    {
      EXPIRED; 
      return TRUE;
    }
  else if (mon < EXPIRY_MONTH)
    {
      return FALSE;
    }
  else if (mon > EXPIRY_MONTH)
    {
      EXPIRED;
      return TRUE;
    }
  else if (day <= EXPIRY_DAY)
    {
      return FALSE;
    }
  else 
    {
      EXPIRED;
      return TRUE;
    }
}

void readpass(char password[128]) {

  struct termios oflags, nflags;

  /* disabling echo */
  tcgetattr(fileno(stdin), &oflags);
  nflags = oflags;
  nflags.c_lflag &= ~ECHO;
  nflags.c_lflag |= ECHONL;

  if (tcsetattr(fileno(stdin), TCSANOW, &nflags) != 0) {
    perror("tcsetattr");
    return ;
  }

  printf("password: ");
  fgets(password, 128*sizeof(char), stdin);
  password[strlen(password) - 1] = 0;

  /* restore terminal */
  if (tcsetattr(fileno(stdin), TCSANOW, &oflags) != 0) {
    perror("tcsetattr");
  }
}

/* Converts a string to lower case. */
void lower(char *s) {
  char  *p;

  for (p = s; *p != '\0'; p++) 
    *p = (char) tolower(*p);
}
