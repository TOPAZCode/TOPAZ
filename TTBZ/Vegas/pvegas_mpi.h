#define ALPH 1.5             /* exponent to compress refinement in range     */
#define NDMX 100             /* maximum number of bins per dimension         */
#define MXDIM 25              /* maximum dimension that can be passed by ndim */  /* MARKUS: has to match the one in vegas_common.f !! */
#define FNMX 10              /* maximum number of integrands                 */

/* 
 * Different levels of verbosity.  They substitute the old unflexible nprn
 * integer.  Build up your own printlevel by bitwise or-ing several of these.
 */
#define NPRN_INPUT   0x0001  /* print input parameters                       */
#define NPRN_RESULT  0x0002  /* print results of primary integration         */
#define NPRN_SECRES  0x0004  /* print results of secondary integrations      */
#define NPRN_RESULTS 0x0006  /* print results of all integrations            */
#define NPRN_DEFAULT 0x0007  /* corresponds to old default output mode       */
#define NPRN_GRID_8  0x0010  /* print every eighth line of grid data output  */
#define NPRN_GRID_4  0x0020  /* print every fourth line of grid data output  */
#define NPRN_GRID_2  0x0040  /* print every second line of grid data output  */
#define NPRN_GRID    0x0080  /* print all grid data                          */
#define NPRN_ALL     0xffff  /* gossip monger mode: print everything above   */

/*
 * prototype for vegas
 */
void vegas_mpi_(double regn[], int *ndim_in, void (*fxn)(double x[], double *wgt, double f[]),
           int *init_in, int *ncall_in, int *itmx_in, int *nprn_in,
           int *fcns_in, int *pdim_in, int *wrks_in,
           double tgral[], double sd[], double chi2a[]);

/*
 * The two constants for the GFSR random number generator
 */
#define SR_P 1279
#define SR_Q 418

/*
 * Here is a list of other magic number pairs:
 *   SR_P        SR_Q
 *    127        1, 7, 15, 30, 63
 *    521        32, 48, 158, 168
 *    607        273, 105, 147
 *   1279        418, 216
 *   2281        715, 915, 1029
 *   3217        67, 576
 *   4423        1419, 271, 369, 370, 649, 1393, 2098
 *   9689        1836, 84, 471, 2444, 4187
 *  19937        881, 7083, 9842
 *  23209        1530, 6619, 9739
 *  44497        8575, 21034
 * 110503        25230, 53719
 * 132049        7000, 33912, 41469, 52549, 54454
 * Note that small values of SR_Q and values near SR_P/2 should
 * be avoided. Source: Hamilton, CPC 85 (1995), 127-152
 */
