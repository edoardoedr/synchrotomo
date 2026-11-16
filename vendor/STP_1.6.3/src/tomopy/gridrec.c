//#ifdef __cplusplus
//extern "C"{
//#endif

#define _VERSION_ "2013.04.29"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include "numpy/arrayobject.h"

#define _CRT_SECURE_NO_WARNINGS
#define ANSI

/**** System includes */
#include <math.h>
#include <stddef.h>
#include <time.h>
#include <sys/stat.h>

#ifdef _MSC_VER
typedef unsigned __int8  uint8_t;
typedef unsigned __int16  uint16_t;
typedef unsigned __int32  uint32_t;
typedef unsigned __int64  uint64_t;
#ifdef _WIN64
typedef __int64  ssize_t;
typedef signed __int64  intptr_t;
typedef unsigned __int64  uintptr_t;
#define SSIZE_MAX (9223372036854775808)
#else
typedef int ssize_t;
typedef _W64 signed int  intptr_t;
typedef _W64 unsigned int  uintptr_t;
#define SSIZE_MAX (2147483648)
#endif
#else
/* non MS compilers */
#include <stdint.h>
#include <limits.h>
#endif

/**** Macros and typedefs ****/
#define I(i,j,k,N,M)    ( (j)*(N) + (i) + (k)*(N)*(M) )
#define I2(i,j,N)       ( (j)*(N) + (i) )
#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))
#define free_matrix(A) (free(*(A)),free(A), A=NULL)
#define abs(A) ((A)>0 ?(A):-(A))
#define pi  3.14159265359
#define TOLERANCE 0.1	/* For comparing centers of two sinugrams */
#define LTBL_DEF 512	/* Default lookup table length */

typedef struct {float r;float i;} complex;

typedef struct SG_STRUCT { 		/* Sinugram parameters */

	int n_ang;	/* No. angles in sinugram */
	int n_det;	/* No. elems (detectors) per angle */
	int geom;		/* 0 if array of angles provided;
					1,2 if uniform in half,full circle */ 
	float *angles;	/* Ptr to the array of angles, if used */
	float center;	/* Rotation axis location */

} sg_struct;

typedef struct PSWF_STRUCT{ /*Prolate spheroidal wave fcn (PSWF) data */

	float C;	/* Parameter for particular 0th order pswf being used*/
	int nt; /* Degree of Legendre polynomial expansion */
	float lmbda; 	/* Eigenvalue */
	float coefs[15];	/* Coeffs for Legendre polynomial expansion */

} pswf_struct;

typedef struct GRID_STRUCT { /* Parameters for gridding algorithm */

	pswf_struct *pswf;	/* Ptr to data for PSWF being used  */
	float sampl;	  	/* "Oversampling" ratio */
	float MaxPixSiz; 	/* Max pixel size for reconstruction */
	float R;		/* Region of interest (ROI) relative size */
	float X0;		/* (X0,Y0)=Offset of ROI from rotation axis, */
	float Y0;		/* in units of center-to-edge distance.  */
	char fname[16];		/*  Name of filter function   */		
	float (*filter)(float);	/* Ptr to filter function.  */
	long ltbl;		/* No. elements in convolvent lookup tables. */

} grid_struct;

/** Global variables **/
int verbose = 1;

/**** Function Prototypes ****/

/** Defined in grid_io.c **/
void get_parm(int argc,char **argv, grid_struct *A);
int data_setup(int argc, char *argv[]);
float **get_sgram(int ifile, sg_struct *A);
void rel_sgram(float **S);
void put_image(int ifile,float **image, int size);
void usage (void);
float *malloc_vector_f (long n);
complex *malloc_vector_c (long n) ;
float **malloc_matrix_f (long nr, long nc);
complex **malloc_matrix_c (long nr, long nc);
void prn_gparm(grid_struct *G);
void prn_pswf(pswf_struct *P);
void prn_sgparm(sg_struct *S);

/** Defined in pswf.c **/
void get_pswf(float C, pswf_struct **P);

/** Defined in filters.c  **/
float (*get_filter(char *name))(float);

/** Defined in grid.c  **/
void recon_init(grid_struct *gparm,sg_struct *sgparm);
void do_recon(float **sgram1, float **sgram2, float **image1,
	float **image2);
void recon_free(void);

/** Defined in fft.c (FFT routines from Numerical Recipes) **/
void four1(float data[], unsigned long nn, int isign);
void fourn(float data[], unsigned long nn[], int ndim, int isign);

void rel_sgram(float **S)
{
	free_matrix(S);
}  

static int fexist (char *filename) 
{
	// Check for existence of a file:
	struct stat stbuf;
	int error;
	error=stat (filename, &stbuf);

	if(error) return 0;
	else return 1;
}

float *malloc_vector_f (long n) 
{
	float *v = NULL;

	v = (float *) malloc((Py_ssize_t) (n * sizeof(float)));
	if (!v) {
		fprintf (stderr, "malloc error in malloc_vector_f for length %ld.\n", n);
		v = NULL;
		return v;
	}
	return v;
}

complex *malloc_vector_c (long n) 
{
	complex *v = NULL;

	v = (complex *) malloc((Py_ssize_t) (n * sizeof(complex)));
	if (!v) {
		fprintf (stderr, "malloc error in malloc_vector_c for length %ld.\n", n);
		v = NULL;
		return v;
	}
	return v;
}

float **malloc_matrix_f (long nr, long nc)
{
	float **m = NULL;
	long i;

	// Allocate pointers to rows:
	m = (float **) malloc((Py_ssize_t) (nr * sizeof(float *)));
	if (!m) {
		fprintf (stderr, "malloc error in malloc_matrix_f for %ld row pointers.\n", nr);
		m = NULL;
		return m;
	}
	// Allocate rows and set the pointers to them :
	m[0] = (float *) malloc((Py_ssize_t) (nr * nc * sizeof(float)));
	if (!m[0]) {
		fprintf (stderr, "malloc error in malloc_matrix_f for %ld row with %ld columns.\n", nr, nc);
		m[0] = NULL;
		free (m);
		m = NULL;
		return m;
	}

	for (i = 1; i < nr; i++) m[i] = m[i-1] + nc;

	return m;
}

complex **malloc_matrix_c (long nr, long nc)
{
	complex **m = NULL;
	long i;

	/* Allocate pointers to rows */
	m = (complex **) malloc((Py_ssize_t) (nr * sizeof(complex *)));
	if (!m) {
		fprintf (stderr, "malloc error in malloc_matrix_c for %ld row pointers.\n", nr);
		m = NULL;
		return m;
	}
	/* Allocate rows and set the pointers to them */
	m[0] = (complex *) malloc((Py_ssize_t) (nr * nc * sizeof(complex)));
	if (!m[0]) {
		fprintf (stderr,
			"malloc error in malloc_matrix_c for %ld row with %ld columns.\n",
			nr, nc);
		m[0] = NULL;
		free (m);
		m = NULL;
		return m;
	}
	for (i = 1; i < nr; i++) m[i] = m[i-1] + nc;

	return m;
}

void prn_gparm(grid_struct *G)
{
	//printf("Contents of grid_struct:\n\t");
	prn_pswf(G->pswf);
	/*printf("\n  sampl=%f MaxPixSiz=%fR=%f X0=%f Y0=%f\n  filter=%s ltbl=%d\n",
		G->sampl,G->MaxPixSiz, G->R, G->X0, G->Y0, G->fname, G->ltbl);*/
	return;
}

void prn_pswf(pswf_struct *P)
{
	int i;
	/*printf("Contents of pswf: C=%f nt=%d lmbda=%f\n\tCoefficients:\n",
		P->C, P->nt, P->lmbda);*/
	for(i=0;i<P->nt/2+1;i++)
	{
		//printf("\t\t%f\n",P->coefs[i]);
	}
}

void prn_sgparm(sg_struct *S)
{
	//printf("Contents of sg_struct:\n");
	/*printf("\tn_ang=%d n_det=%d geom=%d center=%f\n",
		S->n_ang,S->n_det,S->geom,S->center);*/
}

static int isconsis(sg_struct *A, sg_struct *B){

	if(A->n_ang!=B->n_ang) return 0;
	if(A->n_det!=B->n_det) return 0;
	if(abs(A->center-B->center) > TOLERANCE) return 0;
	if(A->geom!=B->geom) return 0;
	if(A->geom==0 && A->angles!=B->angles)return 0;
	return 1;
}


static int isnw(sg_struct *A){

	static int n_ang=0, n_det=0, geom=-1;
	static float center=0.0;
	static float* angles=NULL;

	if(
		A->n_ang==n_ang && 
		A->n_det==n_det && 
		abs(A->center-center)< TOLERANCE &&
		A->geom==geom &&
		(geom!=0 || A->angles==angles)
		)
		return 0;
	else{
		n_ang=A->n_ang;
		n_det=A->n_det;		
		center=A->center;
		geom=A->geom;
		if(geom==0) angles=A->angles;
		else angles=NULL;
		return 1;
	}
}   

/**** Macros and typedefs used in this module only ****/

#define Cmult(A,B,C) {(A).r=(B).r*(C).r-(B).i*(C).i;\
	(A).i=(B).r*(C).i+(B).i*(C).r;}
/** A=B*C for complex A, B, C. A must be distinct, and an lvalue */

#ifdef INTERP
#define Cnvlvnt(X) (wtbl[(int)X]+(X-(int)X)*dwtbl[(int)X])
/* Linear interpolation version */
#else
#define Cnvlvnt(X) (wtbl[(int)(X+0.5)])    
/* Nearest nbr version - no interpolation */
#endif

/*** Local variables for this module ***/

static int flag;	
static int n_det, n_ang;
static long pdim, M, M0,M02,ltbl;
static float sampl, scale, L, X0,Y0;
static float *SINE = NULL;
static float *COSE = NULL;
static float *wtbl = NULL;
static float *dwtbl = NULL;
static float *work = NULL;
static float *winv = NULL;
static complex *cproj = NULL;
static complex *filphase = NULL; 
static complex **H = NULL;

/*** Static function prototypes ***/

static void filphase_su(long pd,float fac,
	float(*pf)(float),complex *A);
static void trig_su(sg_struct *SG,float **SP,float **CP);
static void pswf_su(pswf_struct *pswf,long ltbl, 
	long linv, float* wtbl,float* dwtbl,float* winv);
static float legendre(int n,float *coefs, float x);



/*** Function definitions ***/

void recon_init(grid_struct *GP,sg_struct *SGP)

	/************************************************************/
	/* Compute various parameters needed for the reconstruction.*/
	/* Allocate storage for various arrays.                     */
	/* Set up various lookup tables.                            */
	/************************************************************/

{
	float center,C,MaxPixSiz,R,D0,D1;  /* 7/7/98 */
	float (*filter)(float);
	long itmp;

	pswf_struct *pswf;

	n_ang=SGP->n_ang;
	n_det=SGP->n_det;
	center=SGP->center;	

	sampl=GP->sampl;
	MaxPixSiz=GP->MaxPixSiz;
	R=GP->R;
	X0=GP->X0;
	Y0=GP->Y0;
	filter=GP->filter;
	ltbl=GP->ltbl;
	pswf=GP->pswf;
	C=pswf->C;

	/** Set flag if ROI offset exists **/
	if(X0!=0.||Y0!=0.)flag=1;  
	else flag=0;

	/*** Compute pdim = next power of 2 >=n_det */ 
	pdim=1;
	itmp=n_det-1;
	while(itmp)
	{
		pdim<<=1;
		itmp>>=1;
	}

	D0=R*n_det;  /* Size of ROI */  /* 7/7/98 */
	D1=sampl*D0; /* Size of "extended" region. 7/7/98 */

	/*** Compute raster size M for the oversampled 2D array: */
	/*	  M = next power of two >= D1/MaxPixSiz */  /*7/7/98 */
	M=1;
	itmp = (long) (D1/MaxPixSiz-1);  /* 7/7/98 */
	while(itmp)
	{
		M<<=1;
		itmp>>=1;
	}

	/*** Compute M0 = raster size for the ROI =
	= largest ODD integer <= M/sampl */

	M02 = (long) floor(M/2/sampl-0.5);
	M0=2*M02+1;

	/****/

	GP->sampl=sampl=(float)M/M0; /* Adjust value of sampl */
	D1=sampl*D0;		     /*   and of D1 */  /* 7/7/98 */ 	
	L=(float) (2*C*sampl/pi);	  /* Size of convolution support square */
	scale=D1/pdim;   /* Converts freqs to grid units*/ /* 7/7/98 */

	/*** Allocate storage for various arrays */

	cproj=malloc_vector_c(pdim); 
	filphase=malloc_vector_c(pdim/2);	
	H=malloc_matrix_c(M,M);
	wtbl=malloc_vector_f(ltbl+1);	
#ifdef INTERP
	dwtbl=malloc_vector_f(ltbl+1);	
#endif
	winv=malloc_vector_f(M0);
	work=malloc_vector_f((int)L+1);

	/*** Set up table of sines and cosines ***/

	trig_su(SGP,&SINE,&COSE);

	/*** Set up table of combined filter-phase factors */

	filphase_su(pdim,center,filter,filphase);           

	/*** Set up PSWF lookup tables */

	pswf_su(pswf,ltbl,M02,wtbl,dwtbl,winv);

}   /*** End recon_init() ***/




void do_recon(float** G1,float** G2,float** S1,float** S2)

	/************************************************************/
	/* Reconstruct two real slice images from their sinugrams   */
	/* using the "gridding" algorithm applied to complex data.  */
	/************************************************************/

{
	//*imgsiz=M0;

	/*if(verbose)printf(
		"do_recon(): M0=%d M= %d pdim=%d L=%f scale=%f\n",
		M0,M,pdim,L,scale);*/

	{	/*** First clear the array H ***/
		long iu,iv;
		//#pragma omp parallel for private(iv)
		for(iu=0;iu<M;iu++) 
			for(iv=0;iv<M;iv++)
				H[iu][iv].r=H[iu][iv].i=0.0;
	}

	//if(verbose)printf("Start Phase 1\n");

	{	/***Phase 1 ***************************************

		Loop over the n_ang projection angles. For each angle, do the following:

		1. Copy the real projection data from the two slices into the
		real and imaginary parts of the first n_det elements of the 
		complex array, cproj[].  Set the remaining pdim-n_det elements
		to zero (zero-padding).

		2. Carry out a (1D) Fourier transform on the complex data.
		This results in transform data that is arranged in 
		"wrap-around" order, with non-negative spatial frequencies 
		occupying the first half, and negative frequencies the second 
		half, of the array, cproj[].

		3. Multiply each element of the 1-D transform by a complex,
		frequency dependent factor, filphase[].  These factors were
		precomputed as part of recon_init() and combine the 
		tomographic filtering with a phase factor which shifts the 
		origin in configuration space to the projection of the 
		rotation axis as defined by the parameter, "center".  If a 
		region of interest (ROI) centered on a different origin has 
		been specified [(X0,Y0)!=(0,0)], multiplication by an 
		additional phase factor, dependent on angle as well as 
		frequency, is required.

		4. For each data element, find the Cartesian coordinates, 
		<U,V>, of the corresponding point in the 2D frequency plane, 
		in  units of the spacing in the MxM rectangular grid placed 
		thereon; then calculate the upper and lower limits in each 
		coordinate direction of the integer coordinates for the 
		grid points contained in an LxL box centered on <U,V>.  
		Using a precomputed table of the (1-D) convolving function, 
		W, calculate the contribution of this data element to the
		(2-D) convolvent (the 2_D convolvent is the product of
		1_D convolvents in the X and Y directions) at each of these
		grid points, and update the complex 2D array H accordingly.  


		At the end of Phase 1, the array H[][] contains data arranged in 
		"natural", rather than wrap-around order -- that is, the origin in 
		the spatial frequency plane is situated in the middle, rather than 
		at the beginning, of the array, H[][].  This simplifies the code 
		for carrying out the convolution (step 4 above), but necessitates 
		an additional correction -- See Phase 3 below.
		***********************************************************************/

		complex Cdata1,Cdata2,Ctmp;
		float U,V,rtmp;
		float L2 = (float) (L/2);
		float convolv;
		float tblspcg = (float) (2*ltbl/L);

		long pdim2=pdim>>1,M2=M>>1,
			iul,iuh,iu,ivl,ivh,iv,n;

		/* Following are to handle offset ROI case */
		float offset;
		complex phfac;
		int j,k;


		//#pragma omp parallel for private(j,k,offset,Cdata1,Cdata2,Ctmp,phfac,U,V,rtmp,iul,iuh,iu,ivl,ivh,iv,n)
		for(n=0;n<n_ang;n++)     /*** Start loop on angles */
		{

			if(flag) 
				offset= (float) ((X0*COSE[n]+Y0*SINE[n])*pi);


			j=0;
			while(j<n_det)	
			{     
				cproj[j].r=G1[n][j];
				cproj[j].i=G2[n][j];
				j++;
			}

			while(j<pdim)	/** Zero fill the rest of array **/
			{
				cproj[j].r=cproj[j].i=0.0;
				j++;
			}

			four1((float*)cproj-1,pdim,1);   
			/* FFT -- cf. Numerical Recipes */   	

			for(j=1;j<pdim2;j++)
			{  	/* Start loop on transform data */			

				if(!flag)
				{
					Ctmp.r=filphase[j].r;
					Ctmp.i=filphase[j].i;
				}
				else
				{
					phfac.r = cos(j*offset);
					phfac.i = -sin(j*offset);
					Cmult(Ctmp,filphase[j],phfac);
				}

				Cmult(Cdata1,Ctmp,cproj[j])
					Ctmp.i=-Ctmp.i;
				Cmult(Cdata2,Ctmp,cproj[pdim-j])

					U=(rtmp=scale*j)*COSE[n]+M2; /* X direction*/
				V=rtmp*SINE[n]+M2;	      /* Y direction*/	

				/* Note freq space origin is at (M2,M2), but we
				offset the indices U, V, etc. to range from 0 to M-1 */

				iul = (long) (ceil(U-L2));
				iuh = (long) (floor(U+L2));
				ivl = (long) (ceil(V-L2));
				ivh = (long) (floor(V+L2));
				if(iul<1)iul=1;
				if(iuh>=M)iuh=M-1; 
				if(ivl<1)ivl=1;
				if(ivh>=M)ivh=M-1; 

				/* Note aliasing value (at index=0) is forced to zero */	

				for(iv=ivl,k=0;iv<=ivh;iv++,k++)
					work[k]=Cnvlvnt(abs(V-iv)*tblspcg);
				for(iu=iul;iu<=iuh;iu++)
				{
					rtmp=Cnvlvnt(abs(U-iu)*tblspcg);
					for(iv=ivl,k=0;iv<=ivh;iv++,k++)
					{
						convolv = rtmp*work[k];
						H[iu][iv].r += convolv*Cdata1.r;
						H[iu][iv].i += convolv*Cdata1.i;
						H[M-iu][M-iv].r += convolv*Cdata2.r;
						H[M-iu][M-iv].i += convolv*Cdata2.i;
					}
				}
			} /*** End loop on transform data */
		} /*** End loop on angles */

	}  /*** End phase 1 ************************************************/	


	//if(verbose)printf("Start Phase 2\n");

	{	/*** Phase 2 ********************************************

		Carry out a 2D inverse FFT on the array H.

		At the conclusion of this phase, the configuration 
		space data is arranged in wrap-around order with the origin
		(center of reconstructed images) situated at the start of the 
		array.  The first (resp. second) half of the array contains the lower,
		Y<0 (resp, upper Y>0) part of the image, and within each row of the 
		array, the first (resp. second) half contains data for the right [X>0]
		(resp. left [X<0]) half of the image.

		********************************************************************/

		unsigned long H_size[2];
		H_size[0]=H_size[1]=M;
		fourn((float*)(*H)-1,H_size-1,2,-1);  
		/* Inverse FFT- Numer Recipes */

	}  /*** End phase 2 ************************************************/

	/*** Release input buffers, then allocate output buffers ***/

	//if(G1!=G2)rel_sgram(G2);
	//rel_sgram(G1);

	//*S1=malloc_matrix_f(M0,M0);
	//*S2=malloc_matrix_f(M0,M0);


	//if(verbose)printf("Start Phase 3\n");

	{ /*** Phase 3 ******************************************************

	  Copy the real and imaginary parts of the complex data from H[][],
	  into the output buffers for the two reconstructed real images, 
	  simultaneously carrying out a final multiplicative correction.  
	  The correction factors are taken from the array, winv[], previously 
	  computed in pswf_su(), and consist logically of three parts, namely:

	  1. A positive real factor, corresponding to the reciprocal
	  of the inverse Fourier transform, of the convolving
	  function, W, and

	  2. Multiplication by the cell size, (1/D1)^2, in 2D frequency
	  space.  This correctly normalizes the 2D inverse FFT carried
	  out in Phase 2.  (Note that all quantities are ewxpressed in
	  units in which the detector spacing is one.)

	  3. A sign change for the "odd-numbered" elements (in a 
	  checkerboard pattern) of the array.  This compensates
	  for the fact that the 2-D Fourier transform (Phase 2) 
	  started with a frequency array in which the zero frequency 
	  point appears in the middle of the array instead of at 
	  its start.

	  Only the elements in the square M0xM0 subarray of H[][], centered 
	  about the origin, are utilized.  The other elements are not part of the
	  actual region being reconstructed and are discarded.  Because of the 
	  wrap-around ordering, the subarray must actually be taken from the four
	  corners" of the 2D array, H[][] -- See Phase 2 description, above.

	  The final data correponds physically to the linear X-ray absorption
	  coefficient expressed in units of the inverse detector spacing -- to 
	  convert to inverse cm (say), one must divide the data by the detector 
	  spacing in cm.

	  *********************************************************************/

		long iu,iv,j,k,ustart,vstart,ufin,vfin;
		float corrn_u,corrn,**T1,**T2;
		//T1=*S1; T2=*S2;
		T1=S1; T2=S2;

		j=0;
		ustart=(M-M02);
		ufin=M;
		while(j<M0)
		{
			for(iu=ustart;iu<ufin;j++,iu++)
			{
				corrn_u=winv[j];
				k=0;
				vstart=(M-M02);
				vfin=M;
				while(k<M0)
				{
					for(iv=vstart;iv<vfin;k++,iv++)
					{
						corrn=corrn_u*winv[k]; 
						T1[j][k]=corrn*H[iu][iv].r;
						T2[j][k]=corrn*H[iu][iv].i;
					}
					if(k<M0) (vstart=0,vfin=M02+1);
				}
			}
			if(j<M0) (ustart=0,ufin=M02+1);
		}

	}  /*** End phase 3 *******************************************************/

	return;

}  /*** End do_recon() ***/


void recon_free(void)

	/***** Release memory allocated by recon_init() ***/

{
	if (SINE != NULL) {free(SINE);SINE = NULL;}
	if (COSE != NULL) {free(COSE);COSE = NULL;}
	if (cproj != NULL) {free(cproj);cproj = NULL;}
	if (filphase != NULL) {free(filphase);filphase = NULL;}
	if (wtbl != NULL) {free(wtbl);wtbl = NULL;}
#ifdef INTERP
	if (dwtbl != NULL) {free(dwtbl);dwtbl = NULL;}
#endif
	if (winv != NULL) {free(winv);winv = NULL;}
	if (work != NULL) {free(work);work = NULL;}
	if (H != NULL) {free_matrix(H);H = NULL;}
}  /*** End recon-free() ***/


static void filphase_su(long pd, float center,
	float(*pf)(float), complex *A)

	/******************************************************************/
	/* Set up the complex array, filphase[], each element of which    */
	/* consists of a real filter factor [obtained from the function,  */
	/* (*pf)()], multiplying a complex phase factor (derived from the */
	/* parameter, center}.  See Phase 1 comments in do_recon(), above.*/
	/******************************************************************/

{ 
	long j,pd2=pd>>1;
	float x;
	float rtmp1 = (float) (2*pi*center/pd);
	float rtmp2;
	float norm = (float) (pi/pd/n_ang);	/* Normalization factor for
							back transform  7/7/98  */

	for(j=0;j<pd2;j++)
	{
		x=j*rtmp1;
		rtmp2=(*pf)((float)j/pd)*norm;	/* 7/7/98 */
		A[j].r=rtmp2*cos(x);
		A[j].i=-rtmp2*sin(x);
	}

}  /*** End filphase_su() ***/

static void pswf_su(pswf_struct *pswf,long ltbl, long linv, 
	float* wtbl,float* dwtbl,float* winv)

	/*************************************************************/
	/* Set up lookup tables for convolvent (used in Phase 1 of   */
	/* do_recon()), and for the final correction factor (used in */
	/* Phase 3).                                                  */
	/*************************************************************/

{
	float C,*coefs, lmbda, polyz,norm,fac;
	long i;
	int nt;

	C=pswf->C;
	nt=pswf->nt;
	coefs=pswf->coefs;		
	lmbda=pswf->lmbda;	
	polyz=legendre(nt,coefs,0.);


	wtbl[0]=1.0;
	for(i=1;i<=ltbl;i++) 
	{	wtbl[i]=legendre(nt,coefs,(float)i/ltbl)/polyz;
#ifdef INTERP
	dwtbl[i]=wtbl[i]-wtbl[i-1];
#endif
	}

	fac  = (float) (ltbl/(linv+0.5));
	norm = (float) (sqrt(pi/2/C/lmbda)/sampl);	/* 7/7/98 */

	/* Note the final result at end of Phase 3 contains the factor, 
	norm^2.  This incorporates the normalization of the 2D
	inverse FFT in Phase 2 as well as scale factors involved
	in the inverse Fourier transform of the convolvent.
	7/7/98 			*/

	winv[linv]=norm/Cnvlvnt(0.);
	for(i=1;i<=linv;i++)
	{
		norm=-norm; 
		/* Minus sign for alternate entries
		corrects for "natural" data layout
		in array H at end of Phase 1.  */

		winv[linv+i]=winv[linv-i]=norm/Cnvlvnt(i*fac);	
	}

}   /*** End pswf_su() ***/


static void trig_su(sg_struct *SG, float **SP, float **CP)

	/*********** Set up tables of sines and cosines. ***********/

{
	int j, geom=SG->geom, n_ang=SG->n_ang;
	float *S,*C;

	*SP=S=malloc_vector_f(n_ang);
	*CP=C=malloc_vector_f(n_ang);
	switch (geom)
	{
	case 0:
		{
			float theta;
			//float degtorad = (float) (pi/180.0);
			float *angle=SG->angles;

			for(j=0;j<n_ang;j++)
			{
				//theta=degtorad*angle[j];
				theta=angle[j];
				S[j]=sin(theta);
				C[j]=cos(theta);
				//printf("S[%d], C[%d]=[%0.3f ; %0.3f]\n", j, j, S[j], C[j]  );
			}
		}
		break;

	case 1:
	case 2:
		{
			float dtheta = (float) (geom*pi/n_ang);
			float dcos;
			float dsin;

			dcos=cos(dtheta);
			dsin=sin(dtheta);
			S[0]=0.0; C[0]=1.0;
			for(j=1;j<n_ang;j++)
			{
				S[j]=dcos*S[j-1]+dsin*C[j-1];
				C[j]=dcos*C[j-1]-dsin*S[j-1];
				
				//printf("S[%d], C[%d]=[%0.3f ; %0.3f]\n", j, j, S[j], C[j]  );
			}
		}
		break;

	default: 
		fprintf(stderr,
			"Illegal value for angle geometry indicator.\n");
		exit(2);
	}  /** End switch **/

}	/*** End trig_su ***/


static 
	float legendre(int n,float *coefs, float x)

	/***************************************************
	*                                                  *
	*    Compute SUM(coefs(k)*P(2*k,x), for k=0,n/2)  *
	*                                                  *
	*    where P(j,x) is the jth Legendre polynomial   *
	*                                                  *
	***************************************************/
{
	float penult,last,nw,y;
	int j,k,even;
	if(x>1||x<-1)
	{
		fprintf(stderr, 
			"\nInvalid argument to legendre()");
		exit(2);
	}
	y=coefs[0];
	penult=1.;
	last=x;
	even=1;
	k=1;
	for(j=2;j<=n;j++)
	{
		nw=(x*(2*j-1)*last-(j-1)*penult)/j;
		if(even)
		{
			y+=nw*coefs[k];
			even=0;
			k++;
		}
		else	even=1;

		penult=last;
		last=nw;
	}
	return y;

}   /*** End legendre() ***/


/******************************************************************************/
//
// Part for pswf:
//

#define NO_PSWFS 5

pswf_struct pswf_db[NO_PSWFS]={

	{4.0F, 16, 0.99588549F,
	{0.5239891E+01F,	-0.5308499E+01F,	 0.1184591E+01F, 
	-0.1230763E-00F,	 0.7371623E-02F,	-0.2864074E-03F,
	0.7789983E-05F,	-0.1564700E-06F,	 0.2414647E-08F}},

	{4.2F, 16, 0.99657887F,
	{0.6062942E+01F,	-0.6450252E+01F,	 0.1551875E+01F,
	-0.1755960E-01F,	 0.1150712E-01F,	-0.4903653E-03F,
	0.1464986E-04F,	-0.3235110E-06F,	 0.5492141E-08F}},

	{5.0F, 18, 0.99935241F,
	{0.1115509E+02F,	-0.1384861E+02F,	 0.4289811E+01F,
	-0.6514303E-00F,	 0.5844993E-01F,	-0.3447736E-02F,
	0.1435066E-03F,	-0.4433680E-05F,	 0.1056040E-06F,
	-0.1997173E-08F}},

	{6.0F, 18, 0.9990188F,
	{0.2495593E+02F,	-0.3531124E+02F,	 0.1383722E+02F,
	-0.2799028E+01F,	 0.3437217E-00F,	-0.2818024E-01F,
	0.1645842E-02F,	-0.7179160E-04F,	 0.2424510E-05F,
	-0.6520875E-07F}},

	{7.0F, 20, 0.99998546F,
	{0.5767616E+02F,	-0.8931343E+02F,	 0.4167596E+02F,
	-0.1053599E+02F,	 0.1662374E+01F,	-0.1780527E-00F,
	0.1372983E-01F,	-0.7963169E-03F,	 0.3593372E-04F,
	-0.1295941E-05F,	 0.3817796E-07F}}
};

void get_pswf(float C, pswf_struct **P){
	int i=0;
	while(i<NO_PSWFS && abs(C-pswf_db[i].C)>0.01) i++;
	if(i>=NO_PSWFS) 
	{
		fprintf(stderr,
			"Prolate parameter, C = %f not in data base\n",C);
		exit(2);
	}
	*P=&pswf_db[i];
	return;
}

/******************************************************************************/
//
// Part for filters:
//

#define LEN_FTBL 7	/* Number of entries in ftbl[] below */

float shlo(float x){	/* Shepp-Logan filter */
	//return (float) (abs(sin(pi*x)/pi));
	return (float) (x);
}

float hann(float x){	/* Hann filter */
	return (float) (abs(x)*0.5*(1.+cos(2*pi*x)));
}

float hamm(float x){	/* Hamming filter */
	return (float) (abs(x)*(0.54+0.46*cos(2*pi*x)));
}

float ramp(float x){	/* Ramp filter */
	return (float) (abs(x));
}

struct {char* name; float (*fp)(float);} fltbl[]= {
	{"shlo",shlo},		/* The default choice */
	{"shepp",shlo},
	{"hann",hann},
	{"hamm",hamm},
	{"hamming",hamm},
	{"ramp",ramp},
	{"ramlak",ramp}
};	

float (*get_filter(char *name))(float) {
	int i;
	for(i=0;i<LEN_FTBL;i++){
		if(!strcmp(name,fltbl[i].name)){
			return fltbl[i].fp;
		}
	}
	/* If no match, set default choice */
	//printf(" Using default filter: %s.\n",fltbl[0].name);
	strcpy(name,fltbl[0].name);
	return fltbl[0].fp;   
}  /*** End get_filter() ***/


/******************************************************************************/
//
// Part for FFT:
//

/** 1-D and n-D routines from Numerical Recipes. **/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], unsigned long nn, int isign)
{

	/*fftw_init_threads();
	fftw_plan_with_nthreads(2);

	fftw_execute(plan);
	fftw_destroy_plan(plan);
	fftw_cleanup_threads();*/

	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=(float) (wr*data[j]-wi*data[j+1]);
				tempi=(float) (wr*data[j+1]+wi*data[j]);
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

/* (C) Copr. 1986-92 Numerical Recipes Software 9='0='. */

void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
#undef SWAP


/******************************************************************************/
//
// Part for Python wrapper
//
float **_ptrvector(long n) { 
	float **v = malloc( n* sizeof(*v) );
	if (!v)   {
		printf("In **ptrvector. Allocation of memory for array failed.");
		exit(0); }
	return v;
}

float **pymatrix_to_carray(PyArrayObject *arrayin)  {
	float **c, *a;
	int i,n,m;

	n = PyArray_DIM(arrayin, 0);//arrayin->dimensions[0];
	m = PyArray_DIM(arrayin, 1);//m=arrayin->dimensions[1];
	c=_ptrvector(n);
	a = (float*) PyArray_DATA(arrayin); //	a=(float *) arrayin->data; 
	for ( i=0; i<n; i++)  {
		c[i]=a+i*m;  }	
	return c;
}

float **pymatrix_to_carray1D(PyArrayObject *arrayin)  {
	float **c, *a;
	int i,n,m;

	n = PyArray_DIM(arrayin, 0);//arrayin->dimensions[0];
	m = 0;
	c=_ptrvector(n);
	a = (float*) PyArray_DATA(arrayin); //	a=(float *) arrayin->data; 
	for ( i=0; i<n; i++)  {
		c[i]=a+i*m;  }	
	return c;
}

void free_carrayptrs(float **v)  {
	free((char*) v);
}

static PyObject *py_paralrec(PyObject *self, PyObject *args, PyObject *keywds)
{        
	PyArrayObject *in_array, *out_array;
	PyArrayObject *in_array2, *out_array2;
	PyArrayObject *angles;
	float         **cin_array, **cout_array;
	float         **cin_array2, **cout_array2;
	float         *c_angles;

	float center = -1.0; // Default means half of image width
	char *filter = "shlo";
	char *geom   = "half";
	float zoom   = 1.0;
	float zoom_x0 = 0.0;
	float zoom_y0 = 0.0;
	float oversampling = 1.2;
	int dims[2] = {0, 0};
	sg_struct  sg_parm;
	grid_struct grid_parm;
		
	/* Operate on the matrices  */
	int i;
	float C = 6.0;	
    
	static char *kwlist[] = {"angles", "oversampling", "center", "filter", "zoom", "zoom_origin",  NULL};
	

	if (!PyArg_ParseTupleAndKeywords(args, keywds, "O!O!O!|ffsf(ff)", kwlist,
		&PyArray_Type, &in_array, &PyArray_Type, &in_array2, &PyArray_Type, &angles,
		&oversampling, &center, &filter, &zoom, &zoom_x0, &zoom_y0 ))
        return NULL;

	Py_XINCREF(in_array);
	Py_XINCREF(in_array2);
	Py_XINCREF(angles);

	// Check if input array is a bi-dimensional image:
	if ( (PyArray_NDIM(in_array) != 2) ||  (PyArray_NDIM(in_array2) != 2) ) {	
		PyErr_Format(PyExc_ValueError, "Input array is not bi-dimensional.");		
		return NULL; 
	}

	// Check if input array is of type float:
	/*if ( (PyArray_TYPE(in_array) != PyArray_FLOAT) || (PyArray_TYPE(in_array2) != PyArray_FLOAT) ) 		{	
		PyErr_Format(PyExc_TypeError, "Input array is not of type float32.");
		return NULL;
	}*/	

	// Change contiguous arrays into C ** arrays (Memory is Allocated!) */
	cin_array = pymatrix_to_carray(in_array);    
	cin_array2 = pymatrix_to_carray(in_array2); 
	c_angles = (float*) PyArray_DATA(angles);
	

	M0 = PyArray_DIM(in_array, 1);	
	sg_parm.n_det = PyArray_DIM(in_array, 1);//in_array->dimensions[1];
	sg_parm.n_ang = PyArray_DIM(in_array, 0);//in_array->dimensions[0];

	if ( abs(center + 1.0) < 0.001 )
		sg_parm.center = (float) (sg_parm.n_det / 2.0);
	else
		sg_parm.center = center;
	/*sg_parm.angles = NULL;        	
	if (strcmp (geom,"full") == 0 )
		sg_parm.geom   = 2;
	else
		sg_parm.geom   = 1;*/
	sg_parm.geom = 0;
	sg_parm.angles =  c_angles;
	
	grid_parm.sampl = oversampling; //e.g. 1.2      /* "Oversampling" ratio */
	grid_parm.MaxPixSiz = 1.0F; // e,g. 1.0 /* Max pixel size for reconstruction */
	grid_parm.R  = zoom;  // e.g. 1.0	   /* ROI Size (Zoom...) */
	grid_parm.X0 = zoom_x0; // e,g, 0.0  	   /* (X0,Y0)=Offset of ROI from rotation axis, */
	grid_parm.Y0 = zoom_y0; // e,g, 0.0 	   /* in units of center-to-edge distance.  */
	grid_parm.ltbl=LTBL_DEF;	           /* No. elements in convolvent lookup tables. */       
	for (i = 0; i < 16; i++)        
		grid_parm.fname[i] = (filter)[i];   /* {"shlo",shlo},
											{"shepp",shlo},
											{"hann",hann},
											{"hamm",hamm},
											{"hamming",hamm},
											{"ramp",ramp},
											{"ramlak",ramp} */


	grid_parm.filter=get_filter(grid_parm.fname);
	get_pswf(C, &grid_parm.pswf);	


	// Init reconstruction engine:
	recon_init(&grid_parm,&sg_parm);	        


	// Change contiguous arrays into C ** arrays (Memory is Allocated!) */
	dims[0] = M0;
	dims[1] = M0;
	out_array   = (PyArrayObject *) PyArray_FromDims(2,dims,NPY_FLOAT);      
	out_array2  = (PyArrayObject *) PyArray_FromDims(2,dims,NPY_FLOAT);      
	cout_array  = pymatrix_to_carray(out_array);
	cout_array2 = pymatrix_to_carray(out_array2);

	// Do actual reconstruction:    
	do_recon(cin_array,cin_array2,cout_array,cout_array2);     

	
	// Release resources:
	recon_free();
	//free_carrayptrs(c_angles);
	free_carrayptrs(cin_array);
	free_carrayptrs(cin_array2);
	free_carrayptrs(cout_array);
	free_carrayptrs(cout_array2);
	
	return Py_BuildValue("OO", out_array, out_array2);

//fail:
	// Are the following 4 lines really necessary? Boh....
	Py_XDECREF(in_array);
	Py_XDECREF(in_array2);
	Py_XDECREF(out_array);
	Py_XDECREF(out_array2);
	
	Py_INCREF(Py_None);
	return Py_None;
}


char py_gridrec_doc[] = "Do Reconstruction.";


/*****************************************************************************/
/* Create Python module */

char module_doc[] = "A Python C extension module for calling GridRec.\n";


static PyMethodDef module_methods[] = {
	{"paralrecon", (PyCFunction) py_paralrec, METH_VARARGS|METH_KEYWORDS,
	py_gridrec_doc},
	/*{"decodelzw", (PyCFunction)py_decodelzw, METH_VARARGS,
	py_decodelzw_doc},
	{"decodepackbits", (PyCFunction)py_decodepackbits, METH_VARARGS,
	py_decodepackbits_doc},*/
	{NULL, NULL, 0, NULL} /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3

struct module_state {
	PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static int module_traverse(PyObject *m, visitproc visit, void *arg) {
	Py_VISIT(GETSTATE(m)->error);
	return 0;
}

static int module_clear(PyObject *m) {
	Py_CLEAR(GETSTATE(m)->error);
	return 0;
}

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"_gridrec",
	NULL,
	sizeof(struct module_state),
	module_methods,
	NULL,
	module_traverse,
	module_clear,
	NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
	PyInit__gridrec(void)

#else

#define INITERROR return

PyMODINIT_FUNC
	init_gridrec(void)

#endif
{
	PyObject *module;

	char *doc = (char *)PyMem_Malloc(sizeof(module_doc) + sizeof(_VERSION_));
	sprintf(doc, module_doc, _VERSION_);

#if PY_MAJOR_VERSION >= 3
	moduledef.m_doc = doc;
	module = PyModule_Create(&moduledef);
#else
	module = Py_InitModule3("_gridrec", module_methods, doc);
#endif

	PyMem_Free(doc);

	if (module == NULL)
		INITERROR;

	if (_import_array() < 0) {
		Py_DECREF(module);
		INITERROR;
	}

	{
#if PY_MAJOR_VERSION < 3
		PyObject *s = PyString_FromString(_VERSION_);
#else
		PyObject *s = PyUnicode_FromString(_VERSION_);
#endif
		PyObject *dict = PyModule_GetDict(module);
		PyDict_SetItemString(dict, "__version__", s);
		Py_DECREF(s);
	}

#if PY_MAJOR_VERSION >= 3
	return module;
#endif
}

//#ifdef __cplusplus
//} //extern "C"
//#endif
