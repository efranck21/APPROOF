//******************************************************************************
//
//    Program for solving a problem of type Ax = b using the SVD method
//    in the least squares manner.
//
//  Source:
//
//    This code was originally taken from:
//       http://people.sc.fsu.edu/~jburkardt/cpp_src/qr_solve/qr_solve.html
//    We refer to this site or to the author (John Bukardt) for more details.
//    The program was simplified and detached from adjacent libraries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2014
//
//  Original Author:
//
//    John Burkardt
//
//  Modified by:
//
//    Laura S. Mendoza
//    Emmanuel Franck
//
//******************************************************************************


# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <algorithm>

using namespace std;

# include "svd_solver.hpp"

const double epsilon = 2.220446049250313E-016;
//****************************************************************************

void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )

//****************************************************************************
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
//    Input, int INCY, the increment between successive entries of DY.
//
{

  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
    {
      return;
    }

  if ( da == 0.0 )
    {
      return;
    }
  //
  //  Code for unequal increments or equal increments
  //  not equal to 1.
  //
  if ( incx != 1 || incy != 1 )
    {
      if ( 0 <= incx )
       {
        ix = 0;
       }
      else
       {
        ix = ( - n + 1 ) * incx;
       }

      if ( 0 <= incy )
       {
         iy = 0;
       }
      else
       {
         iy = ( - n + 1 ) * incy;
       }

      for ( i = 0; i < n; i++ )
       {
        dy[iy] = dy[iy] + da * dx[ix];
        ix = ix + incx;
        iy = iy + incy;
       }
    }
  //
  //  Code for both increments equal to 1.
  //
  else
    {
      m = n % 4;

      for ( i = 0; i < m; i++ )
	{
	  dy[i] = dy[i] + da * dx[i];
	}

      for ( i = m; i < n; i = i + 4 )
	{
	  dy[i  ] = dy[i  ] + da * dx[i  ];
	  dy[i+1] = dy[i+1] + da * dx[i+1];
	  dy[i+2] = dy[i+2] + da * dx[i+2];
	  dy[i+3] = dy[i+3] + da * dx[i+3];
	}

    }

  return;
}
//****************************************************************************80

double ddot ( int n, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DDOT forms the dot product of two vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries in DX.
//
//    Input, double DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries in DY.
//
//    Output, double DDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
  double dtemp;
  int i;
  int ix;
  int iy;
  int m;

  dtemp = 0.0;

  if ( n <= 0 )
    {
      return dtemp;
    }
  //
  //  Code for unequal increments or equal increments
  //  not equal to 1.
  //
  if ( incx != 1 || incy != 1 )
    {
      if ( 0 <= incx )
	{
	  ix = 0;
	}
      else
	{
	  ix = ( - n + 1 ) * incx;
	}

      if ( 0 <= incy )
	{
	  iy = 0;
	}
      else
	{
	  iy = ( - n + 1 ) * incy;
	}

      for ( i = 0; i < n; i++ )
	{
	  dtemp = dtemp + dx[ix] * dy[iy];
	  ix = ix + incx;
	  iy = iy + incy;
	}
    }
  //
  //  Code for both increments equal to 1.
  //
  else
    {
      m = n % 5;

      for ( i = 0; i < m; i++ )
	{
	  dtemp = dtemp + dx[i] * dy[i];
	}

      for ( i = m; i < n; i = i + 5 )
	{
      dtemp = dtemp + dx[i  ] * dy[i  ]
                    + dx[i+1] * dy[i+1]
                    + dx[i+2] * dy[i+2]
                    + dx[i+3] * dy[i+3]
	+ dx[i+4] * dy[i+4];
	}

    }

  return dtemp;
}
//****************************************************************************80

double dnrm2 ( int n, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DNRM2 returns the euclidean norm of a vector.
//
//  Discussion:
//
//     DNRM2 ( X ) = sqrt ( X' * X )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[*], the vector whose norm is to be computed.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Output, double DNRM2, the Euclidean norm of X.
//
{
  double absxi;
  int i;
  int ix;
  double norm;
  double scale;
  double ssq;
  double value;

  if ( n < 1 || incx < 1 )
    {
      norm = 0.0;
    }
  else if ( n == 1 )
    {
      norm = fabs ( x[0] );
    }
  else
    {
      scale = 0.0;
      ssq = 1.0;
      ix = 0;

      for ( i = 0; i < n; i++ )
	{
	  if ( x[ix] != 0.0 )
	    {
	      absxi = fabs ( x[ix] );
	      if ( scale < absxi )
		{
		  ssq = 1.0 + ssq * ( scale / absxi ) * ( scale / absxi );
		  scale = absxi;
		}
	      else
		{
		  ssq = ssq + ( absxi / scale ) * ( absxi / scale );
		}
	    }
	  ix = ix + incx;
	}

      norm  = scale * sqrt ( ssq );
    }

  return norm;
}
//****************************************************************************80


void drot ( int n, double x[], int incx, double y[], int incy, double c,
	    double s )

//****************************************************************************80
//
//  Purpose:
//
//    DROT applies a plane rotation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to be rotated.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to be rotated.
//
//    Input, int INCY, the increment between successive elements of Y.
//
//    Input, double C, S, parameters (presumably the cosine and
//    sine of some angle) that define a plane rotation.
//
{
  int i;
  int ix;
  int iy;
  double stemp;

  if ( n <= 0 )
    {
    }
  else if ( incx == 1 && incy == 1 )
    {
      for ( i = 0; i < n; i++ )
	{
	  stemp = c * x[i] + s * y[i];
	  y[i]  = c * y[i] - s * x[i];
	  x[i]  = stemp;
	}
    }
  else
    {
      if ( 0 <= incx )
	{
	  ix = 0;
	}
      else
	{
	  ix = ( - n + 1 ) * incx;
	}

      if ( 0 <= incy )
	{
	  iy = 0;
	}
      else
	{
	  iy = ( - n + 1 ) * incy;
	}

      for ( i = 0; i < n; i++ )
	{
	  stemp = c * x[ix] + s * y[iy];
	  y[iy] = c * y[iy] - s * x[ix];
	  x[ix] = stemp;
	  ix = ix + incx;
	  iy = iy + incy;
	}

    }

  return;
}
//****************************************************************************80

void drotg ( double *sa, double *sb, double *c, double *s )

//****************************************************************************80
//
//  Purpose:
//
//    DROTG constructs a Givens plane rotation.
//
//  Discussion:
//
//    Given values A and B, this routine computes
//
//    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
//          = sign ( B ) if abs ( A ) <= abs ( B );
//
//    R     = SIGMA * ( A * A + B * B );
//
//    C = A / R if R is not 0
//      = 1     if R is 0;
//
//    S = B / R if R is not 0,
//        0     if R is 0.
//
//    The computed numbers then satisfy the equation
//
//    (  C  S ) ( A ) = ( R )
//    ( -S  C ) ( B ) = ( 0 )
//
//    The routine also computes
//
//    Z = S     if abs ( A ) > abs ( B ),
//      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
//      = 1     if C is 0.
//
//    The single value Z encodes C and S, and hence the rotation:
//
//    If Z = 1, set C = 0 and S = 1;
//    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
//    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input/output, double *SA, *SB,  On input, SA and SB are the values
//    A and B.  On output, SA is overwritten with R, and SB is
//    overwritten with Z.
//
//    Output, double *C, *S, the cosine and sine of the Givens rotation.
//
{
  double r;
  double roe;
  double scale;
  double z;

  if ( fabs ( *sb ) < fabs ( *sa ) )
    {
      roe = *sa;
    }
  else
    {
      roe = *sb;
    }

  scale = fabs ( *sa ) + fabs ( *sb );

  if ( scale == 0.0 )
    {
      *c = 1.0;
      *s = 0.0;
      r = 0.0;
    }
  else
    {
      r = scale * sqrt ( ( *sa / scale ) * ( *sa / scale )
           + ( *sb / scale ) * ( *sb / scale ) );
      r = copysign(1.0, roe ) * r;
      *c = *sa / r;
      *s = *sb / r;
    }

  if ( 0.0 < fabs ( *c ) && fabs ( *c ) <= *s )
    {
      z = 1.0 / *c;
    }
  else
    {
      z = *s;
    }

  *sa = r;
  *sb = z;

  return;
}
//****************************************************************************80

void dscal ( int n, double sa, double x[], int incx )

//****************************************************************************80
//
//  Purpose:
//
//    DSCAL scales a vector by a constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double SA, the multiplier.
//
//    Input/output, double X[*], the vector to be scaled.
//
//    Input, int INCX, the increment between successive entries of X.
//
{
  int i;
  int ix;
  int m;

  if ( n <= 0 )
    {
    }
  else if ( incx == 1 )
    {
      m = n % 5;

      for ( i = 0; i < m; i++ )
	{
	  x[i] = sa * x[i];
	}

      for ( i = m; i < n; i = i + 5 )
	{
	  x[i]   = sa * x[i];
	  x[i+1] = sa * x[i+1];
	  x[i+2] = sa * x[i+2];
	  x[i+3] = sa * x[i+3];
	  x[i+4] = sa * x[i+4];
	}
    }
  else
    {
      if ( 0 <= incx )
	{
	  ix = 0;
	}
      else
	{
	  ix = ( - n + 1 ) * incx;
	}

      for ( i = 0; i < n; i++ )
	{
	  x[ix] = sa * x[ix];
	  ix = ix + incx;
	}

    }

  return;
}
//****************************************************************************80

int dsvdc ( double a[], int lda, int m, int n, double s[], double e[], 
	    double u[], int ldu, double v[], int ldv, double work[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DSVDC computes the singular value decomposition of a real rectangular matrix.
//
//  Discussion:
//
//    This routine reduces an M by N matrix A to diagonal form by orthogonal
//    transformations U and V.  The diagonal elements S(I) are the singular
//    values of A.  The columns of U are the corresponding left singular
//    vectors, and the columns of V the right singular vectors.
//
//    The form of the singular value decomposition is then
//
//      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 May 2007
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the M by N matrix whose
//    singular value decomposition is to be computed.  On output, the matrix
//    has been destroyed.  Depending on the user's requests, the matrix may 
//    contain other useful information.
//
//    Input, int LDA, the leading dimension of the array A.
//    LDA must be at least M.
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, int N, the number of columns of the matrix A.
//
//    Output, double S[MM], where MM = min(M+1,N).  The first
//    min(M,N) entries of S contain the singular values of A arranged in
//    descending order of magnitude.
//
//    Output, double E[MM], where MM = min(M+1,N), ordinarily contains zeros.
//    However see the discussion of INFO for exceptions.
//
//    Output, double U[LDU*K].  If JOBA = 1 then K = M;
//    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of left singular
//    vectors.  U is not referenced if JOBA = 0.  If M <= N or if JOBA = 2, then
//    U may be identified with A in the subroutine call.
//
//    Input, int LDU, the leading dimension of the array U.
//    LDU must be at least M.
//
//    Output, double V[LDV*N], the N by N matrix of right singular vectors.
//    V is not referenced if JOB is 0.  If N <= M, then V may be identified
//    with A in the subroutine call.
//
//    Input, int LDV, the leading dimension of the array V.
//    LDV must be at least N.
//
//    Workspace, double WORK[M].
//
//    Input, int JOB, controls the computation of the singular
//    vectors.  It has the decimal expansion AB with the following meaning:
//      A =  0, do not compute the left singular vectors.
//      A =  1, return the M left singular vectors in U.
//      A >= 2, return the first min(M,N) singular vectors in U.
//      B =  0, do not compute the right singular vectors.
//      B =  1, return the right singular vectors in V.
//
//    Output, int *DSVDC, status indicator INFO.
//    The singular values (and their corresponding singular vectors)
//    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
//    Thus if *INFO is 0, all the singular values and their vectors are
//    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
//    matrix with the elements of S on its diagonal and the elements of E on
//    its superdiagonal.  Thus the singular values of A and B are the same.
//
{
  double b;
  double c;
  double cs;
  double el;
  double emm1;
  double f;
  double g;
  int i;
  int info;
  int iter;
  int j;
  int jobu;
  int k;
  int kase;
  int kk;
  int l;
  int ll;
  int lls;
  int ls;
  int lu;
  int maxit = 30;
  int mm;
  int mm1;
  int mn;
  int mp1;
  int nct;
  int nctp1;
  int ncu;
  int nrt;
  int nrtp1;
  double scale;
  double shift;
  double sl;
  double sm;
  double smm1;
  double sn;
  double t;
  double t1;
  double test;
  bool wantu;
  bool wantv;
  double ztest;
  //
  //  Determine what is to be computed.
  //
  info = 0;
  wantu = false;
  wantv = false;
  jobu = ( job % 100 ) / 10;

  if ( 1 < jobu )
    {
      ncu = std::min ( m, n );
    }
  else
    {
      ncu = m;
    }

  if ( jobu != 0 )
    {
      wantu = true;
    }

  if ( ( job % 10 ) != 0 )
    {
      wantv = true;
    }
  //
  //  Reduce A to bidiagonal form, storing the diagonal elements
  //  in S and the super-diagonal elements in E.
  //
  nct = std::min ( m-1, n );
  nrt = std::max ( 0, std::min ( m, n-2 ) );
  lu = std::max ( nct, nrt );

  for ( l = 1; l <= lu; l++ )
    {
      //
      //  Compute the transformation for the L-th column and
      //  place the L-th diagonal in S(L).
      //
      if ( l <= nct )
	{
	  s[l-1] = dnrm2 ( m-l+1, a+l-1+(l-1)*lda, 1 );

	  if ( s[l-1] != 0.0 )
	    {
	      if ( a[l-1+(l-1)*lda] != 0.0 )
		{
		  s[l-1] = copysign(1.0, a[l-1+(l-1)*lda] ) * fabs ( s[l-1] );
		}
	      dscal ( m-l+1, 1.0 / s[l-1], a+l-1+(l-1)*lda, 1 );
	      a[l-1+(l-1)*lda] = 1.0 + a[l-1+(l-1)*lda];
	    }
	  s[l-1] = -s[l-1];
	}

      for ( j = l+1; j <= n; j++ )
	{
	  //
	  //  Apply the transformation.
	  //
	  if ( l <= nct && s[l-1] != 0.0 )
	    {
	      t = - ddot ( m-l+1, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 ) 
		/ a[l-1+(l-1)*lda];
	      daxpy ( m-l+1, t, a+l-1+(l-1)*lda, 1, a+l-1+(j-1)*lda, 1 );
	    }
	  //
	  //  Place the L-th row of A into E for the
	  //  subsequent calculation of the row transformation.
	  //
	  e[j-1] = a[l-1+(j-1)*lda];
	}
      //
      //  Place the transformation in U for subsequent back multiplication.
      //
      if ( wantu && l <= nct )
	{
	  for ( i = l; i <= m; i++ )
	    {
	      u[i-1+(l-1)*ldu] = a[i-1+(l-1)*lda];
	    }
	}

      if ( l <= nrt )
	{
	  //
	  //  Compute the L-th row transformation and place the
	  //  L-th superdiagonal in E(L).
	  //
	  e[l-1] = dnrm2 ( n-l, e+l, 1 );

	  if ( e[l-1] != 0.0 )
	    {
	      if ( e[l] != 0.0 )
		{
		  e[l-1] = copysign(1.0, e[l] ) * fabs ( e[l-1] );
		}
	      dscal ( n-l, 1.0 / e[l-1], e+l, 1 );
	      e[l] = 1.0 + e[l];
	    }

	  e[l-1] = -e[l-1];
	  //
	  //  Apply the transformation.
	  //
	  if ( l+1 <= m && e[l-1] != 0.0 )
	    {
	      for ( j = l+1; j <= m; j++ )
		{
		  work[j-1] = 0.0;
		}

	      for ( j = l+1; j <= n; j++ )
		{
		  daxpy ( m-l, e[j-1], a+l+(j-1)*lda, 1, work+l, 1 );
		}

	      for ( j = l+1; j <= n; j++ )
		{
		  daxpy ( m-l, -e[j-1]/e[l], work+l, 1, a+l+(j-1)*lda, 1 );
		}
	    }
	  //
	  //  Place the transformation in V for subsequent back multiplication.
	  //
	  if ( wantv )
	    {
	      for ( j = l+1; j <= n; j++ )
		{
		  v[j-1+(l-1)*ldv] = e[j-1];
		}
	    }
	}
    }
  //
  //  Set up the final bidiagonal matrix of order MN.
  //
  mn = std::min ( m + 1, n );
  nctp1 = nct + 1;
  nrtp1 = nrt + 1;

  if ( nct < n )
    {
      s[nctp1-1] = a[nctp1-1+(nctp1-1)*lda];
    }

  if ( m < mn )
    {
      s[mn-1] = 0.0;
    }

  if ( nrtp1 < mn )
    {
      e[nrtp1-1] = a[nrtp1-1+(mn-1)*lda];
    }

  e[mn-1] = 0.0;
  //
  //  If required, generate U.
  //
  if ( wantu )
    {
      for ( i = 1; i <= m; i++ )
	{
	  for ( j = nctp1; j <= ncu; j++ )
	    {
	      u[(i-1)+(j-1)*ldu] = 0.0;
	    }
	}

      for ( j = nctp1; j <= ncu; j++ )
	{
	  u[j-1+(j-1)*ldu] = 1.0;
	}

      for ( ll = 1; ll <= nct; ll++ )
	{
	  l = nct - ll + 1;

	  if ( s[l-1] != 0.0 )
	    {
	      for ( j = l+1; j <= ncu; j++ )
		{
		  t = - ddot ( m-l+1, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 ) 
		    / u[l-1+(l-1)*ldu];
		  daxpy ( m-l+1, t, u+(l-1)+(l-1)*ldu, 1, u+(l-1)+(j-1)*ldu, 1 );
		}

	      dscal ( m-l+1, -1.0, u+(l-1)+(l-1)*ldu, 1 );
	      u[l-1+(l-1)*ldu] = 1.0 + u[l-1+(l-1)*ldu];
	      for ( i = 1; i <= l-1; i++ )
		{
		  u[i-1+(l-1)*ldu] = 0.0;
		}
	    }
	  else
	    {
	      for ( i = 1; i <= m; i++ )
		{
		  u[i-1+(l-1)*ldu] = 0.0;
		}
	      u[l-1+(l-1)*ldu] = 1.0;
	    }
	}
    }
  //
  //  If it is required, generate V.
  //
  if ( wantv )
    {
      for ( ll = 1; ll <= n; ll++ )
	{
	  l = n - ll + 1;

	  if ( l <= nrt && e[l-1] != 0.0 )
	    {
	      for ( j = l+1; j <= n; j++ )
		{
		  t = - ddot ( n-l, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 ) 
		    / v[l+(l-1)*ldv];
		  daxpy ( n-l, t, v+l+(l-1)*ldv, 1, v+l+(j-1)*ldv, 1 );
		}

	    }
	  for ( i = 1; i <= n; i++ )
	    {
	      v[i-1+(l-1)*ldv] = 0.0;
	    }
	  v[l-1+(l-1)*ldv] = 1.0;
	}
    }
  //
  //  Main iteration loop for the singular values.
  //
  mm = mn;
  iter = 0;

  while ( 0 < mn )
    {
      //
      //  If too many iterations have been performed, set flag and return.
      //
      if ( maxit <= iter )
	{
	  info = mn;
	  return info;
	}
      //
      //  This section of the program inspects for
      //  negligible elements in the S and E arrays.
      //
      //  On completion the variables KASE and L are set as follows:
      //
      //  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
      //  KASE = 2     if S(L) is negligible and L < MN
      //  KASE = 3     if E(L-1) is negligible, L < MN, and
      //               S(L), ..., S(MN) are not negligible (QR step).
      //  KASE = 4     if E(MN-1) is negligible (convergence).
      //
      for ( ll = 1; ll <= mn; ll++ )
	{
	  l = mn - ll;

	  if ( l == 0 )
	    {
	      break;
	    }

	  test = fabs ( s[l-1] ) + fabs ( s[l] );
	  ztest = test + fabs ( e[l-1] );

	  if ( ztest == test )
	    {
	      e[l-1] = 0.0;
	      break;
	    }
	}

      if ( l == mn - 1 )
	{
	  kase = 4;
	}
      else
	{
	  mp1 = mn + 1;

	  for ( lls = l+1; lls <= mn+1; lls++ )
	    {
	      ls = mn - lls + l + 1;

	      if ( ls == l )
		{
		  break;
		}

	      test = 0.0;
	      if ( ls != mn )
		{
		  test = test + fabs ( e[ls-1] );
		}

	      if ( ls != l + 1 )
		{
		  test = test + fabs ( e[ls-2] );
		}

	      ztest = test + fabs ( s[ls-1] );

	      if ( ztest == test )
		{
		  s[ls-1] = 0.0;
		  break;
		}

	    }

	  if ( ls == l )
	    {
	      kase = 3;
	    }
	  else if ( ls == mn )
	    {
	      kase = 1;
	    }
	  else
	    {
	      kase = 2;
	      l = ls;
	    }
	}

      l = l + 1;
      //
      //  Deflate negligible S(MN).
      //
      if ( kase == 1 )
	{
	  mm1 = mn - 1;
	  f = e[mn-2];
	  e[mn-2] = 0.0;

	  for ( kk = 1; kk <= mm1; kk++ )
	    {
	      k = mm1 - kk + l;
	      t1 = s[k-1];
	      drotg ( &t1, &f, &cs, &sn );
	      s[k-1] = t1;

	      if ( k != l )
		{
		  f = -sn * e[k-2];
		  e[k-2] = cs * e[k-2];
		}

	      if ( wantv )
		{
		  drot ( n, v+0+(k-1)*ldv, 1, v+0+(mn-1)*ldv, 1, cs, sn );
		}
	    }
	}
      //
      //  Split at negligible S(L).
      //
      else if ( kase == 2 )
	{
	  f = e[l-2];
	  e[l-2] = 0.0;

	  for ( k = l; k <= mn; k++ )
	    {
	      t1 = s[k-1];
	      drotg ( &t1, &f, &cs, &sn );
	      s[k-1] = t1;
	      f = - sn * e[k-1];
	      e[k-1] = cs * e[k-1];
	      if ( wantu )
		{
		  drot ( m, u+0+(k-1)*ldu, 1, u+0+(l-2)*ldu, 1, cs, sn );
		}
	    }
	}
      //
      //  Perform one QR step.
      //
      else if ( kase == 3 )
	{
	  //
	  //  Calculate the shift.
	  //
	  scale = std::max ( fabs ( s[mn-1] ), 
			   std::max ( fabs ( s[mn-2] ), 
				    std::max ( fabs ( e[mn-2] ), 
					     std::max ( fabs ( s[l-1] ), fabs ( e[l-1] ) ) ) ) );

	  sm = s[mn-1] / scale;
	  smm1 = s[mn-2] / scale;
	  emm1 = e[mn-2] / scale;
	  sl = s[l-1] / scale;
	  el = e[l-1] / scale;
	  b = ( ( smm1 + sm ) * ( smm1 - sm ) + emm1 * emm1 ) / 2.0;
	  c = ( sm * emm1 ) * ( sm * emm1 );
	  shift = 0.0;

	  if ( b != 0.0 || c != 0.0 )
	    {
	      shift = sqrt ( b * b + c );
	      if ( b < 0.0 )
		{
		  shift = -shift;
		}
	      shift = c / ( b + shift );
	    }

	  f = ( sl + sm ) * ( sl - sm ) - shift;
	  g = sl * el;
	  //
	  //  Chase zeros.
	  //
	  mm1 = mn - 1;

	  for ( k = l; k <= mm1; k++ )
	    {
	      drotg ( &f, &g, &cs, &sn );

	      if ( k != l )
		{
		  e[k-2] = f;
		}

	      f = cs * s[k-1] + sn * e[k-1];
	      e[k-1] = cs * e[k-1] - sn * s[k-1];
	      g = sn * s[k];
	      s[k] = cs * s[k];

	      if ( wantv )
		{
		  drot ( n, v+0+(k-1)*ldv, 1, v+0+k*ldv, 1, cs, sn );
		}

	      drotg ( &f, &g, &cs, &sn );
	      s[k-1] = f;
	      f = cs * e[k-1] + sn * s[k];
	      s[k] = -sn * e[k-1] + cs * s[k];
	      g = sn * e[k];
	      e[k] = cs * e[k];

	      if ( wantu && k < m )
		{
		  drot ( m, u+0+(k-1)*ldu, 1, u+0+k*ldu, 1, cs, sn );
		}
	    }
	  e[mn-2] = f;
	  iter = iter + 1;
	}
      //
      //  Convergence.
      //
      else if ( kase == 4 )
	{
	  //
	  //  Make the singular value nonnegative.
	  //
	  if ( s[l-1] < 0.0 )
	    {
	      s[l-1] = -s[l-1];
	      if ( wantv )
		{
		  dscal ( n, -1.0, v+0+(l-1)*ldv, 1 );
		}
	    }
	  //
	  //  Order the singular value.
	  //
	  for ( ; ; )
	    {
	      if ( l == mm )
		{
		  break;
		}

	      if ( s[l] <= s[l-1] )
		{
		  break;
		}

	      t = s[l-1];
	      s[l-1] = s[l];
	      s[l] = t;

	      if ( wantv && l < n )
		{
		  dswap ( n, v+0+(l-1)*ldv, 1, v+0+l*ldv, 1 );
		}

	      if ( wantu && l < m )
		{
		  dswap ( m, u+0+(l-1)*ldu, 1, u+0+l*ldu, 1 );
		}

	      l = l + 1;
	    }
	  iter = 0;
	  mn = mn - 1;
	}
    }

  return info;
}
//****************************************************************************80

void dswap ( int n, double x[], int incx, double y[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DSWAP interchanges two vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979,
//    ISBN13: 978-0-898711-72-1,
//    LC: QA214.L56.
//
//    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input/output, double X[*], one of the vectors to swap.
//
//    Input, int INCX, the increment between successive entries of X.
//
//    Input/output, double Y[*], one of the vectors to swap.
//
//    Input, int INCY, the increment between successive elements of Y.
//
{

  int i;
  int ix;
  int iy;
  int m;
  double temp;

  if ( n <= 0 )
    {
    }
  else if ( incx == 1 && incy == 1 )
    {
      m = n % 3;

      for ( i = 0; i < m; i++ )
	{
	  temp = x[i];
	  x[i] = y[i];
	  y[i] = temp;
	}

      for ( i = m; i < n; i = i + 3 )
	{
	  temp = x[i];
	  x[i] = y[i];
	  y[i] = temp;

	  temp = x[i+1];
	  x[i+1] = y[i+1];
	  y[i+1] = temp;

	  temp = x[i+2];
	  x[i+2] = y[i+2];
	  y[i+2] = temp;
	}
    }
  else
    {
      if ( 0 <= incx )
	{
	  ix = 0;
	}
      else
	{
	  ix = ( - n + 1 ) * incx;
	}

      if ( 0 <= incy )
	{
	  iy = 0;
	}
      else
	{
	  iy = ( - n + 1 ) * incy;
	}

      for ( i = 0; i < n; i++ )
	{
	  temp = x[ix];
	  x[ix] = y[iy];
	  y[iy] = temp;
	  ix = ix + incx;
	  iy = iy + incy;
	}

    }

  return;
}
//****************************************************************************80


void svd_solve ( int m, int n, double *&a, double  *&b, double *& x )

//****************************************************************************80
//
//  Purpose:
//
//    SVD_SOLVE solves a linear system in the least squares sense.
//
//  Discussion:
//
//    The vector X returned by this routine should always minimize the 
//    Euclidean norm of the residual ||A*x-b||.
//
//    If the matrix A does not have full column rank, then there are multiple
//    vectors that attain the minimum residual.  In that case, the vector
//    X returned by this routine is the unique such minimizer that has the 
//    the minimum possible Euclidean norm, that is, ||A*x-b|| and ||x||
//    are both minimized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input, double A[M*N], the matrix.
//
//    Input, double B[M], the right hand side.
//
//    Output, double SVD_SOLVE[N], the least squares solution.
//
{
  
  double *a_copy=NULL;
  double *e=NULL;
  int i=0;
  int info=0;
  int lda=0;
  int ldu=0;
  int ldv=0;
  int job=0;
  double *sdiag=NULL;
  double smax=0;
  double stol=0;
  double *sub=NULL;
  double *u=NULL;
  double *ub=NULL;
  double *v=NULL;
  double *x_temp=NULL;
  
  double *work=NULL;
  //
  //  Get the SVD.
  //
  a_copy = new double[m*n];
  std::copy(a, a+ m*n, a_copy);
  lda = m;
  sdiag = new double[std::max ( m + 1, n )];
  e = new double[std::max ( m + 1, n )];
  u = new double[m*m];
  ldu = m;
  v = new double[n*n];
  ldv = n;
  work = new double[m];
  x_temp = new double[n];
  job = 11;

  info = dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, job );

  delete [] work;
  

  if ( info != 0 )
    {
      cerr << "\n";
      cerr << "SVD_SOLVE - Failure!\n";
      cerr << "  The SVD could not be calculated.\n";
      cerr << "  LINPACK routine DSVDC returned a nonzero\n";
      cerr << "  value of the error flag, INFO = " << info << "\n";
      exit ( 1 );
    }

  ub = mult_mtv ( m, m, u, b );
  
  //
  //  For singular problems, there may be tiny but nonzero singular values
  //  that should be ignored.  This is a reasonable attempt to avoid such 
  //  problems, although in general, the user might wish to control the tolerance.
  //
  smax = *std::max_element(sdiag, sdiag+sizeof(sdiag));
  if ( smax <= epsilon )
    {
      smax = 1.0;
    }

  stol = epsilon * smax;
  
  sub = new double[n];
  
  for ( i = 0; i < n; i++ )
    {
      sub[i] = 0.0;
      if ( i < m )
	{
	  // if ( stol <= sdiag[i])
	  // {
	      sub[i] = ub[i] / sdiag[i];
	      // }
	}
    }

  x_temp = mult_mv ( n, n, v, sub );
  for ( i = 0; i < n; i++ )
    {
      x[i]=x_temp[i]; 
    }
 
  
  delete [] sdiag;
  delete [] sub;
  delete [] u;
  delete [] ub;
  delete [] v;
  delete [] a_copy;
  delete [] e;

 
}

double vec_norm2 ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_NORM2 returns the L2 norm of an R8VEC.
//
//  Discussion:
//
//    Original function name : r8vec_norm (only name modified)
//
//    The vector L2 norm is defined as:
//
//      VEC_NORM2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 November 2014
//
//  Original Author:
//
//    John Burkardt
//
//  Modified by :
//
//     Laura S. Mendoza
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], the vector whose L2 norm is desired.
//
//    Output, double VEC_NORM2, the L2 norm of A.
//
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}

double vec_norm2_aff ( int n, double v0[], double v1[] )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_NORM2_AFF returns the affine L2 norm of an VEC.
//
//  Discussion:
//
//    Original Name : R8VEC_NORM_AFFINE (only name modified)
//
//    The affine vector L2 norm is defined as:
//
//      VEC_NORM2_AFF(V0,V1)
//        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Modified by :
//
//     Laura S. Mendoza
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, double V0[N], the base vector.
//
//    Input, double V1[N], the vector.
//
//    Output, double VEC_NORM2_AFF, the affine L2 norm.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80


double *mult_mv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULT_MV multiplies a matrix times a vector.
//
//  Discussion:
//
//    Original function name : r8MAT_MV_NEW (only named modified)
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 November 2014
//
//  Author:
//
//    John Burkardt
//
//  Modifier:
//
//    Laura S. Mendoza
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double MULT_MV[M], the product A*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
//****************************************************************************80


double *mult_mtv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MULT_MTV multiplies a transposed matrix times a vector.
//
//  Discussion:
//
//    Original function name: R8MAT_MTV_NEW
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2014
//
//  Modified by :
//
//    Laura
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[M], the vector to be multiplied by A.
//
//    Output, double R8MAT_MTV_NEW[N], the product A'*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[n];

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  return y;
}
//****************************************************************************80


void time_stamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
