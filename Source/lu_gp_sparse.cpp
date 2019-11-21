# include <stdio.h>
# include <stdlib.h>
# include "../Include/lu.h"
# include <sys/time.h>
# include <float.h>
# define MICRO_IN_SEC 1000000.00

/* Time Stamp */
double microtime()
{
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv,&tz);

        return tv.tv_sec+tv.tv_usec/MICRO_IN_SEC;
}

double* lu_gp_sparse(double *a, int *asub, int *xa, int n, int nzl, int nzu, int *perm_c, int *perm_r, int *asub_L, int *xa_L, int *asub_U, int *xa_U)
{
  int sum_l = 0, sum_u = 0, row, sum_for = 0;
  double *L, *U, *xx;
  int *row_index, row_column;
  double U_diag;
  double start, end;
  int j, k, current_column;
  int i;
  L = ( double * )_mm_malloc(sizeof(double) * nzl, 64);
  U = ( double * )_mm_malloc(sizeof(double) * nzu, 64);
  xx = ( double *)_mm_malloc(sizeof(double) * n, 64 );
  
  /* Array xx initialization*/
  for ( i = 0; i < n; i++ )
  {
     xx[i] = 0;
    L[xa_L[i]] = 1.0;
  }

  /* column-oriented G/P algorithm without partial pivoting */
  for ( k = 0; k < n; k++ )
  {
     current_column = perm_c[k];
    
    /* xx[] = A[:,current_column] */
    for ( j = xa[current_column]; j < xa[current_column+1]; j++ )
    {
      xx[perm_r[asub[j]]] = a[j];
    }
    row_column = xa_U[k+1] - xa_U[k] - 1;
    row_index = (int *)malloc( sizeof(int) * row_column);
    for ( j = 0; j < row_column; j++ )
    {
        row_index[j] = asub_U[j+xa_U[k]];
    }
    
    /* L[:,0~k]*xx = A[:,current_column], solve for xx*/
    for ( j = 0; j < row_column; j++ )
    {
        row = row_index[j];
        for ( i = xa_L[row]+1; i < xa_L[row+1]; i++ )
        {
            xx[asub_L[i]] -=  xx[row]*L[i];
        }
    }

    /* solve for U[:,k]*/
    for ( i = xa_U[k]; i < xa_U[k+1]; i++ )
    {
      U[i] = xx[asub_U[i]];
      xx[asub_U[i]] = 0;
    }

    /* solve for L[:,k] */
    U_diag = U[i-1];
    for ( i = xa_L[k]+1; i < xa_L[k+1]; i++ )
    {
      L[i] = xx[asub_L[i]] / U_diag;
      xx[asub_L[i]] = 0;
    } 
  }

   /* solve for Ly = b and Ux = y */
   double *y, *x;
   y = ( double *)malloc( sizeof( double ) * n );
   x = ( double *)malloc( sizeof( double ) * n );


  for ( i = 0; i < n; i++ )
  {
    y[i] = 1.0;
  }

  for ( i = 0; i < n; i++ )
  {
     for ( j = xa_L[i]+1; j < xa_L[i+1]; j++ )
    {
        y[asub_L[j]] -= y[i] * L[j];
    }
  }

  //x[n-1] = y[n-1];
  for ( i = 0; i < n; i++ )
  {
    x[i] = y[i];
  }
  x[n-1] = y[n-1]/U[nzu-1];
  for ( i = n-1; i > 0; i-- )
  {
    for ( j = xa_U[i]; j < xa_U[i+1]-1; j++ )
    {
        x[asub_U[j]] -= x[i] *U[j];
    }
    x[i-1] = x[i-1]/U[xa_U[i]-1];
  }
  
  return x;
}

