#include "Include/lu.h"
#include <stdio.h>
#include <stdlib.h>
# include <cs_demo.h>

bool equal( double a, double b )
{
  if ( Abs(a-b) < 0.001  )
  {
    return true;
  }
  else
  {
    return false;
  }
}
double Abs(double x)
{
  return x < 0 ? -x : x;
} 

int main( int argc, char *argv[] )
{
     /* Read the file1 and store  it by CSC: <a, asub, xa> */
    FILE *file1, *file2, *file3;
    FILE *file4, *file5;
    int *pinv, *q;
    int m, n, i, j = 0;
    int l_i, l_j, sum = 0, u_i, u_j, a_i, a_j;
    int nzl, nzu;
    int *Ai, *Ap, *Li, *Ui, *Lp, *Up;
    int *p, *p_u, *p_a;
    double *Ax, a_x;
    double nnz;
    double start, end, start_GP, end_GP, start_GP_sn, end_GP_sn, start_GP_sn_u, end_GP_sn_u;
    int *sn_record, *sn_record_u;
    file1 = fopen(argv[1], "r"); 
    file2 = fopen(argv[2], "r");
    file3 = fopen(argv[3], "r");
    file4 = fopen(argv[4], "r");
    file5 = fopen(argv[5], "r");
    nzl = atoi(argv[6]);
    nzu = atoi(argv[7]);
    
    printf("nzl = %d nzu = %d\n", nzl, nzu );
    fscanf ( file1, "%d %d %lf\n", &m, &n, &nnz);

    Ap = ( int * )malloc( sizeof( int ) * n+1 ) ;
    Ai = ( int * )malloc( sizeof( int ) * nnz ) ;
    Ax = ( double * )malloc( sizeof( double ) * nnz ) ;
    p_a = ( int * )malloc( sizeof( int ) * n+1 ) ;
    Lp = ( int * )malloc( sizeof( int ) * n+1 ) ;
    Up = ( int * )malloc( sizeof( int ) * n+1 ) ;
    p = ( int * )malloc( sizeof( int ) * n ) ;
    p_u = ( int * )malloc( sizeof( int ) * n ) ;
    Li = ( int * )malloc( sizeof( int ) * nzl ) ;
    Ui = ( int * )malloc( sizeof( int ) * nzu ) ;
    memset(Lp, 0, sizeof(int)*(n+1));
    memset(Up, 0, sizeof(int)*(n+1));
    memset(p, 0, sizeof(int)*n);
    memset(p_u, 0, sizeof(int)*n);
    memset(Ap, 0, sizeof(int)*(n+1));
    memset(p_a, 0, sizeof(int)*n);

    while ( fscanf ( file1, "%d %d %lf\n", &a_i, &a_j, &a_x) == 3 )
    {
         Ai[sum] = a_i-1; 
        Ax[sum] = a_x;
        p_a[a_j-1]++;
        sum++ ;
    }
    for ( i = 1; i <=n; i++ )
    {
        Ap[i] = Ap[i-1]+p_a[i-1];
    }
    fclose(file1);

    sum = 0;
    int sum_l_greater_4 = 0;
    while ( fscanf ( file2, "%d %d\n", &l_i, &l_j) == 2 )
    {
        Li[sum] = l_i; 
        p[l_j]++;
          sum++;
    }
    for ( i = 1; i <=n; i++ )
    {
        Lp[i] = Lp[i-1]+p[i-1];
    }
    fclose(file2);

    sum = 0;
    while ( fscanf ( file3, "%d %d\n", &u_i, &u_j) == 2 )
    {
          Ui[sum] = u_i;
        p_u[u_j]++;
         sum++ ;
    }
    for ( i = 1; i <=n; i++ )
    {
        Up[i] = Up[i-1]+p_u[i-1];
    }
    fclose(file3);

    pinv = ( int * )malloc( sizeof( int ) * n ) ;
    while ( fscanf(file4, "%d\n", &i) == 1 )
    {
         pinv[j] = i;  
         j++; 
    }

    q = ( int * )malloc( sizeof( int ) * n ) ;
    j = 0;
    while ( fscanf(file5, "%d\n", &i) == 1 )
    {
          q[j] = i;  
        j++; 
    }
    
     double *l_gp, *l_gp_sn, *x, *l_gp_sn_u;
     double sum_error = 0;
    
    start_GP = microtime();
    l_gp = lu_gp_sparse(Ax, Ai, Ap, n, nzl, nzu, q, pinv, Li, Lp, Ui, Up);
    end_GP = microtime() - start_GP;
    printf("Time of LU_GP: %lf\n", end_GP);
 
    start_GP_sn = microtime();
    l_gp_sn = lu_gp_sparse_avx2(Ax, Ai, Ap, n, nzl, nzu, q, pinv, Li, Lp, Ui, Up);
    end_GP_sn = microtime() - start_GP_sn;
    printf("Time of LU_GP_avx2: %lf\n", end_GP_sn);
     
    for ( i = 0; i < n; i++ )
    {
        if ( !equal ( l_gp[i], l_gp_sn[i] ) ) 
        {
            sum_error++;
        }
    }
    printf("error results are: %lf\n", sum_error);
    
    return 0;
}







