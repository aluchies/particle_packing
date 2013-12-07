////////////////////////////////////////////////////////////////////////////////
//  void matrix_multiply(double *C, double *A, double *B, int nrows,          //
//                                                  int ncols, int mcols)     //
//                                                                            //
//  Description:                                                              //
//     Post multiply the nrows x ncols matrix A by the ncols x mcols matrix   //
//     B to form the nrows x mcols matrix C, i.e. C = A B.                    //
//     The matrix C should be declared as double C[nrows][mcols] in the       //
//     calling routine.  The memory allocated to C should not include any     //
//     memory allocated to A or B.                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *C    Pointer to the first element of the matrix C.             //
//     double *A    Pointer to the first element of the matrix A.             //
//     double *B    Pointer to the first element of the matrix B.             //
//     int    nrows The number of rows of the matrices A and C.               //
//     int    ncols The number of columns of the matrices A and the           //
//                   number of rows of the matrix B.                          //
//     int    mcols The number of columns of the matrices B and C.            //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     #define NB                                                             //
//     double A[M][N],  B[N][NB], C[M][NB];                                   //
//                                                                            //
//     (your code to initialize the matrices A and B)                         //
//                                                                            //
//     matrix_multiply(&C[0][0], &A[0][0], &B[0][0], M, N, NB);               //
//     printf("The matrix C is \n"); ...                                      //
//                                                                            //
//  Source: http://www.mymathlib.com/matrices/arithmetic/mul_mat.html         //
//  Last visisted on: 2013-12-05                                              //
////////////////////////////////////////////////////////////////////////////////
void matrix_multiply(double *C, double *A, double *B, int nrows, int ncols, int mcols) 
{
   double *pB;
   double *p_B;
   int i,j,k;

   for (i = 0; i < nrows; A += ncols, i++) 
      for (p_B = B, j = 0; j < mcols; C++, p_B++, j++) {
         pB = p_B;
         *C = 0.0; 
         for (k = 0; k < ncols; pB += mcols, k++) 
            *C += *(A+k) * *pB;
      }
}



////////////////////////////////////////////////////////////////////////////////
//  void matrix_transpose(double *At, double *A, int nrows, int ncols)        //
//                                                                            //
//  Description:                                                              //
//     Take the transpose of A and store in At, i.e. At = A'.                 //
//     The matrix At should be declared as double At[ncols][nrows] in the     //
//     calling routine, and the matrix A declared as double A[nrows[ncols].   //
//     In general, At and A should be disjoint i.e. their memory locations    //
//     should be distinct.                                                    //
//                                                                            //
//  Arguments:                                                                //
//     double *At   Pointer to the first element of the matrix At.            //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of the matrix A and number of columns  //
//                  of the matrix At.                                         //
//     int    ncols The number of columns of the matrix A and the number of   //
//                  rows of the matrix At.                                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  At[N][M];                                             //
//                                                                            //
//     (your code to initialize the matrix A)                                 //
//                                                                            //
//     matrix_transpose(&At[0][0], &A[0][0], M, N);                           //
//     printf("The transpose of A is the matrix At \n"); ...                  //
//                                                                            //
//  http://www.mymathlib.com/c_source/matrices/datamovement/transpose_matrix.c//
//  Last visisted on: 2013-12-05                                              //
////////////////////////////////////////////////////////////////////////////////
void matrix_transpose(double *At, double *A, int nrows, int ncols) 
{
   double *pA;
   double *pAt;
   int i,j;

   for (i = 0; i < nrows; At += 1, A += ncols, i++) {
      pAt = At;
      pA = A;
      for (j = 0; j < ncols; pAt += nrows, j++) *pAt = pA[j];
   }
}





////////////////////////////////////////////////////////////////////////////////
//  void set_diagonal(double *A, double v[], int nrows, int ncols)            //
//                                                                            //
//  Description:                                                              //
//     Copy the vector v to the diagonal of the matrix A, i.e.                //
//     A[i][i] = v[i], where 0 <= i < min( nrows, ncols ).                    //
//     Note that v should be declared "double v[N]", N >= min( nrows, ncols ) //
//     in the calling routine.                                                //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the destination matrix A. //
//     double v[]   Source of the new diagonal for the matrix A.              //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to initialize the matrix A and the vector v)                //
//                                                                            //
//     set_diagonal(&A[0][0], v, M, N);                                       //
//     printf("The matrix A is \n"); ...                                      //
//                                                                            //
//  http://www.mymathlib.com/c_source/matrices/datamovement/set_diagonal.c    //
//  Last visisted on: 2013-12-05                                              //
////////////////////////////////////////////////////////////////////////////////
void set_diagonal(double *A, double v[], int nrows, int ncols)
{
   int n;

   if (nrows < ncols) n = nrows; else n = ncols;

   for (; n > 0 ; A += (ncols + 1), n--)  *A = *v++;
}