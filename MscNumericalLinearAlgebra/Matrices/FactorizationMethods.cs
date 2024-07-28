using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.Matrices
{
    public static class FactorizationMethods
    {
        /// <summary>
        /// Performs LU decomposition in-place on a given square matrix, storing the results in the original matrix.
        /// </summary>
        /// <param name="matrixA">The input square matrix to be decomposed.</param>
        /// <returns>The input matrix, transformed into the combined form of lower and upper triangular matrices (LU decomposition).</returns>
        public static double[,] FactorizationLUStoredInMatrixA(double[,] matrixA)
        {
            int n = matrixA.GetLength(1);  // Length of Columns (n)
            int m = matrixA.GetLength(0);  // Length of Rows (m)

            if (n != m)  // Check if the matrix is square
            {
                throw new Exception("The matrix should be square!");
            }

            for (int i = 0; i < n; i++)             // i = 0, 1, 2,..., n
            {
                for (int j = i + 1; j < n; j++)     // j = 1, 2, 3,..., n 
                {
                    if (matrixA[i, i] == 0)  // If the element is 0, do not divide
                    {
                        throw new Exception("Cannot divide with 0, Singular Matrix!");
                    }

                    double factor = matrixA[j, i] / matrixA[i, i];  // i.e. A10 / A00

                    for (int k = i; k < n; k++)  // Run all the row elements to subtract the quantity (factor * above element)   
                    {
                        matrixA[j, k] = matrixA[j, k] - factor * matrixA[i, k];  // i.e. A10 = A10 - factor * A00
                    }

                    matrixA[j, i] = factor;  // Assign the Aji as the factor of the L Matrix
                }
            }

/*            Console.WriteLine("-----------------------");
            Console.WriteLine("A into LU = ");
            MatrixMethods.PrintMatrix(matrixA);
            Console.WriteLine("-----------------------");*/

            return matrixA;
        }

        /// <summary>
        /// Performs Cholesky Factorization on a given positive definite matrix to produce a lower triangular matrix L.
        /// </summary>
        /// <param name="matrixA">The input positive definite matrix to be decomposed.</param>
        /// <returns>The lower triangular matrix resulting from the Cholesky decomposition.</returns>
        public static double[,] FactorizationCholesky(double[,] matrixA)
        {
            int n = matrixA.GetLength(0);
            double[,] matrixL = new double[n, n]; // Initialize the lower triangular matrix L

            // Validate if the input matrix is square
            if (matrixA.GetLength(0) != matrixA.GetLength(1))
            {
                throw new Exception("The Matrix A must be a square matrix!");
            }

            for (int i = 0; i < n; i++)       // Iterate over rows (i-th Row)
            {
                for (int j = 0; j <= i; j++)  // Iterate over columns (j-th Column)
                {
                    double sum = 0.0;

                    // Calculate the sum for the elements of the lower triangular matrix
                    for (int k = 0; k < j; k++)  // k < j
                    {
                        sum = sum + matrixL[i, k] * matrixL[j, k];
                    }

                    if (i == j)
                    {
                        // Diagonal elements calculation: L[i,i] = sqrt(A[i,i] - Σ L[i,k]*L[i,k])
                        matrixL[i, j] = Math.Sqrt(matrixA[i, i] - sum);
                    }
                    else
                    {
                        // Off-diagonal elements calculation: L[i,j] = (1.0 / L[j,j]) * (A[i,j] - Σ L[i,k]*L[j,k])
                        matrixL[i, j] = (1.0 / matrixL[j, j]) * (matrixA[i, j] - sum);
                    }
                }
            }

            MatrixMethods.PrintMatrix(matrixL); // Print Matrix L transpose 

            return matrixL; // Return Matrix L
        }
    }
}
