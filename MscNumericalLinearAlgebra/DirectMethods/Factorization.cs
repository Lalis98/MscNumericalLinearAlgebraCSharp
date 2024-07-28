using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.DirectMethods
{
    public static class Factorization
    {
        public static (double[,] matrixL, double[,] matrixU) FactorizationLU(double[,] matrixA)
        {
            int n = matrixA.GetLength(1);  // Length of Columns (n)
            int m = matrixA.GetLength(0);  // Length of Rows (m)

            double[,] matrixU = matrixA;  // Create matrix U
            double[,] matrixL = new double[m, n];  // Create matrix L
            double factor = 0;  // Create factor variable

            if (n != m)  // Check if the matrix is square
            {
                throw new Exception("The matrix should be square!");
            }

            matrixL = MatrixMethods.CreateIdentityMatrix(n);  // Create Identity Matrix of size n

            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    if (matrixU[i, i] == 0)  // If the element is 0, do not divide
                    {
                        throw new Exception("Cannot divide with 0, Singular Matrix!");
                    }
                    factor = matrixU[j, i] / matrixU[i, i];  // i.e. A10 / A00
                    matrixL[j, i] = factor;  // Assign the Aji as the factor

                    for (int k = i; k < n; k++)  // Run all the row elements to subtract the quantity (factor * above element)   
                    {
                        matrixU[j, k] = matrixU[j, k] - factor * matrixU[i, k];  // i.e. A10 = A10 - factor * A00
                    }
                }
            }

            Console.WriteLine("L = ");
            MatrixMethods.PrintMatrix(matrixL);

            Console.WriteLine("U = ");
            MatrixMethods.PrintMatrix(matrixU);

            return (matrixL, matrixU);
        }


        public static double[,] FactorizationLUStoredInMatrixA(double[,] matrixA)
        {
            int n = matrixA.GetLength(1);  // Size of Matrix

            for (int i = 0; i < n; i++)
            {
                // Update upper triangular matrix U (row i)
                for (int j = i; j < n; j++)
                {
                    matrixA[i, j] = matrixA[i, j];  // U[i,j] = A[i,j]
                }

                // Update lower triangular matrix L (column i)
                for (int j = i + 1; j < n; j++)
                {
                    matrixA[j, i] = matrixA[j, i] / matrixA[i, i];  // L[j,i] = A[j,i] / U[i,i]

                    for (int k = i + 1; k < n; k++)
                    {
                        matrixA[j, k] = matrixA[j, k] - matrixA[j, i] * matrixA[i, k];  // U[j,k] = A[j,k] - L[j,i] * U[i,k]
                    }
                }
            }
            Console.WriteLine("A (into LU) = ");
            MatrixMethods.PrintMatrix(matrixA);

            return matrixA;
        }

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

            Console.WriteLine("L = ");
            MatrixMethods.PrintMatrix(matrixL);

            return matrixL;
        }

        public static double[,] FactorizationCholeskyStoredInMatrixA(double[,] matrixA)
        {
            if (matrixA.GetLength(0) != matrixA.GetLength(1))
            {
                throw new Exception("Matrix should be square!");
            }

            int n = matrixA.GetLength(0); // Size of A

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double sum = 0.0;

                    for (int k = 0; k < j; k++)  // Ensure k starts from 0
                    {
                        sum = sum + matrixA[i, k] * matrixA[j, k];
                    }

                    if (i == j)
                    {
                        matrixA[i, j] = Math.Sqrt(matrixA[i, i] - sum);
                    }
                    else
                    {
                        matrixA[i, j] = (1.0 / matrixA[j, j]) * (matrixA[i, j] - sum);
                        matrixA[j, i] = 0.0; // Simultaneously assign zeros into L^T the upper triangular part
                    }
                }
            }

            Console.WriteLine("A (into Cholesky) =");
            MatrixMethods.PrintMatrix(matrixA);

            return matrixA;
        }

        public static (double[], int[]) FactorizationCholeskySkyline(double[] valuesMatrixA, int[] diagonalOffsetsMatrixA)
        {

            int n = diagonalOffsetsMatrixA.Length - 1;  // Size of Positive Definite Matrix (nxn)

            double[] valuesMatrixL = new double[valuesMatrixA.Length];  // Values of L
            int[] diagonalOffsetsMatrixL = diagonalOffsetsMatrixA;  // Diagonal Offsets of L

            for (int i = 0; i < n; i++)
            {
                int hi = diagonalOffsetsMatrixA[i + 1] - diagonalOffsetsMatrixA[i] - 1;  // Height of Column i
                int mi = i - hi;  // Free Height of Column i

                for (int j = mi; j < i; j++)
                {

                    int hj = diagonalOffsetsMatrixA[j + 1] - diagonalOffsetsMatrixA[j] - 1;  // Height of Column j
                    int mj = j - hj;   // Free Height of Column j

                    double sum1 = 0.0;

                    // Calculate the sum for the elements of the lower triangular matrix
                    for (int k = Math.Max(mi, mj); k < j; k++)  // k < j
                    {
                        int L_ik_index = diagonalOffsetsMatrixA[i] + (i - k);
                        int L_jk_index = diagonalOffsetsMatrixA[j] + (j - k);

                        sum1 = sum1 + valuesMatrixL[L_ik_index] * valuesMatrixL[L_jk_index];
                    }

                    int L_ij_index = diagonalOffsetsMatrixA[i] + (i - j);
                    int L_jj_index = diagonalOffsetsMatrixA[j];

                    valuesMatrixL[L_ij_index] = (valuesMatrixA[L_ij_index] - sum1) / valuesMatrixL[L_jj_index];

                }

                double sum2 = 0;
                for (int k = mi; k < i; k++)  // k < j
                {
                    int L_ik_index = diagonalOffsetsMatrixA[i] + (i - k);

                    sum2 = sum2 + valuesMatrixL[L_ik_index] * valuesMatrixL[L_ik_index];
                }

                int L_ii_index = diagonalOffsetsMatrixA[i];
                valuesMatrixL[L_ii_index] = Math.Sqrt(valuesMatrixA[L_ii_index] - sum2);
            }

            Console.Write("L = ");
            VectorMethods.PrintVector(valuesMatrixL);

            return (valuesMatrixL, diagonalOffsetsMatrixL);
        }
    }
}
