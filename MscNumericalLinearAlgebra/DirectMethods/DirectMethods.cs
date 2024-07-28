using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.DirectMethods
{
    public static class DirectMethods
    {
        public static double[,] SolveLinearSystemGaussElimination(double[,] matrixA, double[,] matrixB)
        {

            int lenColMatrixA = matrixA.GetLength(1);
            int lenRowMatrixA = matrixA.GetLength(0);
            int lenColMatrixB = matrixB.GetLength(1);
            int lenRowMatrixB = matrixB.GetLength(0);

            double[,] vectorX = new double[lenRowMatrixA, lenColMatrixB];

            int maxRow = 0;

            for (int k = 0; k < lenRowMatrixA - 1; k++) // k is the column we are making pivoting
            {
                // Find row with the largest pivot value in the current column (k)
                double maxPivot = Math.Abs(matrixA[k, k]);  // Set maxPivot as the first Point 
                maxRow = k;  // Set the k-th row as the maxPivot row

                for (int i = k + 1; i < lenRowMatrixA; i++) // Run the rest of the columns to compare maxPivot with the rest of the rows
                {
                    if (Math.Abs(matrixA[i, k]) > maxPivot)  // If the current maxPivot is less than the next row
                    {
                        maxPivot = Math.Abs(matrixA[i, k]);  // Assign the new maxPivot
                        maxRow = i;  // Assign the row index of the maxPivot
                    }
                }

                // Swap the current row with the row containing the largest pivot
                if (maxRow != k)  // If the maxRow is not the initial one
                {
                    for (int j = 0; j < lenColMatrixA; j++)  // run all the elements of the maxRow
                    {
                        // Switch elements
                        // A[k, j] <-- Change elements --> A[maxRow, j]
                        double tempA = matrixA[k, j];
                        matrixA[k, j] = matrixA[maxRow, j];
                        matrixA[maxRow, j] = tempA;  // save it temporarily
                    }

                    for (int j = 0; j < lenColMatrixB; j++)  // run all the elements of the maxRow
                    {   
                        // Swicth element also in B matrix
                        // B[k, j] < --Change elements-- > B[maxRow, j]
                        double tempB = matrixB[k, j];
                        matrixB[k, j] = matrixB[maxRow, j];
                        matrixB[maxRow, j] = tempB;  // save it temporarily
                    }
                }

                // Gaussian Elimination
                for (int i = k + 1; i < lenRowMatrixA; i++)  // Run all rows under the 1st one
                {
                    double factor = matrixA[i, k] / matrixA[k, k];  // find the coefficient factor (Aik / Akk)
                    for (int j = k; j < lenColMatrixA; j++)  // Run all the columns from k-th element of the i-th row
                    {
                        matrixA[i, j] = -matrixA[i, j] + factor * matrixA[k, j]; // multiply with this so we get 0 in the current one
                    }

                    for (int j = 0; j < lenColMatrixB; j++)  // Run all the columns from k-th elemnt of the ith-row
                    {
                        matrixB[i, j] = -matrixB[i, j] + factor * matrixB[k, j]; // make the same multiplication for Matrix B
                    }
                }
            }

            // Back substitution
            for (int i = lenRowMatrixA - 1; i >= 0; i--) // from (m-1) to 0 for matrix A
            {
                for (int j = 0; j < lenColMatrixB; j++)  // j running all the colums of matrix B
                {
                    double sum = 0.0;  // Initialize sum
                    for (int k = i + 1; k < lenRowMatrixA; k++)  // 
                    {
                        sum += matrixA[i, k] * vectorX[k, j];
                    }
                    vectorX[i, j] = (matrixB[i, j] - sum) / matrixA[i, i];
                }
            }

            Console.WriteLine("x = ");
            MatrixMethods.PrintMatrix(vectorX);

            return vectorX;
        }

        public static double[] SolveFactorizationLULinearSystem(double[,] matrixL, double[,] matrixU, double[] vectorB)
        {

            //  A * x = b
            //  1. L * y = b
            //  2. U * x = y
            //
            // x = U^-1 (L^-1 * b)
            int n = matrixL.GetLength(0);  // Length matrix L & U
            double[] y = new double[n];  // Create y vector
            double[] x = new double[n];  // Create x vector

            if (matrixL.GetLength(0) != matrixL.GetLength(1))  // Check if the matrix is square
            {
                throw new Exception("The matrix should be square!");
            }

            // Forward Substitution (L * y = b), with b known vector
            for (int i = 0; i < n; i++)  // i = 0, 1, 2, ..., n
            {
                y[i] = vectorB[i];

                for (int j = 0; j < i; j++)
                {
                    y[i] = y[i] - matrixL[i, j] * y[j];
                }
            }

            // Backward Substitution (U * x = y), with y known vector
            for (int i = n - 1; i >= 0; i--)  // i = n-1, ..., 1, 0
            {
                x[i] = y[i];  // for a 3x3, x[3] = y[3]
                for (int j = i + 1; j < n; j++)  // for j = 1, 2, ..., n
                {
                    x[i] = x[i] - matrixU[i, j] * x[j];  //  U[1,1]*x1 = y[i](=x[i]) - U[1,2]*x2 - U[1,3]*x3
                }
                x[i] = x[i] / matrixU[i, i];  // x1 = (...) / U[1,1]
            }

            Console.Write("x = ");
            VectorMethods.PrintVector(x);

            return x;
        }

        public static double[] SolveFactorizationCholeskyLinearSystem(double[,] matrixL, double[] vectorB)
        {
            int n = matrixL.GetLength(0); // Number of Elements in a row
            double[] x = new double[n];   // Initilize result vector x
            double[] y = new double[n];   // Initilize vector y

            double[,] matrixLT = Matrices.MatrixMethods.TransposeMatrix(matrixL);

            if (n != vectorB.Length)
            {
                throw new Exception("Size of matrix L and vecotr b should be equal!");
            }

            // Forward Substitution (L * y = b), with b known vector
            for (int i = 0; i < n; i++)    // i = 0, 1, 2,..., n
            {

                y[i] = vectorB[i];

                for (int j = 0; j < i; j++)   // j = 0, 1,..., i - 1
                {

                    y[i] = y[i] - matrixL[i, j] * y[j];  // Subtract the known mebmers of the equation

                }

                y[i] = y[i] / matrixL[i, i];  // Calculate the final value of y
            }

            // Backward Substitution (U * x = y), with y known vector
            for (int i = n - 1; i >= 0; i--)  // i = n-1, ..., 0
            {
                x[i] = y[i];

                for (int j = i + 1; j < n; j++)   // j = i+1, ..., n
                {
                    x[i] = x[i] - matrixLT[i, j] * x[j];  // Subtract the known mebmers of the equation
                }

                x[i] = x[i] / matrixLT[i, i];  // Calculate the final value of x
            }

            Console.WriteLine("Vector x:");
            Matrices.VectorMethods.PrintVector(x);


            return x;
        }


        public static double[] SolveLinearSystemWithGaussElimination(double[,] matrixA, double[] vectorB)
        {
            //    A    *    X    =    b
            // (n x m) * (m x 1) = (n x 1)
            int rowsMatrixA = matrixA.GetLength(0);  // n
            int colsMatrixA = matrixA.GetLength(1);  // m

            int rowsVectorB = vectorB.Length;        // n

            double[] result = new double[colsMatrixA];

            if (rowsMatrixA != rowsVectorB)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b!");
            }

            int maxRow = 0;

            for ( int k = 0; k < rowsMatrixA; k++)
            {

                maxRow = k;
                double maxPivot = Math.Abs(matrixA[k, k]);

                for ( int i = k + 1; i <  rowsMatrixA; i++)  // run the rest rows
                {
                    if (Math.Abs(matrixA[i, k]) > maxPivot)  // Check if the below elements in the column are greater or not
                    {

                        maxRow = i;  // Assign new max row index
                        maxPivot = Math.Abs(matrixA[i, k]);  // Assign new max pivot value
                    }
                }

                // Swap the rows

                if (maxRow != k)
                {
                    for (int j = 0; j < colsMatrixA; j++)  // run all the columns of the swapping rows
                    {
                        // k <------> maxRow, for MatrixA
                        double tempValue1 = matrixA[k, j];     // save temporarily
                        matrixA[k, j] = matrixA[maxRow, j];    // A[k,j] = A[maxRow, j]
                        matrixA[maxRow, j] = tempValue1;       // A[maxRow,j] = A[k, j]

                        // k <------> maxRow, for vectorB
                        double tempValue2 = vectorB[k];     // save temporarily
                        vectorB[k] = vectorB[maxRow];       // b[k] = b[maxRow]
                        vectorB[maxRow] = tempValue2;       // b[maxRow] = b[k]
                    }
                }

                for ( int i = k + 1; i < rowsMatrixA; i++)
                {

                    double factor = matrixA[i, k] / matrixA[k, k];

                    for ( int j = k; j < colsMatrixA; j++)
                    {
                        matrixA[i, j] = - matrixA[i, j] + factor * matrixA[k, j];
                    }

                    vectorB[i] = -vectorB[i] + factor * vectorB[k];
                }

            }

            // Back substitution
            for (int i = rowsMatrixA - 1; i >= 0; i--) // from (m-1) to 0 for matrix A
            {

                double sum = 0.0;

                for (int j = i + 1; j < rowsMatrixA; j++)  // 
                {
                    sum = sum + matrixA[i, j] * result[j];
                }

                result[i] = (vectorB[i] - sum) / matrixA[i, i];
                
            }

            Console.Write("x = ");
            Matrices.VectorMethods.PrintVector(result);

            return result;
        }
    }
}
