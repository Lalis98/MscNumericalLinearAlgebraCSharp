using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.Intrinsics;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.Matrices
{
    public static class MatrixMethods
    {

        public static double[] MatrixVectorMultiplicationWithCOO(double[] valuesMatrixA, int[] rowsMatrixA, int[] columnsMatrixA, double[] vectorB)
        {
            int n = vectorB.Length;  // Size of A and b
            double[] result = new double[n]; // Result Vector

            for (int j = 0; j < valuesMatrixA.Length; j++) // Run all values
            {
                int row = rowsMatrixA[j];        // Row Index
                int column = columnsMatrixA[j];  // Column Index

                result[row] = result[row] + valuesMatrixA[j] * vectorB[column];  // result = A * b
            }

            return result;
        }

        public static double[] TransposeMatrixVectorMultiplicationWithCOO(double[] valuesMatrixA, int[] rowsMatrixA, int[] columnsMatrixA, double[] vectorB)
        {
            int numCowsMatrixA = columnsMatrixA.Max() + 1; 
            double[] result = new double[numCowsMatrixA]; // Result vector

            for (int j = 0; j < valuesMatrixA.Length; j++) // Run through all values
            {
                int row = rowsMatrixA[j];       // Row Index
                int column = columnsMatrixA[j]; // Column Index

                result[column] = result[column] + valuesMatrixA[j] * vectorB[row]; // result = A^T * b
            }

            return result;
        }



        public static double[] MatrixVectorMultiplicationWithCSR(double[] valuesArrayMatrixA, int[] colIndicesMatrixA, int[] rowOffsetsMatrixA, double[] vectorB)
        {
            int numRowsMatrixA = rowOffsetsMatrixA.Length - 1;  // Rows Number of Matrix A
            double[] resultVector = new double[numRowsMatrixA]; // Size of result vector

            for (int i = 0; i < rowOffsetsMatrixA.Length - 1; i++)  // run through all rows
            {
                resultVector[i] = 0.0;

                for (int j = rowOffsetsMatrixA[i]; j < rowOffsetsMatrixA[i + 1]; j++)  // run the elements of the row
                {
                    resultVector[i] = resultVector[i] + valuesArrayMatrixA[j] * vectorB[colIndicesMatrixA[j]];  // result = A * b
                }

            }

            return resultVector;
        }

        public static double[] TransposeMatrixVectorMultiplicationWithCSR(double[] valuesArrayMatrixA, int[] colIndicesMatrixA, int[] rowOffsetsMatrixA, double[] vectorB)
        {
            int numColsMatrixA = colIndicesMatrixA.Max() + 1; // Number of columns in Matrix A
            double[] result = new double[numColsMatrixA];

            for (int i = 0; i < rowOffsetsMatrixA.Length - 1; i++) // Loop over each row
            {
                for (int j = rowOffsetsMatrixA[i]; j < rowOffsetsMatrixA[i + 1]; j++) // Loop over non-zero elements in the row
                {
                    int col = colIndicesMatrixA[j];

                    result[col] = result[col] + valuesArrayMatrixA[j] * vectorB[i]; // Accumulate in result
                }
            }

            return result;
        }


        public static double[] MatrixVectorMultiplicationWithCSC(double[] valuesArrayMatrixA, int[] rowIndicesMatrixA, int[] colOffsetsMatrixA, double[] vectorB)
        {
            int numRowsMatrixA = rowIndicesMatrixA.Max() + 1; // Number of Rows of Matrix A
            double[] result = new double[numRowsMatrixA];     // Size of result vector

            for (int col = 0; col < vectorB.Length; col++)  // Run through all columns
            {

                for (int i = colOffsetsMatrixA[col]; i < colOffsetsMatrixA[col + 1]; i++)  // Run all elements of each column
                {
                    int row = rowIndicesMatrixA[i];  // find the roe-Index of each element

                    result[row] = result[row] + valuesArrayMatrixA[i] * vectorB[col];  // result = A * x
                }
            }

            return result;
        }

        public static double[] TransposeMatrixVectorMultiplicationWithCSC(double[] valuesArrayMatrixA, int[] rowIndicesMatrixA, int[] colOffsetsMatrixA, double[] vectorB)
        {
            int numColsMatrixA = colOffsetsMatrixA.Length - 1; // Number of Columns of Matrix A
            double[] result = new double[numColsMatrixA];      // Size of result vector

            for (int i = 0; i < valuesArrayMatrixA.Length; i++)  // Run all non-zero elements of Matrix A
            {
                int col = 0;  // column-Index of current element

                for (col = 0; col < numColsMatrixA; col++)
                {
                    if (i < colOffsetsMatrixA[col + 1])
                    {
                        break;
                    }
                }

                int row = rowIndicesMatrixA[i]; //  row - Index

                result[col] = result[col] + valuesArrayMatrixA[i] * vectorB[row];  // result = A^T * b
            }

            return result;
        }


        public static double[] MatrixVectorMultiplicationWithSkyline(double[] valuesMatrixA, int[] diagonalOffsetsMatrixA, double[] vectorB)
        {

            int n = vectorB.Length;           // Size of Linear System
            double[] result = new double[n];  // Size of result vector

            for (int i = 0; i < n; i++)
            {
                int diagonalIndex = diagonalOffsetsMatrixA[i];
                result[i] = result[i] + valuesMatrixA[diagonalIndex] * vectorB[i];

                //valuesMatrixA = [21, 1, 22, 2, 0, 23, 1, 3, 4, 1, 24, 2, 0, 3, 2, 25]
                //diagonalOffsetsMatrixA = [0, 1, 3, 6, 10, 15]
                //
                // i.e. diagonalIndex = 3
                //      endIndex = 6

                int endIndex; 
                if (i + 1 < n)
                {
                    endIndex = diagonalOffsetsMatrixA[i + 1]; // if it is not the last diagonal element
                }
                else
                {
                    endIndex = valuesMatrixA.Length;  // the last diagonal element
                }

                for (int j = diagonalIndex + 1; j < endIndex; j++)  // run elements of the specific diagonal element
                {
                    int row = i;                            // i-row
                    int column = i - (j - diagonalIndex);   // j-column
                    result[row] = result[row] + valuesMatrixA[j] * vectorB[column];  // find result for the elements above the diagonal
                    
                    if (row != column)  // if not a diagonal element
                    {
                        result[column] = result[column] + valuesMatrixA[j] * vectorB[row]; // // find result for the transpose elements (left to the diagonal)
                    }
                }
            }

            return result;
        }

        public static double[] MatrixVectorMultiplicationWithBanded(double[] valuesMatrixA, int size, int bandwidth, double[] vectorB)
        {

            double[] result = new double[size];

            if (size != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b!");
            }

            int index = 0;  // index of elements

            for (int i = 0; i < size; i++)  // run all rows
            {
                int start = Math.Max(0, i - bandwidth);         // Start index of row
                int end = Math.Min(size, i + bandwidth + 1);    // End index of row

                for (int j = start; j < end; j++)   // run all the elements of row
                {
                    result[i] = result[i] + valuesMatrixA[index] * vectorB[j];  // A * b
                    index = index + 1;  // incease index by 1
                }
            }

            return vectorB;
        }




        public static double[,] CreateIdentityMatrix(int length)
        {
            double[,] identityMatrix = new double[length, length];

            for (int i = 0; i < length; i++)
            {
                identityMatrix[i, i] = 1;
            }

            return identityMatrix;
        }
        


        // Function to print a matrix with enhanced formatting and rounded to 2 decimals
        public static void PrintMatrix(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            for (int i = 0; i < rows; i++)
            {
                Console.Write("[ ");
                for (int j = 0; j < cols; j++)
                {
                    Console.Write(matrix[i, j].ToString("N3")); // Round to 2 decimals
                    if (j < cols - 1)
                    {
                        Console.Write(", ");
                    }
                }
                Console.WriteLine(" ]");
            }
            Console.WriteLine();
        }


        public static double[,] TransposeMatrix(double[,] matrix)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);

            double[,] transpose = new double[cols, rows];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    transpose[j, i] = matrix[i, j]; // Swap row and column indices to transpose
                }
            }

            return transpose; // Return the transposed matrix
        }

    }
}
