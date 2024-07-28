using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Metadata;
using System.Text;
using System.Threading.Tasks;
using System.Transactions;

namespace MscNumericalLinearAlgebra.ExcerciseSeries1
{
    public static class Ex6
    {
        /// <summary>
        /// Converts a given matrix in the form of a 2D array into Compressed Sparse Row (CSR) format.
        /// The CSR format represents the matrix as three arrays: valuesArray, colIndices, and rowOffsets.
        /// </summary>
        /// <param name="matrixA">The input 2D matrix to be converted to CSR format.</param>
        /// <returns>
        /// A tuple containing the three arrays representing the matrix in CSR format:
        /// - valuesArray: An array containing the non-zero values of the matrix.
        /// - colIndices: An array containing the column indices of the non-zero values.
        /// - rowOffsets: An array indicating the starting position of each row in valuesArray and colIndices.
        /// </returns>
        public static (double[], int[], int[]) SaveMatrixCSR(double [,] matrixA)
        {
            int lenA = matrixA.GetLength(0); // Number of rows
            List<double> valuesList = new List<double>();
            List<int> colIndicesList = new List<int>();
            List<int> rowOffsetsList = new List<int>();


            int totalValues = 0;  // Total number of values
            int index = 0;  // Index indicating the position of the first non-zero element of valuesArray in each row

            for (int i = 0; i < lenA; i++)
            {
                bool nonZeroValue = false;

                for (int j = 0; j < lenA; j++)
                {
                    if (matrixA[i, j] != 0.0)
                    {
                        valuesList.Add(matrixA[i, j]);  // Add the cell value to a list
                        colIndicesList.Add(j);  // Add the values column index to a list
                        
                        if (!nonZeroValue)
                        {
                            rowOffsetsList.Add(index);  // Add the index of the value
                            nonZeroValue = true;  // Make nonZeroValue so in this line there will be no other added value to rowOffsetsList
                        }

                        index++;
                        totalValues++;  // Add to the total non-zero values number
                    }
                    
                }
            }
            rowOffsetsList.Add(valuesList.Count);  // Add at the last place the total number of element values

            // Convert the lists to arrays
            double[] valuesArray = valuesList.ToArray();
            int[] colIndices = colIndicesList.ToArray();
            int[] rowOffsets = rowOffsetsList.ToArray();

            Console.WriteLine("valuesArray = {" + string.Join(", ", valuesArray) + "}");
            Console.WriteLine("colIndices = {" + string.Join(", ", colIndices) + "}");
            Console.WriteLine("rowOffsets = {" + string.Join(", ", rowOffsets) + "}");

            return (valuesArray, colIndices, rowOffsets);
        }

        /// <summary>
        /// Multiplies a given matrix in Compressed Sparse Row (CSR) format with a vector.
        /// The CSR format is assumed to be represented as three arrays: valuesArray, colIndices, and rowOffsets.
        /// </summary>
        /// <param name="matrixA">The matrix in CSR format to be multiplied with the vector.</param>
        /// <param name="vectorX">The vector to be multiplied with the matrix.</param>
        /// <returns>
        /// A new vector resulting from the matrix-vector multiplication.
        /// </returns>
        /// <exception cref="Exception">Thrown if the dimensions of the matrix and vector are not compatible for multiplication.</exception>
        public static double[] MultiplyMatrixWithVector(double[,] matrixA, double[] vectorX)
        {
            if (matrixA.GetLength(0) != vectorX.Length)
            {
                    throw new Exception("The length of the two vectors is not equal!");
            }

            (double[] valuesArray, int[] colIndices, int[] rowOffsets) = SaveMatrixCSR(matrixA);

            double[] vectorB = new double[rowOffsets.Length - 1];

            for (int i = 0; i < (rowOffsets.Length - 1); i++)
            {
                vectorB[i] = 0;

                for (int k = rowOffsets[i]; k < rowOffsets[i + 1]; k++)
                {
                    vectorB[i] = vectorB[i] + valuesArray[k] * vectorX[colIndices[k]];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("b = {" + string.Join(", ", vectorB) + "}");

            return vectorB;
        }


        /// <summary>
        /// Converts a given matrix in the form of a 2D array into Compressed Sparse Column (CSC) format.
        /// The CSC format represents the matrix as three arrays: valuesArray, rowIndices, and colOffsets.
        /// </summary>
        /// <param name="matrixA">The input 2D matrix to be converted to CSC format.</param>
        /// <returns>
        /// A tuple containing the three arrays representing the matrix in CSC format:
        /// - valuesArray: An array containing the non-zero values of the matrix.
        /// - rowIndices: An array containing the row indices of the non-zero values.
        /// - colOffsets: An array indicating the starting position of each column in valuesArray and rowIndices.
        /// </returns>
        public static (double[], int[], int[]) SaveMatrixCSC(double[,] matrixA)
        {
            // Get the number of rows (lenA) and columns (lenA) in the input matrix.
            int lenA = matrixA.GetLength(0);

            // Initialize lists to store the values, row indices, and column offsets.
            List<double> valuesList = new List<double>();
            List<int> rowIndicesList = new List<int>();
            List<int> colOffsetsList = new List<int>();

            int totalValues = 0;  // Total number of non-zero values
            int index = 0;  // Index indicating the position of the first non-zero element of valuesArray in each column

            // Iterate through the matrix in column-major order to populate the CSC format arrays.
            for (int i = 0; i < lenA; i++)
            {
                bool nonZeroValue = false;

                for (int j = 0; j < lenA; j++)
                {
                    if (matrixA[j, i] != 0.0)
                    {
                        valuesList.Add(matrixA[j, i]);  // Add the cell value to a list
                        rowIndicesList.Add(j);  // Add the row index of the non-zero value to a list

                        if (!nonZeroValue)
                        {
                            colOffsetsList.Add(index);  // Add the index of the value
                            nonZeroValue = true;  // Ensure only one col offset is added per column
                        }

                        index++;
                        totalValues++;  // Increment the total non-zero values count
                    }
                }
            }

            colOffsetsList.Add(valuesList.Count);  // Add the total number of element values as the last col offset.

            // Convert the lists to arrays
            double[] valuesArray = valuesList.ToArray();
            int[] rowIndices = rowIndicesList.ToArray();
            int[] colOffsets = colOffsetsList.ToArray();

            // Print the converted CSC format arrays for debugging or analysis purposes.
            Console.WriteLine("valuesArray = {" + string.Join(", ", valuesArray) + "}");
            Console.WriteLine("rowIndices = {" + string.Join(", ", rowIndices) + "}");
            Console.WriteLine("colOffsets = {" + string.Join(", ", colOffsets) + "}");

            return (valuesArray, rowIndices, colOffsets);
        }

        /// <summary>
        /// Multiplies the transpose of a given matrix in Compressed Sparse Column (CSC) format with a vector.
        /// The CSC format is assumed to be represented as three arrays: valuesArray, rowIndices, and colOffsets.
        /// </summary>
        /// <param name="matrixA">The matrix in CSC format to be transposed and multiplied with the vector.</param>
        /// <param name="vectorX">The vector to be multiplied with the transpose of the matrix.</param>
        /// <returns>
        /// A new vector resulting from the multiplication of the transpose of the matrix with the vector.
        /// </returns>
        /// <exception cref="Exception">Thrown if the dimensions of the matrix and vector are not compatible for multiplication.</exception>
        public static double[] MultiplyTransposeMatrixWithVector(double[,] matrixA, double[] vectorX)
        {
            // Check if the dimensions of the matrix and vector are compatible for multiplication.
            if (matrixA.GetLength(0) != vectorX.Length)
            {
                throw new Exception("The length of the matrix and vector is not compatible for multiplication.");
            }

            // Retrieve the CSC format arrays representing the transpose of matrixA.
            (double[] valuesArray, int[] rowIndices, int[] colOffsets) = SaveMatrixCSC(matrixA);

            // Initialize a vector to store the result of the multiplication.
            double[] vectorB = new double[colOffsets.Length - 1];

            // Perform the transpose matrix-vector multiplication.
            for (int i = 0; i < (colOffsets.Length - 1); i++)
            {
                vectorB[i] = 0;

                for (int k = colOffsets[i]; k < colOffsets[i + 1]; k++)
                {
                    vectorB[i] = vectorB[i] + valuesArray[k] * vectorX[rowIndices[k]];
                }
            }

            // Print the result vector for debugging or analysis purposes.
            Console.WriteLine("");
            Console.WriteLine("Result vector (b) = {" + string.Join(", ", vectorB) + "}");

            return vectorB;
        }





    }
}
