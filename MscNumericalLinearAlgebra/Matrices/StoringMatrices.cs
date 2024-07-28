using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.Matrices
{
    public static class StoringMatrices
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
        public static (double[], int[], int[]) SaveMatrixCSR(double[,] matrixA)
        {
            int lenA = matrixA.GetLength(0); // Number of rows
            List<double> valuesList = new List<double>();  // List with values of matrix A
            List<int> colIndicesList = new List<int>();    // List with Column Indices of the Values
            List<int> rowOffsetsList = new List<int>();    // List with the Row Offsets o


            int totalValues = 0;  // Total number of values
            int index = 0;  // Index indicating the position of the first non-zero element of valuesArray in each row

            for (int i = 0; i < lenA; i++)
            {
                bool nonZeroValue = false;

                for (int j = 0; j < matrixA.GetLength(1); j++)
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
        /// Saves a matrix in COO format with specified storing method.
        /// </summary>
        /// <param name="matrixA">The 2D matrix to save in COO format.</param>
        /// <param name="storingMethodMajor">Storing method: "Row" (default) or "Column".</param>
        /// <returns>Tuple of arrays representing COO format: values, row indices, and column indices.</returns>
        /// <remarks>Saves non-zero elements with their values, row, and column indices.</remarks>
        public static (double[], int[], int[]) SaveMatrixCOO(double[,] matrixA, string storingMethodMajor = "Row")
        {
            if (storingMethodMajor == "Row")
            {
                int numRowsA = matrixA.GetLength(0); // Number of rows
                int numColsA = matrixA.GetLength(1); // Number of rows

                List<double> valuesList = new List<double>();  // List with values of matrix A
                List<int> rowsList = new List<int>();    // List with Column Indices of the Values
                List<int> colList = new List<int>();    // List with the Row Offsets

                for (int i = 0; i < numRowsA; i++)
                {
                    for (int j = 0; j < numColsA; j++)
                    {
                        if (matrixA[i, j] != 0)
                        {
                            valuesList.Add(matrixA[i, j]);
                            rowsList.Add(i);
                            colList.Add(j);
                        }
                    }
                }

                double[] valuesArray = valuesList.ToArray();
                int[] rowsArray = rowsList.ToArray();
                int[] colsArray = colList.ToArray();

                Console.WriteLine("Row Major COO Method:");
                Console.WriteLine("valuesArray = {" + string.Join(", ", valuesArray) + "}");
                Console.WriteLine("rowsArray = {" + string.Join(", ", rowsArray) + "}");
                Console.WriteLine("colsArray = {" + string.Join(", ", colsArray) + "}");

                return (valuesArray, rowsArray, colsArray);
            }
            else
            {
                int numRowsA = matrixA.GetLength(0); // Number of rows
                int numColsA = matrixA.GetLength(1); // Number of rows

                List<double> valuesList = new List<double>();  // List with values of matrix A
                List<int> rowsList = new List<int>();    // List with Column Indices of the Values
                List<int> colList = new List<int>();    // List with the Row Offsets

                for (int i = 0; i < numColsA; i++)
                {
                    for (int j = 0; j < numRowsA; j++)
                    {
                        if (matrixA[j, i] != 0)
                        {
                            valuesList.Add(matrixA[j, i]);
                            colList.Add(i);
                            rowsList.Add(j);
                        }
                    }
                }

                double[] valuesArray = valuesList.ToArray();
                int[] rowsArray = rowsList.ToArray();
                int[] colsArray = colList.ToArray();

                Console.WriteLine("Column Major COO Method:");
                Console.WriteLine("valuesArray = {" + string.Join(", ", valuesArray) + "}");
                Console.WriteLine("colsArray= {" + string.Join(", ", colsArray) + "}");
                Console.WriteLine("rowsArray = {" + string.Join(", ", rowsArray) + "}");

                return (valuesArray, colsArray, rowsArray);
            }
        }
    }
}
