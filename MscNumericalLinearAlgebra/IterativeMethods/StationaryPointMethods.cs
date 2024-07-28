using MscNumericalLinearAlgebra.Matrices;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.IterativeMethods
{
    public static class StationaryPointMethods
    {
        public static double[] Jacobi(double[,] matrixA, double[] vectorB, double epsilon, int maxIterations = 100)
        {

            int current_t = maxIterations;  // In case the for loop makes all the iterations
            int n = matrixA.GetLength(0);  // dimension of the matrix A (also for vectorB)

            // Condition to check the dimensions of the matrices
            if (n != matrixA.GetLength(1) || n != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b!");  // Print Message to fix the dimensions
            }

            // Initilization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_norm_new = 0;

            for (int t = 0; t < maxIterations; t++)  // t: Current Number of Iteration
            {

                for (int i = 0; i < n; i++)  // run all the rows of Aij
                {
                    // x(t+1) = [D^-1 * L * x(t)] + [D^-1 * U * x(t)] + [D^-1 * b]
                    // x(t+1) =        [s1]       +        [s2]       + [D^-1 * b]
                    double s1 = 0;
                    double s2 = 0;

                    for (int j = 0; j < i; j++)  // run all the elements of a row of Aij
                    {
                        // x(i) = D^-1 * L * x(j)
                        s1 = s1 + 1 / matrixA[i, i] * (-matrixA[i, j] * xt[j]);  // Sum of xi for L
                    }

                    for (int j = i + 1; j < n; j++)  // run all the elements of a row of Aij
                    {
                        // x(i) = D^-1 * U * x(j)
                        s2 = s2 + 1 / matrixA[i, i] * (-matrixA[i, j] * xt[j]);  // Sum of xi for D
                    }

                    // x = D^-1 * (L + U) x + D^-1 * b
                    xt_new[i] = s1 + s2 + 1 / matrixA[i, i] * vectorB[i];  // xi = (s1 + s2) + 1 / Aii * bi

                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                // Calculate the r(t) = b - A * x(t)
                for (int row = 0; row < n; row++)
                {
                    double sum_row = 0;

                    for (int col = 0; col < n; col++)
                    {
                        sum_row = sum_row + matrixA[row, col] * xt_new[col];
                    }
                    rt_new[row] = vectorB[row] - sum_row;
                }

                rt_norm_new = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_norm_new / rt_norm <= epsilon)  // If ||r(t)|| / ||r0|| < ε ---> Break Loop
                {
                    current_t = t;  // Save current Iteration
                    break;
                }

                xt = xt_new; // Assign x(t) = x(t+1), for the next loop
            }

            Console.WriteLine("Residual ||rt|| = " + rt_norm_new);
            Console.WriteLine("Number of Iterations N = " + (current_t + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new; // Return Solution x
        }

        public static double[] GaussSeidel(double[,] matrixA, double[] vectorB, double epsilon, int maxIterations = 100)
        {

            int currentIteration = maxIterations;  // In case the for loop makes all the iterations
            int n = matrixA.GetLength(0);  // dimension of the matrix A (also for vectorB)

            // Condition to check the dimensions of the matrices
            if (n != matrixA.GetLength(1) || n != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b");  // Print Message to fix the dimensions
            }

            // Initilization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double[] C = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_new_norm = 0;

            for (int t = 0; t < maxIterations; t++)
            {
                //  Solve (D - L) * x[t+1] = U * x[t] + b
                for (int i = 0; i < n; i++)
                {

                    //Calculate C = U * x[t] + b
                    C[i] = 0;
                    double sum = 0;

                    for (int j = i + 1; j < n; j++)
                    {
                        sum = sum - matrixA[i, j] * xt[j];
                    }
                    C[i] = vectorB[i] + sum;  // C = U * x[t] + b


                    // Calculate with Forward Substitution x[t+1] = (D - L)^-1 * C
                    xt_new[i] = 0;
                    sum = 0;
                    for (int j = 0; j < i; j++)
                    {
                        sum = sum - matrixA[i, j] * xt_new[j];
                    }
                    xt_new[i] = (C[i] + sum) / matrixA[i, i];
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                // Calculate the r(t) = b - A * x(t)
                for (int row = 0; row < n; row++)
                {
                    double sum_row = 0;

                    for (int col = 0; col < n; col++)
                    {
                        sum_row = sum_row + matrixA[row, col] * xt_new[col];
                    }
                    rt_new[row] = vectorB[row] - sum_row;
                }

                rt_new_norm = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_new_norm / rt_norm <= epsilon)  // If ||r(t)|| / ||r0|| < ε ---> Break Loop
                {
                    currentIteration = t;  // Save current Break Iteration
                    break;
                }

                xt = xt_new;
            }

            Console.WriteLine("Residual ||rt|| = " + rt_new_norm);
            Console.WriteLine("Number of Iterations N = " + (currentIteration + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new;
        }

        public static double[] SOR(double[,] matrixA, double[] vectorB, double epsilon, double weight, int maxIterations = 100)
        {
            int currentIteration = maxIterations;
            int n = matrixA.GetLength(0);

            // Condition to check the dimensions of the matrices
            if (n != matrixA.GetLength(1) || n != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b");  // Print Message to fix the dimensions
            }

            // Initilization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_new_norm = 0;

            for (int t = 0; t < maxIterations; t++)
            {
                for (int i = 0; i < n; i++)
                {
                    double sum1 = 0;
                    double sum2 = 0;

                    for (int j = 0; j < i; j++)
                    {
                        sum1 = sum1 + matrixA[i, j] * xt_new[j];
                    }

                    for (int j = i + 1; j < n; j++)
                    {
                        sum2 = sum2 + matrixA[i, j] * xt[j];
                    }

                    xt_new[i] = (1 - weight) * xt[i] + (weight / matrixA[i, i]) * (vectorB[i] - sum1 - sum2);
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                // Calculate the r(t) = b - A * x(t)
                for (int row = 0; row < n; row++)
                {
                    double sum_row = 0;

                    for (int col = 0; col < n; col++)
                    {
                        sum_row += matrixA[row, col] * xt_new[col];
                    }

                    rt_new[row] = vectorB[row] - sum_row;
                }

                rt_new_norm = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_new_norm / rt_norm <= epsilon)  // If ||r(t)|| / ||r0|| < ε ---> Break Loop
                {
                    currentIteration = t;  // Save current Break Iteration
                    break;
                }

                xt = xt_new;
            }

            Console.WriteLine("Residual ||rt|| = " + rt_new_norm);
            Console.WriteLine("Number of Iterations N = " + (currentIteration + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new;
        }


        public static double[] JacobiWithCOO(double[] valuesMatrixA, int[] rowsMatrixA, int[] colsMatrixA, double[] vectorB, double epsilon, int maxIterations = 100)
        {
            int current_t = maxIterations;              // In case the for loop makes all the iterations
            int n = rowsMatrixA.Max() + 1;              // dimension of the matrix A (also for vectorB)
            int m = colsMatrixA.Max() + 1;              // dimension of the matrix A (also for vectorB)
            int numberValues = valuesMatrixA.Length;    // Number of all elements 

            Console.WriteLine("m = " + m);
            Console.WriteLine("m = " + n);
            Console.WriteLine("m = " + vectorB.Length);
            if (n != m || n != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b!");  // Print Message to fix the dimensions
            }

            // Initialization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_norm_new = 0;

            double[] diagonalValues = new double[n];  // Array to store diagonal elements

            for (int i = 0; i < valuesMatrixA.Length; i++) // Find diagonal elements
            {
                int row = rowsMatrixA[i];
                int col = colsMatrixA[i];

                if (row == col)
                {
                    diagonalValues[row] = valuesMatrixA[i];
                }
            }

            for (int t = 0; t < maxIterations; t++)  // t: Current Number of Iteration
            {
                // Initialize total s1 and s2 for each row
                double[] totalS1 = new double[n];
                double[] totalS2 = new double[n];

                for (int i = 0; i < numberValues; i++)  // run all elements
                {
                    int row = rowsMatrixA[i];
                    int col = colsMatrixA[i];

                    
                    if (row != col)  // Exclude diagonals
                    {
                        if (row < col)
                        {
                            
                            totalS2[row] = totalS2[row] + valuesMatrixA[i] * xt[col]; // Upper triangular (U)
                        }
                        else if (row > col)
                        {
                            
                            totalS1[row] = totalS1[row] + valuesMatrixA[i] * xt[col]; // Lower triangular
                        }
                    }
                }

           
                for (int row = 0; row < n; row++)  // Calculate x
                {
                    xt_new[row] = (vectorB[row] - totalS1[row] - totalS2[row]) / diagonalValues[row];
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                double[] Ax = Matrices.MatrixMethods.MatrixVectorMultiplicationWithCOO(valuesMatrixA, rowsMatrixA, colsMatrixA, xt_new);
                // Calculate the r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    rt_new[i] = vectorB[i] - Ax[i];
                }

                rt_norm_new = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_norm_new <= epsilon)  // If the residual norm is less than or equal to epsilon ---> Break Loop
                {
                    current_t = t;  // Save current Iteration
                    break;
                }

                xt = xt_new; // Assign x(t) = x(t+1), for the next loop
            }

            Console.WriteLine("Residual ||rt|| = " + rt_norm_new);
            Console.WriteLine("Number of Iterations N = " + (current_t + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new;
        }


        public static double[] JacobiWithCSR(double[] values, int[] colIndices, int[] rowOffsets, double[] vectorB, double epsilon, int maxIterations = 100)
        {
            int current_t = maxIterations;              // In case the for loop makes all the iterations
            int n = rowOffsets.Length - 1;              // dimension of the matrix A (also for vectorB)

            if (n != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b!");  // Print Message to fix the dimensions
            }

            // Initialization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_norm_new = 0;

            for (int t = 0; t < maxIterations; t++)
            {
                // Initialize total s1 and s2 for each row
                double[] totalS1 = new double[n];
                double[] totalS2 = new double[n];

                for (int row = 0; row < n; row++)  // run all rows of Matrix A
                {
                    double diagonalValue = 0;

                    for (int i = rowOffsets[row]; i < rowOffsets[row + 1]; i++)  // run all elements in the row
                    {
                        int col = colIndices[i];

                        if (row != col)  // Exclude diagonals
                        {
                            if (row < col)
                            {
                                totalS2[row] = totalS2[row] + values[i] * xt[col]; // Upper triangular (U)
                            }
                            else if (row > col)
                            {
                                totalS1[row] = totalS1[row] + values[i] * xt[col]; // Lower triangular (L)
                            }
                        }
                        else  // If it is diagonal
                        {
                            diagonalValue = values[i];  // diagonal element value
                        }
                    }
   
                    xt_new[row] = (vectorB[row] - totalS1[row] - totalS2[row]) / diagonalValue; // x[i] = (b[i] - s1 - s2) / A[i,i]
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                double[] Ax = Matrices.MatrixMethods.MatrixVectorMultiplicationWithCSR(values, colIndices, rowOffsets, xt_new);
                
                for (int row = 0; row < n; row++)
                {
                    rt_new[row] = vectorB[row] - Ax[row]; // r(t) = b - A * x(t)
                }

                rt_norm_new = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_norm_new <= epsilon)  // If the residual norm is less than or equal to epsilon ---> Break Loop
                {
                    current_t = t;  // Save current Iteration
                    break;
                }

                xt = xt_new; // Assign x(t) = x(t+1), for the next loop
            }

            Console.WriteLine("Residual ||rt|| = " + rt_norm_new);
            Console.WriteLine("Number of Iterations N = " + (current_t + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new;
        }


        public static double[] JacobiWithCSC(double[] valuesMatrixA, int[] rowIndicesMatrixA, int[] colOffsetsMatrixA, double[] vectorB, double epsilon, int maxIterations = 100)
        {
            int current_t = maxIterations;  // In case the for loop makes all the iterations
            int n = colOffsetsMatrixA.Length - 1;  // dimension of the matrix A (also for vectorB)
            int numberValues = valuesMatrixA.Length;  // Number of all elements 

            if (n != vectorB.Length)
            {
                throw new Exception("Wrong Dimensions of the Linear System Ax=b!");
            }

            // Initialization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_norm_new = 0;

            for (int t = 0; t < maxIterations; t++)
            {
                // Initialize total s1 and s2 for each column so we can calculate them for each row
                double[] totalS1 = new double[n];         // Array for L matrix
                double[] totalS2 = new double[n];         // Array for U matrix
                double[] diagonalValues = new double[n];  // Array to save diagonal elements

                for (int col = 0; col < n; col++)  // run all columns of Aij 
                {
                    // run all elements of each column
                    for (int index = colOffsetsMatrixA[col]; index < colOffsetsMatrixA[col + 1]; index++)  // index: local system indexing
                    {
                        int row = rowIndicesMatrixA[index];  // row index of each element


                        if (row != col) // Exclude diagonals
                        {
                            if (row < col)  // i < j
                            {
                                totalS1[row] += valuesMatrixA[index] * xt[col];  // Lower triangular (L)
                            }
                            else if (row > col)  // i > j
                            {
                                totalS2[row] += valuesMatrixA[index] * xt[col];  // Upper triangular (U)
                            }
                        }
                        else  // Diagonal element
                        {
                            diagonalValues[row] = valuesMatrixA[index];
                        }
                    }
                }

                for (int row = 0; row < n; row++)  // run all columns of Matrix A
                {
                    xt_new[row] = (vectorB[row] - totalS1[row] - totalS2[row]) / diagonalValues[row]; // x[i] = (b[i] - s1[i] - s2[i]) / A[i,i]
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                double[] Ax = Matrices.MatrixMethods.MatrixVectorMultiplicationWithCSC(valuesMatrixA, rowIndicesMatrixA, colOffsetsMatrixA, xt_new);
                // Calculate the r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    rt_new[i] = vectorB[i] - Ax[i];
                }

                rt_norm_new = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_norm_new <= epsilon)  // If the residual norm is less than or equal to epsilon ---> Break Loop
                {
                    current_t = t;  // Save current Iteration
                    break;
                }

                xt = xt_new; // Assign x(t) = x(t+1), for the next loop
            }

            Console.WriteLine("Residual ||rt|| = " + rt_norm_new);
            Console.WriteLine("Number of Iterations N = " + (current_t + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new; // Return Solution x
        }





        public static double[] GaussSeidelWithCOO(double[] valuesMatrixA, int[] rowsMatrixA, int[] colsMatrixA, double[] vectorB, double epsilon, int maxIterations = 100)
        {
            int currentIteration = maxIterations;  // In case the for loop makes all the iterations
            int n = vectorB.Length;  // dimension of the matrix A (also for vectorB)
            int numberValues = valuesMatrixA.Length;  // Number of all elements 

            // Initialization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double[] C = new double[n];
            double[] diagonalValues = new double[n];  // Array to save diagonal elements
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_new_norm = 0;

            // Calculate diagonal values for each row
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < numberValues; j++)
                {
                    if (rowsMatrixA[j] == i && colsMatrixA[j] == i)
                    {
                        diagonalValues[i] = valuesMatrixA[j];
                        break;
                    }
                }
            }

            for (int t = 0; t < maxIterations; t++)
            {
                // Solve (D - L) * x[t+1] = U * x[t] + b
                for (int i = 0; i < n; i++)
                {
                    // Calculate C = U * x[t] + b
                    C[i] = 0;
                    double sum = 0;

                    for (int j = 0; j < numberValues; j++)
                    {
                        if (rowsMatrixA[j] == i)
                        {
                            if (colsMatrixA[j] == i)
                            {
                                C[i] = vectorB[i];
                            }
                            else
                            {
                                sum += valuesMatrixA[j] * xt[colsMatrixA[j]];
                            }
                        }
                    }
                    C[i] += sum;  // C = U * x[t] + b

                    // Calculate with Forward Substitution x[t+1] = (D - L)^-1 * C
                    xt_new[i] = 0;
                    sum = 0;
                    for (int j = 0; j < numberValues; j++)
                    {
                        if (rowsMatrixA[j] == i && colsMatrixA[j] < i)
                        {
                            sum += valuesMatrixA[j] * xt_new[colsMatrixA[j]];
                        }
                    }
                    xt_new[i] = (C[i] + sum) / diagonalValues[i];
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                // Calculate the r(t) = b - A * x(t)
                double[] Ax = Matrices.MatrixMethods.MatrixVectorMultiplicationWithCOO(valuesMatrixA, rowsMatrixA, colsMatrixA, xt_new);
                // Calculate the r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    rt_new[i] = vectorB[i] - Ax[i];
                }

                rt_new_norm = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_new_norm / rt_norm <= epsilon)  // If ||r(t)|| / ||r0|| < ε ---> Break Loop
                {
                    currentIteration = t;  // Save current Break Iteration
                    break;
                }

                xt = xt_new;
            }

            Console.WriteLine("Residual ||rt|| = " + rt_new_norm);
            Console.WriteLine("Number of Iterations N = " + (currentIteration + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new;
        }



        public static double[] GaussSeidelCOO(double[] valuesMatrixA, int[] rowsMatrixA, int[] colsMatrixA, double[] vectorB, double epsilon, int maxIterations = 100)
        {
            int currentIteration = maxIterations;  // In case the for loop makes all the iterations
            int n = vectorB.Length;  // dimension of the matrix A (also for vectorB)

            // Initialization of arrays and values needed for the iterations
            double[] xt = new double[n];
            double[] xt_new = new double[n];
            double rt_norm = VectorMethods.VectorEuclideanNorm(vectorB);
            double rt_new_norm = 0;

            // Find diagonal values for each row
            double[] diagonalValues = new double[n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < valuesMatrixA.Length; j++)
                {
                    if (rowsMatrixA[j] == i && colsMatrixA[j] == i)
                    {
                        diagonalValues[i] = valuesMatrixA[j];
                        break;
                    }
                }
            }

            for (int t = 0; t < maxIterations; t++)
            {

                double[] sum2 = new double[n];  

                for (int j = 0; j < valuesMatrixA.Length; j++)  // run all values
                {
                    int row = rowsMatrixA[j];  // row index
                    int col = colsMatrixA[j];  // column index

                    if (col < row)  // lower triangular elements
                    {
                        sum2[row] = sum2[row] + valuesMatrixA[j] * xt[col];  // save the sum for each row
                    }
                }

                for (int i = 0; i < n; i++)  // run thorugh rows
                {
                    double sum = 0;

                    for (int j = 0; j < valuesMatrixA.Length; j++)  //  run all values
                    {
                        int row = rowsMatrixA[j];  // row index
                        int col = colsMatrixA[j];  // column index

                        if (row == i && col > i) // Being on the upper triangular
                        {
                            // C = U * x[t] + b
                            sum = sum - valuesMatrixA[j] * xt[col];
                        }
                    }

                    // x[t+1] = (D - L)^-1 * C + (D - L)^-1 * b
                    xt_new[i] = (vectorB[i] + sum - sum2[i]) / diagonalValues[i];
                }

                double[] rt_new = new double[n];  // Residual ||r(t)||

                // Calculate the r(t) = b - A * x(t)
                double[] Ax = Matrices.MatrixMethods.MatrixVectorMultiplicationWithCOO(valuesMatrixA, rowsMatrixA, colsMatrixA, xt_new);

                // Calculate the r(t) = b - A * x(t)
                for (int i = 0; i < n; i++)
                {
                    rt_new[i] = vectorB[i] - Ax[i];
                }

                rt_new_norm = VectorMethods.VectorEuclideanNorm(rt_new); // Calculate ||r(t)|| = √(r1^2 + r2^2 + ... + rn^2)

                if (rt_new_norm / rt_norm <= epsilon)  // If ||r(t)|| / ||r0|| < ε ---> Break Loop
                {
                    currentIteration = t;  // Save current Break Iteration
                    break;
                }

                xt = xt_new;
            }

            Console.WriteLine("Residual ||rt|| = " + rt_new_norm);
            Console.WriteLine("Number of Iterations N = " + (currentIteration + 1));
            Console.Write("Solution x = ");
            VectorMethods.PrintVector(xt_new);

            return xt_new;
        }





    }
}
