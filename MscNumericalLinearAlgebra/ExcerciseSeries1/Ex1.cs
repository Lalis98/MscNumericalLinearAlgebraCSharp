using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MscNumericalLinearAlgebra.Matrices;

namespace MscNumericalLinearAlgebra.ExcerciseSeries1
{
    public static class Ex1
    {
        /// <summary>
        /// Adds two vectors and returns the result.
        /// </summary>
        /// <param name="vector1">The first input vector.</param>
        /// <param name="vector2">The second input vector.</param>
        /// <returns>
        /// A new vector containing the element-wise addition of the input vectors.
        /// </returns>
        /// <exception cref="Exception">
        /// Thrown when the input vectors have different lengths.
        /// </exception>
        public static double[] AddVectors(double[] vector1, double[] vector2) 
        {
            if (vector1.Length != vector2.Length)
            {
                throw new Exception("The vectors must have the same length!");
            }

            int vectorLength = vector1.Length;
            double[] resultVector = new double[vectorLength];

            for (int i = 0; i < vectorLength; i++)
            {
                resultVector[i] = vector1[i] + vector2[i];
                /*Console.Write(resultVector[i] + " ");*/
            }

            return resultVector;
        }

        /// <summary>
        /// Subtracts two vectors and returns the result.
        /// </summary>
        /// <param name="vector1">The first input vector.</param>
        /// <param name="vector2">The second input vector.</param>
        /// <returns>
        /// A new vector containing the subtraction of the input vectors.
        /// </returns>
        /// <exception cref="Exception">
        /// Thrown when the input vectors have different lengths.
        /// </exception>
        public static double[] SubtractVectors(double[] vector1, double[] vector2)
        {
            if (vector1.Length != vector2.Length)
            {
                throw new Exception("The vectors must have the same length!");
            }

            int vectorLength = vector1.Length;
            double[] resultVector = new double[vectorLength];

            for (int i = 0; i < vectorLength; i++)
            {
                resultVector[i] = vector1[i] - vector2[i];
                /*Console.Write(resultVector[i] + " ");*/
            }

            return resultVector;
        }

        /// <summary>
        /// Test the AddVectors method by providing sample input and comparing the computed result with the expected result.
        /// </summary>
        public static void TestAddVectors()
        {
            double[] x = { 1.0, 2.0, 5.7 };
            double[] y = { 1.0, 4.0, -3.0 };
            double[] expected = { 2.0, 6.0, 2.7 };

            double[] z = AddVectors(x, y);

            VectorMethods.CheckVector(expected, z);
        }

        /// <summary>
        /// Test the SubtractVectors method by providing sample input and comparing the computed result with the expected result.
        /// </summary>
        public static void TestSubtractVectors()
        {
            double[] x = { 1.0, 2.0, 5.7 };
            double[] y = { 1.0, 4.0, -3.0 };
            double[] expected = { 0.0, -2.0, 8.7 };

            double[] z = SubtractVectors(x, y);

            VectorMethods.CheckVector(expected, z);
        }
    }
}
