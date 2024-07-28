using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries1
{
    public class Ex2
    {
        /// <summary>
        /// Multiplies a vector with a scalar and returns the result.
        /// </summary>
        /// <param name="vector">The first input vector.</param>
        /// <param name="scalar">The second input vector.</param>
        /// <returns>
        /// A new vector resulted from the input vector and scalar.
        /// </returns>
        public static double[] MultiplyVectorWithScalar(double[] vector, double scalar)
        {
            double[] resultVector = new double[vector.Length];

            for (int i = 0; i < vector.Length; i++)
            {
                resultVector[i] = vector[i] * scalar;
                Console.WriteLine(resultVector[i] + " ");
            }
            return resultVector;

        }

        /// <summary>
        /// Test the MultiplyVectorWithScalar method by providing sample input and comparing the computed result with the expected result.
        /// </summary>
        public static void TestMultiplyVectorWithScalar()
        {
            double[] vector = { 1.0, 2.0, 3.0 };
            double scalar = 3.0;

            double[] expected = { 3.0, 6.0, 9.0 };
            double[] computed = MultiplyVectorWithScalar(vector, scalar);

            Matrices.VectorMethods.CheckVector(expected, computed);

        }
    }
}
