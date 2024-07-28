using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries1
{
    public static class Ex4
    {
        /// <summary>
        /// Calculates the vector product between 2 vectors and returns the result.
        /// </summary>
        /// <param name="vector1">The first input vector.</param>
        /// <param name="vector2">The second input vector.</param>
        /// <returns>
        /// The second dot product of two vectors.
        /// </returns>
        /// <exception cref="Exception">
        /// Thrown when the input vectors have different lengths.
        /// </exception>
        public static double VectorDotProduct(double[] vector1, double[] vector2) 
        {
            if (vector1.Length != vector2.Length)
            {
                throw new Exception("The length of the two vectors is not equal!");
            }

            int vectorLength = vector1.Length;
            double scalarResult = 0.0;

            for (int i = 0; i < vectorLength; i++)
            {
                scalarResult += vector1[i] * vector2[i];
            }

            return scalarResult;
        }

        /// <summary>
        /// Test the VectorDotProduct method by providing sample input and comparing the computed result with the expected result.
        /// </summary>
        public static void TestVectorDotProduct()
        {
            double[] vector1 = { 1.0, 2.0, 3.0 };
            double[] vector2 = { 1.0, 2.0, 3.0 };

            double computed = VectorDotProduct(vector1, vector2);
            double expected = 14.0;

            Matrices.ScalarMethods.CheckScalar(expected, computed);
            
        }
    }


}
