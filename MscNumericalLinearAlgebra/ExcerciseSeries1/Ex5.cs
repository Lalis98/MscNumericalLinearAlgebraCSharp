using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries1
{
    public static class Ex5
    {
        /// <summary>
        /// Calculates the second norm and returns the result.
        /// </summary>
        /// <param name="vector">The input vector.</param>
        /// <returns>
        /// The second norm of a vector.
        /// </returns>
        /// <exception cref="Exception">
        /// Thrown when the input vector is not a vector.
        /// </exception>
        public static double Norm2Vector(double[] vector)
        {
            if (!(vector is double[]))
            {
                throw new Exception("The input is not a vector");
            }

            double norm = 0.0;

            for (int i = 0; i < vector.Length; i++)
            {
                norm += vector[i] * vector[i];
            }


            return Math.Sqrt(norm);

        }

        /// <summary>
        /// Test the Norm2Vector method by providing sample input and comparing the computed result with the expected result.
        /// </summary>
        public static void TestNorm2Vector()
        {
            double[] vector = { 6.0, 8.0 };
            double computed = Norm2Vector(vector);
            double expected = 10.0;

            Matrices.ScalarMethods.CheckScalar(computed, expected);

        }
    }
}
