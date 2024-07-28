using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.ExcerciseSeries1
{
    public static class Ex3
    {
        /// <summary>
        /// Multiplies two vectors with a scalars and returns the result.
        /// </summary>
        /// <param name="vector1">The first input vector.</param>
        /// <param name="scalar1">The second input vector.</param>
        /// <param name="vector2">The first input vector.</param>
        /// <param name="scalar2">The second input vector.</param>
        /// <returns>
        /// A new vector resulted from the input vectors and scalars.
        /// </returns>
        public static double[] AddLinearVectors(double[] vector1, double scalar1, double[] vector2, double scalar2)
        {

            if (vector1.Length != vector2.Length)
            {
                throw new Exception("The two vectors length is not equal");
            }

            int vectorLength = vector1.Length;
            double[] resultVector = new double[vectorLength];

            for (int i = 0; i < vectorLength; i++)
            {
                resultVector[i] = vector1[i] * scalar1 + vector2[i] * scalar2;
            }

            return resultVector;

        }

        /// <summary>
        /// Test the AddLinearVectors method by providing sample input and comparing the computed result with the expected result.
        /// </summary>
        public static void TestAddLinearVectors()
        {
            double[] vector1 = { 5.0, 5.0, 5.0 };
            double scalar1 = 1.0;
            double[] vector2 = { 1.0, 1.0, 1.0 };
            double scalar2 = 5.0;

            double[] computed = AddLinearVectors(vector1, scalar1, vector2, scalar2);
            double[] expected = { 10.0, 10.0, 10.0 };

            Matrices.VectorMethods.CheckVector(computed, expected);

        }
    }
}
