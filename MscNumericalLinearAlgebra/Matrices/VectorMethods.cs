using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.Matrices
{
    public static class VectorMethods
    {
        // Function to print a vector with enhanced formatting and rounded to 2 decimals
        public static void PrintVector(double[] vector)
        {
            Console.Write("[ ");
            for (int i = 0; i < vector.Length; i++)
            {
                Console.Write(vector[i].ToString("N2")); // Round to 2 decimals
                if (i < vector.Length - 1)
                {
                    Console.Write(", ");
                }
            }
            Console.WriteLine(" ]\n");
        }

        public static void PrintVectorInt(int[] vector)
        {
            Console.Write("[ ");
            for (int i = 0; i < vector.Length; i++)
            {
                Console.Write(vector[i].ToString("N2")); // Round to 2 decimals
                if (i < vector.Length - 1)
                {
                    Console.Write(", ");
                }
            }
            Console.WriteLine(" ]\n");
        }

        public static void CheckVector(double[] expected, double[] computed)
        {
            for (int i = 0; i < computed.Length; i++)
            {
                if (computed[i] != expected[i])
                {
                    Console.WriteLine($"The vectors are not equal: i = {i}" +
                        $", computed = {computed[i]} expected = {expected[i]} !");
                }

            }
        }

        public static double VectorEuclideanNorm(double[] vector)
        {
            double internal_sum = 0;

            for (int i = 0; i < vector.Length; i++)
            { 
                internal_sum = internal_sum + vector[i] * vector[i];
            }

            return Math.Sqrt(internal_sum);
        }

        public static double VectorDotProduct(double[] vectorA, double[] vectorB) 
        {
            if (vectorA.Length != vectorB.Length) 
            {
                throw new Exception("Wrong Dimensions of the two vectors dot products!");
            }

            double sum = 0;
            for (int i = 0; (i < vectorA.Length); i++)
            {
                sum = sum + vectorA[i] * vectorB[i];
            }

            return sum;
        }

    }
}
