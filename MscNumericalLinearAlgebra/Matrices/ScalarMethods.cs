using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra.Matrices
{
    public static class ScalarMethods
    {
        public static void CheckScalar(double expected, double computed)
        {
            if (expected != computed)
            {
                Console.WriteLine($"The computed and expected are not equal," +
                    $" expected = {expected}, computed = {computed}");
            }
        }
    }
}
