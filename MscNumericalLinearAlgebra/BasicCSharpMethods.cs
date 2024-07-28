using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MscNumericalLinearAlgebra
{
    internal class BasicCSharpMethods
    {
        public static double Add(double x, double y)
        {
            double z = x + y;
            return z;
        }
        
        public static int Sign(double x)
        {
            if (x < 0)
            {
                return -1;
            }
            else if (x == 0) 
            { 
                return 0;
            }
            else
            {
                return +1;
            }
        }

        public static void Print(double[] myArray)
        {
            for (int i = 0; i < myArray.Length; i++)
            {
                Console.Write(myArray[i] + " ");
            }
        }

        public static void Print2D(double[,] myArray)
        {
            Console.WriteLine();
            for (int i = 0; i < myArray.GetLength(0); i++)
            {
                

                for (int j = 0; j < myArray.GetLength(1); j++)
                {
                    Console.Write(myArray[i, j] + " ");
                }
                Console.WriteLine();
            }
        }

        public static void WhileExample(int[] myArray)
        {
            int i = 0;
            while(i < myArray.Length)
            {
                Console.WriteLine(myArray[i] + " ");
                i++;
            }
           
        }



        }
}
