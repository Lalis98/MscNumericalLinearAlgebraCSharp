// See https://aka.ms/new-console-template for more information
using MscNumericalLinearAlgebra;
using MscNumericalLinearAlgebra.ExcerciseSeries1;
using MscNumericalLinearAlgebra.ExcerciseSeries2;
using MscNumericalLinearAlgebra.Matrices;
using System.Drawing;
using System.Numerics;
using System.Runtime.Intrinsics;
using System.Xml.Linq;


double[,] matrixA =
{
    { 21.0, 1.0, 0.0, 4.0, 0.0 },
    { 1.0, 22.0, 2.0, 0.0, 0.0 },
    { 0.0, 2.0, 23.0, 1.0, 3.0 },
    { 4.0, 0.0, 1.0, 24.0, 2.0 },
    { 0.0, 0.0, 3.0, 2.0, 25.0 }
};

double[] vectorB = { 1.0, 1.0, 2.0, 3.0, 1.0 };
double epsilon = 0.001;
int n = 20;

/*double[] values = { 21, 1, 4, 1, 22, 2, 2, 23, 1, 3, 4, 1, 24, 2, 3, 2, 25 };
int[] colIndices = { 0, 1, 3, 0, 1, 2, 1, 2, 3, 4, 0, 2, 3, 4, 2, 3, 4 };
int[] rowOffsets = { 0, 3, 6, 10, 14, 17 };*/

double[] values = { 21, 1, 4, 1, 22, 2, 2, 23, 1, 3, 4, 1, 24, 2, 3, 2, 25 };
int[] rowIndices = { 0, 1, 3, 0, 1, 2, 1, 2, 3, 4, 0, 2, 3, 4, 2, 3, 4 };
int[] colOffsets = { 0, 3, 6, 10, 14, 17 };

/*double[] values = { 21, 1, 4, 1, 22, 2, 2, 23, 1, 3, 4, 1, 24, 2, 3, 2, 25 };
int[] rows = { 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4 };
int[] cols = { 0, 1, 3, 0, 1, 2, 1, 2, 3, 4, 0, 2, 3, 4, 2, 3, 4};*/

double[,] matrixA1 =
{
    { 21.0, 1.0, 0.0, 4.0, 0.0 },
    { 1.0, 22.0, 2.0, 0.0, 0.0 },
    { 0.0, 2.0, 23.0, 1.0, 3.0 },
    { 4.0, 0.0, 1.0, 24.0, 2.0 },
    { 0.0, 0.0, 3.0, 2.0, 25.0 }
};

double[,] matrixB1 =
{
    { 1.0, 1.0, 0.0 },
    { 1.0, 2.0, 2.0 },
    { 0.0, 2.0, 3.0 },
    { 1.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 }
};


double[,] arrayA2 = { 
    { 5, 1, 3, 0 },
    { 0, 7, 0, 0 },
    { 0, 0, 3, 0 },
    { 1, 0, 0, 2 } 
};
double[] arrayB2 = { 1, 2, 3, 4 };


double[] val = { 5, 1, 3, 7, 3, 1, 2 };
int[] ro = { 0, 0, 0, 1, 2, 3, 3 };
int[] co = { 0, 1, 2, 1, 2, 0, 3 };


double[] vector_b = { 1, 2, 3, 4, 5, 6 };
int size = 6;
int bandwidth = 2;
double[] elements = { 1, 2, 5, 1, 2, 3, 0, 3, 2, 3, 4, -1, 1, 3, 4, 5, 1, 2, 4, 5, 0, 4, 5, 1 };

/*MscNumericalLinearAlgebra.DirectMethods.DirectMethods.SolveLinearSystemWithGaussElimination(matrixA, vectorB);
MscNumericalLinearAlgebra.IterativeMethods.StationaryPointMethods.Jacobi(matrixA, vectorB, epsilon);*/