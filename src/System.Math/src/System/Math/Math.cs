using System;

using System.Collections.Generic;

using System.Text;





/// <summary>

/// Author: Steven Ferrell,

/// Date: June 2009

/// </summary>

namespace MAGE

{

    /// <summary>

    /// Templated math class containing math functions and math structures

    /// </summary>

    public static class Math

    {

        #region Math Functions

        public static class Functions

        {

            public static T Max<T>(params T[] p)

            {

                GenericMathSolution.ICalculator<T> Calculator = GenericMathSolution.CalculatorFactory<T>.Calculator;

                T max = p[0];

                foreach (T val in p)

                    if (Calculator.Less(max, val))

                        max = val;

                return max;

            }



            public static T Min<T>(params T[] p)

            {

                GenericMathSolution.ICalculator<T> Calculator = GenericMathSolution.CalculatorFactory<T>.Calculator;

                T min = p[0];

                foreach (T val in p)

                    if (Calculator.Greater(min, val))

                        min = val;

                return min;

            }



            public static T DegreeToRadian<T>(T d)

            {

                GenericMathSolution.ICalculator<T> Calculator = GenericMathSolution.CalculatorFactory<T>.Calculator;

                return Calculator.Div(d, Calculator.Div(Calculator.NUM(180.0), Calculator.PI));

            }



            public static T RadianToDegree<T>(T r)

            {

                GenericMathSolution.ICalculator<T> Calculator = GenericMathSolution.CalculatorFactory<T>.Calculator;

                return Calculator.Mult(r, Calculator.Div(Calculator.NUM(180.0), Calculator.PI));

            }

        }

        #endregion



        #region Math Structures

        public static class Structures

        {

            /// <summary>

            /// Templated Vector Class (Row Vector)

            /// </summary>

            /// <typeparam name="T">The type of numbers in the Vector</typeparam>

            public class Vector<T>

            {

                #region Data Members

                private T[] _values;

                private Vector<T> ZERO_VECTOR;

                private static GenericMathSolution.ICalculator<T> Calculator = GenericMathSolution.CalculatorFactory<T>.Calculator;

                #endregion



                #region Constructors

                /// <summary>

                /// Templated Vector Constructor

                /// </summary>

                /// <param name="p">Parameterized Arguments defining the numbers in the Vector</param>

                public Vector(params T[] p)

                {

                    _values = new T[p.Length];

                    for (int i = 0; i < p.Length; i++)

                        _values[i] = (T)p.GetValue(i);



                    ZERO_VECTOR = new Vector<T>(p.Length);

                }



                /// <summary>

                /// Templated Vector Copy Constructor

                /// </summary>

                /// <param name="v">Vector to copy</param>

                public Vector(Vector<T> v)

                {

                    _values = new T[v.Dimensions];

                    for (int i = 0; i < v.Dimensions; i++)

                        _values[i] = v[i];



                    ZERO_VECTOR = new Vector<T>(v.Dimensions);

                }



                public Vector(Matrix<T> m)

                {

                    _values = new T[m.Columns];

                    for (int i = 0; i < this.Dimensions; i++)

                        _values[i] = m[0, i];



                    ZERO_VECTOR = new Vector<T>(this.Dimensions);

                }



                public Vector(int dimensions)

                {

                    _values = new T[dimensions];

                    for (int i = 0; i < dimensions; i++)

                        _values[i] = Calculator.NUM(0.0);

                }



                /// <summary>

                /// Templated Vector Copy Constructor Facade

                /// </summary>

                /// <returns>Returns Deep Copy of this Vector</returns>

                public Vector<T> Copy()

                {

                    return new Vector<T>(this);

                }

                #endregion



                #region Vector Operations

                /// <summary>

                /// Dot Product Function

                /// </summary>

                /// <param name="rhs">Right Hand Side Vector</param>

                /// <returns>Returns the dot product of this Vector against the right hand side vector</returns>

                public T Dot(Vector<T> rhs)

                {

                    if (this.Dimensions != rhs.Dimensions)

                        throw new Exception("Invalid Dimensions To Compute Dot Product: Must provide Vectors of equal dimensions");



                    T product = Calculator.NUM(0.0);

                    for (int i = 0; i < this.Dimensions; i++)

                        product = Calculator.Add(product, Calculator.Mult(this[i], rhs[i]));



                    return product;

                }



                /// <summary>

                /// Cross Product Function (supports 3D and 7D Vectors)

                /// </summary>

                /// <param name="rhs">Right Hand Side Vector</param>

                /// <returns>Returns the cross product of this Vector against the right hand side vectors</returns>

                public Vector<T> Cross(Vector<T> rhs)

                {

                    if (this.Dimensions == 3 && rhs.Dimensions == 3)

                    {

                        return new Vector<T>(Calculator.Sub(Calculator.Mult(this[1], rhs[2]), Calculator.Mult(this[2], rhs[1])),

                                             Calculator.Sub(Calculator.Mult(this[2], rhs[0]), Calculator.Mult(this[0], rhs[2])),

                                             Calculator.Sub(Calculator.Mult(this[0], rhs[1]), Calculator.Mult(this[1], rhs[0])));

                    }

                    else if (this.Dimensions == 7 && rhs.Dimensions == 7)

                    {

                        return new Vector<T>(

                            Calculator.Add(

                                Calculator.Add(Calculator.Sub(Calculator.Mult(this[1], rhs[3]), Calculator.Mult(this[3], rhs[1])),

                                               Calculator.Sub(Calculator.Mult(this[2], rhs[6]), Calculator.Mult(this[6], rhs[2]))),

                                               Calculator.Sub(Calculator.Mult(this[4], rhs[5]), Calculator.Mult(this[5], rhs[4]))),

                            Calculator.Add(

                                Calculator.Add(Calculator.Sub(Calculator.Mult(this[2], rhs[4]), Calculator.Mult(this[4], rhs[2])),

                                               Calculator.Sub(Calculator.Mult(this[3], rhs[0]), Calculator.Mult(this[0], rhs[3]))),

                                               Calculator.Sub(Calculator.Mult(this[5], rhs[6]), Calculator.Mult(this[6], rhs[5]))),

                            Calculator.Add(

                                Calculator.Add(Calculator.Sub(Calculator.Mult(this[3], rhs[5]), Calculator.Mult(this[5], rhs[3])),

                                               Calculator.Sub(Calculator.Mult(this[4], rhs[1]), Calculator.Mult(this[1], rhs[4]))),

                                               Calculator.Sub(Calculator.Mult(this[6], rhs[0]), Calculator.Mult(this[0], rhs[6]))),

                            Calculator.Add(

                                Calculator.Add(Calculator.Sub(Calculator.Mult(this[4], rhs[6]), Calculator.Mult(this[6], rhs[4])),

                                               Calculator.Sub(Calculator.Mult(this[5], rhs[2]), Calculator.Mult(this[2], rhs[5]))),

                                               Calculator.Sub(Calculator.Mult(this[0], rhs[1]), Calculator.Mult(this[1], rhs[0]))),

                            Calculator.Add(

                                Calculator.Add(Calculator.Sub(Calculator.Mult(this[5], rhs[0]), Calculator.Mult(this[0], rhs[5])),

                                               Calculator.Sub(Calculator.Mult(this[6], rhs[3]), Calculator.Mult(this[3], rhs[6]))),

                                               Calculator.Sub(Calculator.Mult(this[1], rhs[2]), Calculator.Mult(this[2], rhs[1]))),

                            Calculator.Add(

                                Calculator.Add(Calculator.Sub(Calculator.Mult(this[6], rhs[1]), Calculator.Mult(this[1], rhs[6])),

                                               Calculator.Sub(Calculator.Mult(this[0], rhs[4]), Calculator.Mult(this[4], rhs[0]))),

                                               Calculator.Sub(Calculator.Mult(this[2], rhs[3]), Calculator.Mult(this[3], rhs[2]))),

                            Calculator.Add(

                                Calculator.Add(Calculator.Sub(Calculator.Mult(this[0], rhs[2]), Calculator.Mult(this[2], rhs[0])),

                                               Calculator.Sub(Calculator.Mult(this[1], rhs[5]), Calculator.Mult(this[5], rhs[1]))),

                                               Calculator.Sub(Calculator.Mult(this[3], rhs[4]), Calculator.Mult(this[4], rhs[3]))));

                    }



                    throw new Exception("Invalid Dimensions To Compute Cross Product: Must provide 3D or 7D Vectors");

                }



                /// <summary>

                /// Scalar Triple Function

                /// </summary>

                /// <param name="b">Second Vector</param>

                /// <param name="c">Third Vector</param>

                /// <returns>Returns the scalar triple product of three vectors</returns>

                public T ScalarTriple(Vector<T> b, Vector<T> c)

                {

                    return this * (b % c);

                }



                /// <summary>

                /// Angle Between Vectors (supports 2D and 3D Vectors)

                /// </summary>

                /// <param name="rhs">Right Hand Side Vector</param>

                /// <returns>Returns the angle between this vector and the right hand side vector</returns>

                public T Angle(Vector<T> rhs)

                {

                    if (this.Dimensions == 2 && rhs.Dimensions == 2 || this.Dimensions == 3 && rhs.Dimensions == 3)

                        return Functions.RadianToDegree(Calculator.NUM(Calculator.Acos(this.UnitVector * rhs.UnitVector)));

                    else

                        throw new Exception("Invalid Dimensions To Compute Angle: Must provide 2D or 3D Vectors");

                }



                /// <summary>

                /// Distance Between Vectors

                /// </summary>

                /// <param name="rhs">Right Hand Side Vector</param>

                /// <returns>Returns the distance between this vector and the right hand side vector</returns>

                public T Distance(Vector<T> rhs)

                {

                    if (this.Dimensions != rhs.Dimensions)

                        throw new Exception("Invalid Dimensions To Compute Distance: Must provide Vectors of equal dimensions");



                    T sum = Calculator.NUM(0.0);

                    for (int i = 0; i < this.Dimensions; i++)

                        sum = Calculator.Add(sum, Calculator.Pow(Calculator.Sub(this[i], rhs[i]), 2.0));



                    return Calculator.Pow(sum, 0.5);

                }



                /// <summary>

                /// Magnitude of this Vector

                /// </summary>

                public T Magnitude

                {

                    get

                    {

                        T product = Calculator.NUM(0.0);

                        for (int i = 0; i < this.Dimensions; i++)

                            product = Calculator.Add(product, Calculator.Pow(this[i], 2.0));



                        return Calculator.Pow(product, 0.5);

                    }

                }



                /// <summary>

                /// Unit Vector of this Vector

                /// </summary>

                public Vector<T> UnitVector

                {

                    get { return this / this.Magnitude; }

                }



                /// <summary>

                /// Number of Dimensions in this Vector

                /// </summary>

                public int Dimensions

                {

                    get { return _values.Length; }

                }

                #endregion



                #region Operator Overloads

                public static Vector<T>

                    operator +(Vector<T> lhs, Vector<T> rhs)

                {

                    Vector<T> sum = new Vector<T>(lhs._values);

                    for (int i = 0; i < lhs.Dimensions; i++)

                        sum[i] = Calculator.Add(lhs[i], rhs[i]);

                    return sum;

                }



                public static Vector<T>

                    operator -(Vector<T> lhs, Vector<T> rhs)

                {

                    Vector<T> diff = new Vector<T>(lhs._values);

                    for (int i = 0; i < lhs.Dimensions; i++)

                        diff[i] = Calculator.Sub(lhs[i], rhs[i]);

                    return diff;

                }



                public static Vector<T>

                    operator -(Vector<T> v)

                {

                    Vector<T> negate = new Vector<T>(v._values);

                    for (int i = 0; i < v.Dimensions; i++)

                        negate[i] = Calculator.Negate(v[i]);



                    return negate;

                }



                public static Vector<T>

                    operator *(Vector<T> lhs, T scalar)

                {

                    Vector<T> val = new Vector<T>(lhs._values);

                    for (int i = 0; i < lhs.Dimensions; i++)

                        val[i] = Calculator.Mult(lhs[i], scalar);

                    return val;

                }



                public static Vector<T>

                    operator *(T scalar, Vector<T> rhs)

                {

                    Vector<T> val = new Vector<T>(rhs._values);

                    for (int i = 0; i < rhs.Dimensions; i++)

                        val[i] = Calculator.Mult(rhs[i], scalar);

                    return val;

                }



                public static T

                    operator *(Vector<T> lhs, Vector<T> rhs)

                {

                    return lhs.Dot(rhs);

                }



                public static Vector<T>

                    operator *(Vector<T> lhs, Matrix<T> rhs)

                {

                    if (lhs.Dimensions == rhs.Rows)

                        return new Vector<T>(new Matrix<T>(lhs) * rhs);

                    else

                        return null;

                }



                public static Vector<T>

                    operator %(Vector<T> lhs, Vector<T> rhs)

                {

                    return lhs.Cross(rhs);

                }



                public static Vector<T>

                    operator /(Vector<T> lhs, T scalar)

                {

                    Vector<T> val = new Vector<T>(lhs._values);

                    for (int i = 0; i < lhs.Dimensions; i++)

                        val[i] = Calculator.Div(lhs[i], scalar);

                    return val;

                }



                public static bool

                    operator ==(Vector<T> lhs, Vector<T> rhs)

                {

                    if (lhs.Dimensions != rhs.Dimensions)

                        return false;



                    for (int i = 0; i < lhs.Dimensions; i++)

                        if (Calculator.NotEqual(lhs[i], rhs[i]))

                            return false;



                    return true;

                }



                public static bool

                    operator !=(Vector<T> lhs, Vector<T> rhs)

                {

                    if (lhs.Dimensions != rhs.Dimensions)

                        return true;



                    for (int i = 0; i < lhs.Dimensions; i++)

                        if (Calculator.NotEqual(lhs[i], rhs[i]))

                            return true;



                    return false;

                }



                #region Implicit Conversions

                public static implicit operator Vector<T>(Matrix<T> m)

                {

                    return new Vector<T>(m);

                }



                public static implicit operator Vector<T>(Vector<double> v)

                {

                    Vector<T> t = new Vector<T>(v.Dimensions);

                    for (int i = 0; i < v.Dimensions; i++)

                        t[i] = Calculator.NUM(v[i]);



                    return t;

                }



                public static implicit operator Vector<T>(Vector<float> v)

                {

                    Vector<T> t = new Vector<T>(v.Dimensions);

                    for (int i = 0; i < v.Dimensions; i++)

                        t[i] = Calculator.NUM(v[i]);



                    return t;

                }



                public static implicit operator Vector<T>(Vector<int> v)

                {

                    Vector<T> t = new Vector<T>(v.Dimensions);

                    for (int i = 0; i < v.Dimensions; i++)

                        t[i] = Calculator.NUM(v[i]);



                    return t;

                }



                public static implicit operator Vector<T>(Vector<Int64> v)

                {

                    Vector<T> t = new Vector<T>(v.Dimensions);

                    for (int i = 0; i < v.Dimensions; i++)

                        t[i] = Calculator.NUM(v[i]);



                    return t;

                }

                #endregion

                #endregion



                #region Object Overloads

                public override int GetHashCode()

                {

                    return base.GetHashCode();

                }



                public override bool Equals(object obj)

                {

                    if (obj is Vector<T>)

                    {

                        Vector<T> v = (Vector<T>)obj;



                        for (int i = 0; i < this.Dimensions; i++)

                            if (Calculator.NotEqual(v[i], this[i]))

                                return false;



                        return obj.GetType().Equals(this.GetType());

                    }

                    return false;

                }



                public override string ToString()

                {

                    string s = "";

                    for (int i = 0; i < this.Dimensions; i++)

                        s += this[i] + ",";

                    s = s.Remove(s.LastIndexOf(','));



                    return s;

                }

                #endregion



                #region Indexer

                private T Find(int key)

                {

                    return _values[key];

                }



                public T this[int key]

                {

                    get { return Find(key); }

                    set { SetValue(key, value); }

                }



                private void SetValue(int key, T value)

                {

                    _values[key] = value;

                }

                #endregion

            }



            public class Matrix<T>

            {

                #region Data Members

                private int _rows;

                private int _columns;

                private T[,] _values;

                private static GenericMathSolution.ICalculator<T> Calculator = GenericMathSolution.CalculatorFactory<T>.Calculator;

                #endregion



                #region Enumerations

                public enum MatrixType

                {

                    Homogeneous,

                    Regular

                };

                #endregion



                #region Constructors

                public Matrix(int rows, int columns, params T[] p)

                {

                    _rows = rows;

                    _columns = columns;

                    _values = new T[_rows, _columns];

                    for (int i = 0; i < _rows; i++)

                        for (int j = 0; j < _columns; j++)

                            _values[i, j] = (T)p.GetValue(i * _columns + j);

                }



                public Matrix(Matrix<T> m)

                {

                    _rows = m._rows;

                    _columns = m._columns;

                    _values = new T[_rows, _columns];

                    for (int i = 0; i < _rows; i++)

                        for (int j = 0; j < _columns; j++)

                            _values[i, j] = m[i, j];

                }



                public Matrix(Vector<T> v)

                {

                    _rows = 1;

                    _columns = v.Dimensions;

                    _values = new T[_rows, _columns];

                    for (int i = 0; i < _columns; i++)

                            _values[0, i] = v[i];

                }



                public Matrix(int rows, int columns)

                {

                    _rows = rows;

                    _columns = columns;

                    _values = new T[_rows, _columns];

                    for (int i = 0; i < _rows; i++)

                        for (int j = 0; j < _columns; j++)

                            _values[i, j] = Calculator.NUM(0.0);

                }



                public Matrix<T> Copy()

                {

                    return new Matrix<T>(this);

                }

                #endregion



                #region Matrix Operations

                public int Rows

                {

                    get { return _rows; }

                }



                public int Columns

                {

                    get { return _columns; }

                }



                public Matrix<T> IdentityMatrix

                {

                    get

                    {

                        if (this.Rows == this.Columns)

                        {

                            Matrix<T> m = this.Copy();

                            for (int i = 0; i < _rows; i++)

                                for (int j = 0; j < _columns; j++)

                                {

                                    if (i == j)

                                        m[i, j] = Calculator.NUM(1.0);

                                    else

                                        m[i, j] = Calculator.NUM(0.0);

                                }

                            return m;

                        }

                        else

                            return null;

                    }

                }



                public Matrix<T> Transpose

                {

                    get

                    {

                        Matrix<T> t = new Matrix<T>(_columns, _rows);

                        for (int i = 0; i < _columns; i++)

                            for (int j = 0; j < _rows; j++)

                                t[i, j] = this[j, i];

                        return t;

                    }

                }



                public void Clear()

                {

                    for (int i = 0; i < this.Rows; i++)

                        for (int j = 0; j < this.Columns; j++)

                            this[i, j] = Calculator.NUM(0);

                }



                public bool IsEmpty

                {

                    get

                    {

                        for (int i = 0; i < this.Rows; i++)

                            for (int j = 0; j < this.Columns; j++)

                                if (Calculator.NotEqual(this[i, j], Calculator.NUM(0)))

                                    return false;



                        return true;

                    }

                }



                public Matrix<T> Diagonal

                {

                    get

                    {

                        if (this.Rows == this.Columns)

                        {

                            Matrix<T> m = new Matrix<T>(this.Rows, this.Columns);

                            for (int i = 0; i < this.Rows; i++)

                                for (int j = 0; j < this.Columns; j++)

                                    if (i == j)

                                        m[i, j] = this[i, j];

                                    else

                                        m[i, j] = Calculator.NUM(0);



                            return m;

                        }

                        else

                            return null;

                    }

                }



                public Matrix<T> Power(double power)

                {

                    Matrix<T> m = new Matrix<T>(this);

                    for (int i = 0; i < this.Rows; i++)

                        for (int j = 0; j < this.Columns; j++)

                            m[i, j] = Calculator.Pow(m[i, j], power);

                    return m;

                }



                public Matrix<T> Sum

                {

                    get

                    {

                        if (this.Rows == 1)

                            return new Matrix<T>(this);

                        else if (this.Rows == 0)

                            return null;



                        Matrix<T> m = new Matrix<T>(1, this.Columns);

                        for (int i = 0; i < this.Columns; i++)

                            for (int j = 0; j < this.Rows; j++)

                                m[0, i] = Calculator.Add(m[0, i], this[j, i]);

                        return m;

                    }

                }



                /// <summary>

                /// Rotate a regular or homogeneous matrix in 3D around an arbitrary axis using a provided angle.

                /// </summary>

                /// <param name="angle">The amount of rotation about the axis.</param>

                /// <param name="axis">The 3D Axis Vector to rotate about.</param>

                /// <param name="matrixType">The type of matrix used (Homogeneous or Regular)</param>

                /// <returns>Returns the rotated matrix.</returns>

                public Matrix<T> Rotate(double angle, Vector<double> axis, MatrixType matrixType)

                {

                    Matrix<T> m = new Matrix<T>(this);

                    double c = System.Math.Cos(Functions.DegreeToRadian<double>(angle));

                    double s = System.Math.Sin(Functions.DegreeToRadian<double>(angle));

                    Vector<double> n = axis.UnitVector;



                    if (matrixType == MatrixType.Regular)

                    {

                        if (this.Rows == 3 && this.Columns == 3)

                            return m * new Matrix<T>(3, 3,

                                Calculator.NUM((n[0] * n[0]) * (1.0 - c) + c),

                                    Calculator.NUM((n[0] * n[1]) * (1.0 - c) + (n[2] * s)),

                                    Calculator.NUM((n[0] * n[2]) * (1.0 - c) - (n[0] * s)),

                                Calculator.NUM((n[0] * n[1]) * (1.0 - c) - (n[2] * s)),

                                    Calculator.NUM((n[1] * n[1]) * (1.0 - c) + c),

                                    Calculator.NUM((n[1] * n[0]) * (1.0 - c) + (n[0] * s)),

                                Calculator.NUM((n[0] * n[2]) * (1.0 - c) + (n[1] * s)),

                                    Calculator.NUM((n[1] * n[2]) * (1.0 - c) - (n[0] * s)),

                                    Calculator.NUM((n[2] * n[2]) * (1.0 - c) + c));

                        else

                            return null;

                    }

                    else if(matrixType == MatrixType.Homogeneous)

                    {

                        if (this.Rows == 4 && this.Columns == 4)

                            return m * new Matrix<T>(4, 4,

                                Calculator.NUM((n[0] * n[0]) * (1.0 - c) + c),

                                    Calculator.NUM((n[0] * n[1]) * (1.0 - c) + (n[2] * s)),

                                    Calculator.NUM((n[0] * n[2]) * (1.0 - c) - (n[0] * s)),

                                    Calculator.NUM(0),

                                Calculator.NUM((n[0] * n[1]) * (1.0 - c) - (n[2] * s)),

                                    Calculator.NUM((n[1] * n[1]) * (1.0 - c) + c),

                                    Calculator.NUM((n[1] * n[0]) * (1.0 - c) + (n[0] * s)),

                                    Calculator.NUM(0),

                                Calculator.NUM((n[0] * n[2]) * (1.0 - c) + (n[1] * s)),

                                    Calculator.NUM((n[1] * n[2]) * (1.0 - c) - (n[0] * s)),

                                    Calculator.NUM((n[2] * n[2]) * (1.0 - c) + c),

                                    Calculator.NUM(0),

                                Calculator.NUM(0),

                                    Calculator.NUM(0),

                                    Calculator.NUM(0),

                                    Calculator.NUM(1));

                        else

                            return null;

                    }

                    else

                        return null;

                }



                /// <summary>

                /// Rotate a regular or homogeneous matrix counter-clockwise in 2D using a provided angle.

                /// </summary>

                /// <param name="angle">The amount of rotation about the axis.</param>

                /// <param name="matrixType">The type of matrix used (Homogeneous or Regular)</param>

                /// <returns>Returns the rotated matrix.</returns>

                public Matrix<T> Rotate(double angle, MatrixType matrixType)

                {

                    Matrix<T> m = new Matrix<T>(this);

                    double c = System.Math.Cos(Functions.DegreeToRadian<double>(angle));

                    double s = System.Math.Sin(Functions.DegreeToRadian<double>(angle));



                    if (matrixType == MatrixType.Regular)

                    {

                        if (this.Rows == 2 && this.Columns == 2)

                            return m * new Matrix<T>(2, 2,

                                Calculator.NUM(c), Calculator.NUM(-s),

                                Calculator.NUM(s), Calculator.NUM(c));

                        else

                            return null;

                    }

                    else if (matrixType == MatrixType.Homogeneous)

                    {

                        if (this.Rows == 3 && this.Columns == 3)

                            return m * new Matrix<T>(3, 3,

                                Calculator.NUM(c), Calculator.NUM(-s), Calculator.NUM(0),

                                Calculator.NUM(s), Calculator.NUM(c), Calculator.NUM(0),

                                Calculator.NUM(0), Calculator.NUM(0), Calculator.NUM(1));

                        else

                            return null;

                    }

                    else

                        return null;

                }



                /// <summary>

                /// Translate a homogenous 2D or 3D matrix using the provided translation vector.

                /// </summary>

                /// <param name="t">2D or 3D translation vector.</param>

                /// <returns>Returns the translated matrix.</returns>

                public Matrix<T> Translate(Vector<T> t)

                {

                    //must be homogeneous

                    Matrix<T> m = new Matrix<T>(this);



                    if (this.Rows == 3 && this.Columns == 3 && t.Dimensions == 2)

                    {

                        return m * new Matrix<T>(3, 3, Calculator.NUM(1), Calculator.NUM(0), Calculator.NUM(0),

                                                       Calculator.NUM(0), Calculator.NUM(1), Calculator.NUM(0),

                                                       Calculator.NUM(t[0]), Calculator.NUM(t[1]), Calculator.NUM(1));

                    }

                    else if (this.Rows == 4 && this.Columns == 4 && t.Dimensions == 3)

                    {

                        return m * new Matrix<T>(4, 4, Calculator.NUM(1), Calculator.NUM(0), Calculator.NUM(0), Calculator.NUM(0),

                                                       Calculator.NUM(0), Calculator.NUM(1), Calculator.NUM(0), Calculator.NUM(0),

                                                       Calculator.NUM(0), Calculator.NUM(0), Calculator.NUM(1), Calculator.NUM(0),

                                                       Calculator.NUM(t[0]), Calculator.NUM(t[1]), Calculator.NUM(t[2]), Calculator.NUM(1));

                    }

                    else

                        return null;

                }



                /// <summary>

                /// Scales a regular or homogeneous 2D or 3D matrix using the provided scale vector.

                /// </summary>

                /// <param name="s">2D or 3D scale vector</param>

                /// <param name="matrixType">Type of matrix (Regular or Homogeneous)</param>

                /// <returns>Returns the scaled matrix.</returns>

                public Matrix<T> Scale(Vector<T> s, MatrixType matrixType)

                {

                    Matrix<T> m = new Matrix<T>(this);



                    if (matrixType == MatrixType.Regular)

                    {

                        if (this.Rows == 2 && this.Columns == 2)

                            return m * new Matrix<T>(2, 2, s[0], Calculator.NUM(0),

                                                           Calculator.NUM(0), s[1]);

                        else if (this.Rows == 3 && this.Columns == 3)

                            return m * new Matrix<T>(3, 3, s[0], Calculator.NUM(0), Calculator.NUM(0),

                                                           Calculator.NUM(0), s[1], Calculator.NUM(0),

                                                           Calculator.NUM(0), Calculator.NUM(0), s[2]);

                        else

                            return null;

                    }

                    else if (matrixType == MatrixType.Homogeneous)

                    {

                        if (this.Rows == 3 && this.Columns == 3)

                        {

                            return m * new Matrix<T>(3, 3, s[0], Calculator.NUM(0), Calculator.NUM(0),

                                                           Calculator.NUM(0), s[1], Calculator.NUM(0),

                                                           Calculator.NUM(0), Calculator.NUM(0), Calculator.NUM(1));

                        }

                        else if (this.Rows == 4 && this.Columns == 4)

                        {

                            return m * new Matrix<T>(4, 4, s[0], Calculator.NUM(0), Calculator.NUM(0), Calculator.NUM(0),

                                                           Calculator.NUM(0), s[1], Calculator.NUM(0), Calculator.NUM(0),

                                                           Calculator.NUM(0), Calculator.NUM(0), s[2], Calculator.NUM(0),

                                                           Calculator.NUM(0), Calculator.NUM(0), Calculator.NUM(0), Calculator.NUM(1));

                        }

                        else

                            return null;

                    }

                    else

                        return null;

                }

                #endregion



                #region Operator Overloads

                public static Matrix<T>

                    operator +(Matrix<T> lhs, Matrix<T> rhs)

                {

                    if (lhs.Rows == rhs.Rows && lhs.Columns == rhs.Columns)

                    {

                        Matrix<T> m = lhs.Copy();

                        for (int i = 0; i < lhs.Rows; i++)

                            for (int j = 0; j < lhs.Columns; j++)

                                m[i, j] = Calculator.Add(lhs[i, j], rhs[i, j]);

                        return m;

                    }

                    else

                        return null;

                }



                public static Matrix<T>

                    operator -(Matrix<T> lhs, Matrix<T> rhs)

                {

                    if (lhs.Rows == rhs.Rows && lhs.Columns == rhs.Columns)

                    {

                        Matrix<T> m = lhs.Copy();

                        for (int i = 0; i < lhs.Rows; i++)

                            for (int j = 0; j < lhs.Columns; j++)

                                m[i, j] = Calculator.Sub(lhs[i, j], rhs[i, j]);

                        return m;

                    }

                    else

                        return null;

                }



                public static Matrix<T>

                    operator *(Matrix<T> lhs, T scalar)

                {

                    Matrix<T> m = lhs.Copy();

                    for (int i = 0; i < lhs.Rows; i++)

                        for (int j = 0; j < lhs.Columns; j++)

                            m[i, j] = Calculator.Mult(lhs[i, j], scalar);

                    return m;

                }



                public static Matrix<T>

                    operator *(T scalar, Matrix<T> rhs)

                {

                    Matrix<T> m = rhs.Copy();

                    for (int i = 0; i < rhs.Rows; i++)

                        for (int j = 0; j < rhs.Columns; j++)

                            m[i, j] = Calculator.Mult(rhs[i, j], scalar);

                    return m;

                }



                public static Matrix<T>

                    operator *(Matrix<T> lhs, Matrix<T> rhs)

                {

                    if (lhs.Columns == rhs.Rows)

                    {

                        Matrix<T> m = new Matrix<T>(lhs.Rows, rhs.Columns);

                        for (int i = 0; i < lhs.Rows; i++)

                            for (int j = 0; j < rhs.Columns; j++)

                                for (int k = 0; k < lhs.Columns; k++)

                                    m[i, j] = Calculator.Add(m[i, j], Calculator.Mult(lhs[i, k], rhs[k, j]));

                        return m;

                    }

                    else

                        return null;

                }



                public static bool

                    operator ==(Matrix<T> lhs, Matrix<T> rhs)

                {

                    if (lhs.Rows != rhs.Rows || lhs.Columns != rhs.Columns)

                        return false;



                    for (int i = 0; i < lhs.Columns; i++)

                        for (int j = 0; j < lhs.Rows; j++)

                            if (Calculator.NotEqual(lhs[i, j], rhs[i, j]))

                                return false;



                    return true;

                }



                public static bool

                    operator !=(Matrix<T> lhs, Matrix<T> rhs)

                {

                    if (lhs.Rows != rhs.Rows || lhs.Columns != rhs.Columns)

                        return true;



                    for (int i = 0; i < lhs.Columns; i++)

                        for (int j = 0; j < lhs.Rows; j++)

                            if (Calculator.NotEqual(lhs[i, j], rhs[i, j]))

                                return true;



                    return false;

                }



                #region Implicit Conversions

                public static implicit operator Matrix<T>(Matrix<float> m)

                {

                    Matrix<T> t = new Matrix<T>(m.Rows, m.Columns);

                    for (int i = 0; i < t.Rows; i++)

                        for (int j = 0; j < t.Columns; j++)

                            t[i, j] = Calculator.NUM(m[i, j]);

                    return t;

                }



                public static implicit operator Matrix<T>(Matrix<double> m)

                {

                    Matrix<T> t = new Matrix<T>(m.Rows, m.Columns);

                    for (int i = 0; i < t.Rows; i++)

                        for (int j = 0; j < t.Columns; j++)

                            t[i, j] = Calculator.NUM(m[i, j]);

                    return t;

                }



                public static implicit operator Matrix<T>(Matrix<int> m)

                {

                    Matrix<T> t = new Matrix<T>(m.Rows, m.Columns);

                    for (int i = 0; i < t.Rows; i++)

                        for (int j = 0; j < t.Columns; j++)

                            t[i, j] = Calculator.NUM(m[i, j]);

                    return t;

                }



                public static implicit operator Matrix<T>(Matrix<Int64> m)

                {

                    Matrix<T> t = new Matrix<T>(m.Rows, m.Columns);

                    for (int i = 0; i < t.Rows; i++)

                        for (int j = 0; j < t.Columns; j++)

                            t[i, j] = Calculator.NUM(m[i, j]);

                    return t;

                }



                public static implicit operator Matrix<T>(Vector<T> v)

                {

                    return new Matrix<T>(v);

                }

                #endregion

                #endregion



                #region Object Overloads

                public override string ToString()

                {

                    return "Dimensions: [" + Rows + "," + Columns + "]";

                }



                public override int GetHashCode()

                {

                    return base.GetHashCode();

                }



                public override bool Equals(object obj)

                {

                    if (obj is Matrix<T>)

                    {

                        Matrix<T> m = (Matrix<T>)obj;



                        if (this == m)

                            return true;

                    }

                    return false;

                }

                #endregion



                #region Indexer

                private T Find(int row, int col)

                {

                    return _values[row, col];

                }



                public T this[int row, int col]

                {

                    get { return Find(row, col); }

                    set { SetValue(row, col, value); }

                }



                private void SetValue(int row, int col, T value)

                {

                    _values[row, col] = value;

                }

                #endregion

            }



            public class Quaternion<T>

            {

            }



            public class Complex<T>

            {

                #region Data Members

                private T _real;

                private T _imaginary;

                private static GenericMathSolution.ICalculator<T> Calculator = GenericMathSolution.CalculatorFactory<T>.Calculator;

                #endregion



                #region Constructors

                public Complex(T real, T imaginary)

                {

                    _real = real;

                    _imaginary = imaginary;

                }



                public Complex(T real)

                {

                    _real = real;

                    _imaginary = Calculator.NUM(0.0);

                }



                public Complex()

                {

                    _real = Calculator.NUM(0.0);

                    _imaginary = Calculator.NUM(0.0);

                }



                public Complex(Complex<T> c)

                {

                    _real = c._real;

                    _imaginary = c._imaginary;

                }



                public Complex<T> Copy()

                {

                    return new Complex<T>(this);

                }

                #endregion



                #region Complex Operations

                public T Real

                {

                    get { return _real; }

                    set { _real = value; }

                }



                public T Imaginary

                {

                    get { return _imaginary; }

                    set { _imaginary = value; }

                }



                public Complex<T> Conjugate

                {

                    get { return new Complex<T>(this.Real, Calculator.Negate(this.Imaginary)); }

                }



                public Complex<T> Reciprocal

                {

                    get

                    {

                        Complex<T> t = this.Conjugate;

                        T div = Calculator.Add(Calculator.Mult(this.Real, this.Real), Calculator.Mult(this.Imaginary, this.Imaginary));

                        t.Real = Calculator.Div(t.Real, div);

                        t.Imaginary = Calculator.Div(t.Imaginary, div);

                        return t;

                    }

                }



                public T Modulus

                {

                    get

                    {

                        T mod = Calculator.Add(Calculator.Mult(this.Real, this.Real), Calculator.Mult(this.Imaginary, this.Imaginary));

                        return Calculator.Pow(mod, 0.5);

                    }

                }

                #endregion



                #region Operator Overloads

                public static Complex<T>

                    operator +(Complex<T> lhs, Complex<T> rhs)

                {

                    return new Complex<T>(Calculator.Add(lhs.Real, rhs.Real), Calculator.Add(lhs.Imaginary, lhs.Imaginary));

                }



                public static Complex<T>

                    operator -(Complex<T> lhs, Complex<T> rhs)

                {

                    return new Complex<T>(Calculator.Sub(lhs.Real, rhs.Real), Calculator.Sub(lhs.Imaginary, lhs.Imaginary));

                }



                public static Complex<T>

                    operator *(Complex<T> lhs, Complex<T> rhs)

                {

                    return new Complex<T>(

                        Calculator.Sub(Calculator.Mult(lhs.Real, rhs.Real), Calculator.Mult(lhs.Imaginary, rhs.Imaginary)),

                        Calculator.Add(Calculator.Mult(lhs.Real, rhs.Imaginary), Calculator.Mult(lhs.Imaginary, rhs.Real)));

                }



                public static Complex<T>

                    operator /(Complex<T> lhs, Complex<T> rhs)

                {

                    T div = Calculator.Add(Calculator.Mult(rhs.Real, rhs.Real), Calculator.Mult(rhs.Imaginary, rhs.Imaginary));

                    Complex<T> tmp = null;

                    tmp.Real = Calculator.Add(Calculator.Mult(lhs.Real, rhs.Real), Calculator.Mult(lhs.Imaginary, rhs.Imaginary));

                    tmp.Real = Calculator.Div(tmp.Real, div);

                    tmp.Imaginary = Calculator.Sub(Calculator.Mult(lhs.Imaginary, rhs.Real), Calculator.Mult(lhs.Real, rhs.Imaginary));

                    tmp.Imaginary = Calculator.Div(tmp.Imaginary, div);

                    return tmp;

                }



                public static bool

                    operator ==(Complex<T> lhs, Complex<T> rhs)

                {

                    return Calculator.Equal(lhs.Real, rhs.Real) && Calculator.Equal(lhs.Imaginary, rhs.Imaginary) ? true : false;

                }



                public static bool

                    operator !=(Complex<T> lhs, Complex<T> rhs)

                {

                    return Calculator.Equal(lhs.Real, rhs.Real) && Calculator.Equal(lhs.Imaginary, rhs.Imaginary) ? false : true;

                }

                #endregion



                #region Object Overloads

                public override int GetHashCode()

                {

                    return base.GetHashCode();

                }



                public override bool Equals(object obj)

                {

                    if (obj is Complex<T>)

                    {

                        Complex<T> c = (Complex<T>)obj;

                        return Calculator.Equal(this.Real, c.Real) && Calculator.Equal(this.Imaginary, c.Imaginary) ? false : true;

                    }

                    return false;

                }



                public override string ToString()

                {

                    return "Real: " + this.Real + ", Imaginary: " + this.Imaginary;

                }

                #endregion

            }

        }

        #endregion



        private static class GenericMathSolution

        {

            #region Generic Math Solution

            public interface ICalculator<T>

            {

                T Add(T a, T b);

                T Sub(T a, T b);

                T Mult(T a, T b);

                T Div(T a, T b);

                T Negate(T a);



                T Pow(T a, double b);

                T Atan2(T y, T x);

                T Acos(T a);

                T Round(T a);

                T PI { get; }

                T NUM(double num);

                T NUM(T num);



                bool Equal(T a, T b);

                bool NotEqual(T a, T b);

                bool Greater(T a, T b);

                bool Less(T a, T b);

                bool GreaterEq(T a, T b);

                bool LessEq(T a, T b);

            }



            public abstract class CalculatorFactory<T>

            {

                private static ICalculator<T> _calculator = null;



                #region Calculator<T> Interfacer

                private static Type GetCalculatorType()

                {

                    Type tType = typeof(T);

                    Type calculatorType = null;

                    if (tType == typeof(Int32))

                        calculatorType = typeof(Int32Calculator);

                    else if (tType == typeof(Int64))

                        calculatorType = typeof(Int64Calculator);

                    else if (tType == typeof(float))

                        calculatorType = typeof(FloatCalculator);

                    else if (tType == typeof(Double))

                        calculatorType = typeof(DoubleCalculator);

                    else if (tType == typeof(Structures.Complex<double>))

                        calculatorType = typeof(ComplexDoubleCalculator);

                    else if (tType == typeof(Structures.Complex<float>))

                        calculatorType = typeof(ComplexFloatCalculator);

                    else

                    {

                        throw new InvalidCastException(String.Format("Unsupported Type- Type {0}" +

                              " does not have a partner implementation of interface " +

                              "ICalculator<T> and cannot be used in generic " +

                              "arithmetic using type Number<T>", tType.Name));

                    }

                    return calculatorType;

                }



                private static void MakeCalculator()

                {

                    Type calculatorType = GetCalculatorType();

                    _calculator = Activator.CreateInstance(calculatorType) as ICalculator<T>;

                }



                public static ICalculator<T> Calculator

                {

                    get

                    {

                        if (_calculator == null)

                        {

                            MakeCalculator();

                        }

                        return _calculator;

                    }

                }

                #endregion

            }



            #region Int32

            private class Int32Calculator : ICalculator<int>

            {

                public int Add(int a, int b)

                {

                    return a + b;

                }



                public int Sub(int a, int b)

                {

                    return a - b;

                }



                public int Mult(int a, int b)

                {

                    return a * b;

                }



                public int Div(int a, int b)

                {

                    return a / b;

                }



                public int Negate(int a)

                {

                    return -a;

                }



                public int Pow(int a, double b)

                {

                    return (int)System.Math.Pow(Convert.ToDouble(a), b);

                }



                public int Atan2(int y, int x)

                {

                    return (int)System.Math.Atan2(Convert.ToDouble(y), Convert.ToDouble(x));

                }



                public int Acos(int a)

                {

                    return (int)System.Math.Acos(Convert.ToDouble(a));

                }



                public int Round(int a)

                {

                    return (int)System.Math.Round(Convert.ToDouble(a));

                }



                public int PI

                {

                    get { return (int)System.Math.PI; }

                }



                public int NUM(double num)

                {

                    return int.Parse(num.ToString());

                }



                public int NUM(int num)

                {

                    return num;

                }



                public bool Equal(int a, int b)

                {

                    return a == b ? true : false;

                }



                public bool NotEqual(int a, int b)

                {

                    return a != b ? true : false;

                }



                public bool Greater(int a, int b)

                {

                    return a > b ? true : false;

                }



                public bool Less(int a, int b)

                {

                    return a < b ? true : false;

                }



                public bool GreaterEq(int a, int b)

                {

                    return a >= b ? true : false;

                }



                public bool LessEq(int a, int b)

                {

                    return a <= b ? true : false;

                }

            }

            #endregion



            #region Int64

            private class Int64Calculator : ICalculator<Int64>

            {

                public Int64 Add(Int64 a, Int64 b)

                {

                    return a + b;

                }



                public Int64 Sub(Int64 a, Int64 b)

                {

                    return a - b;

                }



                public Int64 Mult(Int64 a, Int64 b)

                {

                    return a * b;

                }



                public Int64 Div(Int64 a, Int64 b)

                {

                    return a / b;

                }



                public Int64 Negate(Int64 a)

                {

                    return -a;

                }



                public Int64 Pow(Int64 a, double b)

                {

                    return (Int64)System.Math.Pow(Convert.ToDouble(a), b);

                }



                public Int64 Atan2(Int64 y, Int64 x)

                {

                    return (Int64)System.Math.Atan2(Convert.ToDouble(y), Convert.ToDouble(x));

                }



                public Int64 Acos(Int64 a)

                {

                    return (Int64)System.Math.Acos(Convert.ToDouble(a));

                }



                public Int64 Round(Int64 a)

                {

                    return (Int64)System.Math.Round(Convert.ToDouble(a));

                }



                public Int64 PI

                {

                    get { return (Int64)System.Math.PI; }

                }



                public Int64 NUM(double num)

                {

                    return Int64.Parse(num.ToString());

                }



                public Int64 NUM(Int64 num)

                {

                    return num;

                }



                public bool Equal(Int64 a, Int64 b)

                {

                    return a == b ? true : false;

                }



                public bool NotEqual(Int64 a, Int64 b)

                {

                    return a != b ? true : false;

                }



                public bool Greater(Int64 a, Int64 b)

                {

                    return a > b ? true : false;

                }



                public bool Less(Int64 a, Int64 b)

                {

                    return a < b ? true : false;

                }



                public bool GreaterEq(Int64 a, Int64 b)

                {

                    return a >= b ? true : false;

                }



                public bool LessEq(Int64 a, Int64 b)

                {

                    return a <= b ? true : false;

                }

            }

            #endregion



            #region ComplexDouble

            private class ComplexDoubleCalculator : ICalculator<Structures.Complex<double>>

            {

                public Structures.Complex<double> Add(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    return a + b;

                }



                public Structures.Complex<double> Sub(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    return a - b;

                }



                public Structures.Complex<double> Mult(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    return a * b;

                }



                public Structures.Complex<double> Div(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    return a / b;

                }



                public Structures.Complex<double> Negate(Structures.Complex<double> a)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<double> Pow(Structures.Complex<double> a, double b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<double> Atan2(Structures.Complex<double> y, Structures.Complex<double> x)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<double> Acos(Structures.Complex<double> a)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<double> Round(Structures.Complex<double> a)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<double> PI

                {

                    get { throw new Exception("The method or operation is not implemented."); }

                }



                public Structures.Complex<double> NUM(double num)

                {

                    return new Structures.Complex<double>(num, num);

                }



                public Structures.Complex<double> NUM(Structures.Complex<double> num)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public bool Equal(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    return a == b;

                }



                public bool NotEqual(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    return a != b;

                }



                public bool Greater(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public bool Less(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public bool GreaterEq(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public bool LessEq(Structures.Complex<double> a, Structures.Complex<double> b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }

            }

            #endregion



            #region ComplexFloat

            private class ComplexFloatCalculator : ICalculator<Structures.Complex<float>>

            {

                public Structures.Complex<float> Add(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    return a + b;

                }



                public Structures.Complex<float> Sub(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    return a - b;

                }



                public Structures.Complex<float> Mult(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    return a * b;

                }



                public Structures.Complex<float> Div(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    return a / b;

                }



                public Structures.Complex<float> Negate(Structures.Complex<float> a)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<float> Pow(Structures.Complex<float> a, double b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<float> Atan2(Structures.Complex<float> y, Structures.Complex<float> x)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<float> Acos(Structures.Complex<float> a)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<float> Round(Structures.Complex<float> a)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public Structures.Complex<float> PI

                {

                    get { throw new Exception("The method or operation is not implemented."); }

                }



                public Structures.Complex<float> NUM(double num)

                {

                    return new Structures.Complex<float>(float.Parse(num.ToString()), float.Parse(num.ToString()));

                }



                public Structures.Complex<float> NUM(Structures.Complex<float> num)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public bool Equal(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    return a == b;

                }



                public bool NotEqual(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    return a != b;

                }



                public bool Greater(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public bool Less(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public bool GreaterEq(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }



                public bool LessEq(Structures.Complex<float> a, Structures.Complex<float> b)

                {

                    throw new Exception("The method or operation is not implemented.");

                }

            }

            #endregion



            #region Double

            private class DoubleCalculator : ICalculator<double>

            {

                public double Add(double a, double b)

                {

                    return a + b;

                }



                public double Sub(double a, double b)

                {

                    return a - b;

                }



                public double Mult(double a, double b)

                {

                    return a * b;

                }



                public double Div(double a, double b)

                {

                    return a / b;

                }



                public double Negate(double a)

                {

                    return -a;

                }



                public double Pow(double a, double b)

                {

                    return System.Math.Pow(a, b);

                }



                public double Atan2(double y, double x)

                {

                    return System.Math.Atan2(y, x);

                }



                public double Acos(double a)

                {

                    return System.Math.Acos(a);

                }



                public double Round(double a)

                {

                    return System.Math.Round(a);

                }



                public double PI

                {

                    get { return System.Math.PI; }

                }



                public double NUM(double num)

                {

                    return num;

                }



                public bool Equal(double a, double b)

                {

                    return a == b ? true : false;

                }



                public bool NotEqual(double a, double b)

                {

                    return a != b ? true : false;

                }



                public bool Greater(double a, double b)

                {

                    return a > b ? true : false;

                }



                public bool Less(double a, double b)

                {

                    return a < b ? true : false;

                }



                public bool GreaterEq(double a, double b)

                {

                    return a >= b ? true : false;

                }



                public bool LessEq(double a, double b)

                {

                    return a <= b ? true : false;

                }

            }

            #endregion



            #region Float

            private class FloatCalculator : ICalculator<float>

            {

                public float Add(float a, float b)

                {

                    return a + b;

                }



                public float Sub(float a, float b)

                {

                    return a - b;

                }



                public float Mult(float a, float b)

                {

                    return a * b;

                }



                public float Div(float a, float b)

                {

                    return a / b;

                }



                public float Negate(float a)

                {

                    return -a;

                }



                public float Pow(float a, double b)

                {

                    return (float)System.Math.Pow(Convert.ToDouble(a), b);

                }



                public float Atan2(float y, float x)

                {

                    return (float)System.Math.Atan2(Convert.ToDouble(y), Convert.ToDouble(x));

                }



                public float Acos(float a)

                {

                    return (float)System.Math.Acos(Convert.ToDouble(a));

                }



                public float Round(float a)

                {

                    return (float)System.Math.Round(Convert.ToDouble(a));

                }



                public float PI

                {

                    get { return (float)System.Math.PI; }

                }



                public float NUM(double num)

                {

                    return float.Parse(num.ToString());

                }



                public float NUM(float num)

                {

                    return float.Parse(num.ToString());

                }



                public bool Equal(float a, float b)

                {

                    return a == b ? true : false;

                }



                public bool NotEqual(float a, float b)

                {

                    return a != b ? true : false;

                }



                public bool Greater(float a, float b)

                {

                    return a > b ? true : false;

                }



                public bool Less(float a, float b)

                {

                    return a < b ? true : false;

                }



                public bool GreaterEq(float a, float b)

                {

                    return a >= b ? true : false;

                }



                public bool LessEq(float a, float b)

                {

                    return a <= b ? true : false;

                }

            }

            #endregion



            #endregion

        }

    }

}
