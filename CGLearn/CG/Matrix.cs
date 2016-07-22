using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CGLearn.CG
{
    class Matrix
    {
        private double[] matrix_;
        private int width_;
        private int height_;

        static private Matrix pointMatrix = new Matrix(1, 4);
        static private Matrix transformedPointMatrix = new Matrix(1, 4);

        bool DoubleEquals(double d1, double d2)
        {
            return Math.Abs(d1 - d2) < 0.0001;
        }

        public Matrix(int width, int height)
        {
            height_ = height;
            width_ = width;
            matrix_ = new double[width * height];
            for (uint i = 0; i < width_ * height_; i++)
                matrix_[i] = 0;
        }

        public Matrix(Matrix other) :this(other.width_, other.height_)
        {
            for (int i = 0; i < height_; i++)
                for (int j = 0; j < width_; j++)
                    SetElement(i, j, other.GetElement(i, j));
        }

        ~Matrix()
        {
            
        }

        public double GetElement(int row, int col)
        {
            return matrix_[row * width_ + col];
        }

        public void SetElement(int row, int col, double val)
        {
            matrix_[row * width_ + col] = val;
        }
        static public void Multiply(Matrix left, Matrix right, Matrix result)
        {
            int height = left.height_;
            int width = right.width_;
            Array.Clear(result.matrix_, 0x0, width * height);

            for( int i = 0; i < left.height_; i++ ) {
                for( int j = 0; j < right.width_; j++ )
                {
                    for(int k = 0; k < left.width_; k++)
                    {
                        result.matrix_[i*right.width_+j] += left.matrix_[i*left.width_+k]*right.matrix_[k * right.width_+j];
                    }
                }
            }
        }
        static public Vector operator*(Matrix left, Vector v)
        {
            pointMatrix.SetElement(0, 0, v.x_);
            pointMatrix.SetElement(1, 0, v.y_);
            pointMatrix.SetElement(2, 0, v.z_);
            pointMatrix.SetElement(3, 0, 1);

            Matrix.Multiply(left, pointMatrix, transformedPointMatrix);

            double w = transformedPointMatrix.GetElement(3, 0);
            return new Vector(transformedPointMatrix.GetElement(0, 0) / w,
                transformedPointMatrix.GetElement(1, 0) / w,
                transformedPointMatrix.GetElement(2, 0) / w);
        }

        static public StVector operator *(Matrix left, StVector v)
        {
            pointMatrix.SetElement(0, 0, v.x_);
            pointMatrix.SetElement(1, 0, v.y_);
            pointMatrix.SetElement(2, 0, v.z_);
            pointMatrix.SetElement(3, 0, 1);

            Matrix.Multiply(left, pointMatrix, transformedPointMatrix);

            double w = transformedPointMatrix.GetElement(3, 0);

            v.x_ = transformedPointMatrix.GetElement(0, 0) / w;
            v.y_ = transformedPointMatrix.GetElement(1, 0) / w;
            v.z_ = transformedPointMatrix.GetElement(2, 0) / w;
            return v;
        }
            

        static public Matrix operator*(Matrix l, double number)
        {
            Matrix result = new Matrix(l.width_, l.height_);

            for (int row = 0; row < l.height_; row++)
                for (int col = 0; col < l.width_; col++)
                    result.SetElement(row, col, l.GetElement(row, col) * number);

            return result;
        }

        static public Matrix operator *(Matrix l, Matrix right)
        {
            Matrix result = new Matrix(l.height_, right.width_);

            Multiply(l, right, result);

            return result;
        }

        static public Matrix operator -(Matrix l, Matrix right)
        {
            Matrix result = new Matrix(l.width_, l.height_);

            for (int row = 0; row < l.height_; row++)
                for (int col = 0; col < l.width_; col++)
                    result.SetElement(row, col, l.GetElement(row, col) - right.GetElement(row, col));

            return result;
        }

        static public Matrix CreateIdentityMatrix(int size)
        {
            Matrix result = new Matrix(size, size);

            for (int i = 0; i < size; i++)
            {
                result.SetElement(i, i, 1);
            }

            return result;
        }


        static public Matrix Create4x4Matrix(double m11, double m12, double m13, double m14,
                double m21, double m22, double m23, double m24,
                double m31, double m32, double m33, double m34,
                double m41, double m42, double m43, double m44)
        {
            Matrix m = new Matrix(4, 4);
            m.matrix_[0] =  m11;
            m.matrix_[1] =  m12;
            m.matrix_[2] =  m13;
            m.matrix_[3] =  m14;

            m.matrix_[4] =  m21;
            m.matrix_[5] =  m22;
            m.matrix_[6] =  m23;
            m.matrix_[7] =  m24;

            m.matrix_[8] =  m31;
            m.matrix_[9] =  m32;
            m.matrix_[10] =  m33;
            m.matrix_[11] =  m34;

            m.matrix_[12] =  m41;
            m.matrix_[13] =  m42;
            m.matrix_[14] =  m43;
            m.matrix_[15] =  m44;

            return m;
        }

        static public Matrix Create3x3Matrix(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33)
        {
            Matrix m = new Matrix(3, 3);

            m.matrix_[0] = m11;
            m.matrix_[1] = m12;
            m.matrix_[2] =  m13;

            m.matrix_[3] =  m21;
            m.matrix_[4] =  m22;
            m.matrix_[5] =  m23;

            m.matrix_[6] =  m31;
            m.matrix_[7] =  m32;
            m.matrix_[8] =  m33;

            return m;
        }

        static public void Set3x3Matrix(Matrix m, double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33)
        {
            m.matrix_[0] = m11;
            m.matrix_[1] = m12;
            m.matrix_[2] =  m13;

            m.matrix_[3] =  m21;
            m.matrix_[4] =  m22;
            m.matrix_[5] =  m23;

            m.matrix_[6] =  m31;
            m.matrix_[7] =  m32;
            m.matrix_[8] =  m33;
        }

        static public Matrix Create3x1Matrix(double m11, double m21, double m31)
        {
            Matrix m = new Matrix(1, 3);

            m.matrix_[0] =  m11;
            m.matrix_[1] =  m21;
            m.matrix_[2] =  m31;

            return m;
        }

        static public Matrix CreateScaleMatrix(double xFactor, double yFactor, double zFactor)
        {
            return Create4x4Matrix(xFactor, 0, 0, 0,
                                      0, yFactor, 0, 0,
                                      0, 0, zFactor, 0,
                                      0, 0, 0, 1);
        }

        static public Matrix CreateTranslationMatrix(double xMove, double yMove, double zMove)
        {
            return Create4x4Matrix(1, 0, 0, xMove,
                                  0, 1, 0, yMove,
                                  0, 0, 1, zMove,
                                  0, 0, 0, 1);
        }

        static public Matrix CreateTranslationMatrix(Vector v)
        {
            return CreateTranslationMatrix(v.x_, v.y_, v.z_);
        }

        static public Matrix CreateXAxisRotationMatrix(double radians)
        {
            return Create4x4Matrix(1, 0, 0, 0,
                  0, Math.Cos(radians), -Math.Sin(radians), 0,
                  0, Math.Sin(radians), Math.Cos(radians), 0,
                  0, 0, 0, 1);
        }

        static public Matrix CreateXAxisRotationMatrixAroundPoint(double angleInRadians, Vector p)
        {
            return CreateTranslationMatrix(p) * CreateXAxisRotationMatrix(angleInRadians) * CreateTranslationMatrix(-p);
        }

        static public Matrix CreateYAxisRotationMatrix(double radians)
        {
            return Create4x4Matrix(Math.Cos(radians), 0, Math.Sin(radians), 0,
                                                          0, 1, 0, 0,
                                                          -Math.Sin(radians), 0, Math.Cos(radians), 0,
                                                          0, 0, 0, 1);
        }

        static public Matrix CreateYAxisRotationMatrixAroundPoint(double angleInRadians, Vector p)
        {
            return CreateTranslationMatrix(p) * CreateYAxisRotationMatrix(angleInRadians) * CreateTranslationMatrix(-p);
        }

        static public Matrix CreateZAxisRotationMatrix(double radians)
        {
            return Create4x4Matrix(Math.Cos(radians), -Math.Sin(radians), 0, 0,
                    Math.Sin(radians), Math.Cos(radians), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);

        }

        static public Matrix CreateZAxisRotationMatrixAroundPoint(double angleInRadians, Vector p)
        {
            return CreateTranslationMatrix(p) * CreateZAxisRotationMatrix(angleInRadians) * CreateTranslationMatrix(-p);
        }

        static public Matrix CreateOrthogonalProjectionMatrix(int x, int y, int viewportWidth, int viewportHeight)
        {
            return 
                CreateTranslationMatrix(x, y, 0) *
                CreateScaleMatrix(viewportWidth / 2, viewportHeight / 2, 1) * 
                Create4x4Matrix(
                    1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);
        }

        static public Matrix CreatePerspectiveProjectionMatrix2(int x, int y, int viewportWidth, int viewportHeight, double viewAngleRad,
                double aspect, double znear, double zfar)
        {
            Matrix result = new Matrix(4, 4);
            double yScale = 1.0 / Math.Tan(viewAngleRad * 0.5);
            double xScale = yScale / aspect;

            double halfWidth = znear / xScale;
            double halfHeight = znear / yScale;
            double left = -halfWidth;
            double right = halfWidth;
            double bottom = -halfHeight;
            double top = halfHeight;

            double zRange = zfar / (zfar - znear);

            result.SetElement(0, 0, 2 * znear / (right - left));
            result.SetElement(1, 1, 2 * znear / (top - bottom));
            result.SetElement(2, 2, -(zfar + znear) / (zfar - znear));
            result.SetElement(2, 3, -2 * znear * zfar / (zfar - znear));
            result.SetElement(3, 2, -1);

            return CreateTranslationMatrix(x, y, 0) * CreateScaleMatrix(viewportWidth / 2, viewportHeight / 2, 1)* result;
        }

        static public Matrix CreateViewMatrix(Vector cameraPosition, Vector observedPoint, Vector upDirection)
        {
            Vector cameraVector = new Vector(0, 0, 0) - cameraPosition;
            Vector zAxis = (observedPoint - cameraVector).Normalize();
            Vector xAxis = upDirection.Cross(zAxis).Normalize();
            Vector yAxis = zAxis.Cross(xAxis).Normalize();

            Matrix result =  CreateIdentityMatrix(4);

            result.SetElement(0, 0, xAxis.x_);
            result.SetElement(0, 1, xAxis.y_);
            result.SetElement(0, 2, xAxis.z_);

            result.SetElement(1, 0, yAxis.x_);
            result.SetElement(1, 1, yAxis.y_);
            result.SetElement(1, 2, yAxis.z_);

            result.SetElement(2, 0, zAxis.x_);
            result.SetElement(2, 1, zAxis.y_);
            result.SetElement(2, 2, zAxis.z_);

            result.SetElement(0, 3, -xAxis.Dot(cameraVector));
            result.SetElement(1, 3, -yAxis.Dot(cameraVector));
            result.SetElement(2, 3, -zAxis.Dot(cameraVector));

            return result;
        }

        Matrix Clone()
        {
            Matrix result = new Matrix(width_, height_);
            for(int i = 0; i < width_*height_; i++)
            {
                result.matrix_[i] = matrix_[i];
            }
            return result;
        }

        //四阶
        double Determinant()
        {
            double a00 = matrix_[0], a01 = matrix_[1], a02 = matrix_[2], a03 = matrix_[3],
                   a10 = matrix_[4], a11 = matrix_[5], a12 = matrix_[6], a13 = matrix_[7],
                   a20 = matrix_[8], a21 = matrix_[9], a22 = matrix_[10], a23 = matrix_[11],
                   a30 = matrix_[12], a31 = matrix_[13], a32 = matrix_[14], a33 = matrix_[15],

            b00 = a00 * a11 - a01 * a10,
            b01 = a00 * a12 - a02 * a10,
            b02 = a00 * a13 - a03 * a10,
            b03 = a01 * a12 - a02 * a11,
            b04 = a01 * a13 - a03 * a11,
            b05 = a02 * a13 - a03 * a12,
            b06 = a20 * a31 - a21 * a30,
            b07 = a20 * a32 - a22 * a30,
            b08 = a20 * a33 - a23 * a30,
            b09 = a21 * a32 - a22 * a31,
            b10 = a21 * a33 - a23 * a31,
            b11 = a22 * a33 - a23 * a32;

            // Calculate the determinant
            return b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
        }

        public Matrix Invert4x4Matrix()
        {
            double a00 = matrix_[0], a01 = matrix_[1], a02 = matrix_[2], a03 = matrix_[3],
            a10 = matrix_[4], a11 = matrix_[5], a12 = matrix_[6], a13 = matrix_[7],
            a20 = matrix_[8], a21 = matrix_[9], a22 = matrix_[10], a23 = matrix_[11],
            a30 = matrix_[12], a31 = matrix_[13], a32 = matrix_[14], a33 = matrix_[15],

            b00 = a00 * a11 - a01 * a10,
            b01 = a00 * a12 - a02 * a10,
            b02 = a00 * a13 - a03 * a10,
            b03 = a01 * a12 - a02 * a11,
            b04 = a01 * a13 - a03 * a11,
            b05 = a02 * a13 - a03 * a12,
            b06 = a20 * a31 - a21 * a30,
            b07 = a20 * a32 - a22 * a30,
            b08 = a20 * a33 - a23 * a30,
            b09 = a21 * a32 - a22 * a31,
            b10 = a21 * a33 - a23 * a31,
            b11 = a22 * a33 - a23 * a32,

            // Calculate the determinant
            det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

            if (Math.Abs(det) < 0.001)
            {
                return null;
            }
            det = 1.0 / det;

            Matrix result = this.Clone();

            result.matrix_[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
            result.matrix_[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
            result.matrix_[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
            result.matrix_[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
            result.matrix_[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
            result.matrix_[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
            result.matrix_[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
            result.matrix_[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
            result.matrix_[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
            result.matrix_[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
            result.matrix_[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
            result.matrix_[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
            result.matrix_[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
            result.matrix_[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
            result.matrix_[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
            result.matrix_[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;

            return result;
        }

    }
}
