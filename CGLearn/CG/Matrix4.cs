using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CGLearn.CG
{
    class Matrix4
    {
        public double[] matrix_;

        bool DoubleEquals(double d1, double d2)
        {
            return Math.Abs(d1 - d2) < 0.0001;
        }

        public Matrix4()
        {
            matrix_ = new double[16];
            Array.Clear(matrix_, 0, 16);
        }

        public Matrix4(Matrix4 other) :this()
        {
            Array.Copy(other.matrix_, matrix_, 16);
        }

        ~Matrix4()
        {
            
        }

        public double GetElement(int row, int col)
        {
            return matrix_[row * 4 + col];
        }

        public void SetElement(int row, int col, double val)
        {
            matrix_[row * 4 + col] = val;
        }

        static public void Multiply(Matrix4 left, Matrix4 right, Matrix4 result)
        {
            double a00 = right.matrix_[0], a01 = right.matrix_[1], a02 = right.matrix_[2], a03 = right.matrix_[3],
                a10 = right.matrix_[4], a11 = right.matrix_[5], a12 = right.matrix_[6], a13 = right.matrix_[7],
                a20 = right.matrix_[8], a21 = right.matrix_[9], a22 = right.matrix_[10], a23 = right.matrix_[11],
                a30 = right.matrix_[12], a31 = right.matrix_[13], a32 = right.matrix_[14], a33 = right.matrix_[15];

            // Cache only the current line of the second matrix
            double b0  = left.matrix_[0], b1 = left.matrix_[1], b2 = left.matrix_[2], b3 = left.matrix_[3];
            result.matrix_[0] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
            result.matrix_[1] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
            result.matrix_[2] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
            result.matrix_[3] = b0*a03 + b1*a13 + b2*a23 + b3*a33;

            b0 = left.matrix_[4]; b1 = left.matrix_[5]; b2 = left.matrix_[6]; b3 = left.matrix_[7];
            result.matrix_[4] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
            result.matrix_[5] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
            result.matrix_[6] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
            result.matrix_[7] = b0*a03 + b1*a13 + b2*a23 + b3*a33;

            b0 = left.matrix_[8]; b1 = left.matrix_[9]; b2 = left.matrix_[10]; b3 = left.matrix_[11];
            result.matrix_[8] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
            result.matrix_[9] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
            result.matrix_[10] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
            result.matrix_[11] = b0*a03 + b1*a13 + b2*a23 + b3*a33;

            b0 = left.matrix_[12]; b1 = left.matrix_[13]; b2 = left.matrix_[14]; b3 = left.matrix_[15];
            result.matrix_[12] = b0*a00 + b1*a10 + b2*a20 + b3*a30;
            result.matrix_[13] = b0*a01 + b1*a11 + b2*a21 + b3*a31;
            result.matrix_[14] = b0*a02 + b1*a12 + b2*a22 + b3*a32;
            result.matrix_[15] = b0*a03 + b1*a13 + b2*a23 + b3*a33;
        }

        static public Vector3 operator*(Matrix4 left, Vector3 v)
        {
            Vector3 re = new Vector3();
            re.x_ = left.matrix_[0] * v.x_ + left.matrix_[1] * v.y_ + left.matrix_[2] * v.z_ + left.matrix_[3];
            re.y_ = left.matrix_[4] * v.x_ + left.matrix_[5] * v.y_ + left.matrix_[6] * v.z_ + left.matrix_[7];
            re.z_ = left.matrix_[8] * v.x_ + left.matrix_[9] * v.y_ + left.matrix_[10] * v.z_ + left.matrix_[11];
            double w = 1 / (left.matrix_[12] * v.x_ + left.matrix_[13] * v.y_ + left.matrix_[14] * v.z_ + left.matrix_[15]);

            re.x_ = re.x_ * w;
            re.y_ = re.y_ * w;
            re.z_ = re.z_ * w;
            return re;
        }

        static public StVector operator *(Matrix4 left, StVector v)
        {
            StVector re;
            re.x_ = left.matrix_[0] * v.x_ + left.matrix_[1] * v.y_ + left.matrix_[2] * v.z_ + left.matrix_[3];
            re.y_ = left.matrix_[4] * v.x_ + left.matrix_[5] * v.y_ + left.matrix_[6] * v.z_ + left.matrix_[7];
            re.z_ = left.matrix_[8] * v.x_ + left.matrix_[9] * v.y_ + left.matrix_[10] * v.z_ + left.matrix_[11];
            double w = 1/(left.matrix_[12] * v.x_ + left.matrix_[13] * v.y_ + left.matrix_[14] * v.z_ + left.matrix_[15]);

            re.x_ = re.x_ * w;
            re.y_ = re.y_ * w;
            re.z_ = re.z_ * w;
            return re;
        }
            
        static public Matrix4 operator*(Matrix4 left, double number)
        {
            Matrix4 result = new Matrix4();

            for (int i = 0; i < 16; i++)
                result.matrix_[i] = left.matrix_[i] *= number;

            return result;
        }

        static public Matrix4 operator *(Matrix4 left, Matrix4 right)
        {
            Matrix4 result = new Matrix4();

            Multiply(left, right, result);

            return result;
        }

        static public Matrix4 operator -(Matrix4 left, Matrix4 right)
        {
            Matrix4 result = new Matrix4();
            for (int i = 0; i < 16; i++)
                result.matrix_[i] = left.matrix_[i] - right.matrix_[i];

            return result;
        }

        static public Matrix4 CreateIdentityMatrix()
        {
            Matrix4 result = new Matrix4();

            for (int i = 0; i < 4; i++)
            {
                result.SetElement(i, i, 1);
            }
            return result;
        }

        static public Matrix4 Create4x4Matrix(double m11, double m12, double m13, double m14,
                double m21, double m22, double m23, double m24,
                double m31, double m32, double m33, double m34,
                double m41, double m42, double m43, double m44)
        {
            Matrix4 m = new Matrix4();
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

        static public Matrix4 CreateScaleMatrix(double xFactor, double yFactor, double zFactor)
        {
            return Create4x4Matrix(xFactor, 0, 0, 0,
                                      0, yFactor, 0, 0,
                                      0, 0, zFactor, 0,
                                      0, 0, 0, 1);
        }

        static public Matrix4 CreateTranslationMatrix(double xMove, double yMove, double zMove)
        {
            return Create4x4Matrix(1, 0, 0, xMove,
                                  0, 1, 0, yMove,
                                  0, 0, 1, zMove,
                                  0, 0, 0, 1);
        }

        static public Matrix4 CreateTranslationMatrix(Vector3 v)
        {
            return CreateTranslationMatrix(v.x_, v.y_, v.z_);
        }

        static public Matrix4 CreateXAxisRotationMatrix(double radians)
        {
            return Create4x4Matrix(1, 0, 0, 0,
                  0, Math.Cos(radians), -Math.Sin(radians), 0,
                  0, Math.Sin(radians), Math.Cos(radians), 0,
                  0, 0, 0, 1);
        }

        static public Matrix4 CreateYAxisRotationMatrix(double radians)
        {
            return Create4x4Matrix(Math.Cos(radians), 0, Math.Sin(radians), 0,
                                                          0, 1, 0, 0,
                                                          -Math.Sin(radians), 0, Math.Cos(radians), 0,
                                                          0, 0, 0, 1);
        }

        static public Matrix4 CreateZAxisRotationMatrix(double radians)
        {
            return Create4x4Matrix(Math.Cos(radians), -Math.Sin(radians), 0, 0,
                    Math.Sin(radians), Math.Cos(radians), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);

        }

        static public Matrix4 CreatePerspectiveProjectionMatrix(double fovy, double aspect, double znear, double zfar)
        {
            Matrix4 result = new Matrix4();
            double f = 1.0 / Math.Tan(fovy / 2),
            nf = 1 / (znear - zfar);
            result.matrix_[0] = f / aspect;
            result.matrix_[1] = 0;
            result.matrix_[2] = 0;
            result.matrix_[3] = 0;
            result.matrix_[4] = 0;
            result.matrix_[5] = f;
            result.matrix_[6] = 0;
            result.matrix_[7] = 0;
            result.matrix_[8] = 0;
            result.matrix_[9] = 0;
            result.matrix_[10] = (zfar + znear) * nf;
            result.matrix_[11] = (2 * zfar * znear) * nf; 
            result.matrix_[12] = 0;
            result.matrix_[13] = 0;
            result.matrix_[14] = -1;
            result.matrix_[15] = 0;

            return result;
        }

        static public Matrix4 CreateViewportMatrix(int cX, int cY, int viewportWidth, int viewportHeight)
        {
            return CreateTranslationMatrix(cX, cY, 0) * CreateScaleMatrix(viewportWidth / 2, viewportHeight / 2, 1);
        }

        static public Matrix4 CreateViewMatrix(Vector3 cameraPosition, Vector3 observedPoint, Vector3 upDirection)
        {
            Vector3 cameraVector = cameraPosition;
            Vector3 zAxis = (cameraVector - observedPoint ).Normalize();
            Vector3 xAxis = upDirection.Cross(zAxis).Normalize();
            Vector3 yAxis = zAxis.Cross(xAxis).Normalize();

            Matrix4 result =  CreateIdentityMatrix();

            result.matrix_[0] = xAxis.x_;
            result.matrix_[1] = xAxis.y_;
            result.matrix_[2] = xAxis.z_;

            result.matrix_[4] = yAxis.x_;
            result.matrix_[5] = yAxis.y_;
            result.matrix_[6] = yAxis.z_;

            result.matrix_[8]=  zAxis.x_;
            result.matrix_[9] =  zAxis.y_;
            result.matrix_[10] = zAxis.z_;

            result.matrix_[3] = -xAxis.Dot(cameraVector);
            result.matrix_[7] = -yAxis.Dot(cameraVector);
            result.matrix_[11] = -zAxis.Dot(cameraVector);

            return result;
        }

        Matrix4 Clone()
        {
            Matrix4 result = new Matrix4();
            Array.Copy(matrix_, result.matrix_, 16);
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

        public Matrix4 Invert4x4Matrix()
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

            Matrix4 result = this.Clone();

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
