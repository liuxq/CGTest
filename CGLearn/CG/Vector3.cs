using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CGLearn.CG
{
    class Vector3
    {
        public double x_;
        public double y_;
        public double z_;

        public Vector3() : this(0, 0, 0)
        {
        }

        public Vector3(double x, double y, double z)
        {
            x_ = x;
            y_ = y;
            z_ = z;
        }

        static public Vector3 operator+(Vector3 left, Vector3 right)
        {
            return new Vector3(left.x_ + right.x_, left.y_ + right.y_, left.z_ + right.z_);
        }

        static public Vector3 operator -(Vector3 left, Vector3 right)
        {
            return new Vector3(left.x_ - right.x_, left.y_ - right.y_, left.z_ - right.z_);
        }

        static public Vector3 operator *(Vector3 left, double scalar)
        {
            return new Vector3(left.x_ * scalar, left.y_ * scalar, left.z_ * scalar);
        }

        static public Vector3 operator *(Vector3 left, Vector3 right)
        {
            return new Vector3(left.x_ * right.x_, left.y_ * right.y_, left.z_ * right.z_);
        }

        static public Vector3 operator -(Vector3 left)
        {
            return new Vector3(-left.x_, -left.y_,-left.z_);
        }

        public Vector3 Cross(Vector3 right)
        {
            double i = y_ * right.z_ - z_ * right.y_;
            double j = z_ * right.x_ - x_ * right.z_;
            double k = x_ * right.y_ - y_ * right.x_;

            return new Vector3(i, j, k);
        }

        public Vector3 Clone()
        {
            return new Vector3(x_, y_, z_);
        }

        public Vector3 Normalize()
        {
            double length = Length();

            return new Vector3(x_ / length, y_ / length, z_ / length);
        }

        public double Length()
        {
            return Math.Sqrt(x_ * x_ + y_ * y_ + z_ * z_);
        }

        public double Dot(Vector3 right)
        {
            return x_ * right.x_ + y_ * right.y_ + z_ * right.z_;
        }
   

    }
    
    //使用结构体用来防止过多gc
    struct StVector
    {
        public double x_;
        public double y_;
        public double z_;

        public StVector(double x, double y, double z)
        {
            x_ = x;
            y_ = y;
            z_ = z;
        }
    }
}
