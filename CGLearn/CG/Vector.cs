﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CGLearn.CG
{
    class Vector
    {
        public double x_;
        public double y_;
        public double z_;

        static private Matrix pointMatrix = new Matrix(1,4);
        static private Matrix transformedPointMatrix = new Matrix(1,4);

        
        public Vector() : this(0, 0, 0)
        {
        }

        public Vector(double x, double y, double z)
        {
            x_ = x;
            y_ = y;
            z_ = z;
        }

        static public Vector operator+(Vector left, Vector right)
        {
            return new Vector(left.x_ + right.x_, left.y_ + right.y_, left.z_ + right.z_);
        }

        static public Vector operator -(Vector left, Vector right)
        {
            return new Vector(left.x_ - right.x_, left.y_ - right.y_, left.z_ - right.z_);
        }

        static public Vector operator *(Vector left, double scalar)
        {
            return new Vector(left.x_ * scalar, left.y_ * scalar, left.z_ * scalar);
        }

        static public Vector operator *(Vector left, Vector right)
        {
            return new Vector(left.x_ * right.x_, left.y_ * right.y_, left.z_ * right.z_);
        }

        static public Vector operator -(Vector left)
        {
            return new Vector(-left.x_, -left.y_,-left.z_);
        }

        public Vector Cross(Vector right)
        {
            double i = y_ * right.z_ - z_ * right.y_;
            double j = z_ * right.x_ - x_ * right.z_;
            double k = x_ * right.y_ - y_ * right.x_;

            return new Vector(i, j, k);
        }

        public Vector Clone()
        {
            return new Vector(x_, y_, z_);
        }

        public Vector Normalize()
        {
            double length = Length();

            return new Vector(x_ / length, y_ / length, z_ / length);
        }

        public double Length()
        {
            return Math.Sqrt(x_ * x_ + y_ * y_ + z_ * z_);
        }

        public double Dot(Vector right)
        {
            return x_ * right.x_ + y_ * right.y_ + z_ * right.z_;
        }

        public Vector Transform(Matrix transformationMatrix)
        {
            pointMatrix.SetElement(0, 0, x_);
            pointMatrix.SetElement(1, 0, y_);
            pointMatrix.SetElement(2, 0, z_);
            pointMatrix.SetElement(3, 0, 1);

            Matrix.Multiply(transformationMatrix, pointMatrix, transformedPointMatrix);

            double w = transformedPointMatrix.GetElement(3, 0);
            return new Vector(transformedPointMatrix.GetElement(0, 0) / w,
                transformedPointMatrix.GetElement(1, 0) / w,
                transformedPointMatrix.GetElement(2, 0) / w);
        }

        //static public bool operator == (Vector left, Vector right)
        //{
        //    return Math.Abs(left.x_ - right.x_) <= 0.0001 &&
        //        Math.Abs(left.y_ - right.y_) <= 0.0001 &&
        //        Math.Abs(left.z_ - right.z_) <= 0.0001;
        //}

        //static public bool operator !=(Vector left, Vector right)
        //{
        //    return !(Math.Abs(left.x_ - right.x_) <= 0.0001 &&
        //        Math.Abs(left.y_ - right.y_) <= 0.0001 &&
        //        Math.Abs(left.z_ - right.z_) <= 0.0001);
        //}

        

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