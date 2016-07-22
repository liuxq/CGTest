using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CGLearn.CG
{

    struct Triangle2DInterpolator
    {
        public Vector p1_;
        public Vector p1Value_;

        public Vector p2_;
        public Vector p2Value_;

        public Vector p3_;
        public Vector p3Value_;

        public double p2y_p3y;
        public double p3y_p1y;
        public double p3x_p2x;
        public double p1x_p3x;

        public double den;

        //public Matrix transform;

        public void SetVector1(Vector p, Vector value)
        {
            p1_ = p;
            p1Value_ = value;
        }

        public void SetVector2(Vector p, Vector value)
        {
            p2_ = p;
            p2Value_ = value;
        }

        public void SetVector3(Vector p, Vector value)
        {
            p3_ = p;
            p3Value_ = value;
        }

        public void PreCalculate()
        {
            den = 1.0 / ((p2_.y_ - p3_.y_) * (p1_.x_ - p3_.x_) + (p3_.x_ - p2_.x_) * (p1_.y_ - p3_.y_));
            p2y_p3y = p2_.y_ - p3_.y_;
            p3y_p1y = p3_.y_ - p1_.y_;
            p3x_p2x = p3_.x_ - p2_.x_;
            p1x_p3x = p1_.x_ - p3_.x_;   
        }
    }
    struct Triangle3DInterpolator
    {
        public Vector p1_;
        public Vector p1Value_;

        public Vector p2_;
        public Vector p2Value_;

        public Vector p3_;
        public Vector p3Value_;

        public double p2y_p3y;
        public double p3y_p1y;
        public double p3x_p2x;
        public double p1x_p3x;

        public double den;

        public Matrix transform;

        public void SetVector1(Vector p, Vector value)
        {
            p1_ = p;
            p1Value_ = value;
        }

        public void SetVector2(Vector p, Vector value)
        {
            p2_ = p;
            p2Value_ = value;
        }

        public void SetVector3(Vector p, Vector value)
        {
            p3_ = p;
            p3Value_ = value;
        }

        public void PreCalculate()
        {
            //把三角形投影到一个平面上
            Vector zAxis = ((p2_ - p1_).Cross(p3_ - p1_)).Normalize();
            Vector xAxis = (p2_ - p1_).Cross(zAxis).Normalize();
            Vector yAxis = zAxis.Cross(xAxis).Normalize();

            transform = Matrix.CreateIdentityMatrix(4);

            transform.SetElement(0, 0, xAxis.x_);
            transform.SetElement(0, 1, xAxis.y_);
            transform.SetElement(0, 2, xAxis.z_);

            transform.SetElement(1, 0, yAxis.x_);
            transform.SetElement(1, 1, yAxis.y_);
            transform.SetElement(1, 2, yAxis.z_);

            transform.SetElement(2, 0, zAxis.x_);
            transform.SetElement(2, 1, zAxis.y_);
            transform.SetElement(2, 2, zAxis.z_);

            p1_ = transform * p1_;
            p2_ = transform * p2_;
            p3_ = transform * p3_;

            den = 1.0 / ((p2_.y_ - p3_.y_) * (p1_.x_ - p3_.x_) + (p3_.x_ - p2_.x_) * (p1_.y_ - p3_.y_));
            p2y_p3y = p2_.y_ - p3_.y_;
            p3y_p1y = p3_.y_ - p1_.y_;
            p3x_p2x = p3_.x_ - p2_.x_;
            p1x_p3x = p1_.x_ - p3_.x_;
        }
    }
}
