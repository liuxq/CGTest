using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CGLearn.CG
{
    class Triangle3D
    {
        public uint p1_;
        public uint p2_;
        public uint p3_;

        //Vector normal_;

        public Vector3 p1TextureCordinates_;
        public Vector3 p2TextureCordinates_;
        public Vector3 p3TextureCordinates_;

        public Triangle3D(uint p1, uint p2,uint p3)
        {
            p1_ = p1;
            p2_ = p2;
            p3_ = p3;
        }

        public void SetTextureCoordinates(Vector3 t1, Vector3 t2, Vector3 t3)
        {
            p1TextureCordinates_ = t1;
            p2TextureCordinates_ = t2;
            p3TextureCordinates_ = t3;
        }

    }

    class Triangle2D
    {
        public Vector3 p1_;
        public Vector3 p2_;
        public Vector3 p3_;

        public Triangle2D(Vector3 p1, Vector3 p2, Vector3 p3)  
        {
            p1_ = p1;
            p2_ = p2;
            p3_ = p3;
        }
    }

    //三角形全数据
    class Triangle
    {
        public Vector3 p1_;
        public Vector3 p2_;
        public Vector3 p3_;

        public Vector3 p1TextureCordinates_;
        public Vector3 p2TextureCordinates_;
        public Vector3 p3TextureCordinates_;

        public Vector3 p1Color;
        public Vector3 p2Color;
        public Vector3 p3Color;

        public Triangle(Vector3 p1, Vector3 p2, Vector3 p3)
        {
            p1_ = p1;
            p2_ = p2;
            p3_ = p3;
        }
    }
}
