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

        public Vector p1TextureCordinates_;
        public Vector p2TextureCordinates_;
        public Vector p3TextureCordinates_;

        public Triangle3D(uint p1, uint p2,uint p3)
        {
            p1_ = p1;
            p2_ = p2;
            p3_ = p3;
        }

        public void SetTextureCoordinates(Vector t1, Vector t2, Vector t3)
        {
            p1TextureCordinates_ = t1;
            p2TextureCordinates_ = t2;
            p3TextureCordinates_ = t3;
        }

    }

    class Triangle2D
    {
        public Vector p1_;
        public Vector p2_;
        public Vector p3_;

        public Triangle2D(Vector p1, Vector p2, Vector p3)  
        {
            p1_ = p1;
            p2_ = p2;
            p3_ = p3;
        }
    }
}
