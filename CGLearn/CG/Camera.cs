using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CGLearn.CG
{
    class Camera
    {
        public Vector3 position;
        public Vector3 upDirection;
        public Vector3 target;

        double zmin;
        double zmax;

        public Camera(Vector3 pos, Vector3 up, Vector3 tar, double zmin, double zmax)
        {
            position = pos;
            upDirection = up;
            target = tar;
            this.zmin = zmin;
            this.zmax = zmax;
        }

        public Matrix4 GetViewMatrix()
        {
            return Matrix4.CreateViewMatrix(position, target, upDirection);
        }

    }
}
