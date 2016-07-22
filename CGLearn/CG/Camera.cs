using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CGLearn.CG
{
    class Camera
    {
        public Vector position;
        public Vector upDirection;
        public Vector target;

        double zmin;
        double zmax;

        public Camera(Vector pos, Vector up, Vector tar, double zmin, double zmax)
        {
            position = pos;
            upDirection = up;
            target = tar;
            this.zmin = zmin;
            this.zmax = zmax;
        }

        public void Transform(Matrix transformationMatrix)
        {
            position = position.Transform(transformationMatrix);
            upDirection = upDirection.Transform(transformationMatrix);
            target = target.Transform(transformationMatrix);
        }

        public Matrix GetViewMatrix()
        {
            return Matrix.CreateViewMatrix(position, target, upDirection);
        }

    }
}
