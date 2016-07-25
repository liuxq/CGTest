using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;

namespace CGLearn.CG
{
    class Clip
    {
        List<Vector> points = new List<Vector>();
        public void Init()
        {
            points.Add(new Vector(1, 0, 0));//x = 1
            points.Add(new Vector(-1, 0, 0));//x = -1
            points.Add(new Vector(0, 1, 0));//y = 1
            points.Add(new Vector(0, -1, 0));//y = -1
            points.Add(new Vector(0, 0, 1));//z = 1
            points.Add(new Vector(0, 0, -1));//z = -1
        }
        Vector Intersect(Vector s, Vector p, int clipBoundary)
        {
            Vector I = new Vector();
            if(clipBoundary == 0)
            {
                I.x_ = 1;
                I.y_ = s.y_ + (1 - s.x_) * (p.y_ - s.y_) / (p.x_ - s.x_);
                I.z_ = s.z_ + (1 - s.x_) * (p.z_ - s.z_) / (p.x_ - s.x_);
            }
            else if(clipBoundary == 1)
            {
                I.x_ = -1;
                I.y_ = s.y_ + (-1 - s.x_) * (p.y_ - s.y_) / (p.x_ - s.x_);
                I.z_ = s.z_ + (-1 - s.x_) * (p.z_ - s.z_) / (p.x_ - s.x_);
            }
            else if (clipBoundary == 2)
            {
                I.y_ = 1;
                I.x_ = s.x_ + (1 - s.y_) * (p.x_ - s.x_) / (p.y_ - s.y_);
                I.z_ = s.z_ + (1 - s.y_) * (p.z_ - s.z_) / (p.y_ - s.y_);
            }
            else if (clipBoundary == 3)
            {
                I.y_ = -1;
                I.x_ = s.x_ + (-1 - s.y_) * (p.x_ - s.x_) / (p.y_ - s.y_);
                I.z_ = s.z_ + (-1 - s.y_) * (p.z_ - s.z_) / (p.y_ - s.y_);
            }
            else if (clipBoundary == 4)
            {
                I.z_ = 1;
                I.x_ = s.x_ + (1 - s.z_) * (p.x_ - s.x_) / (p.z_ - s.z_);
                I.y_ = s.y_ + (1 - s.z_) * (p.y_ - s.y_) / (p.z_ - s.z_);
            }
            else if (clipBoundary == 5)
            {
                I.z_ = -1;
                I.x_ = s.x_ + (-1 - s.z_) * (p.x_ - s.x_) / (p.z_ - s.z_);
                I.y_ = s.y_ + (-1 - s.z_) * (p.y_ - s.y_) / (p.z_ - s.z_);
            }

            return I;
        }
        bool Inside(Vector testVertex, int clipBoundary)
        {
            return (testVertex - points[clipBoundary]).Dot(points[clipBoundary]) <= 0;
        }
        public void SutherlandHodgmanPolygonClip(List<Vector> inPoints, List<Vector> outPoints, int clipBoundary)
        {
            if (inPoints.Count == 0)
                return;
            Vector s = inPoints[inPoints.Count - 1];
            Vector p;
            for(int i = 0; i < inPoints.Count; i++)
            {
                p = inPoints[i];
                if(Inside(p,clipBoundary))
                {
                    if (Inside(s, clipBoundary))
                        outPoints.Add(p);
                    else
                    {
                        outPoints.Add(Intersect(s, p, clipBoundary));
                        outPoints.Add(p);
                    }
                }
                else if(Inside(s,clipBoundary))
                {
                    outPoints.Add(Intersect(s, p, clipBoundary));
                }
                s = p;
            }
        }
    }
}
