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
        List<Vector3> m_listPoints = new List<Vector3>();
        public void Init()
        {
            m_listPoints.Add(new Vector3(1, 0, 0));//x = 1
            m_listPoints.Add(new Vector3(-1, 0, 0));//x = -1
            m_listPoints.Add(new Vector3(0, 1, 0));//y = 1
            m_listPoints.Add(new Vector3(0, -1, 0));//y = -1
            m_listPoints.Add(new Vector3(0, 0, 1));//z = 1
            m_listPoints.Add(new Vector3(0, 0, -1));//z = -1
        }
        Vector3 Intersect(Vector3 s, Vector3 p, int clipBoundary)
        {
            Vector3 I = new Vector3();
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
        bool Inside(Vector3 testVertex, int clipBoundary)
        {
            return (testVertex - m_listPoints[clipBoundary]).Dot(m_listPoints[clipBoundary]) <= 0;
        }
        public void SutherlandHodgmanPolygonClip(List<Vector3> inPoints, List<Vector3> outPoints, int clipBoundary)
        {
            if (inPoints.Count == 0)
                return;
            Vector3 s = inPoints[inPoints.Count - 1];
            Vector3 p;
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
