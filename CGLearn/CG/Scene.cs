using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using System.Drawing;
using System.Drawing.Imaging;
using System.Diagnostics;
using System.IO;


namespace CGLearn.CG
{
    class Scene3D
    {
        List<Vector3> m_listPoints = new List<Vector3>();
        List<Triangle3D> m_listTriangles = new List<Triangle3D>();
        List<Vector3> m_listPointsNormals = new List<Vector3>();
        List<Vector3> m_listPointsColors = new List<Vector3>();

        Vector3 m_vLightDir = new Vector3(0,-1,0);//向下照射的方向光
        Vector3 m_vDiffuse = new Vector3(1, 0, 0);
        Vector3 m_vSpecular = new Vector3(1, 1, 1);
        double m_dShininess = 2;
        Vector3 m_vAmbientLightColor = new Vector3(0.1, 0.1, 0.1);
        Vector3 m_vLightColor = new Vector3(0, 1, 0);

        byte[] m_listZBuffer;
        byte[] m_listZBufferInit;

        int m_iWidth;
        int m_iHeight;
        int m_iShowMode = 0;//0:line, 1:color, 2:texture
        bool m_bIsLight = false;
        Matrix4 m_mWorldTransformation;

        Bitmap m_bmpMain;
        BitmapData m_dataMain;

        Bitmap m_bmpTexture;
        BitmapData m_dataTex;
        
        Triangle2DInterpolator m_colorInterpolator = new Triangle2DInterpolator();

        Triangle3DInterpolator m_texInterpolator = new Triangle3DInterpolator();

        Clip m_clipManager = new Clip();

        Matrix4 m_mInverseMVPVTransform;
        Matrix4 m_mInverseMVPTransform;

        Stopwatch watch = new Stopwatch();

        public Bitmap GetBmpMain()
        {
            return m_bmpMain;
        }
        
        public void SwitchShowMode()
        {
            m_iShowMode = (m_iShowMode + 1) % 3;
        }
        public void SwitchLight()
        {
            m_bIsLight = !m_bIsLight;
        }

        public Scene3D(int width, int height)
        {
            m_iWidth = width;
            m_iHeight = height;
            m_mWorldTransformation = new Matrix4();
            m_mWorldTransformation = Matrix4.CreateIdentityMatrix();
            m_clipManager.Init();

            m_bmpMain = new Bitmap(m_iWidth, m_iHeight);
            m_bmpTexture = new Bitmap("env2.bmp");

            m_listZBuffer = new byte[width * height];
            m_listZBufferInit = new byte[width * height];
            for (int index = 0; index < m_listZBufferInit.Length; index++)
                m_listZBufferInit[index] = 255;

            m_dataTex = m_bmpTexture.LockBits(new Rectangle(0, 0, m_bmpTexture.Width, m_bmpTexture.Height), ImageLockMode.ReadOnly, PixelFormat.Format24bppRgb);
        }

        ~Scene3D()
        {
            
        }

        public void AddPoint(int x, int y, int z)
        {
            m_listPoints.Add(new Vector3(x, y, z));
        }

        public void AddPoint(Vector3 p)
        {
            m_listPoints.Add(p);
        }

        public void AddNormalPoint(int x, int y, int z)
        {
            this.m_listPointsNormals.Add(new Vector3(x, y, z));
            this.m_listPointsColors.Add(new Vector3());
        }

        public void NormalizeNormalPoints()
        {
            for(int i = 0; i <m_listPointsNormals.Count; i++)
            {
                m_listPointsNormals[i] = m_listPointsNormals[i].Normalize();
            }
        }

        public void AddTriangle(uint p1, uint p2, uint p3)
        {
            Triangle3D triangle = new Triangle3D(p1, p2, p3);
            AddTriangle(triangle);
        }

        public void AddTriangle(Triangle3D triangle)
        {
            m_listTriangles.Add(triangle);
        }

        public Triangle3D GetTriangle(int i)
        {
            return m_listTriangles[i];
        }

        Triangle2D ProjectTriangle(Triangle3D triangle, Matrix4 transformationMatrix)
        {
            return new Triangle2D(transformationMatrix * m_listPoints[(int)triangle.p1_],
                transformationMatrix * m_listPoints[(int)triangle.p2_],
                transformationMatrix * m_listPoints[(int)triangle.p3_]);
        }

        private int double2Int(double d)
        {
            return (int)(d + .5);
        }

        public bool isBack(Triangle2D t2)
        {
            Vector3 dir12 = t2.p1_ - t2.p2_;
            Vector3 dir13 = t2.p1_ - t2.p3_;
            dir12.z_ = 0;
            dir13.z_ = 0;
            return dir12.Cross(dir13).z_ < 0;
        }

        List<Triangle> clip(Triangle tr, Matrix4 mvp)
        {
            List<Triangle> triangles = new List<Triangle>();
            List<Vector3> inPoints = new List<Vector3>();
            List<Vector3> outPoints = new List<Vector3>();
            inPoints.Add(mvp * tr.p1_);
            inPoints.Add(mvp * tr.p2_);
            inPoints.Add(mvp * tr.p3_);

            //设置颜色插值参数
            Triangle3DInterpolator interpolator = new Triangle3DInterpolator();
            interpolator.SetVector1(inPoints[0], tr.p1Color);
            interpolator.SetVector2(inPoints[1], tr.p2Color);
            interpolator.SetVector3(inPoints[2], tr.p3Color);
            interpolator.PreCalculate();

            Triangle3DInterpolator texinterpolator = new Triangle3DInterpolator();
            texinterpolator.SetVector1(tr.p1_, tr.p1TextureCordinates_);
            texinterpolator.SetVector2(tr.p2_, tr.p2TextureCordinates_);
            texinterpolator.SetVector3(tr.p3_, tr.p3TextureCordinates_);
            texinterpolator.PreCalculate();

            for (int i = 0; i < 6; i++)
            {
                m_clipManager.SutherlandHodgmanPolygonClip(inPoints, outPoints, i);
                List<Vector3> t = outPoints;
                outPoints = inPoints;
                inPoints = t;
                outPoints.Clear();
            }

            if (inPoints.Count < 3)
                return triangles;

            List<Vector3> colors = new List<Vector3>();
            List<Vector3> texs = new List<Vector3>();
            for (int i = 0; i < inPoints.Count; i++ )
            {
                Vector3 P = interpolator.transform * inPoints[i];

                double cl1 = (interpolator.p2y_p3y * (P.x_ - interpolator.p3_.x_) + (interpolator.p3x_p2x) * (P.y_ - interpolator.p3_.y_)) * interpolator.den;
                double cl2 = ((interpolator.p3y_p1y) * (P.x_ - interpolator.p3_.x_) + (interpolator.p1x_p3x) * (P.y_ - interpolator.p3_.y_)) * interpolator.den;
                double cl3 = 1 - cl1 - cl2;

                double x = tr.p1Color.x_ * cl1 + tr.p2Color.x_ * cl2 + tr.p3Color.x_ * cl3;
                double y = tr.p1Color.y_ * cl1 + tr.p2Color.y_ * cl2 + tr.p3Color.y_ * cl3;
                double z = tr.p1Color.z_ * cl1 + tr.p2Color.z_ * cl2 + tr.p3Color.z_ * cl3;

                Vector3 rP = texinterpolator.transform * m_mInverseMVPTransform * inPoints[i];
                double tl1 = (texinterpolator.p2y_p3y * (rP.x_ - texinterpolator.p3_.x_) + (texinterpolator.p3x_p2x) * (rP.y_ - texinterpolator.p3_.y_)) * texinterpolator.den;
                double tl2 = ((texinterpolator.p3y_p1y) * (rP.x_ - texinterpolator.p3_.x_) + (texinterpolator.p1x_p3x) * (rP.y_ - texinterpolator.p3_.y_)) * texinterpolator.den;
                double tl3 = 1 - tl1 - tl2;

                double u = tr.p1TextureCordinates_.x_ * tl1 + tr.p2TextureCordinates_.x_ * tl2 + tr.p3TextureCordinates_.x_ * tl3;
                double v = tr.p1TextureCordinates_.y_ * tl1 + tr.p2TextureCordinates_.y_ * tl2 + tr.p3TextureCordinates_.y_ * tl3;
                
                colors.Add(new Vector3(x, y, z));
                texs.Add(new Vector3(u, v, 0));
            }

            for (int i = 2; i < inPoints.Count; i++)
            {
                Triangle tt = new Triangle(inPoints[0], inPoints[i - 1], inPoints[i]);
                tt.p1Color = colors[0];
                tt.p2Color = colors[i-1];
                tt.p3Color = colors[i];
                tt.p1TextureCordinates_ = texs[0];
                tt.p2TextureCordinates_ = texs[i-1];
                tt.p3TextureCordinates_ = texs[i];
                triangles.Add(tt);
            }
            return triangles;
        }

        public void DrawOriginalTriangle(Triangle3D t, Matrix4 viewport, Matrix4 mvp, Camera camera)
        {
            if (isBack(ProjectTriangle(t, viewport*mvp)))//背面消隐
                return;

            Triangle tr = new Triangle(m_listPoints[(int)t.p1_], m_listPoints[(int)t.p2_],m_listPoints[(int)t.p3_]);
            tr.p1Color = m_listPointsColors[(int)t.p1_];
            tr.p2Color = m_listPointsColors[(int)t.p2_];
            tr.p3Color = m_listPointsColors[(int)t.p3_];
            tr.p1TextureCordinates_ = t.p1TextureCordinates_;
            tr.p2TextureCordinates_ = t.p2TextureCordinates_;
            tr.p3TextureCordinates_ = t.p3TextureCordinates_;
            List<Triangle> triangles = clip(tr, mvp);

            //切分好的三角形
            for (int i = 0; i < triangles.Count; i++)
            {
                Vector3 a = viewport * triangles[i].p1_;
                Vector3 b = viewport * triangles[i].p2_;
                Vector3 c = viewport * triangles[i].p3_;
                a.x_ = Math.Floor(a.x_); a.y_ = Math.Floor(a.y_);
                b.x_ = Math.Floor(b.x_); b.y_ = Math.Floor(b.y_);
                c.x_ = Math.Floor(c.x_); c.y_ = Math.Floor(c.y_);

                //去掉面积为0的三角形
                if (a.x_ == b.x_ && a.y_ == b.y_ || a.x_ == c.x_ && a.y_ == c.y_ || b.x_ == c.x_ && b.y_== c.y_)
                    continue;

                //设置颜色插值参数
                m_colorInterpolator.SetVector1(a, triangles[i].p1Color);
                m_colorInterpolator.SetVector2(b, triangles[i].p2Color);
                m_colorInterpolator.SetVector3(c, triangles[i].p3Color);
                m_colorInterpolator.PreCalculate();

                //设置文理坐标插值参数
                m_texInterpolator.SetVector1(m_mInverseMVPVTransform * a, triangles[i].p1TextureCordinates_);
                m_texInterpolator.SetVector2(m_mInverseMVPVTransform * b, triangles[i].p2TextureCordinates_);
                m_texInterpolator.SetVector3(m_mInverseMVPVTransform * c, triangles[i].p3TextureCordinates_);
                m_texInterpolator.PreCalculate();

                if (m_iShowMode == 0)
                {
                    DDA_DrawLine(m_dataMain, double2Int(a.x_), double2Int(a.y_), double2Int(b.x_), double2Int(b.y_), triangles[i].p1Color, triangles[i].p2Color);
                    DDA_DrawLine(m_dataMain, double2Int(b.x_), double2Int(b.y_), double2Int(c.x_), double2Int(c.y_), triangles[i].p2Color, triangles[i].p3Color);
                    DDA_DrawLine(m_dataMain, double2Int(c.x_), double2Int(c.y_), double2Int(a.x_), double2Int(a.y_), triangles[i].p3Color, triangles[i].p1Color); 
                }
                else if (m_iShowMode == 1)
                {
                    Scanline_DrawTriangle(m_dataMain, a, b, c, false);
                }
                else if (m_iShowMode == 2)
                {
                    Scanline_DrawTriangle(m_dataMain, a, b, c, true);
                }
            }
        }

        private void DDA_DrawLine(BitmapData data, int p1x, int p1y, int p2x, int p2y, Vector3 color1, Vector3 color2)
        {
            byte c1R = (byte)(color1.x_ * 255);
            byte c1G = (byte)(color1.y_ * 255);
            byte c1B = (byte)(color1.z_ * 255);
            byte c2R = (byte)(color2.x_ * 255);
            byte c2G = (byte)(color2.y_ * 255);
            byte c2B = (byte)(color2.z_ * 255);

	        long YDis = (p2y - p1y);
	        long XDis = (p2x - p1x);
	        long MaxStep = Math.Max(Math.Abs(XDis),Math.Abs(YDis)); // 步进的步数
	        float fXUnitLen = 1.0f;  // X方向的单位步进
	        float fYUnitLen = 1.0f;  // Y方向的单位步进
	        fYUnitLen = (float)YDis/MaxStep;
	        fXUnitLen = (float)XDis/MaxStep;
	        
            int stride = data.Stride;
            unsafe
            {
                byte* ptr = (byte*)data.Scan0;

                // 设置起点像素颜色
                if(p1x >=0 && p1x < m_iWidth&& p1y >= 0 && p1y < m_iHeight)
                {
                    ptr[(p1x * 3) + p1y * stride] = c1B;
                    ptr[(p1x * 3) + p1y * stride + 1] = c1G;
                    ptr[(p1x * 3) + p1y * stride + 2] = c1R;
                }

                float x = p1x;
                float y = p1y;
                // 循环步进
                for (long i = 1; i <= MaxStep; i++)
                {
                    x = x + fXUnitLen;
                    y = y + fYUnitLen;
                    int ix = (int)x;
                    int iy = (int)y;

                    double beta = (double)i/MaxStep;

                    if (ix >= 0 && ix < m_iWidth && iy >= 0 && iy < m_iHeight)
                    {
                        ptr[(ix * 3) + iy * stride] = (byte)(c1B * (1 - beta) + c2B * beta);
                        ptr[(ix * 3) + iy * stride + 1] = (byte)(c1G * (1 - beta) + c2G * beta);
                        ptr[(ix * 3) + iy * stride + 2] = (byte)(c1R * (1 - beta) + c2R * beta);
                    }
                }
            }
        }
        
        public void DrawScene(Bitmap image, Matrix4 MV, Matrix4 P, Matrix4 View, Camera camera)
        {
            m_mInverseMVPVTransform = (View * P * MV).Invert4x4Matrix();
            m_mInverseMVPTransform = (P * MV).Invert4x4Matrix();
            if(m_bIsLight)
            {
                //计算光照
                for (int i = 0; i < m_listPointsNormals.Count; i++)
                {
                    Vector3 normal = m_listPointsNormals[i];
                    normal = normal.Normalize();
                    double ndot = m_vLightDir.Dot(normal);

                    Vector3 lambert = m_vDiffuse * Math.Max(ndot, 0);

                    Vector3 halfv = (m_vLightDir + camera.position).Normalize();
                    double spec = normal.Dot(halfv);
                    Vector3 specucolor = m_vSpecular * m_vLightColor * Math.Pow(Math.Max(spec, 0), m_dShininess);
                    m_listPointsColors[i] = m_vAmbientLightColor + lambert + specucolor;

                    m_listPointsColors[i].x_ = Math.Min(m_listPointsColors[i].x_, 1.0);
                    m_listPointsColors[i].y_ = Math.Min(m_listPointsColors[i].y_, 1.0);
                    m_listPointsColors[i].z_ = Math.Min(m_listPointsColors[i].z_, 1.0);
                }
            }
            else
            {
                m_listPointsColors[0].x_ = 0; m_listPointsColors[0].y_ = 0; m_listPointsColors[0].z_ = 0;
                m_listPointsColors[1].x_ = 0; m_listPointsColors[1].y_ = 0; m_listPointsColors[1].z_ = 1;
                m_listPointsColors[2].x_ = 0; m_listPointsColors[2].y_ = 1; m_listPointsColors[2].z_ = 1;
                m_listPointsColors[3].x_ = 0; m_listPointsColors[3].y_ = 1; m_listPointsColors[3].z_ = 0;
                m_listPointsColors[4].x_ = 1; m_listPointsColors[4].y_ = 0; m_listPointsColors[4].z_ = 0;
                m_listPointsColors[5].x_ = 1; m_listPointsColors[5].y_ = 0; m_listPointsColors[5].z_ = 1;
                m_listPointsColors[6].x_ = 1; m_listPointsColors[6].y_ = 1; m_listPointsColors[6].z_ = 1;
                m_listPointsColors[7].x_ = 1; m_listPointsColors[7].y_ = 1; m_listPointsColors[7].z_ = 0;
            }

            for (int i = 0; i < m_listTriangles.Count; i++)
            {
                DrawOriginalTriangle(m_listTriangles[i], View, P*MV, camera);
            }
        }

        public void Render(int width, int height, Camera camera, Matrix4 worldTransform)
        {
            Array.Copy(m_listZBufferInit, m_listZBuffer, m_listZBuffer.Length);//清空深度缓存
     
            int m = Math.Min(width, height);
            Matrix4 MV = camera.GetViewMatrix() * worldTransform;
            Matrix4 P = Matrix4.CreatePerspectiveProjectionMatrix(1.1, 1, 1, 15);
            Matrix4 ViewPort = Matrix4.CreateViewportMatrix(width / 2, height / 2, m, m);
            
            m_dataMain = m_bmpMain.LockBits(new Rectangle(0, 0, m_bmpMain.Width, m_bmpMain.Height), ImageLockMode.WriteOnly, PixelFormat.Format24bppRgb);
            
            DrawScene(m_bmpMain, MV, P, ViewPort, camera);

            m_bmpMain.UnlockBits(m_dataMain);
        }

        static void Swap<T>(ref T a, ref T b)
        {
            T t = a;
            a = b;
            b = t;
        }
        void SortVertices(ref Vector3 top, ref Vector3 middle, ref Vector3 bottom)
        {
            if (top.y_ > middle.y_)
                Swap(ref top, ref middle);

            if (middle.y_ > bottom.y_)
                Swap(ref middle, ref bottom);

            if (top.y_ > middle.y_)
                Swap(ref top, ref middle);
        }

        
        void Scanline_DrawTriangle(BitmapData data, Vector3 p1, Vector3 p2, Vector3 p3, bool beTexture)
        {  
            Vector3 top = p1, middle = p2, bottom = p3;
            SortVertices(ref top, ref middle, ref bottom);

            // inv: p1.y <= p2.y <= p3.y
            double y = middle.y_;
            double betay = (top.y_ == bottom.y_) ? 1 : ((double)y - bottom.y_) / (top.y_ - bottom.y_);
            double xr = betay * top.x_ + (1 - betay) * bottom.x_;
            double zr = betay * top.z_ + (1 - betay) * bottom.z_;

            DrawTriangleWithXParellGround(data, top, new Vector3(xr, y, zr), middle, beTexture);
            DrawTriangleWithXParellGround(data, bottom, new Vector3(xr, y, zr), middle, beTexture);

        }

        void DrawTriangleWithXParellGround(BitmapData data, Vector3 p1, Vector3 p2, Vector3 p3, bool beTexture)
        {
            // inv: p2 and p3 are on the same line, p1 is peak of triangle
            Debug.Assert(p2.y_ == p3.y_);
            if ((int)p1.y_ == (int)p2.y_)
                return;
            // we want to draw from left to right
            if (p2.x_ > p3.x_)
                Swap(ref p2, ref p3);

            // in case of problems try changing to double
            int ystart = (int) p1.y_;
            int yend = (int) p2.y_;

            float leftDxDivDy = (float)((p1.x_- p2.x_) / (p1.y_ - p2.y_));
            float rightDxDivDy = (float)((p1.x_ - p3.x_) / (p1.y_ - p3.y_));

            float xLeft = (float)p2.x_;
            float xRight = (float)p3.x_;

            int texWidth = m_bmpTexture.Width;
            int texHeight = m_bmpTexture.Height;

            int direction = ystart < yend ? 1 : -1;
            if(direction == -1)
            {
                leftDxDivDy = -leftDxDivDy;
                rightDxDivDy = -rightDxDivDy;
            }
            StVector rgb;
            int stride = data.Stride;
            int texStride = m_dataTex.Stride;
            Matrix4 texInterpolatorTransform = m_texInterpolator.transform * m_mInverseMVPVTransform;
            double[] matrix_ = texInterpolatorTransform.matrix_;
            unsafe
            {
                byte* ptr = (byte*)data.Scan0;
                byte* texPtr = (byte*)m_dataTex.Scan0;
                for (int y = yend; y != ystart - direction; y -= direction)
                {
                    int xLeft_int = (int)xLeft;
                    int xRight_int = (int)(xRight) + 1;
                    double betay = (ystart == yend) ? 1 : ((double)y - yend) / (ystart - yend);
                    double zl = betay * p1.z_ + (1 - betay) * p2.z_;
                    double zr = betay * p1.z_ + (1 - betay) * p3.z_;
                    //优化代码
                    double xl_xrDiv1 = 1.0 / (xLeft_int - xRight_int);
                    double m1_y_m3 = matrix_[1] * y + matrix_[3];
                    double m5_y_m7 = matrix_[5] * y + matrix_[7];
                    double m13_y_m15 = matrix_[13] * y + matrix_[15];
                    int yStride = y * stride;
                    int ptrIndex = 0;
                    int texPtrIndex = 0;
                    int y_width = y * m_iWidth;

                    for (int x = xLeft_int; x < xRight_int; x++)
                    {
                        double betax = (xLeft_int == xRight_int) ? 1 : (x - xRight_int) * xl_xrDiv1;
                        double z = betax * zl + (1 - betax) * zr;

                        if (x >= 0 && y >= 0 && x < m_iWidth && y < m_iHeight && z * 255 < m_listZBuffer[x + y_width])
                        {
                            m_listZBuffer[x + y_width] = (byte)(z * 255);

                            if (beTexture)
                            {
                                //纹理渲染
                                
                                //矩阵乘法优化
                                double rx, ry;
                                rx = matrix_[0] * x + matrix_[2] * z + m1_y_m3;
                                ry = matrix_[4] * x + matrix_[6] * z + m5_y_m7;
                                double w = 1 / (matrix_[12] * x + matrix_[14] * z + m13_y_m15);
                                rx *= w;
                                ry *= w;

                                double tl1 = (m_texInterpolator.p2y_p3y * (rx - m_texInterpolator.p3_.x_) + (m_texInterpolator.p3x_p2x) * (ry - m_texInterpolator.p3_.y_)) * m_texInterpolator.den;
                                double tl2 = ((m_texInterpolator.p3y_p1y) * (rx - m_texInterpolator.p3_.x_) + (m_texInterpolator.p1x_p3x) * (ry - m_texInterpolator.p3_.y_)) * m_texInterpolator.den;
                                double tl3 = 1 - tl1 - tl2;

                                int u = (int)((m_texInterpolator.p1Value_.x_ * tl1 + m_texInterpolator.p2Value_.x_ * tl2 + m_texInterpolator.p3Value_.x_ * tl3) * texWidth + 0.5);
                                int v = (int)((m_texInterpolator.p1Value_.y_ * tl1 + m_texInterpolator.p2Value_.y_ * tl2 + m_texInterpolator.p3Value_.y_ * tl3) * texHeight + 0.5);
                                u = u < 0 ? 0 : u;
                                u = u >= texWidth ? texWidth-1 : u;
                                v = v < 0 ? 0 : v;
                                v = v >= texHeight ? texHeight-1 : v;

                                ptrIndex = x+x+x + yStride;
                                texPtrIndex = u+u+u + v * texStride;
                                ptr[ptrIndex] = texPtr[texPtrIndex];
                                ptr[ptrIndex + 1] = texPtr[texPtrIndex + 1];
                                ptr[ptrIndex + 2] = texPtr[texPtrIndex + 2];
                            }
                            else
                            {
                                //颜色渲染
                                double l1 = (m_colorInterpolator.p2y_p3y * (x - m_colorInterpolator.p3_.x_) + (m_colorInterpolator.p3x_p2x) * (y - m_colorInterpolator.p3_.y_)) * m_colorInterpolator.den;
                                double l2 = ((m_colorInterpolator.p3y_p1y) * (x - m_colorInterpolator.p3_.x_) + (m_colorInterpolator.p1x_p3x) * (y - m_colorInterpolator.p3_.y_)) * m_colorInterpolator.den;
                                double l3 = 1 - l1 - l2;

                                rgb.x_ = m_colorInterpolator.p1Value_.x_ * l1 + m_colorInterpolator.p2Value_.x_ * l2 + m_colorInterpolator.p3Value_.x_ * l3;
                                rgb.y_ = m_colorInterpolator.p1Value_.y_ * l1 + m_colorInterpolator.p2Value_.y_ * l2 + m_colorInterpolator.p3Value_.y_ * l3;
                                rgb.z_ = m_colorInterpolator.p1Value_.z_ * l1 + m_colorInterpolator.p2Value_.z_ * l2 + m_colorInterpolator.p3Value_.z_ * l3;
                                rgb.x_ = rgb.x_ > 1 ? 1 : rgb.x_;
                                rgb.y_ = rgb.y_ > 1 ? 1 : rgb.y_;
                                rgb.z_ = rgb.z_ > 1 ? 1 : rgb.z_;

                                rgb.x_ = rgb.x_ < 0 ? 0 : rgb.x_;
                                rgb.y_ = rgb.y_ < 0 ? 0 : rgb.y_;
                                rgb.z_ = rgb.z_ < 0 ? 0 : rgb.z_;

                                ptrIndex = x + x + x + yStride;
                                ptr[ptrIndex] = (byte)(rgb.z_ * 255);
                                ptr[ptrIndex + 1] = (byte)(rgb.y_ * 255);
                                ptr[ptrIndex + 2] = (byte)(rgb.x_ * 255);
                            }
                            
                        }
                    }
                    xLeft -= leftDxDivDy;
                    xRight -= rightDxDivDy;
                } 
            }
        }
    }
}
