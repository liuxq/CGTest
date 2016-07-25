using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;
using System.Drawing;
using System.Drawing.Imaging;
using System.Diagnostics;

namespace CGLearn.CG
{
    class Scene3D
    {
        List<Vector> points_ = new List<Vector>();
        List<Triangle3D> triangles_ = new List<Triangle3D>();
        List<Vector> pointsNormals_ = new List<Vector>();
        List<Vector> pointsColors_ = new List<Vector>();

        Vector lightDir = new Vector(0,-1,0);//向下照射的方向光
        Vector diffuse = new Vector(1, 0, 0);
        Vector specular = new Vector(1, 1, 1);
        double shininess = 2;
        Vector ambientLightColor = new Vector(0.1, 0.1, 0.1);
        Vector lightColor = new Vector(0, 1, 0);

        double[,] zBuffer_;

        int width_;
        int height_;
        Matrix worldTransformation_;

        int showMode = 0;//0:line, 1:color, 2:texture
        bool isLight = false;

        Bitmap bmpMain;
        BitmapData dataMain;

        Bitmap bmpTexture;
        BitmapData dataTex;

        Bitmap bmpDepth;
        BitmapData dataDepth;

        Graphics graphMain;
        Graphics graDepth;
        
        Triangle2DInterpolator colorInterpolator = new Triangle2DInterpolator();

        Triangle3DInterpolator texInterpolator = new Triangle3DInterpolator();

        Clip clipManager = new Clip();

        Font drawFont = new Font("宋体", 16);
        SolidBrush drawBrush = new SolidBrush(Color.White);

        Matrix inverseMVPVTransform;
        Matrix inverseMVPTransform;

        Stopwatch watch = new Stopwatch();

        public Bitmap GetBmpMain()
        {
            return bmpMain;
        }
        
        public void SwitchShowMode()
        {
            showMode = (showMode + 1) % 3;
        }
        public void SwitchLight()
        {
            isLight = !isLight;
        }

        public Scene3D(int width, int height)
        {
            width_ = width;
            height_ = height;
            worldTransformation_ = new Matrix(4,4);
            worldTransformation_ = Matrix.CreateIdentityMatrix(4);
            clipManager.Init();
            bmpMain = new Bitmap(width_, height_);
            bmpDepth = new Bitmap(width_, height_);
            bmpTexture = new Bitmap("env2.bmp");
            graphMain = Graphics.FromImage(bmpMain);
            graDepth = Graphics.FromImage(bmpDepth);

            dataTex = bmpTexture.LockBits(new Rectangle(0, 0, bmpTexture.Width, bmpTexture.Height), ImageLockMode.ReadOnly, PixelFormat.Format24bppRgb);

            
        }

        ~Scene3D()
        {
            
            bmpTexture.UnlockBits(dataTex);
            
        }

        public void Init(int width, int height)
        {
            width_ = width;
            height_ = height;
            zBuffer_ = new double[width, height];
        }

        public void ClearBuffer()
        {
            
            int i, j;
            for (i = 0; i < height_ * width_; i++)
                for (j = 0; j < width_; j++)
                    zBuffer_[j, i] = double.MaxValue;
        }

        public void AddPoint(int x, int y, int z)
        {
            points_.Add(new Vector(x, y, z));
        }

        public void AddPoint(Vector p)
        {
            points_.Add(p);
        }

        public void AddNormalPoint(int x, int y, int z)
        {
            this.pointsNormals_.Add(new Vector(x, y, z));
            this.pointsColors_.Add(new Vector());
        }

        public void NormalizeNormalPoints()
        {
            for(int i = 0; i <pointsNormals_.Count; i++)
            {
                pointsNormals_[i] = pointsNormals_[i].Normalize();
            }
        }

        public void AddTriangle(uint p1, uint p2, uint p3)
        {
            Triangle3D triangle = new Triangle3D(p1, p2, p3);
            AddTriangle(triangle);
        }

        public void AddTriangle(Triangle3D triangle)
        {
            triangles_.Add(triangle);
        }
        public Triangle3D GetTriangle(int i)
        {
            return triangles_[i];
        }

        Triangle2D ProjectTriangle(Triangle3D triangle, Matrix transformationMatrix)
        {
            return new Triangle2D(transformationMatrix * points_[(int)triangle.p1_],
                transformationMatrix * points_[(int)triangle.p2_],
                transformationMatrix * points_[(int)triangle.p3_]);
        }

        private int float2Int(float f)
        {
            return (int)(f+.5);
        }
        private int double2Int(double d)
        {
            return (int)(d + .5);
        }

        public bool isBack(Triangle2D t2)
        {
            Vector dir12 = t2.p1_ - t2.p2_;
            Vector dir13 = t2.p1_ - t2.p3_;
            dir12.z_ = 0;
            dir13.z_ = 0;
            return dir12.Cross(dir13).z_ < 0;
        }
        List<Triangle> clip(Triangle tr, Matrix mvp)
        {
            List<Triangle> triangles = new List<Triangle>();
            List<Vector> inPoints = new List<Vector>();
            List<Vector> outPoints = new List<Vector>();
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
                clipManager.SutherlandHodgmanPolygonClip(inPoints, outPoints, i);
                List<Vector> t = outPoints;
                outPoints = inPoints;
                inPoints = t;
                outPoints.Clear();
            }

            if (inPoints.Count < 3)
                return triangles;

            List<Vector> colors = new List<Vector>();
            List<Vector> texs = new List<Vector>();
            for (int i = 0; i < inPoints.Count; i++ )
            {
                Vector P = interpolator.transform * inPoints[i];

                double cl1 = (interpolator.p2y_p3y * (P.x_ - interpolator.p3_.x_) + (interpolator.p3x_p2x) * (P.y_ - interpolator.p3_.y_)) * interpolator.den;
                double cl2 = ((interpolator.p3y_p1y) * (P.x_ - interpolator.p3_.x_) + (interpolator.p1x_p3x) * (P.y_ - interpolator.p3_.y_)) * interpolator.den;
                double cl3 = 1 - cl1 - cl2;

                double x = tr.p1Color.x_ * cl1 + tr.p2Color.x_ * cl2 + tr.p3Color.x_ * cl3;
                double y = tr.p1Color.y_ * cl1 + tr.p2Color.y_ * cl2 + tr.p3Color.y_ * cl3;
                double z = tr.p1Color.z_ * cl1 + tr.p2Color.z_ * cl2 + tr.p3Color.z_ * cl3;

                Vector rP = texinterpolator.transform * inverseMVPTransform * inPoints[i];
                double tl1 = (texinterpolator.p2y_p3y * (rP.x_ - texinterpolator.p3_.x_) + (texinterpolator.p3x_p2x) * (rP.y_ - texinterpolator.p3_.y_)) * texinterpolator.den;
                double tl2 = ((texinterpolator.p3y_p1y) * (rP.x_ - texinterpolator.p3_.x_) + (texinterpolator.p1x_p3x) * (rP.y_ - texinterpolator.p3_.y_)) * texinterpolator.den;
                double tl3 = 1 - tl1 - tl2;

                double u = tr.p1TextureCordinates_.x_ * tl1 + tr.p2TextureCordinates_.x_ * tl2 + tr.p3TextureCordinates_.x_ * tl3;
                double v = tr.p1TextureCordinates_.y_ * tl1 + tr.p2TextureCordinates_.y_ * tl2 + tr.p3TextureCordinates_.y_ * tl3;
                
                colors.Add(new Vector(x, y, z));
                texs.Add(new Vector(u, v, 0));
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

        public void DrawProjectedTriangle(Triangle3D t, Matrix viewport, Matrix mvp, Camera camera)
        {
            Triangle2D t2 = ProjectTriangle(t, viewport*mvp);

            if (isBack(t2))//背面消隐
                return;

            Triangle tr = new Triangle(points_[(int)t.p1_], points_[(int)t.p2_],points_[(int)t.p3_]);
            tr.p1Color = pointsColors_[(int)t.p1_];
            tr.p2Color = pointsColors_[(int)t.p2_];
            tr.p3Color = pointsColors_[(int)t.p3_];
            tr.p1TextureCordinates_ = t.p1TextureCordinates_;
            tr.p2TextureCordinates_ = t.p2TextureCordinates_;
            tr.p3TextureCordinates_ = t.p3TextureCordinates_;
            List<Triangle> triangles = clip(tr, mvp);

            //切分好的三角形
            for (int i = 0; i < triangles.Count; i++)
            {
                Vector a = viewport * triangles[i].p1_;
                Vector b = viewport * triangles[i].p2_;
                Vector c = viewport * triangles[i].p3_;
                a.x_ = Math.Floor(a.x_); a.y_ = Math.Floor(a.y_);
                b.x_ = Math.Floor(b.x_); b.y_ = Math.Floor(b.y_);
                c.x_ = Math.Floor(c.x_); c.y_ = Math.Floor(c.y_);
                //设置颜色插值参数
                colorInterpolator.SetVector1(a, triangles[i].p1Color);
                colorInterpolator.SetVector2(b, triangles[i].p2Color);
                colorInterpolator.SetVector3(c, triangles[i].p3Color);
                colorInterpolator.PreCalculate();

                //设置文理坐标插值参数
                texInterpolator.SetVector1(inverseMVPVTransform * a, triangles[i].p1TextureCordinates_);
                texInterpolator.SetVector2(inverseMVPVTransform * b, triangles[i].p2TextureCordinates_);
                texInterpolator.SetVector3(inverseMVPVTransform * c, triangles[i].p3TextureCordinates_);
                texInterpolator.PreCalculate();

                if (showMode == 0)
                {
                    DDA_DrawLine(dataMain, double2Int(a.x_), double2Int(a.y_), double2Int(b.x_), double2Int(b.y_), triangles[i].p1Color, triangles[i].p2Color);
                    DDA_DrawLine(dataMain, double2Int(b.x_), double2Int(b.y_), double2Int(c.x_), double2Int(c.y_), triangles[i].p2Color, triangles[i].p3Color);
                    DDA_DrawLine(dataMain, double2Int(c.x_), double2Int(c.y_), double2Int(a.x_), double2Int(a.y_), triangles[i].p3Color, triangles[i].p1Color);

                    //DDA_DrawLine(data, double2Int(t2.p1_.x_), double2Int(t2.p1_.y_), double2Int(t2.p2_.x_), double2Int(t2.p2_.y_), pointsColors_[(int)t.p1_], pointsColors_[(int)t.p2_]);
                    //DDA_DrawLine(data, double2Int(t2.p2_.x_), double2Int(t2.p2_.y_), double2Int(t2.p3_.x_), double2Int(t2.p3_.y_), pointsColors_[(int)t.p2_], pointsColors_[(int)t.p3_]);
                    //DDA_DrawLine(data, double2Int(t2.p3_.x_), double2Int(t2.p3_.y_), double2Int(t2.p1_.x_), double2Int(t2.p1_.y_), pointsColors_[(int)t.p3_], pointsColors_[(int)t.p1_]);
                }
                else if (showMode == 1)
                {
                    DrawTriangle(dataMain, a, b, c, false);
                }
                else if (showMode == 2)
                {
                    DrawTriangle(dataMain, a, b, c, true);
                }
            }

        }
        private void DDA_DrawLine(BitmapData data, int p1x, int p1y, int p2x, int p2y, Vector color1, Vector color2)
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

                //graph.SetPixel(p1x, p1y, Color.FromArgb(c1R, c1G, c1B));
                // 设置起点像素颜色
                if(p1x >=0 && p1x < width_&& p1y >= 0 && p1y < height_)
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

                    if (ix >= 0 && ix < width_ && iy >= 0 && iy < height_)
                    {
                        ptr[(ix * 3) + iy * stride] = (byte)(c1B * (1 - beta) + c2B * beta);
                        ptr[(ix * 3) + iy * stride + 1] = (byte)(c1G * (1 - beta) + c2G * beta);
                        ptr[(ix * 3) + iy * stride + 2] = (byte)(c1R * (1 - beta) + c2R * beta);
                    }
                }
            }
        }
        
        public void DrawScene(Bitmap image, Matrix MV, Matrix P, Matrix View, Camera camera)
        {
            
            inverseMVPVTransform = (View * P * MV).Invert4x4Matrix();
            inverseMVPTransform = (P * MV).Invert4x4Matrix();
            if(isLight)
            {
                //计算光照
                for (int i = 0; i < pointsNormals_.Count; i++)
                {
                    Vector normal = pointsNormals_[i];
                    normal = normal.Normalize();
                    double ndot = lightDir.Dot(normal);

                    Vector lambert = diffuse * Math.Max(ndot, 0);

                    Vector halfv = (lightDir + camera.position).Normalize();
                    double spec = normal.Dot(halfv);
                    Vector specucolor = specular * lightColor * Math.Pow(Math.Max(spec, 0), shininess);
                    pointsColors_[i] = ambientLightColor + lambert + specucolor;

                    pointsColors_[i].x_ = Math.Min(pointsColors_[i].x_, 1.0);
                    pointsColors_[i].y_ = Math.Min(pointsColors_[i].y_, 1.0);
                    pointsColors_[i].z_ = Math.Min(pointsColors_[i].z_, 1.0);
                }
            }
            else
            {
                pointsColors_[0].x_ = 0; pointsColors_[0].y_ = 0; pointsColors_[0].z_ = 0;
                pointsColors_[1].x_ = 0; pointsColors_[1].y_ = 0; pointsColors_[1].z_ = 1;
                pointsColors_[2].x_ = 0; pointsColors_[2].y_ = 1; pointsColors_[2].z_ = 1;
                pointsColors_[3].x_ = 0; pointsColors_[3].y_ = 1; pointsColors_[3].z_ = 0;
                pointsColors_[4].x_ = 1; pointsColors_[4].y_ = 0; pointsColors_[4].z_ = 0;
                pointsColors_[5].x_ = 1; pointsColors_[5].y_ = 0; pointsColors_[5].z_ = 1;
                pointsColors_[6].x_ = 1; pointsColors_[6].y_ = 1; pointsColors_[6].z_ = 1;
                pointsColors_[7].x_ = 1; pointsColors_[7].y_ = 1; pointsColors_[7].z_ = 0;
            }

            for (int i = 0; i < triangles_.Count; i++)
            {
                DrawProjectedTriangle(triangles_[i], View, P*MV, camera);
            }
            

           
        }

        public void Render(int width, int height, Camera camera, Matrix worldTransform)
        {
            

            graphMain.Clear(Color.Black);
            graDepth.Clear(Color.White);
            graphMain.DrawString("W 切换显示模式（线框、填充、纹理）\nL 开关光照\n鼠标控制旋转和缩放", drawFont, drawBrush, 0, 0);
            
            int m = Math.Min(width, height);
            Matrix MV = camera.GetViewMatrix() * worldTransform;
            Matrix P = Matrix.CreatePerspectiveProjectionMatrix(1.1, 1, 1, 15);
            Matrix View = Matrix.CreateViewportMatrix(width / 2, height / 2, m, m);
            
            dataMain = bmpMain.LockBits(new Rectangle(0, 0, bmpMain.Width, bmpMain.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            dataDepth = bmpDepth.LockBits(new Rectangle(0, 0, bmpDepth.Width, bmpDepth.Height), ImageLockMode.ReadWrite, PixelFormat.Format24bppRgb);
            
            DrawScene(bmpMain, MV, P, View, camera);
            watch.Reset();
            watch.Start();
            bmpMain.UnlockBits(dataMain);
            bmpDepth.UnlockBits(dataDepth);
            watch.Stop();
           
        }
        static void Swap<T>(ref T a, ref T b)
        {
            T t = a;
            a = b;
            b = t;
        }
        void SortVertices(ref Vector top, ref Vector middle, ref Vector bottom)
        {
            if (top.y_ > middle.y_)
                Swap(ref top, ref middle);

            if (middle.y_ > bottom.y_)
                Swap(ref middle, ref bottom);

            if (top.y_ > middle.y_)
                Swap(ref top, ref middle);
        }

        double CalculateDeltaY(double y, double ystart, double yend)
        {
            return (ystart == yend) ? 1 : (y - yend) / (ystart - yend);
        }
        double zinterpolator(double betay, double start, double end)
        {
            return betay * start + (1 - betay) * end;
        }
        void DrawTriangle(BitmapData data, Vector p1, Vector p2, Vector p3, bool beTexture)
        {  
            Vector top = p1, middle = p2, bottom = p3;
            SortVertices(ref top, ref middle, ref bottom);

            // inv: p1.y <= p2.y <= p3.y
            double y = middle.y_;
            double betay = CalculateDeltaY(y, top.y_, bottom.y_);
            double xr = betay * top.x_ + (1 - betay) * bottom.x_;
            double zr = betay * top.z_ + (1 - betay) * bottom.z_;

            
            DrawTriangleWithXParellGround(data, top, new Vector(xr, y, zr), middle, beTexture);
            DrawTriangleWithXParellGround(data, bottom, new Vector(xr, y, zr), middle, beTexture);
            //watch.Stop();

        }

        void DrawTriangleWithXParellGround(BitmapData data, Vector p1, Vector p2, Vector p3, bool beTexture)
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

            int texWidth = bmpTexture.Width;
            int texHeight = bmpTexture.Height;

            int direction = ystart < yend ? 1 : -1;
            if(direction == -1)
            {
                leftDxDivDy = -leftDxDivDy;
                rightDxDivDy = -rightDxDivDy;
            }
            StVector rgb;
            int stride = data.Stride;
            int texStride = dataTex.Stride;
            int depthStride = dataDepth.Stride;
            Matrix texInterpolatorTransform = texInterpolator.transform * inverseMVPVTransform;
            unsafe
            {
                byte* ptr = (byte*)data.Scan0;
                byte* texPtr = (byte*)dataTex.Scan0;
                byte* depthPtr = (byte*)dataDepth.Scan0;
                for (int y = yend; y != ystart - direction; y -= direction)
                {
                    int xLeft_int = (int)xLeft;
                    int xRight_int = (int)(xRight)+1;
                    double betay = (ystart == yend) ? 1 : ((double)y - yend) / (ystart - yend);
                    double zl = betay * p1.z_ + (1 - betay) * p2.z_;
                    double zr = betay * p1.z_ + (1 - betay) * p3.z_;
                    
                    for (int x = xLeft_int; x < xRight_int; x++)
                    {
                        double betax = (xLeft_int == xRight_int) ? 1 : ((double)x - xRight_int) / (xLeft_int - xRight_int);
                        double z = betax * zl + (1 - betax) * zr;

                        if (x >= 0 && y >= 0 && x < width_ && y < height_ && z * 255 < depthPtr[(x * 3) + y * depthStride])
                        {
                            depthPtr[(x * 3) + y * depthStride] = (byte)(z*255);
                            StVector cur;
                            cur.x_ = x; cur.y_ = y; cur.z_ = 1;

                            if (beTexture)
                            {
                                //纹理渲染
                                StVector projectP;
                                projectP.x_ = x;
                                projectP.y_ = y;
                                projectP.z_ = z;
                                StVector P = texInterpolatorTransform * projectP;

                                double tl1 = (texInterpolator.p2y_p3y * (P.x_ - texInterpolator.p3_.x_) + (texInterpolator.p3x_p2x) * (P.y_ - texInterpolator.p3_.y_)) * texInterpolator.den;
                                double tl2 = ((texInterpolator.p3y_p1y) * (P.x_ - texInterpolator.p3_.x_) + (texInterpolator.p1x_p3x) * (P.y_ - texInterpolator.p3_.y_)) * texInterpolator.den;
                                double tl3 = 1 - tl1 - tl2;

                                int u = (int)((texInterpolator.p1Value_.x_ * tl1 + texInterpolator.p2Value_.x_ * tl2 + texInterpolator.p3Value_.x_ * tl3) * texWidth + 0.5);
                                int v = (int)((texInterpolator.p1Value_.y_ * tl1 + texInterpolator.p2Value_.y_ * tl2 + texInterpolator.p3Value_.y_ * tl3) * texHeight + 0.5);
                                u = u < 0 ? 0 : u;
                                u = u >= texWidth ? texWidth-1 : u;
                                v = v < 0 ? 0 : v;
                                v = v >= texHeight ? texHeight-1 : v;

                                ptr[(x * 3) + y * stride] = texPtr[(u * 3) + v * texStride];
                                ptr[(x * 3) + y * stride + 1] = texPtr[(u * 3) + v * texStride + 1];
                                ptr[(x * 3) + y * stride + 2] = texPtr[(u * 3) + v * texStride + 2];
                            }
                            else
                            {
                                //颜色渲染
                                double l1 = (colorInterpolator.p2y_p3y * (cur.x_ - colorInterpolator.p3_.x_) + (colorInterpolator.p3x_p2x) * (cur.y_ - colorInterpolator.p3_.y_)) * colorInterpolator.den;
                                double l2 = ((colorInterpolator.p3y_p1y) * (cur.x_ - colorInterpolator.p3_.x_) + (colorInterpolator.p1x_p3x) * (cur.y_ - colorInterpolator.p3_.y_)) * colorInterpolator.den;
                                double l3 = 1 - l1 - l2;

                                rgb.x_ = colorInterpolator.p1Value_.x_ * l1 + colorInterpolator.p2Value_.x_ * l2 + colorInterpolator.p3Value_.x_ * l3;
                                rgb.y_ = colorInterpolator.p1Value_.y_ * l1 + colorInterpolator.p2Value_.y_ * l2 + colorInterpolator.p3Value_.y_ * l3;
                                rgb.z_ = colorInterpolator.p1Value_.z_ * l1 + colorInterpolator.p2Value_.z_ * l2 + colorInterpolator.p3Value_.z_ * l3;
                                rgb.x_ = rgb.x_ > 1 ? 1 : rgb.x_;
                                rgb.y_ = rgb.y_ > 1 ? 1 : rgb.y_;
                                rgb.z_ = rgb.z_ > 1 ? 1 : rgb.z_;

                                rgb.x_ = rgb.x_ < 0 ? 0 : rgb.x_;
                                rgb.y_ = rgb.y_ < 0 ? 0 : rgb.y_;
                                rgb.z_ = rgb.z_ < 0 ? 0 : rgb.z_;

                                ptr[(x * 3) + y * stride] = (byte)(rgb.z_ * 255);
                                ptr[(x * 3) + y * stride + 1] = (byte)(rgb.y_ * 255);
                                ptr[(x * 3) + y * stride + 2] = (byte)(rgb.x_ * 255);

                                //ptr[(x * 3) + y * stride] = (byte)(255);
                                //ptr[(x * 3) + y * stride + 1] = (byte)(255);
                                //ptr[(x * 3) + y * stride + 2] = (byte)(255);

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
