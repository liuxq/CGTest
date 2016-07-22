using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using CGLearn.CG;

namespace CGLearn
{
    public partial class View : Form
    {
        private Color backColor = Color.Black;
        private Bitmap theImage;
        private Camera camera = null;
        Scene3D scene = new Scene3D();
        Graphics graph;
        Matrix worldTransform = Matrix.CreateIdentityMatrix(4);

        private bool isPress = false;
        private int preX = 0;
        private int preY = 0;

        public View()
        {
            InitializeComponent();
            this.SetStyle(ControlStyles.OptimizedDoubleBuffer, true);
            this.SetStyle(ControlStyles.AllPaintingInWmPaint, true);
            this.SetStyle(ControlStyles.UserPaint, true);

            theImage = new Bitmap(this.ClientRectangle.Width, this.ClientRectangle.Height);
            graph = Graphics.FromImage(theImage);
            
            //顶点
            scene.AddPoint(1, -1, -1);
            scene.AddPoint(1, -1, 1);
            scene.AddPoint(-1, -1, 1);
            scene.AddPoint(-1, -1, -1);
            scene.AddPoint(1, 1, -1);
            scene.AddPoint(1, 1, 1);
            scene.AddPoint(-1, 1, 1);
            scene.AddPoint(-1, 1, -1);
            //法线
            scene.AddNormalPoint(1, -1, -1);
            scene.AddNormalPoint(1, -1, 1);
            scene.AddNormalPoint(-1, -1, 1);
            scene.AddNormalPoint(-1, -1, -1);
            scene.AddNormalPoint(1, 1, -1);
            scene.AddNormalPoint(1, 1, 1);
            scene.AddNormalPoint(-1, 1, 1);
            scene.AddNormalPoint(-1, 1, -1);
            scene.NormalizeNormalPoints();
            //三角形和文理坐标
            scene.AddTriangle(0, 1, 2); scene.GetTriangle(0).SetTextureCoordinates(new Vector(0, 0, 0), new Vector(0, 1, 0), new Vector(1, 1, 0));
            scene.AddTriangle(0, 2, 3); scene.GetTriangle(1).SetTextureCoordinates(new Vector(0, 0, 0), new Vector(1, 1, 0), new Vector(1, 0, 0));

            scene.AddTriangle(0, 5, 1); scene.GetTriangle(2).SetTextureCoordinates(new Vector(1, 0, 0), new Vector(0, 1, 0), new Vector(1, 1, 0));
            scene.AddTriangle(0, 4, 5); scene.GetTriangle(3).SetTextureCoordinates(new Vector(1, 0, 0), new Vector(0, 0, 0), new Vector(0, 1, 0));

            scene.AddTriangle(4, 6, 5); scene.GetTriangle(4).SetTextureCoordinates(new Vector(1, 0, 0), new Vector(0, 1, 0), new Vector(1, 1, 0));
            scene.AddTriangle(4, 7, 6); scene.GetTriangle(5).SetTextureCoordinates(new Vector(1, 0, 0), new Vector(0, 0, 0), new Vector(0, 1, 0));

            scene.AddTriangle(7, 2, 6); scene.GetTriangle(6).SetTextureCoordinates(new Vector(1, 0, 0), new Vector(0, 1, 0), new Vector(1, 1, 0));
            scene.AddTriangle(7, 3, 2); scene.GetTriangle(7).SetTextureCoordinates(new Vector(1, 0, 0), new Vector(0, 0, 0), new Vector(0, 1, 0));

            scene.AddTriangle(1, 5, 6); scene.GetTriangle(8).SetTextureCoordinates(new Vector(1, 1, 0), new Vector(1, 0, 0), new Vector(0, 0, 0));
            scene.AddTriangle(1, 6, 2); scene.GetTriangle(9).SetTextureCoordinates(new Vector(1, 1, 0), new Vector(0, 0, 0), new Vector(0, 1, 0));

            scene.AddTriangle(0, 7, 4); scene.GetTriangle(10).SetTextureCoordinates(new Vector(1, 0, 0), new Vector(0, 1, 0), new Vector(1, 1, 0));
            scene.AddTriangle(0, 3, 7); scene.GetTriangle(11).SetTextureCoordinates(new Vector(1, 0, 0), new Vector(0, 0, 0), new Vector(0, 1, 0));

            camera = new Camera(new Vector(0, 0, 5), new Vector(0,1,0), new Vector(0,0,0), 0, 100);

            scene.Init(this.ClientRectangle.Width, this.ClientRectangle.Height);
            scene.Render(graph, theImage, this.ClientRectangle.Width, this.ClientRectangle.Height, camera, worldTransform);
            
        }

        private void View_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            g.Clear(backColor);
            //创建一个Bitmap  
            
            g.DrawImage(theImage, this.ClientRectangle);
            
        }

        private void View_MouseMove(object sender, MouseEventArgs e)
        {
            if(isPress)
            {
                int deltaX = e.X - preX;
                int deltaY = e.Y - preY;
                if(Math.Abs(deltaX) > 2 || Math.Abs(deltaY) > 2)
                {
                    worldTransform = Matrix.CreateYAxisRotationMatrix(deltaX / 400.0) * worldTransform;
                    worldTransform = Matrix.CreateXAxisRotationMatrix(-deltaY / 400.0) * worldTransform;

                    scene.Render(graph, theImage, this.ClientRectangle.Width, this.ClientRectangle.Height, camera, worldTransform);
                    Invalidate(false);
                    preX = e.X;
                    preY = e.Y;
                }
                
            }
        }

        private void View_MouseDown(object sender, MouseEventArgs e)
        {
            preX = e.X;
            preY = e.Y;
            isPress = true;
        }

        private void View_MouseUp(object sender, MouseEventArgs e)
        {
            isPress = false;
        }

        private void View_KeyDown(object sender, KeyEventArgs e)
        {
            if(e.KeyCode == Keys.W)
            {
                scene.SwitchShowMode();
                scene.Render(graph, theImage, this.ClientRectangle.Width, this.ClientRectangle.Height, camera, worldTransform);
                Invalidate(false);
            }
            else if(e.KeyCode == Keys.C)
            {
                scene.SwitchCameraMode();
                scene.Render(graph, theImage, this.ClientRectangle.Width, this.ClientRectangle.Height, camera, worldTransform);
                Invalidate(false);
            }
            else if (e.KeyCode == Keys.L)
            {
                scene.SwitchLight();
                scene.Render(graph, theImage, this.ClientRectangle.Width, this.ClientRectangle.Height, camera, worldTransform);
                Invalidate(false);
            }
        }

        private void View_MouseWheel(object sender, MouseEventArgs e)
        {
            if(e.Delta > 0)
            {
                scene.ChangeScale(1.1);
                scene.Render(graph, theImage, this.ClientRectangle.Width, this.ClientRectangle.Height, camera, worldTransform);
                Invalidate(false);
            }
            else
            {
                scene.ChangeScale(0.9);
                scene.Render(graph, theImage, this.ClientRectangle.Width, this.ClientRectangle.Height, camera, worldTransform);
                Invalidate(false);
            }
        }
    }
}
