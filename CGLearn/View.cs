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
using System.Diagnostics;

namespace CGLearn
{
    public partial class View : Form
    {        
        Camera m_camera = null;
        Scene3D m_scene;
        
        Matrix4 m_mWorldTransform = Matrix4.CreateIdentityMatrix();

        bool m_bIsPress = false;
        int m_iPreX = 0;
        int m_iPreY = 0;

        Font m_drawFont = new Font("宋体", 16);
        SolidBrush m_drawBrush = new SolidBrush(Color.White);

        Stopwatch watch = new Stopwatch();
        public View()
        {
            InitializeComponent();
            this.SetStyle(ControlStyles.OptimizedDoubleBuffer, true);
            this.SetStyle(ControlStyles.AllPaintingInWmPaint, true);
            this.SetStyle(ControlStyles.UserPaint, true);

            m_scene = new Scene3D(this.ClientRectangle.Width, this.ClientRectangle.Height);
            
            //顶点
            m_scene.AddPoint(1, -1, -1);
            m_scene.AddPoint(1, -1, 1);
            m_scene.AddPoint(-1, -1, 1);
            m_scene.AddPoint(-1, -1, -1);
            m_scene.AddPoint(1, 1, -1);
            m_scene.AddPoint(1, 1, 1);
            m_scene.AddPoint(-1, 1, 1);
            m_scene.AddPoint(-1, 1, -1);
            //法线
            m_scene.AddNormalPoint(1, -1, -1);
            m_scene.AddNormalPoint(1, -1, 1);
            m_scene.AddNormalPoint(-1, -1, 1);
            m_scene.AddNormalPoint(-1, -1, -1);
            m_scene.AddNormalPoint(1, 1, -1);
            m_scene.AddNormalPoint(1, 1, 1);
            m_scene.AddNormalPoint(-1, 1, 1);
            m_scene.AddNormalPoint(-1, 1, -1);
            m_scene.NormalizeNormalPoints();
            //三角形和文理坐标
            m_scene.AddTriangle(0, 1, 2); m_scene.GetTriangle(0).SetTextureCoordinates(new Vector3(0, 0, 0), new Vector3(0, 1, 0), new Vector3(1, 1, 0));
            m_scene.AddTriangle(0, 2, 3); m_scene.GetTriangle(1).SetTextureCoordinates(new Vector3(0, 0, 0), new Vector3(1, 1, 0), new Vector3(1, 0, 0));

            m_scene.AddTriangle(0, 5, 1); m_scene.GetTriangle(2).SetTextureCoordinates(new Vector3(1, 0, 0), new Vector3(0, 1, 0), new Vector3(1, 1, 0));
            m_scene.AddTriangle(0, 4, 5); m_scene.GetTriangle(3).SetTextureCoordinates(new Vector3(1, 0, 0), new Vector3(0, 0, 0), new Vector3(0, 1, 0));

            m_scene.AddTriangle(4, 6, 5); m_scene.GetTriangle(4).SetTextureCoordinates(new Vector3(1, 0, 0), new Vector3(0, 1, 0), new Vector3(1, 1, 0));
            m_scene.AddTriangle(4, 7, 6); m_scene.GetTriangle(5).SetTextureCoordinates(new Vector3(1, 0, 0), new Vector3(0, 0, 0), new Vector3(0, 1, 0));

            m_scene.AddTriangle(7, 2, 6); m_scene.GetTriangle(6).SetTextureCoordinates(new Vector3(1, 0, 0), new Vector3(0, 1, 0), new Vector3(1, 1, 0));
            m_scene.AddTriangle(7, 3, 2); m_scene.GetTriangle(7).SetTextureCoordinates(new Vector3(1, 0, 0), new Vector3(0, 0, 0), new Vector3(0, 1, 0));

            m_scene.AddTriangle(1, 5, 6); m_scene.GetTriangle(8).SetTextureCoordinates(new Vector3(1, 1, 0), new Vector3(1, 0, 0), new Vector3(0, 0, 0));
            m_scene.AddTriangle(1, 6, 2); m_scene.GetTriangle(9).SetTextureCoordinates(new Vector3(1, 1, 0), new Vector3(0, 0, 0), new Vector3(0, 1, 0));

            m_scene.AddTriangle(0, 7, 4); m_scene.GetTriangle(10).SetTextureCoordinates(new Vector3(1, 0, 0), new Vector3(0, 1, 0), new Vector3(1, 1, 0));
            m_scene.AddTriangle(0, 3, 7); m_scene.GetTriangle(11).SetTextureCoordinates(new Vector3(1, 0, 0), new Vector3(0, 0, 0), new Vector3(0, 1, 0));

            m_camera = new Camera(new Vector3(0, 0, -5), new Vector3(0,1,0), new Vector3(0,0,0), 0, 100);

            m_scene.Render(this.ClientRectangle.Width, this.ClientRectangle.Height, m_camera, m_mWorldTransform);
 
            
        }

        private void View_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            g.DrawImage(m_scene.GetBmpMain(), this.ClientRectangle);
            g.DrawString("W 切换显示模式（线框、填充、纹理）\nL 开关光照\n鼠标控制旋转和缩放", m_drawFont, m_drawBrush, 0, 0);
        }

        private void View_MouseMove(object sender, MouseEventArgs e)
        {
            if(m_bIsPress)
            {
                int deltaX = e.X - m_iPreX;
                int deltaY = e.Y - m_iPreY;
                
                    m_mWorldTransform = Matrix4.CreateYAxisRotationMatrix(deltaX / 400.0) * m_mWorldTransform;
                    m_mWorldTransform = Matrix4.CreateXAxisRotationMatrix(deltaY / 400.0) * m_mWorldTransform;

                    m_scene.Render(this.ClientRectangle.Width, this.ClientRectangle.Height, m_camera, m_mWorldTransform);
                    Invalidate(false);
                    m_iPreX = e.X;
                    m_iPreY = e.Y;
                
                
            }
        }

        private void View_MouseDown(object sender, MouseEventArgs e)
        {
            m_iPreX = e.X;
            m_iPreY = e.Y;
            m_bIsPress = true;
        }

        private void View_MouseUp(object sender, MouseEventArgs e)
        {
            m_bIsPress = false;
        }

        private void View_KeyDown(object sender, KeyEventArgs e)
        {
            if(e.KeyCode == Keys.W)
            {
                m_scene.SwitchShowMode();
                m_scene.Render(this.ClientRectangle.Width, this.ClientRectangle.Height, m_camera, m_mWorldTransform);
                Invalidate(false);
            }
            else if (e.KeyCode == Keys.L)
            {
                m_scene.SwitchLight();
                m_scene.Render(this.ClientRectangle.Width, this.ClientRectangle.Height, m_camera, m_mWorldTransform);
                Invalidate(false);
            }
        }

        private void View_MouseWheel(object sender, MouseEventArgs e)
        {
            if(e.Delta > 0)
            {
                if(m_camera.position.z_ < -2.2)
                    m_camera.position.z_ += .3;
                m_scene.Render(this.ClientRectangle.Width, this.ClientRectangle.Height, m_camera, m_mWorldTransform);
                Invalidate(false);
            }
            else
            {
                if (m_camera.position.z_ > -20)
                    m_camera.position.z_ -= .3;
                m_scene.Render(this.ClientRectangle.Width, this.ClientRectangle.Height, m_camera, m_mWorldTransform);
                Invalidate(false);
            }
        }
    }
}
