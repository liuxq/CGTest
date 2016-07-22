namespace CGLearn
{
    partial class View
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.SuspendLayout();
            // 
            // View
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 12F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(813, 579);
            this.Name = "View";
            this.Text = "View";
            this.Paint += new System.Windows.Forms.PaintEventHandler(this.View_Paint);
            this.KeyDown += new System.Windows.Forms.KeyEventHandler(this.View_KeyDown);
            this.MouseDown += new System.Windows.Forms.MouseEventHandler(this.View_MouseDown);
            this.MouseMove += new System.Windows.Forms.MouseEventHandler(this.View_MouseMove);
            this.MouseUp += new System.Windows.Forms.MouseEventHandler(this.View_MouseUp);
            this.MouseWheel += new System.Windows.Forms.MouseEventHandler(this.View_MouseWheel);
            this.ResumeLayout(false);

        }

        #endregion
    }
}

