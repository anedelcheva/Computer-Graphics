#include <GL/glfw.h>
#include <cmath>
#define M_PI 3.141592

bool running()
{
    return (!glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam(GLFW_OPENED));
}
///////////////////////////////////////////////////////////////////////////////
void init()
{
    int width, height;

    glfwInit();
    if( !glfwOpenWindow( 640, 480, 0, 0, 0, 0, 8, 0, GLFW_WINDOW ) ) return;

    glfwGetWindowSize( &width, &height );
    height = height > 0 ? height : 1;

    glViewport( 0, 0, width, height );
    glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 65.0f, (GLfloat)width/(GLfloat)height, 1.0f, 100.0f );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    gluLookAt(0.0f, -7.0f, 0.0f,
            0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 1.0f );
}
///////////////////////////////////////////////////////////////////////////////
void finit()
{
    glfwTerminate();
}
///////////////////////////////////////////////////////////////////////////////
void drawParallelepiped(float x, float y, float z, float a, float b, float c)
{
    // -Y
    glBegin(GL_POLYGON);
        glNormal3f(0.0, +1.0, 0.0);
        glVertex3f(x-a, y+b, z+c);
        glVertex3f(x+a, y+b, z+c);
        glVertex3f(x+a, y+b, z-c);
        glVertex3f(x-a, y+b, z-c);
    glEnd();

    // Y
    glBegin(GL_POLYGON);
        glNormal3f(0.0, -1.0, 0.0);
        glVertex3f(x-a, y-b, z+c);
        glVertex3f(x+a, y-b, z+c);
        glVertex3f(x+a, y-b, z-c);
        glVertex3f(x-a, y-b, z-c);
    glEnd();

    // X
    glBegin(GL_POLYGON);
        glNormal3f(-1.0, 0.0, 0.0);
        glVertex3f(x-a, y+b, z+c);
        glVertex3f(x-a, y-b, z+c);
        glVertex3f(x-a, y-b, z-c);
        glVertex3f(x-a, y+b, z-c);
    glEnd();

    // -X
    glBegin(GL_POLYGON);
        glNormal3f(+1.0, 0.0, 0.0);
        glVertex3f(x+a, y+b, z+c);
        glVertex3f(x+a, y-b, z+c);
        glVertex3f(x+a, y-b, z-c);
        glVertex3f(x+a, y+b, z-c);
    glEnd();

    // Z
    glBegin(GL_POLYGON);
        glNormal3f(0.0, 0.0, -1.0);
        glVertex3f(x-a, y+b, z-c);
        glVertex3f(x-a, y-b, z-c);
        glVertex3f(x+a, y-b, z-c);
        glVertex3f(x+a, y+b, z-c);
    glEnd();

    // -Z
    glBegin(GL_POLYGON);
        glNormal3f(0.0, 0.0, +1.0);
        glVertex3f(x-a, y+b, z+c);
        glVertex3f(x-a, y-b, z+c);
        glVertex3f(x+a, y-b, z+c);
        glVertex3f(x+a, y+b, z+c);
    glEnd();
}
/////////////////////////////////////////////////////////////////////////////////////////////////
void drawTwoCylinders(float x, float y, float z, float r, float R, float h, float from, float to)
{
    int n = 64;
    //float alpha = 0.0;
    //float alpha = 2*M_PI/2;
    float alpha = from;
    float angle = M_PI/n/to;

    for(int i = 0; i < n; i++)
    {
        float dx1r = r*cos(alpha);
        float dx1R = R*cos(alpha);

        float dz1r = r*sin(alpha);
        float dz1R = R*sin(alpha);

        float dx2r = r*cos(alpha+angle);
        float dx2R = R*cos(alpha+angle);

        float dz2r = r*sin(alpha+angle);
        float dz2R = R*sin(alpha+angle);

        // рисуваме долната основа
        glBegin(GL_POLYGON);
            glNormal3f( 0, -1, 0 );
            glVertex3f(x+dx1r, y, z+dz1r);
            glVertex3f(x+dx1R, y, z+dz1R);
            glVertex3f(x+dx2R, y, z+dz2R);
            glVertex3f(x+dx2r, y, z+dz2r);
        glEnd();

        // рисуваме горната основа
        glBegin(GL_POLYGON);
            glNormal3f( 0, 1, 0 );
            glVertex3f(x+dx1r, y+h, z+dz1r);
            glVertex3f(x+dx1R, y+h, z+dz1R);
            glVertex3f(x+dx2R, y+h, z+dz2R);
            glVertex3f(x+dx2r, y+h, z+dz2r);
        glEnd();

        // рисуваме вътрешната околна стена
        glBegin(GL_POLYGON);
            glNormal3f( -cos(alpha),  0, -sin(alpha));
            glVertex3f(x+dx1r, y, z+dz1r);
            glVertex3f(x+dx1r, y+h, z+dz1r);
            glNormal3f( -cos(alpha+angle),  0, -sin(alpha+angle));
            glVertex3f(x+dx2r, y+h, z+dz2r);
            glVertex3f(x+dx2r, y, z+dz2r);
        glEnd();

        //рисуваме външната околна стена
        glBegin(GL_POLYGON);
            glNormal3f( cos(alpha),  0, sin(alpha));
            glVertex3f(x+dx1R, y, z+dz1R);
            glVertex3f(x+dx1R, y+h, z+dz1R);
            glNormal3f( cos(alpha+angle),  0, sin(alpha+angle));
            glVertex3f(x+dx2R, y+h, z+dz2R);
            glVertex3f(x+dx2R, y, z+dz2R);
        glEnd();
        alpha += angle;
    }
}
///////////////////////////////////////////////////////////////////
int main()
{
    init();

    glClearColor( 0.5f, 0.7f, 1.0f, 1.0f );
    glEnable( GL_DEPTH_TEST );
    glEnable( GL_LIGHTING );
    glEnable( GL_COLOR_MATERIAL );
    glEnable( GL_LIGHT0 );

    glEnable(GL_SMOOTH);
	glShadeModel(GL_SMOOTH);

    while(running())
    {
        glClear(GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT);
        glRotatef( 0.01, 0.4, -0.1, 0.7);

        glColor3ub(255, 255, 0);

        drawParallelepiped(0, 0, 0.1, 0.85, 0.5, 0.15);
        drawParallelepiped(0, 0, 0, 0.15, 0.5, 0.85);
        drawParallelepiped(0, 0, 0.71, 0.8, 0.5, 0.15);
        drawParallelepiped(0, 0, 1.4, 0.8, 0.5, 0.15);
        drawParallelepiped(0, 0, 2.1, 0.8, 0.5, 0.15);
        drawParallelepiped(0, 0, 2.8, 0.8, 0.5, 0.15);
        drawParallelepiped(0.65, 0, 1.1, 0.15, 0.5, 0.45);
        drawParallelepiped(-0.65, 0, 1.1, 0.15, 0.5, 0.45);
        drawParallelepiped(0.65, 0, 1.7, 0.15, 0.5, 0.4);
        drawParallelepiped(-0.65, 0, 1.7, 0.15, 0.5, 0.4);
        drawParallelepiped(0, 0, 2.7, 0.15, 0.5, 0.65);
        drawParallelepiped(2, 0, 0.435, 0.15, 0.5, 1.3);
        drawParallelepiped(2, 0, 0.6, 1.0, 0.5, 0.15);
        drawParallelepiped(2, 0, 1.6, 0.5, 0.5, 0.15);

        drawTwoCylinders(3.5, -0.5, 3.6, 1.5, 1.8, 1, 2*M_PI/1.85, 3);
        drawTwoCylinders(0.3, -0.5, 3.4, 1.5, 1.8, 1, 3.3*M_PI/2, 3);

        glBegin(GL_POLYGON);
            glVertex3f(3.5 + 1.5*cos(2*M_PI/1.85 + 65*(M_PI/65/3)), -0.5, 3.6 + 1.5*sin(2*M_PI/1.85 + 65*(M_PI/65/3)));
            glVertex3f(3.5 + 1.8*cos(2*M_PI/1.85 + 65*(M_PI/65/3)), -0.5, 3.6 + 1.8*sin(2*M_PI/1.85 + 65*(M_PI/65/3)));
            glVertex3f(3.5 + 1.8*cos(2*M_PI/1.85 + 65*(M_PI/65/3)), -0.5 + 1, 3.6 + 1.8*sin(2*M_PI/1.85 + 65*(M_PI/65/3)));
            glVertex3f(3.5 + 1.5*cos(2*M_PI/1.85 + 65*(M_PI/65/3)), -0.5 + 1, 3.6 + 1.5*sin(2*M_PI/1.85 + 65*(M_PI/65/3)));
        glEnd();

        glBegin(GL_POLYGON);
            glNormal3f( 0, -cos(2*M_PI/1.85), 0 );
            glVertex3f(3.5 + 1.5*cos(2*M_PI/1.85), -0.5, 3.6 + 1.5*sin(2*M_PI/1.85));
            glVertex3f(3.5 + 1.8*cos(2*M_PI/1.85), -0.5, 3.6 + 1.8*sin(2*M_PI/1.85));
            glVertex3f(3.5 + 1.8*cos(2*M_PI/1.85), -0.5 + 1, 3.6 + 1.8*sin(2*M_PI/1.85));
            glVertex3f(3.5 + 1.5*cos(2*M_PI/1.85), -0.5 + 1, 3.6 + 1.5*sin(2*M_PI/1.85));
        glEnd();

        glBegin(GL_POLYGON);
            //glNormal3f(0, cos(3.3*M_PI/2 + M_PI), 0);
            glVertex3f(0.3 + 1.5*cos(3.3*M_PI/2), -0.5, 3.4 + 1.5*sin(3.3*M_PI/2));
            glVertex3f(0.3 + 1.8*cos(3.3*M_PI/2), -0.5, 3.4 + 1.8*sin(3.3*M_PI/2));
            glVertex3f(0.3 + 1.8*cos(3.3*M_PI/2), -0.5 + 1, 3.4 + 1.8*sin(3.3*M_PI/2));
            glVertex3f(0.3 + 1.5*cos(3.3*M_PI/2), -0.5 + 1, 3.4 + 1.5*sin(3.3*M_PI/2));
        glEnd();

        glfwSwapBuffers();
    }

    finit();

    return 0;
}
