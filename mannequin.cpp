/* Описание на сцената:

    Сцената описва двама влюбени. Младежът предлага брак на дамата, като коленичи пред нея, държащ в ръката си годежен
    пръстен. Дамата леко кима с глава в знак на съгласие. Тя е протегнала ръката си, за да вземе пръстена.

*/

#include <math.h>
#include <GL/glfw.h>

#define nL 21	// брой точки в крайниците
#define nH 50	// четен брой точки в главата
#define nJ 16	// четен брой точки в ставите

#define PI 3.141592653589793
#define EPS 0.1

#define turnForward(a)	glRotatef((a),1,0,0)
#define turnSideway(a)	glRotatef((a),0,1,0)
#define turnAround(a)	glRotatef((a),0,0,1)

void drawParallelepiped(float x, float y, float z, float a, float b, float c)
{
    // рисуване на стена откъм -Y
    glBegin(GL_POLYGON);
        glNormal3f(0.0, +1.0, 0.0);
        glVertex3f(x-a, y+b, z+c);
        glVertex3f(x+a, y+b, z+c);
        glVertex3f(x+a, y+b, z-c);
        glVertex3f(x-a, y+b, z-c);
    glEnd();

    // рисуване на стена откъм Y
    glBegin(GL_POLYGON);
        glNormal3f(0.0, -1.0, 0.0);
        glVertex3f(x-a, y-b, z+c);
        glVertex3f(x+a, y-b, z+c);
        glVertex3f(x+a, y-b, z-c);
        glVertex3f(x-a, y-b, z-c);
    glEnd();

    // рисуване на стена откъм X
    glBegin(GL_POLYGON);
        glNormal3f(-1.0, 0.0, 0.0);
        glVertex3f(x-a, y+b, z+c);
        glVertex3f(x-a, y-b, z+c);
        glVertex3f(x-a, y-b, z-c);
        glVertex3f(x-a, y+b, z-c);
    glEnd();

    // рисуване на стена откъм -X
    glBegin(GL_POLYGON);
        glNormal3f(+1.0, 0.0, 0.0);
        glVertex3f(x+a, y+b, z+c);
        glVertex3f(x+a, y-b, z+c);
        glVertex3f(x+a, y-b, z-c);
        glVertex3f(x+a, y+b, z-c);
    glEnd();

    // рисуване на стена откъм Z
    glBegin(GL_POLYGON);
        glNormal3f(0.0, 0.0, -1.0);
        glVertex3f(x-a, y+b, z-c);
        glVertex3f(x-a, y-b, z-c);
        glVertex3f(x+a, y-b, z-c);
        glVertex3f(x+a, y+b, z-c);
    glEnd();

    // рисуване на стена откъм -Z
    glBegin(GL_POLYGON);
        glNormal3f(0.0, 0.0, +1.0);
        glVertex3f(x-a, y+b, z+c);
        glVertex3f(x-a, y-b, z+c);
        glVertex3f(x+a, y-b, z+c);
        glVertex3f(x+a, y+b, z+c);
    glEnd();
}

// дефиниция на 3D точка или вектор
struct OKG_POINT
{
    float x;
    float y;
    float z;
} ;

void drawTwoCylinders(float x, float y, float z, float r, float R, float h)
{
    int n = 64;
    float alpha = 0.0;
    float angle = 2*PI/n;

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
            glNormal3f( cos(alpha),  0, sin(alpha));
            glVertex3f(x+dx1r, y, z+dz1r);
            glVertex3f(x+dx1r, y+h, z+dz1r);
            glNormal3f( cos(alpha+angle),  0, sin(alpha+angle));
            glVertex3f(x+dx2r, y+h, z+dz2r);
            glVertex3f(x+dx2r, y, z+dz2r);
        glEnd();

        //рисуваме външната околна стена
        glBegin(GL_POLYGON);
            glNormal3f( -cos(alpha),  0, -sin(alpha));
            glVertex3f(x+dx1R, y, z+dz1R);
            glVertex3f(x+dx1R, y+h, z+dz1R);
            glNormal3f( -cos(alpha+angle),  0, -sin(alpha+angle));
            glVertex3f(x+dx2R, y+h, z+dz2R);
            glVertex3f(x+dx2R, y, z+dz2R);
        glEnd();
        alpha += angle;
    }
}

// векторно произведение на вектори
OKG_POINT vmul(OKG_POINT a, OKG_POINT b, OKG_POINT p)
{
	a = {x:a.x-p.x, y:a.y-p.y, z:a.z-p.z};
	b = {x:b.x-p.x, y:b.y-p.y, z:b.z-p.z};
	return {x:a.y*b.z-a.z*b.y, y:a.z*b.x-a.x*b.z, z:a.x*b.y-a.y*b.x};
}

// пресмята декартови координати от сферични
OKG_POINT spherical( float alpha, float beta, float r )
{
    return {x:r*cosf(alpha)*cosf(beta), y:r*sinf(alpha)*cosf(beta), z:r*sinf(beta)};
}

// пресмята точка по крайник - линейна комбинация между ротационна синусоида и полусфера
OKG_POINT limbal( float alpha, float beta, float from, float to, float width, float sx1, float sx2, float sy1, float sy2, float sz )
{
    float k1 = sinf(beta)/2+0.5;
    float k = pow(fabs(sin(beta)),10);
    float r = 0.5+0.5*width*cosf((from+(to-from)*k1)*PI/180);
    r = r*(1-k+k*cosf(beta));
    return {x:(sx1*(1-k1)+sx2*k1)*r*cosf(alpha), y:(sy1*(1-k1)+sy2*k1)*r*sinf(alpha), z:sz*k1};
}

// изчислява издутина в диапазона от min до max
double cosser( double t, double min, double max )
{
	if( min<=t && t<=max )
		t = t-min;
	else
	{
		t -= 2*PI;
		if( min<=t && t<=max )
			t = t-min;
		else
			return 0;
	}

	return 0.5+0.5*cos( t/(max-min)*2*PI-PI );
}

// прилага n издутини, всяка със свои параметри: сила, диапазон по a, диапазон по b
double smartRad( double a, double b, int n, float rads[][5] )
{
	double r = 1;
	for(int i=0; i<n; i++)
		r = r+rads[i][0]*cosser( a-PI/2,rads[i][1],rads[i][2] )*cosser( b,rads[i][3],rads[i][4] );
	return r;
}

// дефинира връх и нормален вектор на точка от крайник
void vertex(float alpha,float beta,float from,float span,float width, float sx1, float sx2, float sy1, float sy2, float sz )
{
	OKG_POINT p  = limbal(alpha,beta,from,span,width,sx1,sx2,sy1,sy2,sz);
	OKG_POINT n  = vmul( limbal(alpha+EPS,beta,from,span,width,sx1,sx2,sy1,sy2,sz),
						 limbal(alpha,beta+EPS,from,span,width,sx1,sx2,sy1,sy2,sz),p);
	glNormal3f(n.x,n.y,n.z);
	if( fabs(beta)>PI/2-EPS) glNormal3f(0,0,beta);
	glVertex3f(p.x,p.y,p.z);
}

// дефинира връх и нормален вектор на точка от става
void hertex(float alpha,float beta, int m, float rads[][5], float r=1 )
{
	OKG_POINT p  = spherical(alpha,beta,r*smartRad(alpha,beta,m,rads));
	OKG_POINT n  = vmul( spherical(alpha+EPS,beta,r*smartRad(alpha+EPS,beta,m,rads)),
						 spherical(alpha,beta+EPS,r*smartRad(alpha,beta+EPS,m,rads)),p);
	glNormal3f(n.x,n.y,n.z);
	if( fabs(beta)>PI/2-EPS) glNormal3f(0,0,beta);
	glVertex3f(p.x,p.y,p.z);
}

// рисува крайник
void drawLimb( float from, float span, float width, float sx1, float sy1, float sx2, float sy2, float sz, bool useColor=true )
{
	if( useColor ) glColor3f(0,0.4,0.45);

	float alpha = 0;
    float dalpha = 2.0*PI/nL;
    for( int i=0; i<nL; i++, alpha+=dalpha)
    {
        glBegin( GL_QUAD_STRIP );
		float beta = PI/2;
		for( int j=0; j<nL+1; j++, beta-=dalpha/2)
        {
			vertex(alpha,beta,from,span,width,sx1,sx2,sy1,sy2,sz);
			vertex(alpha+dalpha,beta,from,span,width,sx1,sx2,sy1,sy2,sz);
        }
        glEnd( );
    }
    glTranslatef(0,0,sz); // преместваме се в края на крайника
}

// рисува става
void drawJoint( float r )
{
	glColor3f(1,0.5,0);

    float alpha = 0;
    float dalpha = 2.0*PI/nJ;
    for( int i=0; i<nJ; i++, alpha+=dalpha)
    {
        glBegin( GL_QUAD_STRIP );
		float beta = PI/2;
		for( int j=0; j<nJ/2+1; j++, beta-=dalpha)
		{
			hertex(alpha,beta,0,NULL,r);
			hertex(alpha+dalpha,beta,0,NULL,r);
        }
        glEnd( );
    }
}

// рисува оцветена с ленти глава
void drawHead( int n, float rads[][5] )
{
	double alpha = PI-PI/2;
	double dalpha = 2.0*PI/nH;
	for( int i=0; i<nH; i++, alpha+=dalpha )
	{
		glBegin( GL_QUAD_STRIP );
		double beta = PI/2;
		for( int j=1; j<nH; j++, beta-=dalpha/2 )
		{
			if( (i%2==0 && (j<nH/3 || i<nH/2-10 || i>nH/2+10) && j<nH-15)
				|| (j==nH/2-3 && ((i<nH/2-2 && i>nH/2-8) || (i>nH/2+2 && i<nH/2+8)) )
				|| (j==nH/2-2 && ((i<nH/2-2 && i>nH/2-6) || (i>nH/2+2 && i<nH/2+6)) )
				|| (j==nH/2-1 && ((i<nH/2-2 && i>nH/2-4) || (i>nH/2+2 && i<nH/2+4)) )
				|| (i==nH/2 && j>nH/2-5 && j<nH/2+10))
				glColor3f(0.3,0.5,0.55);
			else
				glColor3f(0,0.4,0.45);
			hertex( alpha,beta,n,rads );
			hertex( alpha+dalpha,beta,n,rads );
		}
		glEnd( );
	}
}

// рисува крак в дадена поза
void drawLeg(float angle[])
{
	glPushMatrix();
		drawJoint(0.5);
		turnForward(angle[0]);	// завъртане на крака
		turnSideway(angle[1]);
		drawLimb (-90,230,0.3,	0.9,1, 0.7,1, 3);
		drawJoint(0.3);
		turnForward(angle[2]);	// завъртане на коляното
		drawLimb (-90,300,0.2,	0.5,0.7, 0.4,0.5, 3);
		drawJoint(0.2);
		glTranslatef(0,0,0.3);
		turnForward(90);	// завъртане на ходилото
		turnForward(angle[3]);
		glTranslatef(0,0,-0.4);
		drawLimb (-90,420,0.1,	0.6,0.5, 0.9,0.2, 1.5);
	glPopMatrix();
}

// рисува ръка в дадена поза
void drawHand(float angle[])
{
	glPushMatrix();
		drawJoint(0.3);
		turnSideway(angle[0]);	// завъртане в рамото
		turnForward(angle[1]);
		turnAround(angle[2]);
		drawLimb (-90,200,0.4,	0.5,0.6, 0.9,0.7, 2.5);
		drawJoint(0.25);
		turnForward(angle[3]);	// завъртане в лакътя
		drawLimb (-90,180,0.2,0.5,0.5, 0.4,0.3, 2);
		drawJoint(0.15);
		turnForward(angle[4]);	// завъртане в китката
		turnSideway(angle[5]);
		turnAround(angle[6]);
		drawLimb (-90,180,0.3,	0.4,0.3, 0.6,0.1, 0.9);
	glPopMatrix();
}

// рисува тяло в дадена поза
void drawBody(float pos[], float ori[], float scale[], float waist[], float head[], float leftHand[], float rightHand[], float leftLeg[], float rightLeg[])
{
	glPushMatrix();
		glTranslatef(pos[0],pos[1],pos[2]); 	// позиция на тялото
		turnAround(ori[0]);			// обща ориентация на тялото
		turnForward(ori[1]);

		glScalef(scale[1],scale[3],scale[1]);	// размер на тялото

		glPushMatrix();
			glTranslatef(0,-0.1,-0.3);	// таз
			drawLimb(-120,30,0.9,	2,1.5, 0.7,0.4, 1.5);
			drawJoint(0.2);
			turnForward(waist[0]);		// завъртане в кръста
			turnAround(waist[1]);
			drawLimb(-130,150,0.4,	1.2,1, 1.9,0.9, 3.0, true);

			glPushMatrix();
				glTranslatef(0,0,-0.5);
				turnForward(180);
				glTranslatef(-1.1,0,0);
				drawHand(leftHand);	// лява ръка
				glTranslatef(2.2,0,0);
				drawHand(rightHand);	// дясна ръка
			glPopMatrix();

			turnAround(head[0]+180);	// завъртане във врата и глава
			turnForward(head[1]);
			turnSideway(head[2]);
			glTranslatef(0,-0.3,0.85);
			glScalef(scale[3]*0.8,scale[3],scale[3]);
			float rads[2][5] = {{-0.25,-3,+3,-1.8,0.0}, {-0.20,PI-3,PI+3,-1.5,1.0}};
			drawHead(2,rads);
		glPopMatrix();

		turnForward(180);
		glTranslatef(0.7,0,0);
		drawLeg(rightLeg);			// десен крак
		glTranslatef(-1.4,0,0);
		drawLeg(leftLeg);			// ляв крак
	glPopMatrix();
}

bool running()
{
    int width, height;
    glfwGetWindowSize( &width, &height );
    glViewport( 0, 0, width, height );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 40.0f, (GLfloat)width/(GLfloat)height, 4, 100 );
    glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
    glfwSwapBuffers();
	glClear( GL_COLOR_BUFFER_BIT+GL_DEPTH_BUFFER_BIT );
    return( !glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam( GLFW_OPENED) );
}

int main()
{
    glfwInit();
	glfwOpenWindowHint(GLFW_FSAA_SAMPLES,GL_TRUE);
    if( !glfwOpenWindow( 1000, 600, 8, 8, 8, 0, 16, 0, GLFW_WINDOW ) ) return 1;
	glfwSetWindowTitle("Marriage Proposal");

    glEnable( GL_DEPTH_TEST );
    glEnable( GL_LIGHTING );
    glEnable( GL_LIGHT0 );
	glEnable( GL_COLOR_MATERIAL );
	glEnable( GL_NORMALIZE );
	glShadeModel(GL_SMOOTH);

    glMaterialf(GL_FRONT,GL_SHININESS,10);
    float F[4] = {0.4,0.3,0,1};
    glMaterialfv(GL_FRONT,GL_SPECULAR,F);

    while( running() )
    {
        float t=glfwGetTime()/20;
        float dist=15+5*cos(t);
        gluLookAt(dist*cos(3*t),dist*sin(3*t),dist*0.3, 0,0,0, 0,0,1 );

        // кутийката за пръстена
        glColor3ub(255, 0, 0);
        drawParallelepiped(-1.5, 2.2, -1.35, 0.2, 0.2, 0.1);
        drawParallelepiped(-1.5, 2.4, -1.2, 0.2, 0.07, 0.3);

        // пръстенът
        glColor3f(1,0.5,0);
        drawTwoCylinders(-1.5, 2.1, -1.2, 0.15,0.09,0.1);

		{	// жената
			float position[]={0,0,0};
			float orientation[]={0,0};
			float scale[]={1,1,1,1};
			float waist[]={0,-17};
			float head[]={10+sinf(130*t),10+sinf(130*t),0};
			float leftHand[]={-5,0,-33,25,0,0,100};
			float rightHand[]={0,21,-25,15,33,0,0};
			float leftLeg[]={0,0,0,0};
            float rightLeg[]={0,0,0,0};
			drawBody (position,orientation,scale,waist,head,leftHand,rightHand,leftLeg,rightLeg);
		}

		{	// мъжът
			float position[]={0,6,-2.7};
			float orientation[]={180,-10};
			float scale[]={1.1,1.1,1.1,1.1};
			float waist[]={-350,0};
			float head[]={0,-20,-5};
			float leftHand[]={-5,25,0,75,0,0,0};
			float rightHand[]={5,20,0,60, 15,0,0};
			float leftLeg[]={100,0,-100,0};
			float rightLeg[]={35,0,-120,0};
			drawBody (position,orientation,scale,waist,head,leftHand,rightHand,leftLeg,rightLeg);
		}
    }

    glfwTerminate();
    return 0;
}
