

#include<gl/glut.h>
#include <iostream>
#include "Tick.h"

#include "MassSpring2D.h"


// FPS
unsigned int tick_		= 0;
unsigned int lastTick_	= 0;
float elapsedSeconds_	= 0.0f;
int fps_				= 0;
int	fpsFrames_			= 0;
float fpsDuration_		= 0.0f;

bool paused = true;


MassSpring2D massSpring;
int numParticlesHor = 10;
int numParticlesVer = 10;
float deltaX = 10.0f;
float deltaY = 10.0f;
float f = 50.0f;

//int width = numParticlesHor*deltaX;
//int height = numParticlesVer*deltaY;

int width = (numParticlesHor-1)*deltaX;
int height = (numParticlesVer-1)*deltaY;


//-----------------------------------------------------------------------------                                                   
void initMassSpring()
{
	massSpring.addParticle(10.0f, 0.0f, true);
	massSpring.addParticle(20.0f, 0.0f, true); 
	massSpring.addParticle(0.0f, 10.0f, true);
	massSpring.addParticle(10.0f, 10.0f, false); // false
	massSpring.addParticle(20.0f, 10.0f, false); // false
	massSpring.addParticle(30.0f, 10.0f, true);
	massSpring.addParticle(0.0f, 20.0f, true);
	massSpring.addParticle(10.0f, 20.0f, false); // false
	massSpring.addParticle(20.0f, 20.0f, false); // false
	massSpring.addParticle(30.0f, 20.0f, true);
	massSpring.addParticle(10.0f, 30.0f, true);
	massSpring.addParticle(20.0f, 30.0f, true);

	
	/*
	// inside
	massSpring.addForce(3, f, f);
	massSpring.addForce(4, -f, f);
	massSpring.addForce(8, -f, -f);
	massSpring.addForce(7, f, -f);
	//*/

	/*
	// unica
	massSpring.addForce(3, 1.0f, 1.0f);
	massSpring.addForce(4, 1.0f, 1.0f);
	massSpring.addForce(8, 1.0f, 1.0f);
	massSpring.addForce(7, 1.0f, 1.0f);
	*/

	//*
	//outside
	massSpring.addForce(3, -f, -f);
	massSpring.addForce(4, f, -f);
	massSpring.addForce(8, f, f);
	massSpring.addForce(7, -f, f);
	//*/

	massSpring.addSpring(0, 3, 10.0f, 10.0f);
	massSpring.addSpring(1, 4, 10.0f, 10.0f);
	massSpring.addSpring(2, 3, 10.0f, 10.0f);
	massSpring.addSpring(4, 5, 10.0f, 10.0f);
	massSpring.addSpring(6, 7, 10.0f, 10.0f);
	massSpring.addSpring(8, 9, 10.0f, 10.0f);
	massSpring.addSpring(7, 10, 10.0f, 10.0f);
	massSpring.addSpring(8, 11, 10.0f, 10.0f);
	massSpring.addSpring(3, 4, 25.0f, 10.0f); // esta
	massSpring.addSpring(3, 7, 25.0f, 10.0f); // esta
	massSpring.addSpring(4, 8, 25.0f, 10.0f); // esta
	massSpring.addSpring(7, 8, 25.0f, 10.0f); // esta
	//massSpring.addSpring(3, 8, 25.0f, 10.0f); // esta
	//massSpring.addSpring(4, 7, 25.0f, 10.0f); // esta

	massSpring.setDamping(0.5f);
}


void createMassSpringGrid()
{
	bool fixed = false;

	// creating mass
	for (int i = 0; i < numParticlesHor; ++i)
	{
		for (int j = 0; j < numParticlesVer; ++j)
		{
			// fix the mass at borders
			if ((i == 0) 
				|| (j == 0) 
				|| (i == numParticlesHor-1) 
				|| (j == numParticlesVer-1))
			{
				fixed = true;
			}
			else
			{
				fixed = false;
			}

			massSpring.addParticle(i*deltaX, j*deltaY, fixed);
		}
	}

	const float restLength = 10.0f;
	const float stiffness = 25.0f;

	
	// creating springs
	int step = 0;
	for (int i = 0; i < numParticlesHor; ++i)
	{
		for (int j = 0; j < numParticlesVer-1; ++j)
		{
			massSpring.addSpring(j+step, (j+1)+step, stiffness, restLength);
		}
		step += numParticlesVer;
	}
	step = 0;
	
	for (int j = 0; j < numParticlesVer; ++j)
	{
		for (int i = 0; i < numParticlesHor; ++i)
		{
			massSpring.addSpring(i+step, (i+numParticlesHor)+step, stiffness, restLength);
		}
		step = j*numParticlesHor;
	}

	massSpring.addForce(22, -f,-f);
	massSpring.addForce(23, -f,f);
	massSpring.addForce(33, f,f);
	massSpring.addForce(32, f,-f);

	/*
	// cima
	massSpring.addForce(34, 0,f);
	massSpring.addForce(44, f,f);
	massSpring.addForce(24, 0,f);
	massSpring.addForce(14, -f,f);

	// direita
	massSpring.addForce(43, f,0);
	massSpring.addForce(42, f,0);
	*/
}


//-----------------------------------------------------------------------------
void changeForces()
{
	/*
	// inside
	massSpring.addForce(3, f, f);
	massSpring.addForce(4, -f, f);
	massSpring.addForce(8, -f, -f);
	massSpring.addForce(7, f, -f);
	//*/

	/*
	// unica
	massSpring.addForce(3, 1.0f, 1.0f);
	massSpring.addForce(4, 1.0f, 1.0f);
	massSpring.addForce(8, 1.0f, 1.0f);
	massSpring.addForce(7, 1.0f, 1.0f);
	*/

	//*
	//outside
	massSpring.addForce(22, -f, -f);
	massSpring.addForce(32, f, -f);
	massSpring.addForce(33, f, f);
	massSpring.addForce(23, -f, f);
	//*/
}


//-----------------------------------------------------------------------------
void drawAxis()
{
	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(width, 0.0f, 0.0f);
	glEnd();

	glColor3f(0.0f, 1.0f, 0.0f);
	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, height, 0.0f);
	glEnd();
}



//-----------------------------------------------------------------------------
void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		

	/*
	glPushMatrix();
	glTranslatef(0.1f, 0.1f, 0.0f);
	drawAxis();
	glPopMatrix();
	*/

	massSpring.drawSprings();
	massSpring.drawForces();
	massSpring.drawMass();

	glutSwapBuffers();
}


void idle()
{
	// computing FPS
	tick_ = Tick::getTick();
	elapsedSeconds_ = (tick_ - lastTick_) / 1000.0f;
	fpsDuration_ += elapsedSeconds_;
	++fpsFrames_;

	if (fpsDuration_ > 1.0f)
	{
		fps_			= static_cast<int>(fpsFrames_ / fpsDuration_);
		fpsDuration_	= 0.0f;
		fpsFrames_		= 0;
	}
	lastTick_ = tick_;

	if (!paused)
	{
		massSpring.doStep(0.005f);
	}

	glutPostRedisplay();
}


//-----------------------------------------------------------------------------
void myReshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	gluOrtho2D(0.0f ,width, 0.0f, height);
    glMatrixMode(GL_MODELVIEW);
}



//----------------------------------------------------------------------------
void onKey(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 27:
		case 'q':
		case 'Q':
			exit(0);
		break;

		case 'p':
		case 'P':
			paused = !paused;
			break;

		case '=':
		case '+':
			f += 50.0f;
			changeForces();
			std::cout << f << std::endl;
			break;

		case '-':
		case '_':
			f -= 50.0f;
			changeForces();
			std::cout << f << std::endl;
			break;

		default:
			break;
	}

	glutPostRedisplay();
}


//-----------------------------------------------------------------------------
void main(int argc, char **argv)
{
    glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
    glutCreateWindow("Mass-spring");
   
	glClearColor(0.0, 0.0, 0.0, 1.0);

	glutReshapeFunc(myReshape);
    glutDisplayFunc(display); 
	glutIdleFunc(idle);
	glutKeyboardFunc(onKey);


	//initMassSpring();
	createMassSpringGrid();

    glutMainLoop();
}