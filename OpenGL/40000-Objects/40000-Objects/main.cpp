 /*
  * 40000 objects
  * This program is meant to replicate the webgl demo given at Google IO 2011. 
  */
#include <GL/glut.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <iostream>
//#include "fast.h"

# define PI           3.14159265358979323846  /* pi */

using namespace std;

// global variables
float g_eyeSpeed = 0.1;
int g_eyeHeight = 2;
int g_eyeRadius = 19;
bool updatePositions = true;

float clock = 0.0;
float scale = 0.2;
float time = 0;
float fps = 0;
float frames = 0;
float base_time = 0;


float colors[12][4] = {
	{0.000, 0.282, 0.439, 1},
		{0.157, 0.259, 0.565, 1},
		{0.000, 0.400, 0.200, 1},
		{0.059, 0.671, 1.000, 1},
		{0.933, 0.698, 0.067, 1},
		{1.000, 0.827, 0.098, 1},
		{0.000, 0.659, 0.208, 1},
		{0.627, 0.808, 0.404, 1},
		{0.357, 0.078, 0.000, 1},
		{0.525, 0.063, 0.004, 1},
		{1.000, 0.000, 0.000, 1},
		{1.000, 0.537, 0.514, 1}
};

float instances[40000][4]; 

GLfloat xRotated, yRotated, zRotated;
GLdouble radius = .2;

// pre-allocate a bunch of arrays
vector<float> projection;
vector<float> view;
vector<float> world;
vector<float> worldInverse;
vector<float> worldInverseTranspose;
vector<float> viewProjection;
vector<float> worldViewProjection;
vector<float> viewInverse;
vector<float> viewProjectionInverse;
vector<float> eyePosition;
vector<float> target;
vector<float> up;
vector<float> lightWorldPos;
vector<float> v3t0;
vector<float> v3t1;
vector<float> v3t2;
vector<float> v3t3;
vector<float> m4t0;
vector<float> m4t1;
vector<float> m4t2;
vector<float> m4t3;
vector<float> zero4;
vector<float> one4;


// array[int] models = 40000;



void init()
{
	//glClearColor(0.0, 0.0, 0.0, 0.0);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	//glShadeModel(GL_SMOOTH);

	//projection.resize(16);
	view.resize(16);
	//world.resize(16);
	//worldInverse.resize(16);
	//worldInverseTranspose.resize(16);
	//viewProjection.resize(16);
	//worldViewProjection.resize(16);
	//viewInverse.resize(16);
	//viewProjectionInverse.resize(16);
	eyePosition.resize(3);
	target.resize(3);
	up.resize(3); up = {0, 1, 0};
	//lightWorldPos.resize(3);
	//v3t0.resize(3);
	//v3t1.resize(3);
	//v3t2.resize(3);
	//v3t3.resize(3);
	//m4t0.resize(16);
	//m4t1.resize(16);
	//m4t2.resize(16);
	//m4t3.resize(16);
	//zero4.resize(4);
	//one4.resize(4); one4 = { 1, 1, 1, 1 };
	glLoadIdentity();
	glScalef(.1, .1, .1);
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_MODELVIEW);
	glClearColor(1, 1, 1, 0);
	glColorMask(true, true, true, true);
	glDepthMask(true);
	glDepthRange(10, -10);



	for (int i = 0; i < 34; i++)
	{
		for (int j = 0; j < 34; j++)
		{
			for (int k = 0; k < 34; k++)
			{

			}
		}

	}
	float x, y, z;

	x = -.5;
	y = -.5;
	z = -.5;


	for (int i = 0; i < 40000; i++)
	{

		instances[i][0] = x;
		instances[i][1] = y;
		instances[i][2] = z;

	}

	//glViewport(10, 10, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
	gluLookAt(0, 0, -50, 0, 0, 0, 0, 1, 0);
}


float degToRad(float deg) {
	return (deg * (PI / 180));
}


void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// housekeeping
	frames++; 
	// Make the camera rotate around the scene.
	eyePosition[0] = sin(frames * g_eyeSpeed) * g_eyeRadius; // x 
	eyePosition[1] = g_eyeHeight; // y
	eyePosition[2] = cos(frames * g_eyeSpeed) * g_eyeRadius; // z 

	//glutSolidSphere(radius, 20, 20);

	//float x = -3.0;
	//float y = -3.0;
	//float z = -3.0;
	//glTranslatef(x, y, z);




	for (int k = 0; k < 34; k++)
	{
		for (int j = 0; j < 34; j++)
		{
			for (int i = 0; i < 34; i++)
			{
				glColor3f(colors[i % 12][0], colors[i % 12][1], colors[i % 12][2]);
				glutSolidCube(.1);
				glTranslatef(.17, 0, 0);
			}
			glTranslatef(-(34 * .17), .17, 0);
		}
		glTranslatef(0, -(34 * .17), .17);
	}
	//glColor3f(colors[2][0], colors[2][1], colors[2][2]);
	//glBegin(GL_POLYGON);
	//glVertex3f(-1.0, 1.0, 0.0);
	//glVertex3f(1.0, 1.0, 0.0);
	//glVertex3f(1.0, -1.0, 0.0);
	//glVertex3f(-1.0, -1.0, 0.0);
	//glEnd();
	//glBegin(GL_POLYGON);

	//glVertex3f(-1.0, 1.0, -1.0);
	//glVertex3f(1.0, 1.0, -1.0);
	//glVertex3f(1.0, -1.0, -1.0);
	//glVertex3f(-1.0, -1.0, -1.0);
	//glEnd();
	//glBegin(GL_POLYGON);
	//glVertex3f(-1.0, 1.0, 0.0);
	//glVertex3f(-1.0, 1.0, -1.0);
	//glVertex3f(1.0, 1.0, -1.0);
	//glVertex3f(1.0, 1.0, 0.0);
	//glEnd();
	




	// update fps
	time = glutGet(GLUT_ELAPSED_TIME);
	if ((time - base_time) > 1000.0)
	{
		fps = frames * 1000.0 / (time - base_time);
		base_time = time;
		frames = 0;
	}
	printf(to_string(fps).c_str());
	printf("\n");


	// look at mass
	glLoadIdentity();
	gluLookAt(eyePosition[0], eyePosition[1], eyePosition[2], 0, 0, 0, 0, 1, 0);
	//gluLookAt(0, 0, -50, 0, 0, 0, 0, 1, 0);

	// breaks without flush
	glFlush();

	// call redisplay
	glutPostRedisplay();

}


void reshape(int x, int y)
{
	if (y == 0 || x == 0) return;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(39.0, (GLdouble)x / (GLdouble)y, 0.6, 21.0);
	glMatrixMode(GL_MODELVIEW);
	glViewport(10, 10, x, y);  //Use the whole window for rendering
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 27:
		exit(0);
		break;
	}
	updatePositions = !updatePositions;

}



int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitWindowSize(700, 700);
	glutCreateWindow("40000 Objects");
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMainLoop();

	return 0;
}