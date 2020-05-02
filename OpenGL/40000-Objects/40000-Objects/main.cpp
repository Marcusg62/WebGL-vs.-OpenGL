 /*
  * 40000 objects
  * This program is meant to replicate the webgl demo given at Google IO 2011. 
  */
#include <GL/glut.h>
#include <stdlib.h>
#include <vector>
#include <math.h>       /* sin */

# define PI           3.14159265358979323846  /* pi */

using namespace std;

// global variables
float g_eyeSpeed = 0.5;
int g_eyeHeight = 2;
int g_eyeRadius = 9;
bool updatePositions = true;

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

/**
 * Computes the dot product of two vectors; assumes both vectors have
 * three entries.
 * @param {!tdl.fast.Vector} a Operand vector.
 * @param {!tdl.fast.Vector} b Operand vector.
 * @return {number} dot product
 */
float tdl_fast_dot(vector<float>a, vector<float> b) {
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
};

/**
 * Computes the cross product of two vectors; assumes both vectors have
 * three entries.
 * @param {!tdl.fast.Vector} dst vector.
 * @param {!tdl.fast.Vector} a Operand vector.
 * @param {!tdl.fast.Vector} b Operand vector.
 * @return {!tdl.fast.Vector} The vector a cross b.
 */
vector<float> tdl_fast_cross(vector<float> dst, vector<float>a, vector<float>b) {
	dst[0] = a[1] * b[2] - a[2] * b[1];
	dst[1] = a[2] * b[0] - a[0] * b[2];
	dst[2] = a[0] * b[1] - a[1] * b[0];
	return dst;
};

/**
 * Multiplies two 4-by-4 matrices; assumes that the given matrices are 4-by-4.
 * Note: It is faster to call this than tdl.fast.mul.
 * @param {!tdl.fast.Matrix4} a The matrix on the left.
 * @param {!tdl.fast.Matrix4} b The matrix on the right.
 * @return {!tdl.fast.Matrix4} The matrix product of a and b.
 */
vector<float> tdl_fast_matrix4_mul(vector<float>dst, vector<float>a, vector<float>b) {
	return tdl.fast.mulMatrixMatrix4(dst, a, b);
};

/**
 * Subtracts two vectors.
 * @param {!tdl.fast.Vector} dst vector.
 * @param {!tdl.fast.Vector} a Operand vector.
 * @param {!tdl.fast.Vector} b Operand vector.
 */
vector<float> tdl_fast_subVector(vector<float> dst, vector<float> a, vector<float> b) {
	int aLength = a.size;
	for (int i = 0; i < aLength; ++i)
		dst[i] = a[i] - b[i];
	return dst;
};




void init(void)
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glShadeModel(GL_SMOOTH);

	projection.resize(16);
	view.resize(16);
	world.resize(16);
	worldInverse.resize(16);
	worldInverseTranspose.resize(16);
	viewProjection.resize(16);
	worldViewProjection.resize(16);
	viewInverse.resize(16);
	viewProjectionInverse.resize(16);
	eyePosition.resize(3);
	target.resize(3);
	up.resize(3); up = {0, 1, 0};
	lightWorldPos.resize(3);
	v3t0.resize(3);
	v3t1.resize(3);
	v3t2.resize(3);
	v3t3.resize(3);
	m4t0.resize(16);
	m4t1.resize(16);
	m4t2.resize(16);
	m4t3.resize(16);
	zero4.resize(4);
	one4.resize(4); one4 = { 1, 1, 1, 1 };


}

/**
 * Divides a vector by its Euclidean length and returns the quotient.
 * @param {!tdl.fast.Vector} dst vector.
 * @param {!tdl.fast.Vector} a The vector.
 * @return {!tdl.fast.Vector} The normalized vector.
 */
vector<float> tdl_fast_normalize(vector<float> dst, vector<float> a) {
	float n = 0.0;
	float aLength = a.size;
	for (int i = 0; i < aLength; ++i)
		n += a[i] * a[i];
	n = sqrt(n);
	if (n > 0.00001) {
		for (int i = 0; i < aLength; ++i)
			dst[i] = a[i] / n;
	}
	else {
		for (int i = 0; i < aLength; ++i)
			dst[i] = 0;
	}
	return dst;
};

/**
 * Computes a 4-by-4 look-at transformation.  The transformation generated is
 * an orthogonal rotation matrix with translation component.  The translation
 * component sends the eye to the origin.  The rotation component sends the
 * vector pointing from the eye to the target to a vector pointing in the
 * negative z direction, and also sends the up vector into the upper half of
 * the yz plane.
 * @return {!tdl.fast.Matrix4} dst matrix.
 * @param {(!tdl.fast.Vector3} eye The
 *     position of the eye.
 * @param {(!tdl.fast.Vector3} target The
 *     position meant to be viewed.
 * @param {(!tdl.fast.Vector3} up A vector
 *     pointing up.
 * @return {!tdl.fast.Matrix4} The look-at matrix.
 */
vector<float> tdl_fast_matrix4_lookAt(vector<float> dst, vector<float> eye, vector<float> target, vector<float> up) {
	vector<float> t0;
	t0.push_back(0);
	t0.push_back(0);
	t0.push_back(0);
	vector<float> t1;
	t1.push_back(0);
	t1.push_back(0);
	t1.push_back(0);
	vector<float> t2;
	t2.push_back(0);
	t2.push_back(0);
	t2.push_back(0);

	vector<float> vz = tdl_fast_normalize(t0, tdl_fast_subVector(t0, eye, target));
	vector<float> vx = tdl_fast_normalize(t1, tdl_fast_cross(t1, up, vz));
	vector<float> vy = tdl_fast_cross(t2, vz, vx);

	dst[0] = vx[0];
	dst[1] = vy[0];
	dst[2] = vz[0];
	dst[3] = 0;
	dst[4] = vx[1];
	dst[5] = vy[1];
	dst[6] = vz[1];
	dst[7] = 0;
	dst[8] = vx[2];
	dst[9] = vy[2];
	dst[10] = vz[2];
	dst[11] = 0;
	dst[12] = -tdl_fast_dot(vx, eye);
	dst[13] = -tdl_fast_dot(vy, eye);
	dst[14] = -tdl_fast_dot(vz, eye);
	dst[15] = 1;

	return;
};

/**
 * Computes a 4-by-4 perspective transformation matrix given the angular height
 * of the frustum, the aspect ratio, and the near and far clipping planes.  The
 * arguments define a frustum extending in the negative z direction.  The given
 * angle is the vertical angle of the frustum, and the horizontal angle is
 * determined to produce the given aspect ratio.  The arguments near and far are
 * the distances to the near and far clipping planes.  Note that near and far
 * are not z coordinates, but rather they are distances along the negative
 * z-axis.  The matrix generated sends the viewing frustum to the unit box.
 * We assume a unit box extending from -1 to 1 in the x and y dimensions and
 * from 0 to 1 in the z dimension.
 * @param {!tdl.fast.Matrix4} dst matrix.
 * @param {number} angle The camera angle from top to bottom (in radians).
 * @param {number} aspect The aspect ratio width / height.
 * @param {number} zNear The depth (negative z coordinate)
 *     of the near clipping plane.
 * @param {number} zFar The depth (negative z coordinate)
 *     of the far clipping plane.
 * @return {!tdl.fast.Matrix4} The perspective matrix.
 */
vector<float> tdl_fast_matrix4_perspective(float angle,float aspect,float zNear, float zFar) {\
	float f = tan(PI * 0.5 - 0.5 * angle);
	float rangeInv = 1.0 / (zNear - zFar);

	vector<float> dst;
	dst.push_back(f / aspect);
	dst.push_back(0);
	dst.push_back(0);
	dst.push_back(0);

	dst.push_back(0);
	dst.push_back(f);
	dst.push_back(0);
	dst.push_back(0);

	dst.push_back(0);
	dst.push_back(0);
	dst.push_back((zNear + zFar) * rangeInv);
	dst.push_back(-1);

	dst.push_back(0);
	dst.push_back(0);
	dst.push_back(zNear * zFar * rangeInv * 2);
	dst.push_back(0);

	return dst;
};

float degToRad(float deg) {
	return (deg * (PI / 180));
}


void renderBegin() {
	// clear the screen.
	glColorMask(true, true, true, true);
	glDepthMask(true);
	glClearColor(1, 1, 1, 0);
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glColorMask(true, true, true, true);


	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);

	// Compute a projection and view matrices.
	vector <float> m4 = tdl_fast_matrix4_perspective(
		degToRad(60),
		glutGet(GLUT_WINDOW_WIDTH) / glutGet(GLUT_WINDOW_HEIGHT),
		1,
		5000);

	m4 = tdl_fast_matrix4_lookAt(view, eyePosition, target, up);

	m4 = mul(viewProjection, view, projection);
	m4.inverse(viewInverse, view);
	m4.inverse(viewProjectionInverse, viewProjection);

	// Put the light near the camera
	m4.getAxis(v3t0, viewInverse, 0); // x
	m4.getAxis(v3t1, viewInverse, 1); // y
	m4.getAxis(v3t2, viewInverse, 2); // z
	fast.mulScalarVector(v3t0, 10, v3t0);
	fast.mulScalarVector(v3t1, 10, v3t1);
	fast.mulScalarVector(v3t2, 10, v3t2);
	fast.addVector(lightWorldPos, eyePosition, v3t0);
	fast.addVector(lightWorldPos, lightWorldPos, v3t1);
	fast.addVector(lightWorldPos, lightWorldPos, v3t2);

	// compute shared matrices
	m4.translation(world, [0, 0, 0]);
	m4.mul(worldViewProjection, world, viewProjection);
	m4.inverse(worldInverse, world);
	m4.transpose(worldInverseTranspose, worldInverse);
}



void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	renderBegin();
	renderScene();
	renderEnd();
	glFlush();

}


void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (w <= h)
		gluOrtho2D(0.0, 30.0, 0.0, 30.0 * (GLfloat)h / (GLfloat)w);
	else
		gluOrtho2D(0.0, 30.0 * (GLfloat)w / (GLfloat)h, 0.0, 30.0);
	glMatrixMode(GL_MODELVIEW);
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
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(700, 700);
	glutInitWindowPosition(100, 100);
	glutCreateWindow(argv[0]);
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMainLoop();
	return 0;
}