#include "fast.h"
#include <vector>
#include <math.h>
# define PI           3.14159265358979323846  /* pi */

using namespace std;

/**
 * Adds two vectors; assumes a and b have the same dimension.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 */
vector<float> tdl_fast_addVector(vector<float> dst, vector<float> a, vector<float> b) {
	int aLength = a.size();
	for (int i = 0; i < aLength; ++i)
		dst[i] = a[i] + b[i];
	return dst;
};

/**
 * Subtracts two vectors.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 */
vector<float> tdl_fast_subVector(vector<float> dst, vector<float> a, vector<float> b) {
	int aLength = a.size();
	for (int i = 0; i < aLength; ++i)
		dst[i] = a[i] - b[i];
	return dst;
};

/**
 * Performs linear interpolation on two vectors.
 * Given vectors a and b and interpolation coefficient t, returns
 * (1 - t) * a + t * b.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @param {number} t Interpolation coefficient.
 */
vector<float> tdl_fast_lerpVector(vector<float> dst, vector<float> a, vector<float> b, float t) {
	int aLength = a.size();
	for (int i = 0; i < aLength; ++i)
		dst[i] = (1 - t) * a[i] + t * b[i];
	return dst;
};

/**
 * Divides a vector by a scalar.
 * @param {!tdl_fast_Vector} dst The vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @param {number} k The scalar.
 * @return {!tdl_fast_Vector} dst.
 */
vector<float> tdl_fast_divVectorScalar(vector<float> dst, vector<float> v, float k) {
	int vLength = v.size();
	for (int i = 0; i < vLength; ++i)
		dst[i] = v[i] / k;
	return dst;
};

/**
 * Computes the cross product of two vectors; assumes both vectors have
 * three entries.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @return {!tdl_fast_Vector} The vector a cross b.
 */
vector<float> tdl_fast_cross(vector<float> dst, vector<float> a, vector<float> b) {
	dst[0] = a[1] * b[2] - a[2] * b[1];
	dst[1] = a[2] * b[0] - a[0] * b[2];
	dst[2] = a[0] * b[1] - a[1] * b[0];
	return dst;
};

/**
 * Computes the dot product of two vectors; assumes both vectors have
 * three entries.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @return {number} dot product
 */
float tdl_fast_dot(vector<float> a, vector<float> b) {
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
};

/**
 * Divides a vector by its Euclidean length and returns the quotient.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a The vector.
 * @return {!tdl_fast_Vector} The normalized vector.
 */
vector<float> tdl_fast_normalize(vector<float> dst, vector<float>  a) {
	float n = 0.0;
	int aLength = a.size();
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
 * Negates a vector.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} -v.
 */
vector<float> tdl_fast_negativeVector(vector<float> dst, vector<float> v) {
	int vLength = v.size();
	for (int i = 0; i < vLength; ++i) {
		dst[i] = -v[i];
	}
	return dst;
};

/**
 * Negates a matrix.
 * @param {!tdl_fast_Matrix} dst matrix.
 * @param {!tdl_fast_Matrix} v The matrix.
 * @return {!tdl_fast_Matrix} -v.
 */
vector<float> tdl_fast_negativeMatrix(vector<float> dst, vector<float>  v) {
	int vLength = v.size();
	for (int i = 0; i < vLength; ++i) {
		dst[i] = -v[i];
	}
	return dst;
};

/**
 * Copies a vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} A copy of v.
 */
vector<float> tdl_fast_copyVector(vector<float> dst, vector<float> v) {
	dst.assign(v.begin(), v.end()); 	
	return dst;
};

/**
 * Copies a matrix.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Matrix} A copy of m.
 */
vector<float> tdl_fast_copyMatrix(vector<float> dst, vector<float> m) {
	dst.assign(m.begin(), m.end());
	return dst;
};

/**
 * Multiplies a scalar by a vector.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {number} k The scalar.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} The product of k and v.
 */
vector<float> tdl_fast_mulScalarVector(vector<float> dst, float k, vector<float> v) {
	int vLength = v.size();
	for (int i = 0; i < vLength; ++i) {
		dst[i] = k * v[i];
	}
	return dst;
};

/**
 * Multiplies a vector by a scalar.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @param {number} k The scalar.
 * @return {!tdl_fast_Vector} The product of k and v.
 */
vector<float> tdl_fast_mulVectorScalar(vector<float> dst, vector<float> v, float k) {
	return tdl_fast_mulScalarVector(dst, k, v);
};

/**
 * Multiplies a scalar by a matrix.
 * @param {!tdl_fast_Matrix} dst matrix.
 * @param {number} k The scalar.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Matrix} The product of m and k.
 */
vector<float> tdl_fast_mulScalarMatrix(vector<float> dst, float k, vector<float> m) {
	int mLength = m.size();
	for (int i = 0; i < mLength; ++i) {
		dst[i] = k * m[i];
	}
	return dst;
};

/**
 * Multiplies a matrix by a scalar.
 * @param {!tdl_fast_Matrix} dst matrix.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} k The scalar.
 * @return {!tdl_fast_Matrix} The product of m and k.
 */
vector<float> tdl_fast_mulMatrixScalar(vector<float> dst, vector<float>m, float k) {
	return tdl_fast_mulScalarMatrix(dst, k, m);
};

/**
 * Multiplies a vector by another vector (component-wise); assumes a and
 * b have the same length.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @return {!tdl_fast_Vector} The vector of products of entries of a and
 *     b.
 */
vector<float>  tdl_fast_mulVectorVector(vector<float> dst, vector<float> a, vector<float> b) {
	int aLength = a.size();
	for (int i = 0; i < aLength; ++i)
		dst[i] = a[i] * b[i];
	return dst;
};

/**
 * Divides a vector by another vector (component-wise); assumes a and
 * b have the same length.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @return {!tdl_fast_Vector} The vector of quotients of entries of a and
 *     b.
 */
vector<float> tdl_fast_divVectorVector(vector<float> dst, vector<float> a, vector<float> b) {
	int aLength = a.size();
	for (int i = 0; i < aLength; ++i)
		dst[i] = a[i] / b[i];
	return dst;
};

/**
 * Multiplies a vector by a matrix; treats the vector as a row vector; assumes
 * matrix entries are accessed in [row][column] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Vector} The product of v and m as a row vector.
 */
vector<float> tdl_fast_rowMajor_mulVectorMatrix4(vector<float> dst, vector<float>v, vector<float>m) {
	for (int i = 0; i < 4; ++i) {
		dst[i] = 0.0;
		for (int j = 0; j < 4; ++j)
			dst[i] += v[j] * m[j * 4 + i];
	}
	return dst;
};

/**
 * Multiplies a vector by a matrix; treats the vector as a row vector; assumes
 * matrix entries are accessed in [column][row] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Vector} The product of v and m as a row vector.
 */
vector<float> tdl_fast_columnMajor_mulVectorMatrix4(vector<float> dst, vector<float> v, vector<float> m) {
	int mLength = m.size();
	int vLength = v.size();
	for (int i = 0; i < 4; ++i) {
		dst[i] = 0.0;
		int col = i * 4;
		for (int j = 0; j < 4; ++j)
			dst[i] += v[j] * m[col + j];
	}
	return dst;
};



/**
 * Multiplies a matrix by a vector; treats the vector as a column vector.
 * assumes matrix entries are accessed in [row][column] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} The product of m and v as a column vector.
 */
vector<float> tdl_fast_rowMajor_mulMatrix4Vector(vector<float> dst, vector<float>m, vector<float>v) {
	for (int i = 0; i < 4; ++i) {
		dst[i] = 0.0;
		int row = i * 4;
		for (int j = 0; j < 4; ++j)
			dst[i] += m[row + j] * v[j];
	}
	return dst;
};

/**
 * Multiplies a matrix by a vector; treats the vector as a column vector;
 * assumes matrix entries are accessed in [column][row] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} The product of m and v as a column vector.
 */
vector<float> tdl_fast_columnMajor_mulMatrix4Vector(vector<float> dst, vector<float> m, vector<float> v) {
	for (int i = 0; i < 4; ++i) {
		dst[i] = 0.0;
		for (int j = 0; j < 4; ++j)
			dst[i] += v[j] * m[j * 4 + i];
	}
	return dst;
};

/**
 * Multiplies two 3-by-3 matrices; assumes that the given matrices are 3-by-3;
 * assumes matrix entries are accessed in [row][column] fashion.
 * @param {!tdl_fast_Matrix3} dst matrix.
 * @param {!tdl_fast_Matrix3} a The matrix on the left.
 * @param {!tdl_fast_Matrix3} b The matrix on the right.
 * @return {!tdl_fast_Matrix3} The matrix product of a and b.
 */
vector<float> tdl_fast_rowMajor_mulMatrixMatrix3(vector<float> dst, vector<float> a, vector<float> b) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a10 = a[3 + 0];
	float a11 = a[3 + 1];
	float a12 = a[3 + 2];
	float a20 = a[6 + 0];
	float a21 = a[6 + 1];
	float a22 = a[6 + 2];
	float b00 = b[0];
	float b01 = b[1];
	float b02 = b[2];
	float b10 = b[3 + 0];
	float b11 = b[3 + 1];
	float b12 = b[3 + 2];
	float b20 = b[6 + 0];
	float b21 = b[6 + 1];
	float b22 = b[6 + 2];
	dst[0] = a00 * b00 + a01 * b10 + a02 * b20;
	dst[1] = a00 * b01 + a01 * b11 + a02 * b21;
	dst[2] = a00 * b02 + a01 * b12 + a02 * b22;
	dst[3] = a10 * b00 + a11 * b10 + a12 * b20;
	dst[4] = a10 * b01 + a11 * b11 + a12 * b21;
	dst[5] = a10 * b02 + a11 * b12 + a12 * b22;
	dst[6] = a20 * b00 + a21 * b10 + a22 * b20;
	dst[7] = a20 * b01 + a21 * b11 + a22 * b21;
	dst[8] = a20 * b02 + a21 * b12 + a22 * b22;
	return dst;
};

/**
 * Multiplies two 3-by-3 matrices; assumes that the given matrices are 3-by-3;
 * assumes matrix entries are accessed in [column][row] fashion.
 * @param {!tdl_fast_Matrix3} dst matrix.
 * @param {!tdl_fast_Matrix3} a The matrix on the left.
 * @param {!tdl_fast_Matrix3} b The matrix on the right.
 * @return {!tdl_fast_Matrix3} The matrix product of a and b.
 */
vector<float> tdl_fast_columnMajor_mulMatrixMatrix3(vector<float> dst, vector<float> a, vector<float> b) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a10 = a[3 + 0];
	float a11 = a[3 + 1];
	float a12 = a[3 + 2];
	float a20 = a[6 + 0];
	float a21 = a[6 + 1];
	float a22 = a[6 + 2];
	float b00 = b[0];
	float b01 = b[1];
	float b02 = b[2];
	float b10 = b[3 + 0];
	float b11 = b[3 + 1];
	float b12 = b[3 + 2];
	float b20 = b[6 + 0];
	float b21 = b[6 + 1];
	float b22 = b[6 + 2];
	dst[0] = a00 * b00 + a10 * b01 + a20 * b02;
	dst[1] = a01 * b00 + a11 * b01 + a21 * b02;
	dst[2] = a02 * b00 + a12 * b01 + a22 * b02;
	dst[3] = a00 * b10 + a10 * b11 + a20 * b12;
	dst[4] = a01 * b10 + a11 * b11 + a21 * b12;
	dst[5] = a02 * b10 + a12 * b11 + a22 * b12;
	dst[6] = a00 * b20 + a10 * b21 + a20 * b22;
	dst[7] = a01 * b20 + a11 * b21 + a21 * b22;
	dst[8] = a02 * b20 + a12 * b21 + a22 * b22;
	return dst;
};



/**
 * Multiplies two 4-by-4 matrices; assumes that the given matrices are 4-by-4;
 * assumes matrix entries are accessed in [row][column] fashion.
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @param {!tdl_fast_Matrix4} a The matrix on the left.
 * @param {!tdl_fast_Matrix4} b The matrix on the right.
 * @return {!tdl_fast_Matrix4} The matrix product of a and b.
 */
vector<float> tdl_fast_rowMajor_mulMatrixMatrix4(vector<float> dst, vector<float> a, vector<float> b) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a03 = a[3];
	float a10 = a[4 + 0];
	float a11 = a[4 + 1];
	float a12 = a[4 + 2];
	float a13 = a[4 + 3];
	float a20 = a[8 + 0];
	float a21 = a[8 + 1];
	float a22 = a[8 + 2];
	float a23 = a[8 + 3];
	float a30 = a[12 + 0];
	float a31 = a[12 + 1];
	float a32 = a[12 + 2];
	float a33 = a[12 + 3];
	float b00 = b[0];
	float b01 = b[1];
	float b02 = b[2];
	float b03 = b[3];
	float b10 = b[4 + 0];
	float b11 = b[4 + 1];
	float b12 = b[4 + 2];
	float b13 = b[4 + 3];
	float b20 = b[8 + 0];
	float b21 = b[8 + 1];
	float b22 = b[8 + 2];
	float b23 = b[8 + 3];
	float b30 = b[12 + 0];
	float b31 = b[12 + 1];
	float b32 = b[12 + 2];
	float b33 = b[12 + 3];
	dst[0] = a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30;
	dst[1] = a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31;
	dst[2] = a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32;
	dst[3] = a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33;
	dst[4] = a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30;
	dst[5] = a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31;
	dst[6] = a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32;
	dst[7] = a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33;
	dst[8] = a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30;
	dst[9] = a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31;
	dst[10] = a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32;
	dst[11] = a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33;
	dst[12] = a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30;
	dst[13] = a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31;
	dst[14] = a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32;
	dst[15] = a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33;
	return dst;
};

/**
 * Multiplies two 4-by-4 matrices; assumes that the given matrices are 4-by-4;
 * assumes matrix entries are accessed in [column][row] fashion.
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @param {!tdl_fast_Matrix4} a The matrix on the left.
 * @param {!tdl_fast_Matrix4} b The matrix on the right.
 * @return {!tdl_fast_Matrix4} The matrix product of a and b.
 */
vector<float> tdl_fast_columnMajor_mulMatrixMatrix4(vector<float> dst, vector<float> a, vector<float> b) {
	float a00 = a[0];
	float a01 = a[1];
	float a02 = a[2];
	float a03 = a[3];
	float a10 = a[4 + 0];
	float a11 = a[4 + 1];
	float a12 = a[4 + 2];
	float a13 = a[4 + 3];
	float a20 = a[8 + 0];
	float a21 = a[8 + 1];
	float a22 = a[8 + 2];
	float a23 = a[8 + 3];
	float a30 = a[12 + 0];
	float a31 = a[12 + 1];
	float a32 = a[12 + 2];
	float a33 = a[12 + 3];
	float b00 = b[0];
	float b01 = b[1];
	float b02 = b[2];
	float b03 = b[3];
	float b10 = b[4 + 0];
	float b11 = b[4 + 1];
	float b12 = b[4 + 2];
	float b13 = b[4 + 3];
	float b20 = b[8 + 0];
	float b21 = b[8 + 1];
	float b22 = b[8 + 2];
	float b23 = b[8 + 3];
	float b30 = b[12 + 0];
	float b31 = b[12 + 1];
	float b32 = b[12 + 2];
	float b33 = b[12 + 3];
	dst[0] = a00 * b00 + a10 * b01 + a20 * b02 + a30 * b03;
	dst[1] = a01 * b00 + a11 * b01 + a21 * b02 + a31 * b03;
	dst[2] = a02 * b00 + a12 * b01 + a22 * b02 + a32 * b03;
	dst[3] = a03 * b00 + a13 * b01 + a23 * b02 + a33 * b03;
	dst[4] = a00 * b10 + a10 * b11 + a20 * b12 + a30 * b13;
	dst[5] = a01 * b10 + a11 * b11 + a21 * b12 + a31 * b13;
	dst[6] = a02 * b10 + a12 * b11 + a22 * b12 + a32 * b13;
	dst[7] = a03 * b10 + a13 * b11 + a23 * b12 + a33 * b13;
	dst[8] = a00 * b20 + a10 * b21 + a20 * b22 + a30 * b23;
	dst[9] = a01 * b20 + a11 * b21 + a21 * b22 + a31 * b23;
	dst[10] = a02 * b20 + a12 * b21 + a22 * b22 + a32 * b23;
	dst[11] = a03 * b20 + a13 * b21 + a23 * b22 + a33 * b23;
	dst[12] = a00 * b30 + a10 * b31 + a20 * b32 + a30 * b33;
	dst[13] = a01 * b30 + a11 * b31 + a21 * b32 + a31 * b33;
	dst[14] = a02 * b30 + a12 * b31 + a22 * b32 + a32 * b33;
	dst[15] = a03 * b30 + a13 * b31 + a23 * b32 + a33 * b33;
	return dst;
};


/**
 * Gets the jth column of the given matrix m; assumes matrix entries are
 * accessed in [row][column] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} j The index of the desired column.
 * @return {!tdl_fast_Vector} The jth column of m as a vector.
 */
vector<float> tdl_fast_rowMajor_column4(vector<float> dst, vector<float>  m, int j) {
	for (int i = 0; i < 4; ++i) {
		dst[i] = m[i * 4 + j];
	}
	return dst;
};

/**
 * Gets the jth column of the given matrix m; assumes matrix entries are
 * accessed in [column][row] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} j The index of the desired column.
 * @return {!tdl_fast_Vector} The jth column of m as a vector.
 */
vector<float> tdl_fast_columnMajor_column4(vector<float> dst, vector<float>  m, int j) {
	int off = j * 4;
	dst[0] = m[off + 0];
	dst[1] = m[off + 1];
	dst[2] = m[off + 2];
	dst[3] = m[off + 3];
	return dst;
};


/**
 * Gets the ith row of the given matrix m; assumes matrix entries are
 * accessed in [row][column] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} i The index of the desired row.
 * @return {!tdl_fast_Vector} The ith row of m.
 */
vector<float> tdl_fast_rowMajor_row4(vector<float> dst, vector<float>  m, int i) {
	int off = i * 4;
	dst[0] = m[off + 0];
	dst[1] = m[off + 1];
	dst[2] = m[off + 2];
	dst[3] = m[off + 3];
	return dst;
};

/**
 * Gets the ith row of the given matrix m; assumes matrix entries are
 * accessed in [column][row] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} i The index of the desired row.
 * @return {!tdl_fast_Vector} The ith row of m.
 */
vector<float> tdl_fast_columnMajor_row4(vector<float> dst, vector<float>  m, int i) {
	for (int j = 0; j < 4; ++j) {
		dst[j] = m[j * 4 + i];
	}
	return dst;
};


/**
 * Creates an n-by-n identity matrix.
 *
 * @param {!tdl_fast_Matrix} dst matrix.
 * @return {!tdl_fast_Matrix} An n-by-n identity matrix.
 */
vector<float> tdl_fast_identity4(vector<float> dst) {
	dst[0] = 1;
	dst[1] = 0;
	dst[2] = 0;
	dst[3] = 0;
	dst[4] = 0;
	dst[5] = 1;
	dst[6] = 0;
	dst[7] = 0;
	dst[8] = 0;
	dst[9] = 0;
	dst[10] = 1;
	dst[11] = 0;
	dst[12] = 0;
	dst[13] = 0;
	dst[14] = 0;
	dst[15] = 1;
	return dst;
};

/**
 * Takes the transpose of a matrix.
 * @param {!tdl_fast_Matrix} dst matrix.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Matrix} The transpose of m.
 */
vector<float> tdl_fast_transpose4(vector<float> dst, vector<float>  m) {
	if (dst == m) {
		float t;

		t = m[1];
		m[1] = m[4];
		m[4] = t;

		t = m[2];
		m[2] = m[8];
		m[8] = t;

		t = m[3];
		m[3] = m[12];
		m[12] = t;

		t = m[6];
		m[6] = m[9];
		m[9] = t;

		t = m[7];
		m[7] = m[13];
		m[13] = t;

		t = m[11];
		m[11] = m[14];
		m[14] = t;
		return dst;
	}

	float m00 = m[0 * 4 + 0];
	float m01 = m[0 * 4 + 1];
	float m02 = m[0 * 4 + 2];
	float m03 = m[0 * 4 + 3];
	float m10 = m[1 * 4 + 0];
	float m11 = m[1 * 4 + 1];
	float m12 = m[1 * 4 + 2];
	float m13 = m[1 * 4 + 3];
	float m20 = m[2 * 4 + 0];
	float m21 = m[2 * 4 + 1];
	float m22 = m[2 * 4 + 2];
	float m23 = m[2 * 4 + 3];
	float m30 = m[3 * 4 + 0];
	float m31 = m[3 * 4 + 1];
	float m32 = m[3 * 4 + 2];
	float m33 = m[3 * 4 + 3];

	dst[0] = m00;
	dst[1] = m10;
	dst[2] = m20;
	dst[3] = m30;
	dst[4] = m01;
	dst[5] = m11;
	dst[6] = m21;
	dst[7] = m31;
	dst[8] = m02;
	dst[9] = m12;
	dst[10] = m22;
	dst[11] = m32;
	dst[12] = m03;
	dst[13] = m13;
	dst[14] = m23;
	dst[15] = m33;
	return dst;
};

/**
 * Computes the inverse of a 4-by-4 matrix.
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @return {!tdl_fast_Matrix4} The inverse of m.
 */
vector<float> tdl_fast_inverse4(vector<float> dst, vector<float>  m) {
	float m00 = m[0 * 4 + 0];
	float m01 = m[0 * 4 + 1];
	float m02 = m[0 * 4 + 2];
	float m03 = m[0 * 4 + 3];
	float m10 = m[1 * 4 + 0];
	float m11 = m[1 * 4 + 1];
	float m12 = m[1 * 4 + 2];
	float m13 = m[1 * 4 + 3];
	float m20 = m[2 * 4 + 0];
	float m21 = m[2 * 4 + 1];
	float m22 = m[2 * 4 + 2];
	float m23 = m[2 * 4 + 3];
	float m30 = m[3 * 4 + 0];
	float m31 = m[3 * 4 + 1];
	float m32 = m[3 * 4 + 2];
	float m33 = m[3 * 4 + 3];
	float tmp_0 = m22 * m33;
	float tmp_1 = m32 * m23;
	float tmp_2 = m12 * m33;
	float tmp_3 = m32 * m13;
	float tmp_4 = m12 * m23;
	float tmp_5 = m22 * m13;
	float tmp_6 = m02 * m33;
	float tmp_7 = m32 * m03;
	float tmp_8 = m02 * m23;
	float tmp_9 = m22 * m03;
	float tmp_10 = m02 * m13;
	float tmp_11 = m12 * m03;
	float tmp_12 = m20 * m31;
	float tmp_13 = m30 * m21;
	float tmp_14 = m10 * m31;
	float tmp_15 = m30 * m11;
	float tmp_16 = m10 * m21;
	float tmp_17 = m20 * m11;
	float tmp_18 = m00 * m31;
	float tmp_19 = m30 * m01;
	float tmp_20 = m00 * m21;
	float tmp_21 = m20 * m01;
	float tmp_22 = m00 * m11;
	float tmp_23 = m10 * m01;

	float t0 = (tmp_0 * m11 + tmp_3 * m21 + tmp_4 * m31) -
		(tmp_1 * m11 + tmp_2 * m21 + tmp_5 * m31);
	float t1 = (tmp_1 * m01 + tmp_6 * m21 + tmp_9 * m31) -
		(tmp_0 * m01 + tmp_7 * m21 + tmp_8 * m31);
	float t2 = (tmp_2 * m01 + tmp_7 * m11 + tmp_10 * m31) -
		(tmp_3 * m01 + tmp_6 * m11 + tmp_11 * m31);
	float t3 = (tmp_5 * m01 + tmp_8 * m11 + tmp_11 * m21) -
		(tmp_4 * m01 + tmp_9 * m11 + tmp_10 * m21);

	float d = 1.0 / (m00 * t0 + m10 * t1 + m20 * t2 + m30 * t3);

	dst[0] = d * t0;
	dst[1] = d * t1;
	dst[2] = d * t2;
	dst[3] = d * t3;
	dst[4] = d * ((tmp_1 * m10 + tmp_2 * m20 + tmp_5 * m30) -
		(tmp_0 * m10 + tmp_3 * m20 + tmp_4 * m30));
	dst[5] = d * ((tmp_0 * m00 + tmp_7 * m20 + tmp_8 * m30) -
		(tmp_1 * m00 + tmp_6 * m20 + tmp_9 * m30));
	dst[6] = d * ((tmp_3 * m00 + tmp_6 * m10 + tmp_11 * m30) -
		(tmp_2 * m00 + tmp_7 * m10 + tmp_10 * m30));
	dst[7] = d * ((tmp_4 * m00 + tmp_9 * m10 + tmp_10 * m20) -
		(tmp_5 * m00 + tmp_8 * m10 + tmp_11 * m20));
	dst[8] = d * ((tmp_12 * m13 + tmp_15 * m23 + tmp_16 * m33) -
		(tmp_13 * m13 + tmp_14 * m23 + tmp_17 * m33));
	dst[9] = d * ((tmp_13 * m03 + tmp_18 * m23 + tmp_21 * m33) -
		(tmp_12 * m03 + tmp_19 * m23 + tmp_20 * m33));
	dst[10] = d * ((tmp_14 * m03 + tmp_19 * m13 + tmp_22 * m33) -
		(tmp_15 * m03 + tmp_18 * m13 + tmp_23 * m33));
	dst[11] = d * ((tmp_17 * m03 + tmp_20 * m13 + tmp_23 * m23) -
		(tmp_16 * m03 + tmp_21 * m13 + tmp_22 * m23));
	dst[12] = d * ((tmp_14 * m22 + tmp_17 * m32 + tmp_13 * m12) -
		(tmp_16 * m32 + tmp_12 * m12 + tmp_15 * m22));
	dst[13] = d * ((tmp_20 * m32 + tmp_12 * m02 + tmp_19 * m22) -
		(tmp_18 * m22 + tmp_21 * m32 + tmp_13 * m02));
	dst[14] = d * ((tmp_18 * m12 + tmp_23 * m32 + tmp_15 * m02) -
		(tmp_22 * m32 + tmp_14 * m02 + tmp_19 * m12));
	dst[15] = d * ((tmp_22 * m22 + tmp_16 * m02 + tmp_21 * m12) -
		(tmp_20 * m12 + tmp_23 * m22 + tmp_17 * m02));
	return dst;
};

/**
 * Computes the inverse of a 4-by-4 matrix.
 * Note: It is faster to call this than tdl_fast_inverse.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @return {!tdl_fast_Matrix4} The inverse of m.
 */
vector<float> tdl_fast_matrix4_inverse(vector<float> dst, vector<float>  m) {
	return tdl_fast_inverse4(dst, m);
};

/**
 * Multiplies two 4-by-4 matrices; assumes that the given matrices are 4-by-4.
 * Note: It is faster to call this than tdl_fast_mul.
 * @param {!tdl_fast_Matrix4} a The matrix on the left.
 * @param {!tdl_fast_Matrix4} b The matrix on the right.
 * @return {!tdl_fast_Matrix4} The matrix product of a and b.
 */
vector<float> tdl_fast_matrix4_mul(vector<float> dst, vector<float> a, vector<float> b) {
	return tdl_fast_rowMajor_mulMatrix4Vector(dst, a, b);
};

/**
 * Copies a Matrix4.
 * Note: It is faster to call this than tdl_fast_copy.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @return {!tdl_fast_Matrix4} A copy of m.
 */
vector<float> tdl_fast_matrix4_copy(vector<float> dst, vector<float>  m) {
	return tdl_fast_copyMatrix(dst, m);
};

/**
 * Sets the translation component of a 4-by-4 matrix to the given
 * vector.
 * @param {!tdl_fast_Matrix4} a The matrix.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} v The vector.
 * @return {!tdl_fast_Matrix4} a once modified.
 */
vector<float> tdl_fast_matrix4_setTranslation(vector<float> a, vector<float> v) {
	a[12] = v[0];
	a[13] = v[1];
	a[14] = v[2];
	a[15] = 1;
	return a;
};

/**
 * Returns the translation component of a 4-by-4 matrix as a vector with 3
 * entries.
 * @return {!tdl_fast_Vector3} dst vector..
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @return {!tdl_fast_Vector3} The translation component of m.
 */
vector<float> tdl_fast_matrix4_getTranslation(vector<float> dst, vector<float>  m) {
	dst[0] = m[12];
	dst[1] = m[13];
	dst[2] = m[14];
	return dst;
};

/**
 * Creates a 4-by-4 identity matrix.
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @return {!tdl_fast_Matrix4} The 4-by-4 identity.
 */
vector<float> tdl_fast_matrix4_identity(vector<float> dst) {
	return tdl_fast_identity4(dst);
};

vector<float> tdl_fast_matrix4_getAxis(vector<float> dst, vector<float>  m, int axis) {
	int off = axis * 4;
	dst[0] = m[off + 0];
	dst[1] = m[off + 1];
	dst[2] = m[off + 2];
	return dst;
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
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @param {number} angle The camera angle from top to bottom (in radians).
 * @param {number} aspect The aspect ratio width / height.
 * @param {number} zNear The depth (negative z coordinate)
 *     of the near clipping plane.
 * @param {number} zFar The depth (negative z coordinate)
 *     of the far clipping plane.
 * @return {!tdl_fast_Matrix4} The perspective matrix.
 */
vector<float> tdl_fast_matrix4_perspective(vector<float> dst, float angle, float aspect, float zNear, float zFar) {
	float f = tan(PI * 0.5 - 0.5 * angle);
	float rangeInv = 1.0 / (zNear - zFar);

	dst[0] = f / aspect;
	dst[1] = 0;
	dst[2] = 0;
	dst[3] = 0;

	dst[4] = 0;
	dst[5] = f;
	dst[6] = 0;
	dst[7] = 0;

	dst[8] = 0;
	dst[9] = 0;
	dst[10] = (zNear + zFar) * rangeInv;
	dst[11] = -1;

	dst[12] = 0;
	dst[13] = 0;
	dst[14] = zNear * zFar * rangeInv * 2;
	dst[15] = 0;

	return dst;
};


/**
 * Computes a 4-by-4 othogonal transformation matrix given the left, right,
 * bottom, and top dimensions of the near clipping plane as well as the
 * near and far clipping plane distances.
 * @param {!tdl_fast_Matrix4} dst Output matrix.
 * @param {number} left Left side of the near clipping plane viewport.
 * @param {number} right Right side of the near clipping plane viewport.
 * @param {number} top Top of the near clipping plane viewport.
 * @param {number} bottom Bottom of the near clipping plane viewport.
 * @param {number} near The depth (negative z coordinate)
 *     of the near clipping plane.
 * @param {number} far The depth (negative z coordinate)
 *     of the far clipping plane.
 * @return {!tdl_fast_Matrix4} The perspective matrix.
 */
vector<float> tdl_fast_matrix4_ortho(vector<float> dst, float left, float right, float bottom, float top, float near, float far) {


	dst[0] = 2 / (right - left);
	dst[1] = 0;
	dst[2] = 0;
	dst[3] = 0;

	dst[4] = 0;
	dst[5] = 2 / (top - bottom);
	dst[6] = 0;
	dst[7] = 0;

	dst[8] = 0;
	dst[9] = 0;
	dst[10] = -1 / (far - near);
	dst[11] = 0;

	dst[12] = (right + left) / (left - right);
	dst[13] = (top + bottom) / (bottom - top);
	dst[14] = -near / (near - far);
	dst[15] = 1;

	return dst;
}

/**
 * Computes a 4-by-4 perspective transformation matrix given the left, right,
 * top, bottom, near and far clipping planes. The arguments define a frustum
 * extending in the negative z direction. The arguments near and far are the
 * distances to the near and far clipping planes. Note that near and far are not
 * z coordinates, but rather they are distances along the negative z-axis. The
 * matrix generated sends the viewing frustum to the unit box. We assume a unit
 * box extending from -1 to 1 in the x and y dimensions and from 0 to 1 in the z
 * dimension.
 * @param {number} left The x coordinate of the left plane of the box.
 * @param {number} right The x coordinate of the right plane of the box.
 * @param {number} bottom The y coordinate of the bottom plane of the box.
 * @param {number} top The y coordinate of the right plane of the box.
 * @param {number} near The negative z coordinate of the near plane of the box.
 * @param {number} far The negative z coordinate of the far plane of the box.
 * @return {!tdl_fast_Matrix4} The perspective projection matrix.
 */
vector<float> tdl_fast_matrix4_frustum(vector<float> dst, float left, float right, float bottom, float top, float near, float far) {
	float dx = (right - left);
	float dy = (top - bottom);
	float dz = (near - far);

	dst[0] = 2 * near / dx;
	dst[1] = 0;
	dst[2] = 0;
	dst[3] = 0;
	dst[4] = 0;
	dst[5] = 2 * near / dy;
	dst[6] = 0;
	dst[7] = 0;
	dst[8] = (left + right) / dx;
	dst[9] = (top + bottom) / dy;
	dst[10] = far / dz;
	dst[11] = -1;
	dst[12] = 0;
	dst[13] = 0;
	dst[14] = near * far / dz;
	dst[15] = 0;

	return dst;
};

/**
 * Computes a 4-by-4 look-at transformation.  The transformation generated is
 * an orthogonal rotation matrix with translation component.  The translation
 * component sends the eye to the origin.  The rotation component sends the
 * vector pointing from the eye to the target to a vector pointing in the
 * negative z direction, and also sends the up vector into the upper half of
 * the yz plane.
 * @return {!tdl_fast_Matrix4} dst matrix.
 * @param {(!tdl_fast_Vector3} eye The
 *     position of the eye.
 * @param {(!tdl_fast_Vector3} target The
 *     position meant to be viewed.
 * @param {(!tdl_fast_Vector3} up A vector
 *     pointing up.
 * @return {!tdl_fast_Matrix4} The look-at matrix.
 */
vector<float> tdl_fast_matrix4_lookAt(vector<float> dst, vector<float>eye, vector<float>target, vector<float>up) {
	vector<float> t0;
	t0.resize(3);
	vector<float> t1;
	t0.resize(3);
	vector<float> t2;
	t0.resize(3);

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

	return dst;
};

/**
 * Computes a 4-by-4 camera look-at transformation. This is the
 * inverse of lookAt The transformation generated is an
 * orthogonal rotation matrix with translation component.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} eye The position
 *     of the eye.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} target The
 *     position meant to be viewed.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} up A vector
 *     pointing up.
 * @return {!tdl_fast_Matrix4} The camera look-at matrix.
 */
vector<float> tdl_fast_matrix4_cameraLookAt(vector<float> dst, vector<float>eye, vector<float>target, vector<float>up) {
	vector<float> t0;
	t0.resize(3);
	vector<float> t1;
	t0.resize(3);
	vector<float> t2;
	t0.resize(3);

	vector<float> vz = tdl_fast_normalize(t0, tdl_fast_subVector(t0, eye, target));
	vector<float> vx = tdl_fast_normalize(t1, tdl_fast_cross(t1, up, vz));
	vector<float> vy = tdl_fast_cross(t2, vz, vx);

	dst[0] = vx[0];
	dst[1] = vx[1];
	dst[2] = vx[2];
	dst[3] = 0;
	dst[4] = vy[0];
	dst[5] = vy[1];
	dst[6] = vy[2];
	dst[7] = 0;
	dst[8] = vz[0];
	dst[9] = vz[1];
	dst[10] = vz[2];
	dst[11] = 0;
	dst[12] = eye[0];
	dst[13] = eye[1];
	dst[14] = eye[2];
	dst[15] = 1;

	return dst;
};

/**
 * Creates a 4-by-4 matrix which translates by the given vector v.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} v The vector by
 *     which to translate.
 * @return {!tdl_fast_Matrix4} The translation matrix.
 */
vector<float> tdl_fast_matrix4_translation(vector<float> dst, vector<float>v) {
	dst[0] = 1;
	dst[1] = 0;
	dst[2] = 0;
	dst[3] = 0;
	dst[4] = 0;
	dst[5] = 1;
	dst[6] = 0;
	dst[7] = 0;
	dst[8] = 0;
	dst[9] = 0;
	dst[10] = 1;
	dst[11] = 0;
	dst[12] = v[0];
	dst[13] = v[1];
	dst[14] = v[2];
	dst[15] = 1;
	return dst;
};

/**
 * Modifies the given 4-by-4 matrix by translation by the given vector v.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} v The vector by
 *     which to translate.
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_translate(vector<float>m, vector<float>v) {
	float v0 = v[0];
	float v1 = v[1];
	float v2 = v[2];
	float m00 = m[0];
	float m01 = m[1];
	float m02 = m[2];
	float m03 = m[3];
	float m10 = m[1 * 4 + 0];
	float m11 = m[1 * 4 + 1];
	float m12 = m[1 * 4 + 2];
	float m13 = m[1 * 4 + 3];
	float m20 = m[2 * 4 + 0];
	float m21 = m[2 * 4 + 1];
	float m22 = m[2 * 4 + 2];
	float m23 = m[2 * 4 + 3];
	float m30 = m[3 * 4 + 0];
	float m31 = m[3 * 4 + 1];
	float m32 = m[3 * 4 + 2];
	float m33 = m[3 * 4 + 3];

	m[12] = m00 * v0 + m10 * v1 + m20 * v2 + m30;
	m[13] = m01 * v0 + m11 * v1 + m21 * v2 + m31;
	m[14] = m02 * v0 + m12 * v1 + m22 * v2 + m32;
	m[15] = m03 * v0 + m13 * v1 + m23 * v2 + m33;

	return m;
};


/**
 * Creates a 4-by-4 matrix which rotates around the x-axis by the given angle.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} The rotation matrix.
 */
vector<float> tdl_fast_matrix4_rotationX(vector<float> dst, float angle) {
	float c = cos(angle);
	float s = sin(angle);

	dst[0] = 1;
	dst[1] = 0;
	dst[2] = 0;
	dst[3] = 0;
	dst[4] = 0;
	dst[5] = c;
	dst[6] = s;
	dst[7] = 0;
	dst[8] = 0;
	dst[9] = -s;
	dst[10] = c;
	dst[11] = 0;
	dst[12] = 0;
	dst[13] = 0;
	dst[14] = 0;
	dst[15] = 1;

	return dst;
};

/**
 * Modifies the given 4-by-4 matrix by a rotation around the x-axis by the given
 * angle.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_rotateX(vector<float>m, float angle) {
	float m10 = m[4];
	float m11 = m[5];
	float m12 = m[6];
	float m13 = m[7];
	float m20 = m[8];
	float m21 = m[9];
	float m22 = m[10];
	float m23 = m[11];
	float c = cos(angle);
	float s = sin(angle);

	m[4] = c * m10 + s * m20;
	m[5] = c * m11 + s * m21;
	m[6] = c * m12 + s * m22;
	m[7] = c * m13 + s * m23;
	m[8] = c * m20 - s * m10;
	m[9] = c * m21 - s * m11;
	m[10] = c * m22 - s * m12;
	m[11] = c * m23 - s * m13;

	return m;
};

/**
 * Creates a 4-by-4 matrix which rotates around the y-axis by the given angle.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} The rotation matrix.
 */
vector<float> tdl_fast_matrix4_rotationY(vector<float> dst, float angle) {
	float c = cos(angle);
	float s = sin(angle);

	dst[0] = c;
	dst[1] = 0;
	dst[2] = -s;
	dst[3] = 0;
	dst[4] = 0;
	dst[5] = 1;
	dst[6] = 0;
	dst[7] = 0;
	dst[8] = s;
	dst[9] = 0;
	dst[10] = c;
	dst[11] = 0;
	dst[12] = 0;
	dst[13] = 0;
	dst[14] = 0;
	dst[15] = 1;

	return dst;
};

/**
 * Modifies the given 4-by-4 matrix by a rotation around the y-axis by the given
 * angle.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_rotateY(vector<float>m, float angle) {
	float m00 = m[0 * 4 + 0];
	float m01 = m[0 * 4 + 1];
	float m02 = m[0 * 4 + 2];
	float m03 = m[0 * 4 + 3];
	float m20 = m[2 * 4 + 0];
	float m21 = m[2 * 4 + 1];
	float m22 = m[2 * 4 + 2];
	float m23 = m[2 * 4 + 3];
	float c = cos(angle);
	float s = sin(angle);

	m[0] = c * m00 - s * m20;
	m[1] = c * m01 - s * m21;
	m[2] = c * m02 - s * m22;
	m[3] = c * m03 - s * m23;
	m[8] = c * m20 + s * m00;
	m[9] = c * m21 + s * m01;
	m[10] = c * m22 + s * m02;
	m[11] = c * m23 + s * m03;

	return m;
};

/**
 * Creates a 4-by-4 matrix which rotates around the z-axis by the given angle.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} The rotation matrix.
 */
vector<float> tdl_fast_matrix4_rotationZ(vector<float> dst, float angle) {
	float c = cos(angle);
	float s = sin(angle);

	dst[0] = c;
	dst[1] = s;
	dst[2] = 0;
	dst[3] = 0;
	dst[4] = -s;
	dst[5] = c;
	dst[6] = 0;
	dst[7] = 0;
	dst[8] = 0;
	dst[9] = 0;
	dst[10] = 1;
	dst[11] = 0;
	dst[12] = 0;
	dst[13] = 0;
	dst[14] = 0;
	dst[15] = 1;

	return dst;
};

/**
 * Modifies the given 4-by-4 matrix by a rotation around the z-axis by the given
 * angle.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_rotateZ(vector<float>m, float angle) {
	float m00 = m[0 * 4 + 0];
	float m01 = m[0 * 4 + 1];
	float m02 = m[0 * 4 + 2];
	float m03 = m[0 * 4 + 3];
	float m10 = m[1 * 4 + 0];
	float m11 = m[1 * 4 + 1];
	float m12 = m[1 * 4 + 2];
	float m13 = m[1 * 4 + 3];
	float c = cos(angle);
	float s = sin(angle);

	m[0] = c * m00 + s * m10;
	m[1] = c * m01 + s * m11;
	m[2] = c * m02 + s * m12;
	m[3] = c * m03 + s * m13;
	m[4] = c * m10 - s * m00;
	m[5] = c * m11 - s * m01;
	m[6] = c * m12 - s * m02;
	m[7] = c * m13 - s * m03;

	return m;
};

/**
 * Creates a 4-by-4 matrix which rotates around the given axis by the given
 * angle.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} axis The axis
 *     about which to rotate.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} A matrix which rotates angle radians
 *     around the axis.
 */
vector<float> tdl_fast_matrix4_axisRotation(vector<float> dst, vector<float>axis, float angle) {
	float x = axis[0];
	float y = axis[1];
	float z = axis[2];
	float n = sqrt(x * x + y * y + z * z);
	x /= n;
	y /= n;
	z /= n;
	float xx = x * x;
	float yy = y * y;
	float zz = z * z;
	float c = cos(angle);
	float s = sin(angle);
	float oneMinusCosine = 1 - c;

	dst[0] = xx + (1 - xx) * c;
	dst[1] = x * y * oneMinusCosine + z * s;
	dst[2] = x * z * oneMinusCosine - y * s;
	dst[3] = 0;
	dst[4] = x * y * oneMinusCosine - z * s;
	dst[5] = yy + (1 - yy) * c;
	dst[6] = y * z * oneMinusCosine + x * s;
	dst[7] = 0;
	dst[8] = x * z * oneMinusCosine + y * s;
	dst[9] = y * z * oneMinusCosine - x * s;
	dst[10] = zz + (1 - zz) * c;
	dst[11] = 0;
	dst[12] = 0;
	dst[13] = 0;
	dst[14] = 0;
	dst[15] = 1;

	return dst;
};

/**
 * Modifies the given 4-by-4 matrix by rotation around the given axis by the
 * given angle.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} axis The axis
 *     about which to rotate.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_axisRotate(vector<float>m, vector<float>axis, float angle) {
	float x = axis[0];
	float y = axis[1];
	float z = axis[2];
	float n = sqrt(x * x + y * y + z * z);
	x /= n;
	y /= n;
	z /= n;
	float xx = x * x;
	float yy = y * y;
	float zz = z * z;
	float c = cos(angle);
	float s = sin(angle);
	float oneMinusCosine = 1 - c;

	float r00 = xx + (1 - xx) * c;
	float r01 = x * y * oneMinusCosine + z * s;
	float r02 = x * z * oneMinusCosine - y * s;
	float r10 = x * y * oneMinusCosine - z * s;
	float r11 = yy + (1 - yy) * c;
	float r12 = y * z * oneMinusCosine + x * s;
	float r20 = x * z * oneMinusCosine + y * s;
	float r21 = y * z * oneMinusCosine - x * s;
	float r22 = zz + (1 - zz) * c;

	float m00 = m[0];
	float m01 = m[1];
	float m02 = m[2];
	float m03 = m[3];
	float m10 = m[4];
	float m11 = m[5];
	float m12 = m[6];
	float m13 = m[7];
	float m20 = m[8];
	float m21 = m[9];
	float m22 = m[10];
	float m23 = m[11];
	float m30 = m[12];
	float m31 = m[13];
	float m32 = m[14];
	float m33 = m[15];

	m[0] = r00 * m00 + r01 * m10 + r02 * m20;
	m[1] = r00 * m01 + r01 * m11 + r02 * m21;
	m[2] = r00 * m02 + r01 * m12 + r02 * m22;
	m[3] = r00 * m03 + r01 * m13 + r02 * m23;
	m[4] = r10 * m00 + r11 * m10 + r12 * m20;
	m[5] = r10 * m01 + r11 * m11 + r12 * m21;
	m[6] = r10 * m02 + r11 * m12 + r12 * m22;
	m[7] = r10 * m03 + r11 * m13 + r12 * m23;
	m[8] = r20 * m00 + r21 * m10 + r22 * m20;
	m[9] = r20 * m01 + r21 * m11 + r22 * m21;
	m[10] = r20 * m02 + r21 * m12 + r22 * m22;
	m[11] = r20 * m03 + r21 * m13 + r22 * m23;

	return m;
};

/**
 * Creates a 4-by-4 matrix which scales in each dimension by an amount given by
 * the corresponding entry in the given vector; assumes the vector has three
 * entries.
 * @param {!tdl_fast_Vector3} v A vector of
 *     three entries specifying the factor by which to scale in each dimension.
 * @return {!tdl_fast_Matrix4} The scaling matrix.
 */
vector<float> tdl_fast_matrix4_scaling(vector<float> dst, vector<float> v) {
	dst[0] = v[0];
	dst[1] = 0;
	dst[2] = 0;
	dst[3] = 0;
	dst[4] = 0;
	dst[5] = v[1];
	dst[6] = 0;
	dst[7] = 0;
	dst[8] = 0;
	dst[9] = 0;
	dst[10] = v[2];
	dst[11] = 0;
	dst[12] = 0;
	dst[13] = 0;
	dst[14] = 0;
	dst[15] = 1;
	return dst;
};

/**
 * Modifies the given 4-by-4 matrix, scaling in each dimension by an amount
 * given by the corresponding entry in the given vector; assumes the vector has
 * three entries.
 * @param {!tdl_fast_Matrix4} m The matrix to be modified.
 * @param {!tdl_fast_Vector3} v A vector of three entries specifying the
 *     factor by which to scale in each dimension.
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_scale(vector<float>m, vector<float>v) {
	float v0 = v[0];
	float v1 = v[1];
	float v2 = v[2];

	m[0] = v0 * m[0 * 4 + 0];
	m[1] = v0 * m[0 * 4 + 1];
	m[2] = v0 * m[0 * 4 + 2];
	m[3] = v0 * m[0 * 4 + 3];
	m[4] = v1 * m[1 * 4 + 0];
	m[5] = v1 * m[1 * 4 + 1];
	m[6] = v1 * m[1 * 4 + 2];
	m[7] = v1 * m[1 * 4 + 3];
	m[8] = v2 * m[2 * 4 + 0];
	m[9] = v2 * m[2 * 4 + 1];
	m[10] = v2 * m[2 * 4 + 2];
	m[11] = v2 * m[2 * 4 + 3];

	return m;
};


