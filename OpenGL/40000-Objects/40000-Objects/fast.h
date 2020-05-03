#ifndef FAST
#define FAST
#include <vector>
using namespace std;

/**
 * Adds two vectors; assumes a and b have the same dimension.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 */
vector<float> tdl_fast_addVector(vector<float> dst, vector<float> a, vector<float> b) ;

/**
 * Subtracts two vectors.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 */
vector<float> tdl_fast_subVector(vector<float> dst, vector<float> a, vector<float> b) ;

/**
 * Performs linear interpolation on two vectors.
 * Given vectors a and b and interpolation coefficient t, returns
 * (1 - t) * a + t * b.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @param {number} t Interpolation coefficient.
 */
vector<float> tdl_fast_lerpVector(vector<float> dst, vector<float> a, vector<float> b, float t) ;

/**
 * Divides a vector by a scalar.
 * @param {!tdl_fast_Vector} dst The vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @param {number} k The scalar.
 * @return {!tdl_fast_Vector} dst.
 */
vector<float> tdl_fast_divVectorScalar(vector<float> dst, vector<float> v, float k) ;

/**
 * Computes the cross product of two vectors; assumes both vectors have
 * three entries.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @return {!tdl_fast_Vector} The vector a cross b.
 */
vector<float> tdl_fast_cross(vector<float> dst, vector<float> a, vector<float> b) ;

/**
 * Computes the dot product of two vectors; assumes both vectors have
 * three entries.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @return {number} dot product
 */
float tdl_fast_dot(vector<float> a, vector<float> b) ;

/**
 * Divides a vector by its Euclidean length and returns the quotient.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a The vector.
 * @return {!tdl_fast_Vector} The normalized vector.
 */
vector<float> tdl_fast_normalize(vector<float> dst, vector<float>  a) ;

/**
 * Negates a vector.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} -v.
 */
vector<float> tdl_fast_negativeVector(vector<float> dst, vector<float> v) ;

/**
 * Negates a matrix.
 * @param {!tdl_fast_Matrix} dst matrix.
 * @param {!tdl_fast_Matrix} v The matrix.
 * @return {!tdl_fast_Matrix} -v.
 */
vector<float> tdl_fast_negativeMatrix(vector<float> dst, vector<float>  v) ;

/**
 * Copies a vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} A copy of v.
 */
vector<float> tdl_fast_copyVector(vector<float> dst, vector<float> v) ;

/**
 * Copies a matrix.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Matrix} A copy of m.
 */
vector<float> tdl_fast_copyMatrix(vector<float> dst, vector<float> m) ;

/**
 * Multiplies a scalar by a vector.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {number} k The scalar.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} The product of k and v.
 */
vector<float> tdl_fast_mulScalarVector(vector<float> dst, float k, vector<float> v) ;

/**
 * Multiplies a vector by a scalar.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @param {number} k The scalar.
 * @return {!tdl_fast_Vector} The product of k and v.
 */
vector<float> tdl_fast_mulVectorScalar(vector<float> dst, vector<float> v, float k) ;

/**
 * Multiplies a scalar by a matrix.
 * @param {!tdl_fast_Matrix} dst matrix.
 * @param {number} k The scalar.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Matrix} The product of m and k.
 */
vector<float> tdl_fast_mulScalarMatrix(vector<float> dst, float k, vector<float> m) ;

/**
 * Multiplies a matrix by a scalar.
 * @param {!tdl_fast_Matrix} dst matrix.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} k The scalar.
 * @return {!tdl_fast_Matrix} The product of m and k.
 */
vector<float> tdl_fast_mulMatrixScalar(vector<float> dst, vector<float>m, float k) ;

/**
 * Multiplies a vector by another vector (component-wise); assumes a and
 * b have the same length.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @return {!tdl_fast_Vector} The vector of products of entries of a and
 *     b.
 */
vector<float>  tdl_fast_mulVectorVector(vector<float> dst, vector<float> a, vector<float> b) ;

/**
 * Divides a vector by another vector (component-wise); assumes a and
 * b have the same length.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} a Operand vector.
 * @param {!tdl_fast_Vector} b Operand vector.
 * @return {!tdl_fast_Vector} The vector of quotients of entries of a and
 *     b.
 */
vector<float> tdl_fast_divVectorVector(vector<float> dst, vector<float> a, vector<float> b) ;

/**
 * Multiplies a vector by a matrix; treats the vector as a row vector; assumes
 * matrix entries are accessed in [row][column] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Vector} The product of v and m as a row vector.
 */
vector<float> tdl_fast_rowMajor_mulVectorMatrix4(vector<float> dst, vector<float>v, vector<float>m) ;

/**
 * Multiplies a vector by a matrix; treats the vector as a row vector; assumes
 * matrix entries are accessed in [column][row] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Vector} v The vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Vector} The product of v and m as a row vector.
 */
vector<float> tdl_fast_columnMajor_mulVectorMatrix4(vector<float> dst, vector<float> v, vector<float> m) ;



/**
 * Multiplies a matrix by a vector; treats the vector as a column vector.
 * assumes matrix entries are accessed in [row][column] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} The product of m and v as a column vector.
 */
vector<float> tdl_fast_rowMajor_mulMatrix4Vector(vector<float> dst, vector<float>m, vector<float>v) ;

/**
 * Multiplies a matrix by a vector; treats the vector as a column vector;
 * assumes matrix entries are accessed in [column][row] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {!tdl_fast_Vector} v The vector.
 * @return {!tdl_fast_Vector} The product of m and v as a column vector.
 */
vector<float> tdl_fast_columnMajor_mulMatrix4Vector(vector<float> dst, vector<float> m, vector<float> v) ;

/**
 * Multiplies two 3-by-3 matrices; assumes that the given matrices are 3-by-3;
 * assumes matrix entries are accessed in [row][column] fashion.
 * @param {!tdl_fast_Matrix3} dst matrix.
 * @param {!tdl_fast_Matrix3} a The matrix on the left.
 * @param {!tdl_fast_Matrix3} b The matrix on the right.
 * @return {!tdl_fast_Matrix3} The matrix product of a and b.
 */
vector<float> tdl_fast_rowMajor_mulMatrixMatrix3(vector<float> dst, vector<float> a, vector<float> b) ;

/**
 * Multiplies two 3-by-3 matrices; assumes that the given matrices are 3-by-3;
 * assumes matrix entries are accessed in [column][row] fashion.
 * @param {!tdl_fast_Matrix3} dst matrix.
 * @param {!tdl_fast_Matrix3} a The matrix on the left.
 * @param {!tdl_fast_Matrix3} b The matrix on the right.
 * @return {!tdl_fast_Matrix3} The matrix product of a and b.
 */
vector<float> tdl_fast_columnMajor_mulMatrixMatrix3(vector<float> dst, vector<float> a, vector<float> b) ;



/**
 * Multiplies two 4-by-4 matrices; assumes that the given matrices are 4-by-4;
 * assumes matrix entries are accessed in [row][column] fashion.
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @param {!tdl_fast_Matrix4} a The matrix on the left.
 * @param {!tdl_fast_Matrix4} b The matrix on the right.
 * @return {!tdl_fast_Matrix4} The matrix product of a and b.
 */
vector<float> tdl_fast_rowMajor_mulMatrixMatrix4(vector<float> dst, vector<float> a, vector<float> b) ;

/**
 * Multiplies two 4-by-4 matrices; assumes that the given matrices are 4-by-4;
 * assumes matrix entries are accessed in [column][row] fashion.
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @param {!tdl_fast_Matrix4} a The matrix on the left.
 * @param {!tdl_fast_Matrix4} b The matrix on the right.
 * @return {!tdl_fast_Matrix4} The matrix product of a and b.
 */
vector<float> tdl_fast_columnMajor_mulMatrixMatrix4(vector<float> dst, vector<float> a, vector<float> b) ;


/**
 * Gets the jth column of the given matrix m; assumes matrix entries are
 * accessed in [row][column] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} j The index of the desired column.
 * @return {!tdl_fast_Vector} The jth column of m as a vector.
 */
vector<float> tdl_fast_rowMajor_column4(vector<float> dst, vector<float>  m, int j) ;

/**
 * Gets the jth column of the given matrix m; assumes matrix entries are
 * accessed in [column][row] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} j The index of the desired column.
 * @return {!tdl_fast_Vector} The jth column of m as a vector.
 */
vector<float> tdl_fast_columnMajor_column4(vector<float> dst, vector<float>  m, int j) ;


/**
 * Gets the ith row of the given matrix m; assumes matrix entries are
 * accessed in [row][column] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} i The index of the desired row.
 * @return {!tdl_fast_Vector} The ith row of m.
 */
vector<float> tdl_fast_rowMajor_row4(vector<float> dst, vector<float>  m, int i) ;

/**
 * Gets the ith row of the given matrix m; assumes matrix entries are
 * accessed in [column][row] fashion.
 * @param {!tdl_fast_Vector} dst vector.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @param {number} i The index of the desired row.
 * @return {!tdl_fast_Vector} The ith row of m.
 */
vector<float> tdl_fast_columnMajor_row4(vector<float> dst, vector<float>  m, int i) ;


/**
 * Creates an n-by-n identity matrix.
 *
 * @param {!tdl_fast_Matrix} dst matrix.
 * @return {!tdl_fast_Matrix} An n-by-n identity matrix.
 */
vector<float> tdl_fast_identity4(vector<float> dst) ;

/**
 * Takes the transpose of a matrix.
 * @param {!tdl_fast_Matrix} dst matrix.
 * @param {!tdl_fast_Matrix} m The matrix.
 * @return {!tdl_fast_Matrix} The transpose of m.
 */
vector<float> tdl_fast_transpose4(vector<float> dst, vector<float>  m) ;

/**
 * Computes the inverse of a 4-by-4 matrix.
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @return {!tdl_fast_Matrix4} The inverse of m.
 */
vector<float> tdl_fast_inverse4(vector<float> dst, vector<float>  m) ;

/**
 * Computes the inverse of a 4-by-4 matrix.
 * Note: It is faster to call this than tdl_fast_inverse.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @return {!tdl_fast_Matrix4} The inverse of m.
 */
vector<float> tdl_fast_matrix4_inverse(vector<float> dst, vector<float>  m) ;

/**
 * Multiplies two 4-by-4 matrices; assumes that the given matrices are 4-by-4.
 * Note: It is faster to call this than tdl_fast_mul.
 * @param {!tdl_fast_Matrix4} a The matrix on the left.
 * @param {!tdl_fast_Matrix4} b The matrix on the right.
 * @return {!tdl_fast_Matrix4} The matrix product of a and b.
 */
vector<float> tdl_fast_matrix4_mul(vector<float> dst, vector<float> a, vector<float> b) ;

/**
 * Copies a Matrix4.
 * Note: It is faster to call this than tdl_fast_copy.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @return {!tdl_fast_Matrix4} A copy of m.
 */
vector<float> tdl_fast_matrix4_copy(vector<float> dst, vector<float>  m) ;

/**
 * Sets the translation component of a 4-by-4 matrix to the given
 * vector.
 * @param {!tdl_fast_Matrix4} a The matrix.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} v The vector.
 * @return {!tdl_fast_Matrix4} a once modified.
 */
vector<float> tdl_fast_matrix4_setTranslation(vector<float> a, vector<float> v) ;

/**
 * Returns the translation component of a 4-by-4 matrix as a vector with 3
 * entries.
 * @return {!tdl_fast_Vector3} dst vector..
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @return {!tdl_fast_Vector3} The translation component of m.
 */
vector<float> tdl_fast_matrix4_getTranslation(vector<float> dst, vector<float>  m);

/**
 * Creates a 4-by-4 identity matrix.
 * @param {!tdl_fast_Matrix4} dst matrix.
 * @return {!tdl_fast_Matrix4} The 4-by-4 identity.
 */
vector<float> tdl_fast_matrix4_identity(vector<float> dst);

vector<float> tdl_fast_matrix4_getAxis(vector<float> dst, vector<float>  m, int axis);

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
vector<float> tdl_fast_matrix4_perspective(vector<float> dst, float angle, float aspect, float zNear, float zFar);


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
vector<float> tdl_fast_matrix4_ortho(vector<float> dst, float left, float right, float bottom, float top, float near, float far);

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
vector<float> tdl_fast_matrix4_frustum(vector<float> dst, float left, float right, float bottom, float top, float near, float far);

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
vector<float> tdl_fast_matrix4_lookAt(vector<float> dst, vector<float>eye, vector<float>target, vector<float>up);

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
vector<float> tdl_fast_matrix4_cameraLookAt(vector<float> dst, vector<float>eye, vector<float>target, vector<float>up);

/**
 * Creates a 4-by-4 matrix which translates by the given vector v.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} v The vector by
 *     which to translate.
 * @return {!tdl_fast_Matrix4} The translation matrix.
 */
vector<float> tdl_fast_matrix4_translation(vector<float> dst, vector<float>v);

/**
 * Modifies the given 4-by-4 matrix by translation by the given vector v.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} v The vector by
 *     which to translate.
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_translate(vector<float>m, vector<float>v);


/**
 * Creates a 4-by-4 matrix which rotates around the x-axis by the given angle.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} The rotation matrix.
 */
vector<float> tdl_fast_matrix4_rotationX(vector<float> dst, float angle);

/**
 * Modifies the given 4-by-4 matrix by a rotation around the x-axis by the given
 * angle.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_rotateX(vector<float>m, float angle);

/**
 * Creates a 4-by-4 matrix which rotates around the y-axis by the given angle.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} The rotation matrix.
 */
vector<float> tdl_fast_matrix4_rotationY(vector<float> dst, float angle);

/**
 * Modifies the given 4-by-4 matrix by a rotation around the y-axis by the given
 * angle.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_rotateY(vector<float>m, float angle);

/**
 * Creates a 4-by-4 matrix which rotates around the z-axis by the given angle.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} The rotation matrix.
 */
vector<float> tdl_fast_matrix4_rotationZ(vector<float> dst, float angle);

/**
 * Modifies the given 4-by-4 matrix by a rotation around the z-axis by the given
 * angle.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_rotateZ(vector<float>m, float angle);

/**
 * Creates a 4-by-4 matrix which rotates around the given axis by the given
 * angle.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} axis The axis
 *     about which to rotate.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} A matrix which rotates angle radians
 *     around the axis.
 */
vector<float> tdl_fast_matrix4_axisRotation(vector<float> dst, vector<float>axis, float angle);

/**
 * Modifies the given 4-by-4 matrix by rotation around the given axis by the
 * given angle.
 * @param {!tdl_fast_Matrix4} m The matrix.
 * @param {(!tdl_fast_Vector3|!tdl_fast_Vector4)} axis The axis
 *     about which to rotate.
 * @param {number} angle The angle by which to rotate (in radians).
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_axisRotate(vector<float>m, vector<float>axis, float angle);

/**
 * Creates a 4-by-4 matrix which scales in each dimension by an amount given by
 * the corresponding entry in the given vector; assumes the vector has three
 * entries.
 * @param {!tdl_fast_Vector3} v A vector of
 *     three entries specifying the factor by which to scale in each dimension.
 * @return {!tdl_fast_Matrix4} The scaling matrix.
 */
vector<float> tdl_fast_matrix4_scaling(vector<float> dst, vector<float>  v);
/**
 * Modifies the given 4-by-4 matrix, scaling in each dimension by an amount
 * given by the corresponding entry in the given vector; assumes the vector has
 * three entries.
 * @param {!tdl_fast_Matrix4} m The matrix to be modified.
 * @param {!tdl_fast_Vector3} v A vector of three entries specifying the
 *     factor by which to scale in each dimension.
 * @return {!tdl_fast_Matrix4} m once modified.
 */
vector<float> tdl_fast_matrix4_scale(vector<float>m, vector<float>v); 


#endif