// #ifndef _GENLIB_H
//#define _GENLIB_H

// Need gl.h for GLfloat. gl.h requires windows.h
//#include <windows.h>
//#include <gl/gl.h>


#define	ROUND(a)		(int)((a) + 0.5f)
#define PI				3.1415926535f
#define ONERAD			0.0174532925f
#define RAD(x)			x * ONERAD

typedef struct _Matrix { GLfloat v[16]; } Matrix;

typedef struct _Point { GLfloat x, y; } Point;

typedef struct _Coord { GLfloat x, y, z; } Coord;

typedef struct Plane {
	GLfloat a, b, c, d;
	Plane() { a = b = c = d = 0.0f; }
	Plane(GLfloat _a, GLfloat _b, GLfloat _c, GLfloat _d) { a = _a; b = _b; c = _c; d = _d; }
} Plane;

typedef struct float4 {
	float v[4];
	float4(float a, float b, float c, float d) { v[0] = a; v[1] = b; v[2] = c; v[3] = d; }
} float4;

inline Coord newCoord() { Coord c; c.x = c.y = c.z = 0.0f; return c; }
inline Coord newCoord(GLfloat x, GLfloat y, GLfloat z) { Coord c; c.x = x; c.y = y; c.z = z; return c; }

inline const Coord &operator += (Coord &c1, const Coord &c2) { c1.x += c2.x; c1.y += c2.y; c1.z += c2.z; return c1; }
inline const Coord &operator -= (Coord &c1, const Coord &c2) { c1.x -= c2.x; c1.y -= c2.y; c1.z -= c2.z; return c1; }
inline const Coord &operator *= (Coord &c1, GLfloat a) { c1.x *= a; c1.y *= a; c1.z *= a; return c1; }
inline const Coord &operator /= (Coord &c1, GLfloat a) { c1.x /= a; c1.y /= a; c1.z /= a; return c1; }
inline Coord operator  + (const Coord &c1, const Coord &c2) { Coord c3 = c1; c3 += c2; return c3; }
inline Coord operator  - (const Coord &c1, const Coord &c2) { Coord c3 = c1; c3 -= c2; return c3; }
inline Coord operator  * (const Coord &c1, GLfloat a) { Coord c3 = c1; c3 *= a; return c3; }
inline Coord operator  * (const Coord &c1, const Coord &c2) { Coord c3 = c1; c3.x *= c2.x; c3.y *= c2.y; c3.z *= c2.z; return c3; }
inline Coord operator  / (const Coord &c1, GLfloat a) { Coord c3 = c1; c3 /= a; return c3; }

inline const Plane &operator /= (Plane &p, GLfloat a) { p.a /= a; p.b /= a; p.c /= a; p.d /= a; return p; }

/* Matrix cross-multiplication */
inline Matrix operator * (Matrix &m1, Matrix &m2) {
	Matrix m;

	m.v[0]  = m1.v[0] * m2.v[0]  + m1.v[1] * m2.v[4]  + m1.v[2] * m2.v[8]   + m1.v[3] * m2.v[12];
	m.v[1]  = m1.v[0] * m2.v[1]  + m1.v[1] * m2.v[5]  + m1.v[2] * m2.v[9]   + m1.v[3] * m2.v[13];
	m.v[2]  = m1.v[0] * m2.v[2]  + m1.v[1] * m2.v[6]  + m1.v[2] * m2.v[10]  + m1.v[3] * m2.v[14];
	m.v[3]  = m1.v[0] * m2.v[3]  + m1.v[1] * m2.v[7]  + m1.v[2] * m2.v[11]  + m1.v[3] * m2.v[15];

	m.v[4]  = m1.v[4] * m2.v[0]  + m1.v[5] * m2.v[4]  + m1.v[6] * m2.v[8]   + m1.v[7] * m2.v[12];
	m.v[5]  = m1.v[4] * m2.v[1]  + m1.v[5] * m2.v[5]  + m1.v[6] * m2.v[9]   + m1.v[7] * m2.v[13];
	m.v[6]  = m1.v[4] * m2.v[2]  + m1.v[5] * m2.v[6]  + m1.v[6] * m2.v[10]  + m1.v[7] * m2.v[14];
	m.v[7]  = m1.v[4] * m2.v[3]  + m1.v[5] * m2.v[7]  + m1.v[6] * m2.v[11]  + m1.v[7] * m2.v[15];

	m.v[8]  = m1.v[8] * m2.v[0]  + m1.v[9] * m2.v[4]  + m1.v[10] * m2.v[8]  + m1.v[11] * m2.v[12];
	m.v[9]  = m1.v[8] * m2.v[1]  + m1.v[9] * m2.v[5]  + m1.v[10] * m2.v[9]  + m1.v[11] * m2.v[13];
	m.v[10] = m1.v[8] * m2.v[2]  + m1.v[9] * m2.v[6]  + m1.v[10] * m2.v[10] + m1.v[11] * m2.v[14];
	m.v[11] = m1.v[8] * m2.v[3]  + m1.v[9] * m2.v[7]  + m1.v[10] * m2.v[11] + m1.v[11] * m2.v[15];

	m.v[12] = m1.v[12] * m2.v[0] + m1.v[13] * m2.v[4] + m1.v[14] * m2.v[8]  + m1.v[15] * m2.v[12];
	m.v[13] = m1.v[12] * m2.v[1] + m1.v[13] * m2.v[5] + m1.v[14] * m2.v[9]  + m1.v[15] * m2.v[13];
	m.v[14] = m1.v[12] * m2.v[2] + m1.v[13] * m2.v[6] + m1.v[14] * m2.v[10] + m1.v[15] * m2.v[14];
	m.v[15] = m1.v[12] * m2.v[3] + m1.v[13] * m2.v[7] + m1.v[14] * m2.v[11] + m1.v[15] * m2.v[15];

	return m;
}

/*
typedef struct { GLfloat w, x, y, z; } Quaternion;

inline Quaternion operator * (Quaternion &q1, Quaternion &q2) {
	Quaternion q3;
	q3.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
	q3.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
	q3.y = q1.w * q2.y + q1.y * q2.w + q1.z * q2.x - q1.x * q2.z;
	q3.z = q1.w * q2.z + q1.z * q2.w + q1.x * q2.y - q1.y * q2.x;
	return q3;
}
*/

typedef struct {
	GLfloat xx, xy, xz, yx, yy, yz, zx, zy, zz, xx_xy, yx_yy, zx_zy;
} RotationInits;

typedef struct Rotation {
	bool enabled;
	Coord angles;
	RotationInits inits;
	Rotation() { angles.x = angles.y = angles.z = 0.0; }
	Rotation(GLfloat _xy, GLfloat _xz, GLfloat _yz) { angles.z = _xy; angles.y = _xz; angles.x = _yz; }
} Rotation;


#define BR	-1	// Break in the tristrip... will be found in fc_vtI
typedef struct {
	// Number of faces, clockwise vertex indices (for collisions), vertex colour
	int fc_n, *fc_vt_n, *fc_edge_n;
	GLubyte **fc_vtI, **fc_edgeI;
//	Colour **fc_clr;
	// Number of vertices, vertices on x-, y- and z-planes
	int vt_n, nv_n;
	Coord *nv;	// Normals and vertices interleaved
	Coord *sn;	// Surface normals
	Point *tx;	// Texture coordinates
	// Central reference point of the object (axis of rotation and scaling)
	Coord cp; // previously .ref, seems to be a reserved keyword

	/////////////Quaternal qt[4];
} Object;

typedef struct {
	GLuint textureId;
	Object obj;								// Object to draw
	Rotation rotation;
	Coord translation;
	GLfloat scale;	// Animation properties
	bool exists, interactive;
	unsigned char Eid;
	int EnumSprings;
	unsigned char Esprings[10];	// Max 10 springs per object
	float Edensity;
} Entity; // E meaning (Physical) Energy (Movement) (Animation)

typedef struct {
	Coord *nv;
	GLushort *indices;
	int numIndices;
} Barrier;

typedef struct Quadrant {
	bool hasChildren;
	struct Quadrant *frontLeft, *frontRight, *backLeft, *backRight;
	// Any quadrant, with or without children, can hold entities
	//int numEntities;
	//Entity *entities;
	// Quadrants that do NOT have children can hold barriers (walls, chunks of terrain)
	int numBarriers;
	Barrier *barriers;

	Coord ltn, rbf, centre;	// Points
} Quadrant;

typedef struct Terrain {
	Coord *nv;
	Point *tx;
	GLuint textureId;
	int nv_n, walls_n;
	int width, depth;
	GLushort *indices, *wallStarts, *wallCounts;
	int boxWidth, boxDepth, numBoxesX, numBoxesZ, contentSize;
} Terrain;


typedef struct {
	Entity *entities;
	Terrain *terrains;
	GLuint *textures;
	Quadrant *quadrants;	// Order is unknown. All we know is quadrants[0] is the one we start the search in
	Coord quadSize;
	int numEntities, numTerrains, numTextures, numQuadrants, numQuadrantLevels;
} World;


/* Return dot product of two vectors */
inline GLfloat dot(Coord &c1, Coord &c2);

/* Return magnitude of a vector */
inline GLfloat mag(Coord &c);

/* Return direction of a vector */
//extern "C" inline Coord dir(Coord &c);
inline Coord dir(Coord &c);


/* Rotate vector by given matrix */
//Coord rotate(GLMatrix &m, Coord &c);

/* *******************************************************************
 *	TranslateObject - Moves an object by the given offset
 *		- Passes the vertices and normals to the Translate() function
 *		- Moves the object's reference point so future rotations and
 *		scales are not affected
 *	Parameters:	Object to move, xyz offset
 *	Return:		None.
 *	Effect:		Moves all vertices, normals and the reference point
 *				by the given offset
 */
void translateObject(Object *o, Coord t);

/* *******************************************************************
 *	Translate - Moves vertices by given offsets
 *		- Iterates through each point of the given vertex array and moves its
 *		position to that which corresponds to the old point plus each of the
 *		three specified offsets
 *	Parameters:	Vertex array, number of vertices, reference point, x offset,
 *				y offset, z offset
 *	Return:		None.
 *	Effect:		Changes all vertices by the given offsets
 */
void translate(Coord *vt, int vt_n, GLfloat x, GLfloat y, GLfloat z);

/* *******************************************************************
 *	RotateObject - Rotates an object by given angles (in radians)
 *		- Passes the vertices and the normals to the Scale() function
 *		- Supports both negative and positive points and angle values
 *	Parameters:	Object containing the list of points, angles in radians
 *	Return:		None.
 *	Effect:		Changes all vertices and normals of the given object
 */
void rotateObject(Object *o, GLfloat angle_xy, GLfloat angle_xz, GLfloat angle_yz);

/* *******************************************************************
 *	Rotate - Rotates vertices by given angles
 *		- Iterates through each point of the given vertex array and moves its
 *		position to that which corresponds to the old point offset by each of
 *		the three specified angles (relative to the given reference point)
 *	Parameters:	Vertex array, number of vertices, reference point, angle on xy
 *				plane, angle on xz plane, angle on yz plane
 *	Return:		None.
 *	Effect:		Changes all vertices by the given angles
 */
void rotate(Coord *vt, int vt_n, Coord *zero, GLfloat angle_xy, GLfloat angle_xz, GLfloat angle_yz);

/* *******************************************************************
 *	CalculateRotationInits - Precalculates multiplications for rotation
 *		- When rotating in 3D, there are many elements of the procedure
 *		that only need to be calculated once and can be reused to speed
 *		it up.
 *		- Takes a Rotation object and uses its xy, xz and yz angle
 *		variables to populate its inits structure.
 *	Parameters:	Rotation object to calculate inits for
 *	Return:		None.
 *	Effect:		Updates the inits structure inside the Rotation object
 */
void calculateRotationInits(Rotation *r);

/* *******************************************************************
 *	QuickRotateObject - Quickly rotate an object
 *		- Works the same way as RotateObject(), but uses the precalculated
 *		rotation inits to speed up the process.
 *	Parameters:	Object to rotate, precalculated variables to use
 *	Return:		None.
 *	Effect:		Changes all vertices by the given angles
 */
void quickRotateObject(Object *o, RotationInits *r);

/* *******************************************************************
 *	QuickRotate - Quickly rotate a vertex
 *		- Uses the precalculated rotation inits to rotate a single vertex.
 *	Parameters:	Vertex to rotate, precalculated rotation inits
 *	Return:		Rotated vertex.
 *	Effect:		None.
 */
Coord quickRotate(Coord vt, RotationInits *r);

/* *******************************************************************
 *	ScaleObject - Changes object size by the given ratio
 *		- Passes the vertices and the normals to the Scale() function
 *	Parameters:	Object containing the list of points, ratio against
 *				:1 (the decimal equivalent to a percentage).
 *	Return:		None.
 *	Effect:		Changes all vertices and normals of the given object
 */
void scaleObject(Object *o, GLfloat scale);

/* *******************************************************************
 *	Scale - Changes object size by the given ratio
 *		- Iterates through each point of the given vertex array and moves its
 *		position to that which, in polar form, would have the same angle
 *		relative to the origin but a distance corresponding to the ratio
 *		- A ratio value of 1.0 equals no change: it will be 1:1 or 100%
 *		of its previous size. 1.1 = 1.1:1 or 110%, 0.5 = 0.5:1 or 50%
 *	Parameters:	Vertex array, number of vertices, reference point, ratio
 *				against :1 (the decimal equivalent to a percentage).
 *	Return:		None.
 *	Effect:		Changes all points by the given ratio
 */
void scale(Coord *vt, int vt_n, Coord *zero, GLfloat scale);

/* *******************************************************************
 *	SkewObject - Changes object size by given ratios
 *		- Passes the vertices and the normals to the Skew() function
 *	Parameters:	Object containing the list of points, x,y,z ratio against
 *				:1 (the decimal equivalent to a percentage).
 *	Return:		None.
 *	Effect:		Changes all vertices and normals of the given object
 */
void skewObject(Object *o, GLfloat ratio_x, GLfloat ratio_y, GLfloat ratio_z);

/* *******************************************************************
 *	Skew - Changes object size by given ratios
 *		- Iterates through each point of the given vertex array and moves its
 *		x, y and z position by individual ratios
 *		- A ratio value of 1.0 equals no change: it will be 1:1 or 100%
 *		of its previous size. 1.1 = 1.1:1 or 110%, 0.5 = 0.5:1 or 50%
 *	Parameters:	Vertex array, x, y, z ratios against :1 (the decimal equivalent
 *				to a percentage).
 *	Return:		None.
 *	Effect:		Changes all points by the given ratios
 */
void skew(Coord *vt, int vt_n, Coord *zero, GLfloat skewx, GLfloat skewy, GLfloat skewz);


/* *******************************************************************
 *	calculateNormal - Calculates the normal for a face
 *		- Given two non-parallel vectors (vt1 - vt2 and vt3 - vt2),
 *		this function calculates the angle which is normal to the
 *		vectors (which belong to a face). If a reference point is given,
 *		the point where a a ray in this direction crosses the face's
 *		plane is the returned coordinate. This is done so that regardless
 *		of any movement or rotation (object or camera), the angle between
 *		the object's reference point and this point will always be
 *		normal to the face.
 *	Parameters: Vertices that make up the vectors v1 - v2 and v3 - v2
 *	Return:		The point at which a ray from the object's reference in
 *				the direction of the normal crosses the relevant face's
 *				plane.
 *	Effect:		None.
 */
Coord calculateNormal(Coord vt1, Coord vt2, Coord vt3, Coord *zero);

/* *******************************************************************
 *	calculatePlane - Calculates the constants for a plane
 *		- Given two non-parallel vectors (vt1 - vt2 and vt3 - vt2),
 *		this function calculates the angle which is normal to the
 *		vectors and then returns the A, B, C and D constants that
 *		correspond to an infinite plane at that surface.
 *	Parameters: Vertices that make up the vectors v1 - v2 and v3 - v2
 *	Return:		Plane constants which correspond to the vector surface
 *	Effect:		None.
 */
Plane calculatePlane(Coord vt1, Coord vt2, Coord vt3);

/* *******************************************************************
 *	calculateVertexNormal -  Calculates the normal for an object vertex
 *		- Given an object and the desired vertex, this function looks
 *		at all faces on the object and finds the average direction of the
 *		surface normals belonging to faces which include this vertex. The
 *		result is an accurate vertex normal.
 *		- Because texture mapping can be difficult with complex objects,
 *		some vertices may have been repeated. The supplied duplicates array
 *		specifies, in GLubyte pairs (0 = isDuplicate?, 1 = indexOfOriginal),
 *		which vertices are duplicates and which indices are the equivalent.
 *		Without this array, some surface normals are excluded when they
 *		shouldn't be.
 *	Parameters:	The object with all vertices/face indices/surface normals,
 *				the index of the vertex to be calculated, and an array
 *				specifying which vertices are duplicates.
 *	Return:		The vertex normal
 *	Effect:		None.
 */
Coord calculateVertexNormal(Object *obj, GLubyte vtIndex, GLubyte *duplicates);

/*
void createQuaternion(Quaternion *q, GLfloat angle, GLfloat x, GLfloat y, GLfloat z);
void createMatrixFromQuaternion(GLMatrix *m, Quaternion q);
*/

//#endif
