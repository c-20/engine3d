//#include <math.h>

//#include "GenLib.h"

inline GLfloat dot(Coord &c1, Coord &c2) {
    return c1.x * c2.x + c1.y * c2.y + c1.z * c2.z;
}

inline GLfloat mag(Coord &c) {
    return sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
}

inline Coord dir(Coord &c) {
	GLfloat m = mag(c);
	if (m == 0.0f) {
		return c;
	} else {
		return c / m;
	}
}

void translateObject(Object *o, Coord t) {
	for (int i = 1; i < o->nv_n; i += 2) {
		o->nv[i] += t;
	}
	for (int i = 0; i < o->fc_n; i++) {
		o->sn[i] += t;
	}
	o->cp += t;
}

void translate(Coord *vt, int vt_n, GLfloat x, GLfloat y, GLfloat z) {
	int inc = 1;
	if (vt_n < 0) {	// Normals in this array... only skew every second element
		vt_n = -vt_n;
		inc = 2;
	}
	for (int i = (inc == 2)?1:0; i < vt_n; i += inc) {
		vt[i].x += x;
		vt[i].y += y;
		vt[i].z += z;
	}
}


void rotateObject(Object *o, GLfloat angle_xy, GLfloat angle_xz, GLfloat angle_yz) {
	for (int i = 0; i < o->nv_n; i += 2) {
		rotate(&o->nv[i], 1, NULL, angle_xy, angle_xz, angle_yz);
	}
	for (int i = 1; i < o->nv_n; i += 2) {
		rotate(&o->nv[i], 1, &o->cp, angle_xy, angle_xz, angle_yz);
	}
	rotate(o->sn, o->fc_n, &o->cp/*NULL*/, angle_xy, angle_xz, angle_yz);
}


void rotate(Coord *vt, int vt_n, Coord *zero, GLfloat angle_xy, GLfloat angle_xz, GLfloat angle_yz) {
	Coord from, to;
	for (int i = 0; i < vt_n; i++) {
		// Make vertices relative to origin and grab values
		from = vt[i];
		if (zero) { from -= (*zero); }
		// Perform the xy rotate
		to.x = cos(angle_xy) * from.x - sin(angle_xy) * from.y;
		to.y = sin(angle_xy) * from.x + cos(angle_xy) * from.y;
		// Update the values for further rotations
		from.x = to.x; from.y = to.y;
		// Perform the xz rotate
		to.x = cos(angle_xz) * from.x - sin(angle_xz) * from.z;
		to.z = sin(angle_xz) * from.x + cos(angle_xz) * from.z;
		// Update value for further rotation
		from.z = to.z;
		// Perform the yz rotate
		to.y = cos(angle_yz) * from.y - sin(angle_yz) * from.z;
		to.z = sin(angle_yz) * from.y + cos(angle_yz) * from.z;
		// Move the vertex relative to the object's true centre
		if (zero) { to += (*zero); }
		vt[i] = to;
	}
}


void calculateRotationInits(Rotation *r) {
	GLfloat sin_angle_xy = sinf(-r->angles.z), cos_angle_xy = cosf(-r->angles.z);
	GLfloat sin_angle_xz = sinf(r->angles.y), cos_angle_xz = cosf(r->angles.y);
	GLfloat sin_angle_yz = sinf(r->angles.x), cos_angle_yz = cosf(r->angles.x);

	RotationInits *i = &r->inits;
	i->xx = cos_angle_xz * cos_angle_xy;
	i->xy = sin_angle_yz * sin_angle_xz * cos_angle_xy - cos_angle_yz * sin_angle_xy;
	i->xz = cos_angle_yz * sin_angle_xz * cos_angle_xy + sin_angle_yz * sin_angle_xy;
	i->yx = cos_angle_xz * sin_angle_xy;
	i->yy = sin_angle_yz * sin_angle_xz * sin_angle_xy + cos_angle_yz * cos_angle_xy;
	i->yz = cos_angle_yz * sin_angle_xz * sin_angle_xy - sin_angle_yz * cos_angle_xy;
	i->zx = -1 * sin_angle_xz;
	i->zy = sin_angle_yz * cos_angle_xz;
	i->zz = cos_angle_yz * cos_angle_xz;
	i->xx_xy = i->xx * i->xy;
	i->yx_yy = i->yx * i->yy;
	i->zx_zy = i->zx * i->zy;
}


void quickRotateObject(Object *o, RotationInits *r) {
	Coord cp = o->cp;
	// Rotate vertices..
	for (int i = 1; i < o->nv_n; i += 2) { o->nv[i] = quickRotate(o->nv[i] - cp, r) + cp; }
	// ..and their normals
	for (int i = 0; i < o->nv_n; i += 2) { o->nv[i] = quickRotate(o->nv[i], r); }
	// Rotate surface normals
	for (int i = 0; i < o->fc_n; i ++) { o->sn[i] = quickRotate(o->sn[i] - cp, r) + cp; }
}


Coord quickRotate(Coord vt, RotationInits *r) {
	Coord out;
	GLfloat x_y = vt.x * vt.y;
	out.x = (r->xx + vt.y) * (r->xy + vt.x) + (vt.z * r->xz) - (r->xx_xy + x_y);
	out.y = (r->yx + vt.y) * (r->yy + vt.x) + (vt.z * r->yz) - (r->yx_yy + x_y);
	out.z = (r->zx + vt.y) * (r->zy + vt.x) + (vt.z * r->zz) - (r->zx_zy + x_y);
	return out;
}


void scaleObject(Object *o, GLfloat scale) {
	for (int i = 1; i < o->nv_n; i += 2) {
		o->nv[i] = ((o->nv[i] - o->cp) * scale) + o->cp;
	}
	for (int i = 0; i < o->fc_n; i++) {
		o->sn[i] = ((o->sn[i] - o->cp) * scale) + o->cp;
	}
}


void scale(Coord *vt, int vt_n, Coord *zero, GLfloat scale) {
	for (int i = 0; i < vt_n; i++) {
		if (zero) { vt[i] -= (*zero); }
		vt[i] *= scale;
		if (zero) { vt[i] += (*zero); }
	}
}


void skewObject(Object *o, GLfloat skewx, GLfloat skewy, GLfloat skewz) {
	for (int i = 1; i < o->nv_n; i += 2) {
		o->nv[i].x = ((o->nv[i].x - o->cp.x) * skewx) + o->cp.x;
		o->nv[i].y = ((o->nv[i].y - o->cp.y) * skewy) + o->cp.y;
		o->nv[i].z = ((o->nv[i].z - o->cp.z) * skewz) + o->cp.z;
	}
	for (int i = 0; i < o->fc_n; i++) {
		o->sn[i].x = ((o->sn[i].x - o->cp.x) * skewx) + o->cp.x;
		o->sn[i].y = ((o->sn[i].y - o->cp.y) * skewy) + o->cp.y;
		o->sn[i].z = ((o->sn[i].z - o->cp.z) * skewz) + o->cp.z;
	}
}


void skew(Coord *vt, int vt_n, Coord *zero, GLfloat skewx, GLfloat skewy, GLfloat skewz) {
	int inc = 1;
	if (vt_n < 0) {	// Normals in this array... only skew every second element
		vt_n = -vt_n;
		inc = 2;
	}
	for (int i = (inc == 2)?1:0; i < vt_n; i += inc) {
		if (zero) { vt[i] -= (*zero); }
		vt[i].x *= skewx;
		vt[i].y *= skewy;
		vt[i].z *= skewz;
		if (zero) { vt[i] += (*zero); }
	}
}

Plane calculatePlane(Coord vt1, Coord vt2, Coord vt3) {
	Plane p;
	p.a = vt1.y * (vt2.z - vt3.z) + vt2.y * (vt3.z - vt1.z) + vt3.y * (vt1.z - vt2.z);
	p.b = vt1.z * (vt2.x - vt3.x) + vt2.z * (vt3.x - vt1.x) + vt3.z * (vt1.x - vt2.x);
	p.c = vt1.x * (vt2.y - vt3.y) + vt2.x * (vt3.y - vt1.y) + vt3.x * (vt1.y - vt2.y);
	p.d = 0.0f - vt1.x * (vt2.y * vt3.z - vt3.y * vt2.z)
		- vt2.x * (vt3.y * vt1.z - vt1.y * vt3.z) - vt3.x * (vt1.y * vt2.z - vt2.y * vt1.z);
	return p;
}

Coord calculateNormal(Coord vt1, Coord vt2, Coord vt3, Coord *zero) {
	// Get the normal directional vector
	Coord normal;
	normal.x = ((vt1.y - vt2.y) * (vt3.z - vt2.z)) - ((vt3.y - vt2.y) * (vt1.z - vt2.z));
	normal.y = ((vt3.x - vt2.x) * (vt1.z - vt2.z)) - ((vt1.x - vt2.x) * (vt3.z - vt2.z));
	normal.z = ((vt1.x - vt2.x) * (vt3.y - vt2.y)) - ((vt3.x - vt2.x) * (vt1.y - vt2.y));
	normal /= mag(normal);

	if (zero == NULL) {	// Pure directional vector... no surface plane
		return normal;
	} else {					// Relative vector that can be moved/rotated like any other
		Coord cp = *zero;
		Plane p = calculatePlane(vt1, vt2, vt3);

		// Find the point at which the ray crosses the plane
		GLfloat t = fabs((-1.0f * (p.a * cp.x + p.b * cp.y + p.c * cp.z + p.d)) /
			(p.a * normal.x + p.b * normal.y + p.c * normal.z));
		cp.x += normal.x * t;
		cp.y += normal.y * t;
		cp.z += normal.z * t;
		return cp;
	}
}

Coord calculateVertexNormal(Object *obj, GLubyte vtIndex, GLubyte *duplicates) {
	Coord avg = newCoord(0.0f, 0.0f, 0.0f);
	GLfloat n = 0.0f;
	// Mapping a duplicate vertex? Use the original instead.
	if (duplicates[vtIndex << 1] == TRUE) {
		vtIndex = duplicates[(vtIndex << 1) + 1];
	}
	// Look for this index in each of the faces
	for (int f = 0; f < obj->fc_n; f++) {
		for (int i = 0; i < obj->fc_vt_n[f]; i++) {
			GLubyte thisIndex = obj->fc_vtI[f][i];
			// Comparing against a duplicate? Compare against the original
			if (duplicates[thisIndex << 1] == TRUE) {
				thisIndex = duplicates[(thisIndex << 1) + 1];
			}
			// Does this vertex exist in the current face? Include its surface normal in the avg
			if (thisIndex == vtIndex) {
				avg += dir(obj->sn[f]/* - obj->ref*/);
				n += 1.0f;
				break;
			}
		}
	}
	avg /= n;				// Compute average
	return avg / mag(avg);	// Renormalise and return
}

/*
void createQuaternion(Quaternion *q, GLfloat angle, GLfloat x, GLfloat y, GLfloat z) {
	GLfloat t = (GLfloat)sin(angle / 2.0f);
	q->w = (GLfloat)cos(angle / 2.0f);
	q->x = x * t;
	q->y = y * y;
	q->z = z * t;
}


void createMatrixFromQuaternion(GLMatrix *m, Quaternion q) {
	// First row
	m->v[0] = 1.0f - 2.0f * (q.y * q.y + q.z * q.z);
	m->v[1] = 2.0f * (q.x * q.y + q.z * q.w);
	m->v[2] = 2.0f * (q.x * q.z - q.y * q.w);
	// Second row
	m->v[4] = 2.0f * (q.x * q.y - q.z * q.w);
	m->v[5] = 1.0f - 2.0f * (q.x * q.x + q.z * q.z);
	m->v[6] = 2.0f * (q.z * q.y + q.x * q.w);
	// Third row
	m->v[8] = 2.0f * ( q.x * q.z + q.y * q.w );
	m->v[9] = 2.0f * ( q.y * q.z - q.x * q.w );
	m->v[10] = 1.0f - 2.0f * (q.x * q.x + q.y * q.y);
	// Fill rest of the matrix
	m->v[3] = m->v[7] = m->v[11] = m->v[12] = m->v[13] = m->v[14] = 0.0f;
	m->v[15] = 1.0f;
}
*/
