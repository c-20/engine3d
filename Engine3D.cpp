
// Standard includes
#include <windows.h>
#include <gl/gl.h>
#include <gl/glu.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;


#define MSGBOX(msg)		MessageBox(HWND_DESKTOP, msg, "Message", MB_OK | MB_ICONEXCLAMATION);

// Search for these libraries when linking
//#pragma comment( lib, "opengl32.lib" )
//#pragma comment( lib, "glu32.lib" )
//#pragma comment( lib, "Physics.lib" )

// CDS_FULLSCREEN is not defined by some compilers. Defining it this way avoids errors.
#ifndef CDS_FULLSCREEN
#define CDS_FULLSCREEN 4
#endif

// Contains functions to create the basic window and switch in and out of fullscreen mode
#include "NeHeGL.h"
// Contains general structures and vector functions
#include "GenLib.h"
#include "FileHandler.h"

#include "GenLib.cpp.h"
#include "FileHandler.cpp.h"
#include "NeHeGL.cpp.h"
// Window/keyboard input information
GL_Window*	g_window;
Keys*		g_keys;

#define F_LEFT		0
#define F_RIGHT		1
#define F_NEAR		2
#define F_FAR		3
#define F_TOP		4
#define F_BOTTOM	5
typedef struct Camera {
	Coord pos, rot;
	GLfloat pan, tilt, roll;
	GLfloat dolly, boom, track;
	Plane frustum[6];
} Camera;
Camera cam;

typedef struct {
    bool useAntialiasing, wireframeMode, useLighting, usePhysics;
} Settings;
Settings settings;

// Game variables
Object *objects;
char **objectNames;
int numObjects;

World world;

int terr_debug = 0;
int debug_level = -1;
bool debug_drawNormals = false;

void updateSetting(bool *setting) {
	if (setting == &settings.useAntialiasing) {
		if (*setting) {
			glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
			glEnable(GL_LINE_SMOOTH);
			OutputDebugString("Enable antialiasing\n");
		} else { glDisable (GL_LINE_SMOOTH); }
	} else if (setting == &settings.wireframeMode) {
		if (*setting) { glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); }
		else { glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); }
	} else if (setting == &settings.useLighting) {
		if (*setting) { glEnable(GL_LIGHTING); }
		else { glDisable(GL_LIGHTING); }
	} else if (setting == &settings.usePhysics) {
		/* if (*setting) { } else { } */
	}
}

/* Called when the program starts */
BOOL init(GL_Window* window, Keys* keys) {
	g_window	= window;
	g_keys		= keys;

	// Set up OpenGL
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glClearColor (0.227f, 0.431f, 0.647f, 0.0f);				// Windows 2000 blue
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_CULL_FACE);		glCullFace(GL_BACK);
	glEnable(GL_DEPTH_TEST);	glDepthFunc(GL_LEQUAL);		glClearDepth(FRUSTUM_ZFAR);
	glShadeModel(GL_SMOOTH);

	// Enable lighting
    glLightfv(GL_LIGHT0, GL_AMBIENT, float4(0.1f, 0.1f, 0.1f, 1.0f).v);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, float4(0.5f, 0.5f, 0.5f, 1.0f).v);
    glLightfv(GL_LIGHT0, GL_POSITION, float4(0.0f, 0.0f, 0.0f, 1.0f).v);
	glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 5.0f);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, 0.0f);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	// Make everything shiny
    glMaterialfv(GL_FRONT, GL_DIFFUSE, float4(0.4f, 0.4f, 0.4f, 1.0f).v);
    glMaterialfv(GL_FRONT, GL_SPECULAR, float4(0.8f, 0.8f, 0.8f, 1.0f).v);
    glMaterialf(GL_FRONT, GL_SHININESS, 20.0f);
	glEnable(GL_COLOR_MATERIAL);

	// Specify that we'll be drawing with vertex arrays
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);

	// Make sure we're on a clean modelview matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Load object templates
	if (!loadObjects("objects.vjs", &objects, &numObjects, &objectNames)) {
		MSGBOX("Failed to load objects file.");
		return FALSE;
	}
	// Load the file which defines terrain and object placement in the world
	if (!loadWorld("world.lvl", &world, objects, numObjects, objectNames)) {
		MSGBOX("Failed to load world file.");
		return FALSE;
	}

	// Set camera
	cam.pos = newCoord(0.0f, 20.0f, -100.0f);
	cam.rot = newCoord(0.0f, 0.0f, 0.0f);
	cam.pan = cam.tilt = cam.roll = 0.0f;
	cam.dolly = cam.track = cam.boom = 0.0f;

	// Define view settings (these can be updated live)
	updateSetting(&(settings.useAntialiasing = true));
	updateSetting(&(settings.wireframeMode = false));
	updateSetting(&(settings.useLighting = true));
	updateSetting(&(settings.usePhysics = false));

	return TRUE;
}

void dealloc() {
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	// Free World Objects
	freeWorld(&world);
	// Free Template Objects
	freeObjects(objects, numObjects, objectNames);
}

/* Called by the window's main loop */
bool update(DWORD milliseconds) {
	if (milliseconds < 15) { return false; }
	if (g_keys->keyDown[VK_ESCAPE]) { TerminateApplication (g_window); }
	if (g_keys->keyDown[VK_F1]) { ToggleFullscreen (g_window); }

	if (g_keys->keyDown['M'] || g_keys->keyDown['m']) {
		updateSetting(&(settings.wireframeMode = !settings.wireframeMode));
		g_keys->keyDown['M'] = g_keys->keyDown['m'] = FALSE;
	}

	if (g_keys->keyDown['N'] || g_keys->keyDown['n']) {
		debug_drawNormals = !debug_drawNormals;
		g_keys->keyDown['N'] = g_keys->keyDown['n'] = FALSE;
	}
	if (g_keys->keyDown['L'] || g_keys->keyDown['l']) {
		debug_level = ((debug_level + 2) % (world.numQuadrantLevels + 1) - 1);
		g_keys->keyDown['L'] = g_keys->keyDown['l'] = FALSE;
	}

	if (g_keys->keyDown[' ']) {
		updateSetting(&(settings.usePhysics = !settings.usePhysics));
		g_keys->keyDown[' '] = FALSE;
	}

	// Camera Look
	if (g_keys->keyDown[VK_UP])		{ cam.tilt -= 1.0f; }
	if (g_keys->keyDown[VK_DOWN])	{ cam.tilt += 1.0f; }
	if (g_keys->keyDown[VK_LEFT])	{ cam.pan -= 1.0f; }
	if (g_keys->keyDown[VK_RIGHT])	{ cam.pan += 1.0f; }
	// Camera Move
	if (g_keys->keyDown['W'] || g_keys->keyDown['w']) { cam.dolly += 0.4f; }
	if (g_keys->keyDown['S'] || g_keys->keyDown['s']) { cam.dolly -= 0.4f; }
	if (g_keys->keyDown['A'] || g_keys->keyDown['a']) { cam.track -= 0.4f; }
	if (g_keys->keyDown['D'] || g_keys->keyDown['d']) { cam.track += 0.4f; }
	if (g_keys->keyDown['R'] || g_keys->keyDown['r']) { cam.boom += 0.4f; }
	if (g_keys->keyDown['F'] || g_keys->keyDown['f']) { cam.boom -= 0.4f; }

	for (int i = 0; i < world.numEntities; i++) {
		if (world.entities[i].rotation.enabled) {
			quickRotateObject(&world.entities[i].obj, &world.entities[i].rotation.inits);
		}
	}
//	if (settings.usePhysics) { PHY_Run(); PHY_Run(); PHY_Run(); }
	//if (modelRotate) { modelAngle += (float) (milliseconds) / 10.0f; }
	return true;
}

/* Moves the camera to its defined position and recalculates the frustum planes */
void setCameraPerspective() {
	Matrix m;

	// Create the equivalent of a yaw * pitch quaternion matrix and...
	// ...retrieve the i and k components of the direction vector
	glLoadIdentity();
	glRotatef(cam.pan, 0.0f, 1.0f, 0.0f);
	glRotatef(cam.tilt, 1.0f, 0.0f, 0.0f);
	glGetFloatv(GL_MODELVIEW_MATRIX, m.v);
	cam.rot.x = m.v[8];
	cam.rot.z = m.v[10];

	// Create the equivalent of a pitch-only quaternion matrix and...
	// ...retrieve the j component of the direction vector
	glLoadIdentity();
	glRotatef(cam.tilt, 1.0f, 0.0f, 0.0f);
	glGetFloatv(GL_MODELVIEW_MATRIX, m.v);
	cam.rot.y = m.v[9];

	glRotatef(cam.pan, 0.0f, 1.0f, 0.0f);

	// Move camera forward/back to new position (minus a bit for frustum culling)
	cam.pos += (cam.rot * (cam.dolly - 1.0f));
	// Add strafing/tracking (rotate vector by 90deg and disregard y component)
	Coord track = cam.rot;
	track.y = 0.0f;
	track = dir(track);
	cam.pos.x += (track.z * cam.track);
	cam.pos.z += -(track.x * cam.track);
	// Add elevating/booming
	cam.pos.y += cam.boom;
	// Tell GL about our fancy new location
	glTranslatef(-cam.pos.x, -cam.pos.y, cam.pos.z);

	cam.dolly = cam.track = cam.boom = 0.0f;

	// Get the frustum planes
	Matrix modelview, projection;
	glGetFloatv(GL_MODELVIEW_MATRIX, modelview.v);
	glGetFloatv(GL_PROJECTION_MATRIX, projection.v);

	// Finish the move
	cam.pos += (cam.rot * 1.0f);
	glLoadIdentity();
	glRotatef(cam.tilt, 1.0f, 0.0f, 0.0f);
	glRotatef(cam.pan, 0.0f, 1.0f, 0.0f);
	glTranslatef(-cam.pos.x, -cam.pos.y, cam.pos.z);

	// Create frustum planes
	m = modelview * projection;
	// Create the RIGHT plane and normalise it
	cam.frustum[F_RIGHT] = Plane(m.v[3] - m.v[0], m.v[7] - m.v[4], m.v[11] - m.v[8], m.v[15] - m.v[12]);
	cam.frustum[F_RIGHT] /= sqrtf(cam.frustum[F_RIGHT].a * cam.frustum[F_RIGHT].a + cam.frustum[F_RIGHT].b * cam.frustum[F_RIGHT].b + cam.frustum[F_RIGHT].c * cam.frustum[F_RIGHT].c);
	// Create the LEFT plane and normalise it
	cam.frustum[F_LEFT] = Plane(m.v[3] + m.v[0], m.v[7] + m.v[4], m.v[11] + m.v[8], m.v[15] + m.v[12]);
	cam.frustum[F_LEFT] /= sqrtf(cam.frustum[F_LEFT].a * cam.frustum[F_LEFT].a + cam.frustum[F_LEFT].b * cam.frustum[F_LEFT].b + cam.frustum[F_LEFT].c * cam.frustum[F_LEFT].c);
	// Create the BOTTOM plane and normalise it
	cam.frustum[F_LEFT] = Plane(m.v[3] + m.v[1], m.v[7] + m.v[5], m.v[11] + m.v[9], m.v[15] + m.v[13]);
	cam.frustum[F_BOTTOM] /= sqrtf(cam.frustum[F_BOTTOM].a * cam.frustum[F_BOTTOM].a + cam.frustum[F_BOTTOM].b * cam.frustum[F_BOTTOM].b + cam.frustum[F_BOTTOM].c * cam.frustum[F_BOTTOM].c);
	// Create the TOP plane and normalise it
	cam.frustum[F_LEFT] = Plane(m.v[3] - m.v[1], m.v[7] - m.v[5], m.v[11] - m.v[9], m.v[15] - m.v[13]);
	cam.frustum[F_TOP] /= sqrtf(cam.frustum[F_TOP].a * cam.frustum[F_TOP].a + cam.frustum[F_TOP].b * cam.frustum[F_TOP].b + cam.frustum[F_TOP].c * cam.frustum[F_TOP].c);
	// Create the FAR plane and normalise it
	cam.frustum[F_LEFT] = Plane(m.v[3] - m.v[2], m.v[7] - m.v[6], m.v[11] - m.v[10], m.v[15] - m.v[14]);
	cam.frustum[F_FAR] /= sqrtf(cam.frustum[F_FAR].a * cam.frustum[F_FAR].a + cam.frustum[F_FAR].b * cam.frustum[F_FAR].b + cam.frustum[F_FAR].c * cam.frustum[F_FAR].c);
	// Create the NEAR plane and normalise it
	cam.frustum[F_LEFT] = Plane(m.v[3] + m.v[2], m.v[7] + m.v[6], m.v[11] + m.v[10], m.v[15] + m.v[14]);
	cam.frustum[F_NEAR] /= sqrtf(cam.frustum[F_NEAR].a * cam.frustum[F_NEAR].a + cam.frustum[F_NEAR].b * cam.frustum[F_NEAR].b + cam.frustum[F_NEAR].c * cam.frustum[F_NEAR].c);
}

/* Returns true if the given point is within the current frustum planes */
bool pointWithinFrustum(Coord p) {
	for (int i = 0; i < 6; i++) {	// Iterate through each plane
		if (cam.frustum[i].a * p.x + cam.frustum[i].a * p.y + cam.frustum[i].c * p.z + cam.frustum[i].d <= 0) {
			return false;
		}
	}
	return true;
}

/*///	Returns 1 if the given square on the xz-plane is partially within the current fristum planes, 2 if it
/*///	is completely within, and 0 if it is not.
int xzSquareInFrustum(Coord c1, Coord c2) {
	int verticesInPlane, inPlanes = 0;
	for (int i = 0; i < 6; i++) {
		if (i == F_TOP || i == F_BOTTOM) { continue; }	// Skip Y planes
		verticesInPlane = 0;
		if (cam.frustum[i].a * c1.x + cam.frustum[i].b * cam.pos.y + cam.frustum[i].c * c1.z + cam.frustum[i].d > 0) { verticesInPlane++; }
		if (cam.frustum[i].a * c2.x + cam.frustum[i].b * cam.pos.y + cam.frustum[i].c * c1.z + cam.frustum[i].d > 0) { verticesInPlane++; }
		if (cam.frustum[i].a * c1.x + cam.frustum[i].b * cam.pos.y + cam.frustum[i].c * c2.z + cam.frustum[i].d > 0) { verticesInPlane++; }
		if (cam.frustum[i].a * c2.z + cam.frustum[i].b * cam.pos.y + cam.frustum[i].c * c2.z + cam.frustum[i].d > 0) { verticesInPlane++; }
		if (!verticesInPlane) { return 0; }
		if (verticesInPlane == 4) { inPlanes++; }
	}
	return (inPlanes == 4)? 2 : 1;
}


void drawQuadrant(Quadrant *q) {
	//if (!xzSquareInFrustum(q->ltn, q->rbf)) { return; }
	for (int i = 0; i < q->numBarriers; i++) {
		glNormalPointer(GL_FLOAT, sizeof(Coord) * 2, &q->barriers[i].nv[0]);
		glVertexPointer(3, GL_FLOAT, sizeof(Coord) * 2, &q->barriers[i].nv[1]);
		glDrawElements(GL_TRIANGLE_STRIP, q->barriers[i].numIndices, GL_UNSIGNED_SHORT, q->barriers[i].indices);
	}
	if (q->hasChildren) {
		drawQuadrant(q->frontLeft);		drawQuadrant(q->frontRight);
		drawQuadrant(q->backLeft);		drawQuadrant(q->backRight);
	}
}

void debug_drawTerrainNormals(); void debug_drawPolygonNormals(); void debug_drawQuadrants();
void drawScene() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear buffers

	// Call the debug functions if necessary
	debug_drawQuadrants(); if (debug_drawNormals) { debug_drawTerrainNormals(); debug_drawPolygonNormals(); }

	// Draw the terrain
	glColor3f(1.0f, 1.0f, 1.0f);
	glBindTexture(GL_TEXTURE_2D, world.terrains[0].textureId);
	glTexCoordPointer(2, GL_FLOAT, 0, world.terrains[0].tx);
	drawQuadrant(&world.quadrants[0]);

	// Draw the objects
	glColor3f(1.0f, 1.0f, 1.0f);
	GLuint currentTexture = 9999;
	Point *currentTexCoords = NULL;
//	Matrix physicsTransform;
	for (int e = 0; e < world.numEntities; e++) {
		const Entity& _ent = world.entities[e];
		if (currentTexture != _ent.textureId) {
			glBindTexture(GL_TEXTURE_2D, (currentTexture = _ent.textureId));
		}
		if (currentTexCoords != _ent.obj.tx) {
			glTexCoordPointer(2, GL_FLOAT, 0, (currentTexCoords = _ent.obj.tx));
		}
//		if (_ent.PHY_id != PHY_NOPHYSICS) {
//			glPushMatrix();
//			glRotatef(-90.0f,1.0f,0.0f,0.0f);				// Someone thought it would be clever to have Z-axis gravity
//			for (int i = 0; i < _ent.PHY_numSprings; i++) {
//				glBegin(GL_LINES);
//					glVertex3fv(PHY_GetSpringVertex3f(_ent.springs[i], 0));
//					glVertex3fv(PHY_GetSpringVertex3f(_ent.springs[i], 1));
//				glEnd();
//			}
//			PHY_GetBodyTransform(_ent.PHY_id, physicsTransform.v);
//
//			// Quick and dirty fix for a bug in the physics engine that allows objects to slide through the floor when motionless
//			// This does NOT fix the slight angular droop! That can't really be fixed without the PhysLib code
//			if (physicsTransform.v[14] < 10.5f) { physicsTransform.v[14] = 10.5f; }
//
//			glMultMatrixf(physicsTransform.v);
//		}
		glNormalPointer(GL_FLOAT, sizeof(Coord) * 2, &_ent.obj.nv[0]);
		glVertexPointer(3, GL_FLOAT, sizeof(Coord) * 2, &_ent.obj.nv[1]);
		for (int i = 0; i < _ent.obj.fc_n; i++) {
			glDrawElements(GL_TRIANGLE_STRIP, _ent.obj.fc_vt_n[i],
								GL_UNSIGNED_BYTE, _ent.obj.fc_vtI[i]);
		}
//		if (_ent.Eid != PHY_NOPHYSICS) {
//			glPopMatrix();
//		}
	}

	setCameraPerspective();	// Since frustum boundaries are being recalced.. move higher^^

	//exit(0);
	SwapBuffers(g_window->hDC);			// Swap buffers (double buffering)
}


























//////DEBUG//////



void debug_drawTerrainNormals() {
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisable(GL_LIGHTING);
	glColor3f(0.0f, 0.0f, 1.0f);
	glBegin(GL_LINES);
	for (int i = 0; i < world.terrains[0].width * world.terrains[0].depth; i++) {
		glVertex3f(world.terrains[0].nv[(i << 1) + 1].x, world.terrains[0].nv[(i << 1) + 1].y, world.terrains[0].nv[(i << 1) + 1].z);
		Coord temp = dir(world.terrains[0].nv[i << 1]) * 1.0f;
		glVertex3f(world.terrains[0].nv[(i << 1) + 1].x + temp.x, world.terrains[0].nv[(i << 1) + 1].y + temp.y, world.terrains[0].nv[(i << 1) + 1].z + temp.z);
	}
	glEnd();
	glEnable(GL_LIGHTING);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
}


void debug_drawPolygonNormals() {
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisable(GL_LIGHTING);
	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_LINES);
	for (int e = 0; e < world.numEntities; e++) {
		// Ignore objects controlled by the physics engine
//		if (world.entities[e].PHY_id != PHY_NOPHYSICS) { continue; }

		const Object& _obj = world.entities[e].obj;
		// Draw vertex normals
		for (int i = 1; i < _obj.nv_n; i += 2) {
			glVertex3f(_obj.nv[i].x, _obj.nv[i].y, _obj.nv[i].z);
			Coord temp = dir(_obj.nv[i - 1]) * 1.0f;
			glVertex3f(_obj.nv[i].x + temp.x, _obj.nv[i].y + temp.y, _obj.nv[i].z + temp.z);
		}
		// Draw surface normals
		for (int i = 0; i < _obj.fc_n; i++) {
			glVertex3f(_obj.sn[i].x, _obj.sn[i].y, _obj.sn[i].z);
			Coord diff = _obj.sn[i] - _obj.cp;
            Coord temp = dir(diff) * 1.0f;
//            Coord temp = dir(_obj.sn[i] - _obj.cp) * 1.0f;
			glVertex3f(_obj.sn[i].x + temp.x, _obj.sn[i].y + temp.y, _obj.sn[i].z + temp.z);
		}
	}
	glEnd();
	glEnable(GL_LIGHTING);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
}



void debug_drawQuadrants() {
	if (debug_level != -1) {
		glColor3f(1.0f, 1.0f, 0.0f);
		glBindTexture(GL_TEXTURE_2D, 0);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisable(GL_CULL_FACE);
		glDisable(GL_LIGHTING);
		Coord lines[8];
		glVertexPointer(3, GL_FLOAT, 0, lines);
		for (int i = (debug_level < 2)?debug_level:(debug_level == 2)?5:21;
			i < ((debug_level == 0)?1:(debug_level == 1)?5:(debug_level == 2)?21:world.numQuadrants); i++) {

			lines[0] = lines[1] = lines[2] = lines[3] = lines[4] = world.quadrants[i].ltn;
			lines[1].x = lines[2].x = world.quadrants[i].rbf.x;
			lines[2].z = lines[3].z = world.quadrants[i].rbf.z;
			glDrawArrays(GL_LINE_STRIP, 0, 5);

			lines[0].y = lines[1].y = lines[2].y = lines[3].y = lines[4].y = world.quadrants[i].rbf.y;
			glDrawArrays(GL_LINE_STRIP, 0, 5);

			lines[0].x = lines[1].x = lines[4].x = lines[5].x = world.quadrants[i].ltn.x;
			lines[2].x = lines[3].x = lines[6].x = lines[7].x = world.quadrants[i].rbf.x;
			lines[0].y = lines[2].y = lines[4].y = lines[6].y = world.quadrants[i].ltn.y;
			lines[1].y = lines[3].y = lines[5].y = lines[7].y = world.quadrants[i].rbf.y;
			lines[0].z = lines[1].z = lines[2].z = lines[3].z = world.quadrants[i].ltn.z;
			lines[4].z = lines[5].z = lines[6].z = lines[7].z = world.quadrants[i].rbf.z;
			glDrawArrays(GL_LINES, 0, 8);

		}
		glEnable(GL_CULL_FACE);
		glEnable(GL_LIGHTING);
		glEnableClientState(GL_NORMAL_ARRAY);
	}
}
