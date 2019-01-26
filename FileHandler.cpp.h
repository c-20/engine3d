
//#include "FileHandler.h"


GLfloat getNextNumber(const char **arr, char delim) {
	GLfloat val = (GLfloat)atof(*arr);
	(*arr) = strchr(*arr, delim);		// *cArr will become NULL if no more values are available
	if ((*arr) != NULL) { (*arr) = &((*arr)[1]); }	// More values? Move past the comma
	return val;		// Return the current value
}


bool loadObjects(const char *filename, Object **objects, int *numObjects, char ***objectNames) {
	ifstream f(filename);
	string s, s2, s3;
	// Find out how many objects we have and declare the memory space
	getline(f, s);
	*numObjects = atoi(s.c_str());
	*objects = (Object *)malloc(sizeof(Object) * (*numObjects));
	*objectNames = (char **)malloc(sizeof(char *) * (*numObjects));
	int objIndex = 0;

	Object *obj = NULL;
	GLubyte temp[128], faceIndex = 0;
	ReadState state = Identity;
	// Read lines until we have as many objects as specified abive
	while (objIndex < (*numObjects) && getline(f, s)) {
		if (s.compare(0, 2, string("//")) != 0) {			// Not a comment? Yay data!
			const char *cArr = s.c_str();
			switch(state) {
				case Identity: {
					// This is a new object... create the structure
					obj = &((*objects)[objIndex]);
					// Find the comma
					int i = 0; while (cArr[i] != ',') { i++; }
					// Assign text left of comma to object name
					(*objectNames)[objIndex] = (char *)malloc(sizeof(char) * (i + 1));
					for (int j = 0; j <= i; j++) { (*objectNames)[objIndex][j] = cArr[j]; }
					(*objectNames)[objIndex][i] = '\0';
					// Assign text right of comma to number of faces
					obj->fc_n = atoi(&cArr[i + 1]);
					// Declare what we now know
					obj->fc_vt_n = (int *) malloc(sizeof(int) * obj->fc_n);
					obj->fc_edge_n = (int *) malloc(sizeof(int) * obj->fc_n);
					obj->fc_vtI = (GLubyte **) malloc(sizeof(GLubyte *) * obj->fc_n);
					obj->fc_edgeI = (GLubyte **) malloc(sizeof(GLubyte *) * obj->fc_n);
					//obj->fc_clr = (Colour **) malloc(sizeof(Colour *) * obj->fc_n);
					obj->vt_n = 0;
					obj->cp = newCoord(0.0f, 0.0f, 0.0f);
					faceIndex = 0;
					state++;
							   } break;
				case Face_Vertices: {
					// Find number of vertices in this face and test for highest vertex index
					int *i = &(obj->fc_vt_n[faceIndex] = 0);
					while (cArr != NULL) {	// While can see numbers
						temp[*i] = (GLubyte)getNextNumber(&cArr, ',');	// Get the point val
						if (temp[*i] > obj->vt_n) { obj->vt_n = temp[*i]; }
						(*i)++;
					}
					obj->fc_vtI[faceIndex] = (GLubyte *) malloc(sizeof(GLubyte) * (*i));
					memcpy(obj->fc_vtI[faceIndex], temp, sizeof(GLubyte) * (*i));
					state++;
									} break;
				case Face_Edges: {
					int *i = &(obj->fc_edge_n[faceIndex] = 0);
					while (cArr != NULL) {	// While can see numbers
						temp[(*i)++] = (GLubyte)getNextNumber(&cArr, ',');	// Get the point val
					}
					obj->fc_edgeI[faceIndex] = (GLubyte *) malloc(sizeof(GLubyte) * (*i));
					memcpy(obj->fc_edgeI[faceIndex], temp, sizeof(GLubyte) * (*i));
					if (++faceIndex >= obj->fc_n) { state++; }	// Move on to x vertices
					else { state = Face_Vertices; }			// Move on to the next face
								   } break;
				case Vertices: {
					obj->vt_n++;	// Currently set to max index... len == max + 1
					obj->nv_n = obj->vt_n << 1;
					obj->nv = (Coord *)malloc(sizeof(Coord) * obj->nv_n);
					obj->sn = (Coord *)malloc(sizeof(Coord) * obj->fc_n);
					obj->tx = (Point *)malloc(sizeof(Point) * obj->vt_n);

					GLubyte *duplicates = (GLubyte *)malloc(sizeof(GLubyte) * obj->nv_n);
					const char *xArr = cArr;
					getline(f, s2); const char *yArr = s2.c_str();
					getline(f, s3); const char *zArr = s3.c_str();

					for (int i = 1; i < obj->nv_n; i += 2) {
						if (xArr[0] == '*') {	// This is a duplicate
							duplicates[i - 1] = TRUE;
							duplicates[i] = (GLubyte)(int)getNextNumber(&yArr, ',');
							getNextNumber(&xArr, ',');	getNextNumber(&zArr, ',');	// Throw away the *s
							obj->nv[i] = obj->nv[(duplicates[i] << 1) + 1];
						} else {
							duplicates[i - 1] = duplicates[i] = FALSE;
							obj->nv[i].x = getNextNumber(&xArr, ',');
							obj->nv[i].y = getNextNumber(&yArr, ',');
							obj->nv[i].z = getNextNumber(&zArr, ',');
						}
					}
					// Calculate surface normals
					for (int f = 0; f < obj->fc_n; f++) {
						GLubyte *fcv = obj->fc_edgeI[f];
						obj->sn[f] = calculateNormal(obj->nv[(fcv[0] << 1) + 1], obj->nv[(fcv[1] << 1) + 1],
							obj->nv[(fcv[2] << 1) + 1], &obj->cp) * -1.0f;
					}
					// Calculate vertex normals
					for (int i = 0; i < obj->vt_n; i++) {
						obj->nv[i << 1] = calculateVertexNormal(obj, i, duplicates);
					}
					free(duplicates);
					state++;
								} break;

				case TexCoords: {
					const char *xArr = cArr;
					getline(f, s2); const char *yArr = s2.c_str();
					for (int i = 0; i < obj->vt_n; i++) {
						if (xArr == NULL) { MSGBOX("ERROR: Not enough texture coords"); return false; }
						obj->tx[i].x = getNextNumber(&xArr, ',');
						obj->tx[i].y = getNextNumber(&yArr, ',');
					}/*
					state++;
								 } break;
				case YTexCoords: {
					for (int i = 0; i < obj->vt_n; i++) {
						if (cArr == NULL) { MSGBOX("ERROR: Not enough Y texture coords"); return false; }
						obj->tx[i].y = getNextNumber(&cArr, ',');
					}*/
					state = Identity;
					objIndex++;
								 } break;
			}
		}
	}
	f.close();

	return true;
}

void freeObjects(Object *objects, int numObjects, char **objectNames) {
	for (int i = 0; i < numObjects; i++) { destroyObject(&objects[i], true); free(objectNames[i]); }
	free(objects);
	free(objectNames);
}

bool loadTerrain(const char *filename, Terrain *t, Coord quadrantSize, int extraSpace) {
	FILE *file = fopen(filename, "rb");

	// Read file header
	BitmapHeader header;
	fread(&header, sizeof(BitmapHeader), 1, file);
	if(header.bfType != 0x4D42) {
		MSGBOX("Terrain file not a bitmap!");
		fclose(file); return false;
	}
	t->width = (int)header.biWidth;
	t->depth = abs((int)header.biHeight);
	if (t->width <= 0 || t->depth <= 0) {
		MSGBOX("Invalid image size");
		fclose(file); return false;
	}

	// Read terrain height data
	fseek(file, header.bfOffBits, SEEK_SET);
	t->nv_n = (t->depth * t->width) << 1;
	t->nv = (Coord *)malloc(sizeof(Coord) * (t->nv_n + extraSpace * 2));
	int bytes = ((int)header.biBitCount) >> 3;	// number of bytes... because 2^3 == 8
	unsigned char *line = (unsigned char *)malloc(sizeof(unsigned char) * t->width * bytes);
	// Set the smaller axis to have points between 0 and 1
	GLfloat size = max(1.0f / (GLfloat)(t->width - 1), 1.0f / (GLfloat)(t->depth - 1));
	// Read each character and define a vertex
	int z, i;
	for (int zv = 0; zv < t->depth; zv++) {	// Because Windows is retarded and stores bitmaps upside-down
		if ((int)header.biHeight < 0) { z = t->depth - 1 - zv; } else { z = zv; }
		fread(line, 1, t->width * bytes, file);

		for (int x = 0; x < t->width; x++) {
			// Set vertices in such a way that all fit within (-.5,-.5,-.5) and (.5,.5,.5)
			i = ((z * t->width + x) << 1) + 1;
			t->nv[i].x = (size * (GLfloat)x) - 0.5f;
			t->nv[i].z = (size * (GLfloat)z) - 0.5f;
			t->nv[i].y = ((GLfloat)line[x * bytes + 1] - 128.0f) / -128.0f;	// Use second component... not sure if RGB/RGBA/ABGR/BGRA so assume second element is always good
		}
	}
	free(line);

	t->numBoxesX = (t->width - 1) / t->boxWidth;
	t->numBoxesZ = (t->depth - 1) / t->boxDepth;

	// Size terrain before calculating normals
	skew(t->nv, -t->nv_n, NULL, quadrantSize.x * t->numBoxesX, quadrantSize.y, quadrantSize.z * t->numBoxesZ);

	// Calculate normals for the vertices
	//t->nm = (Coord *)malloc(sizeof(Coord) * t->depth * t->width);
	Coord zero = newCoord(0.0f, 0.0f, 0.0f);
	int dblWidth = t->width << 1;
	for (int z = 0; z < t->depth; z++) {
		for (int x = 0; x < t->width; x++) {
			int i = (((z * t->width) + x) << 1) + 1;
			Coord avg = zero;

			if (x > 0 || z > 0) {	// Don't add this surface normal if at the top left corner
				avg += calculateNormal((z == 0)?zero:t->nv[i - dblWidth],
					t->nv[i], (x == 0)?zero:t->nv[i - 2], NULL);
			}
			if (x < t->width - 1 || z < t->depth - 1) {	// Don't add this one if at bottom right corner
				avg += calculateNormal((z == t->depth - 1)?zero:t->nv[i + dblWidth],
					t->nv[i], (x == t->width - 1)?zero:t->nv[i + 2], NULL);
			}
			Coord half = avg / 2.0f;
			t->nv[i - 1] = dir(half);
		}
	}

	// Close the file and finish
	fclose(file);
	return true;
}

void freeTerrain(Terrain *t) {
	free(t->nv);
	free(t->indices);
}

// Returns a pointer to the first string in the list, or NULL if there isn't one
// If f is NULL and first char in format is a *, use a given array rather than grabbing a line
char *getLineOfVars(ifstream *f, const char *format, ...) {
	char *ret = NULL;

	va_list varList;
	va_start(varList, format);

	// Get the line as a C array
	string s;
	const char *arr;
	if (f == NULL && format[0] == '*') {	// Array given!
		arr = va_arg(varList, const char *);
		format = &format[1];				// Ditch the *
	} else {								// Retrieve array
		getline(*f, s);
		arr = s.c_str();
	}
	// Find out how many arguments we have... if odd number, will use \0 as last delimeter
	int numChars; for (numChars = 0; format[numChars] != '\0'; numChars++);	// Find length

	// Start processing arguments!
	for (int i = 0; i < numChars; i += 2) {
		if (format[i] == '^') {
			const char **var = va_arg(varList, const char **);
			*var = arr;
		}
		if (format[i] == 'i') {
			int *var = va_arg(varList, int *);
			*var = (int)getNextNumber(&arr, format[i + 1]);
		}
		if (format[i] == 'C') {
			Coord *var = va_arg(varList, Coord *);
			(*var).x = getNextNumber(&arr, ',');
			(*var).y = getNextNumber(&arr, ',');
			(*var).z = getNextNumber(&arr, format[i + 1]);
		}
		if (format[i] == 'f') {
			float *var = va_arg(varList, float *);
			*var = getNextNumber(&arr, format[i + 1]);
		}
		if (format[i] == 's') {
			char *var = va_arg(varList, char *);
			int end; for (end = 0; arr[end] != format[i + 1]; end++); // Find end character
			memcpy(var, arr, sizeof(char) * end);
			var[end] = '\0';
			if (ret == NULL) { ret = var; }
			arr = &arr[end + 1];
		}
	}
	va_end(varList);

	return ret;
}


bool loadWorld(const char *filename, World *w, Object *objects, int numObjects, char **objectNames) {
	ifstream lvl(filename);
	string s;
	char temp[32];

	// Find out how many textures we have and declare the memory space
	getLineOfVars(&lvl, "i", &w->numTextures);
	w->textures = (GLuint *)malloc(sizeof(GLuint) * w->numTextures);
	for (int texIndex = 0; texIndex < w->numTextures; texIndex++) {
		loadTexture(getLineOfVars(&lvl, "s", &temp), &(w->textures[texIndex]));
	}

	// Define the total number of quadrants
	getLineOfVars(&lvl, "i,C", &w->numQuadrantLevels, &w->quadSize);
	w->numQuadrants = 0;
	int *quadrantsAtThisLevel = (int *)malloc(sizeof(int) * w->numQuadrantLevels);
	for (int l = 0; l < w->numQuadrantLevels; l++) {
		w->numQuadrants += (quadrantsAtThisLevel[l] = (int)pow(4.0f, (float)l));	// 1, 4, 16, 64
	}
	w->quadrants = (Quadrant *)malloc(sizeof(Quadrant) * w->numQuadrants);
	// Set the starting corner for the grids (in GL land)
	//GLfloat a = sqrtf(pow(4.0f, (float)w->numQuadrantLevels - 1.0f));	// VAR NAME!
	GLfloat quadsPerAxis = pow(2.0f, (float)w->numQuadrantLevels - 1.0f);
	Coord gridCorner;
	gridCorner.x = -((quadsPerAxis / 2.0f) * w->quadSize.x);
	gridCorner.y = -((w->quadSize.y * 5.0f) / 2.0f);
	gridCorner.z = -((quadsPerAxis / 2.0f) * w->quadSize.z);

	// Define each quadrant
	int start = 0;
	int quad = 1;
	for (int l = 0; l < w->numQuadrantLevels; l++) {
		int numBoxes = (int)sqrtf((GLfloat)quadrantsAtThisLevel[l]);
		GLfloat quadrantSpan = sqrtf((GLfloat)quadrantsAtThisLevel[w->numQuadrantLevels - l - 1]);
		for (int i = start; i < start + quadrantsAtThisLevel[l]; i++) {
			Quadrant& _quad = w->quadrants[i];
			// Define the two opposite quadrant corners (therefore 6 vertices)
			_quad.ltn.x = gridCorner.x + ((i - start) % numBoxes) * w->quadSize.x * quadrantSpan;
			_quad.ltn.y = gridCorner.y;
			_quad.ltn.z = gridCorner.z + ((i - start) / numBoxes) * w->quadSize.z * quadrantSpan;
			_quad.rbf.x = w->quadrants[i].ltn.x + (w->quadSize.x * quadrantSpan);
			_quad.rbf.y = w->quadrants[i].ltn.y + (w->quadSize.y * 5.0f);
			_quad.rbf.z = w->quadrants[i].ltn.z + (w->quadSize.z * quadrantSpan);
			_quad.centre = (_quad.rbf - _quad.ltn) / 2.0f;

			_quad.numBarriers = 0;//_quad.numEntities = 0;

			// Point to the appropriate children
			if ((_quad.hasChildren = (l != w->numQuadrantLevels - 1))) {
				_quad.frontLeft = &w->quadrants[quad++];	_quad.frontRight = &w->quadrants[quad++];
				_quad.backLeft = &w->quadrants[quad++];		_quad.backRight = &w->quadrants[quad++];
			}
		}
		start += quadrantsAtThisLevel[l];
	}

	// Load the terrain that will be inserted into quadrants
	getLineOfVars(&lvl, "i", &w->numTerrains);
	w->terrains = (Terrain *)malloc(sizeof(Terrain) * w->numTerrains);
	for (int i = 0; i < w->numTerrains; i++) {
		Terrain& _ter = w->terrains[i];

		int texId;
		float xRepeat, zRepeat;
		ifstream vec(getLineOfVars(&lvl, "s,i,f,f", &temp, &texId, &xRepeat, &zRepeat));
		_ter.textureId = w->textures[texId - 1];

		int extraSpace;
		Coord offset;
		getLineOfVars(&vec, "i;i,i;C;s ", &extraSpace, &_ter.boxWidth, &_ter.boxDepth, &offset, &temp);
		if (!loadTerrain(temp, &_ter, w->quadSize, extraSpace)) {
			MSGBOX("Failed to load terrain.");
			return false;
		}
		if ((_ter.width != _ter.boxWidth * (int)quadsPerAxis + 1) || (_ter.depth != _ter.boxDepth * (int)quadsPerAxis + 1)) {
			MSGBOX("Terrain size is incompatible with quadrant tree.");
			return false;
		}
		// Set up the texture coordinate array
		GLfloat xInc = xRepeat / (GLfloat)(_ter.width - 1);
		GLfloat zInc = zRepeat / (GLfloat)(_ter.depth - 1);
		_ter.tx = (Point *)malloc(sizeof(Point) * ((_ter.nv_n >> 1) + extraSpace));
		int index = 0;
		for (GLfloat z = 0.0f; z <= zRepeat; z += zInc) {
			for (GLfloat x = 0.0f; x <= xRepeat; x += xInc) {
				_ter.tx[index].x = x;
				_ter.tx[index].y = z;
				index++;
			}
		}
		for (int t = 0; t < extraSpace; t++) {
			_ter.tx[index + t].x = _ter.tx[index + t].y = 0.0f;
		}
		// Move terrain into position
		translate(_ter.nv, -_ter.nv_n, offset.x, offset.y, offset.z);

		// Set up indices array -- split into quadrants of triangle strips!
		int contentWidth = _ter.boxWidth + 1;
		_ter.contentSize = (contentWidth << 1) * _ter.boxDepth;
		int numIndices = _ter.contentSize * _ter.numBoxesX * _ter.numBoxesZ;
		_ter.indices = (GLushort *)malloc(sizeof(GLushort) * (numIndices + extraSpace));
		for (int z_ind = 0; z_ind < _ter.numBoxesZ; z_ind++) {
			for (int x_ind = 0; x_ind < _ter.numBoxesX; x_ind++) {
				int gridIndex = (z_ind * _ter.numBoxesX * _ter.contentSize) + (x_ind * _ter.contentSize);
				for (int z = 0; z < _ter.boxDepth; z++) {
					for (int x = 0; x < contentWidth; x++) {
						int index = ((z * contentWidth) + x) << 1;
						int value = (z_ind * _ter.width * _ter.boxDepth) + (x_ind * _ter.boxWidth) + (z * _ter.width) + x;
						_ter.indices[gridIndex + index] = value;
						_ter.indices[gridIndex + index + 1] = value + _ter.width;
					}
				}
			}
		}

		int wallGridIndex, numWalls;
		getLineOfVars(&vec, "i,i", &wallGridIndex, &numWalls);

		int start = w->numQuadrants - quadrantsAtThisLevel[w->numQuadrantLevels - 1];
		for (int b = start; b < w->numQuadrants; b++) {
			const Quadrant& _quad = w->quadrants[b];
			// THIS WILL CURRENTLY ONLY WORK FOR ONE TERRAIN FILE!
			w->quadrants[b].numBarriers = _ter.boxDepth + ((wallGridIndex == b - start) ? numWalls : 0);	// Remember to allocate extra for walls
			w->quadrants[b].barriers = (Barrier *)malloc(sizeof(Barrier) * _quad.numBarriers);
			int thisBox = ((b - start) % (_ter.numBoxesX * _ter.numBoxesZ)) * _ter.contentSize;
			int n = _ter.contentSize / _ter.boxDepth;
			for (int c = 0; c < _ter.boxDepth; c++) {								// Add the terrain pieces
				_quad.barriers[c].nv = _ter.nv;
				_quad.barriers[c].indices = &(_ter.indices[thisBox + n * c]);
				_quad.barriers[c].numIndices = (_ter.boxWidth + 1) << 1;
			}
			if (wallGridIndex == b - start) {	// This quadrant has walls! Add them
				for (int c = _ter.boxDepth; c < _quad.numBarriers; c++) {	// Walls
					getline(vec, s);
					//cArr = s.c_str();

					w->quadrants[b].barriers[c].nv = _ter.nv;

					// Set up some short variable names for readability
					Coord *_nv = w->quadrants[b].barriers[c].nv;
					const int& _nv_n = _ter.nv_n;
					const Coord& _ltn = w->quadrants[b].ltn;
					const Coord& _rbf = w->quadrants[b].rbf;

					// Add this set of indices onto the end of the existing array
					const char *iterator;
					Coord point;
					getLineOfVars(NULL, "*i;^", s.c_str(), &_quad.barriers[c].numIndices, &iterator);
					//_quad.barriers[c].numIndices = (int)getNextNumber(&cArr, ';');
					_quad.barriers[c].indices = &_ter.indices[numIndices];
					numIndices += _quad.barriers[c].numIndices;

					// Add each of the vertex points in this barrier
					for (int p = 0; p < w->quadrants[b].barriers[c].numIndices; p++) {
						_quad.barriers[c].indices[p] = _nv_n >> 1;
						_ter.nv_n += 2;
						getLineOfVars(NULL, "*C;^", iterator, &point, &iterator);
						_nv[_nv_n - 1] = _ltn + (_rbf - _ltn) * point;
						// Third vertex? Calculate the surface normal
						if (p == 2) {
							_nv[_nv_n - 2] = calculateNormal(_nv[_nv_n - 5], _nv[_nv_n - 3], _nv[_nv_n - 1], NULL);
							_nv[_nv_n - 6] = _nv[_nv_n - 4] = _nv[_nv_n - 2] *= -1.0f;	// Invert normal and copy (quick, dirty fix)
						}
						// > third? Copy the previous normal (assume a flat surface)
						if (p > 2) { _nv[_nv_n - 2] = _nv[_nv_n - 4]; }
					}
				}
				// Ok! Done adding walls... reset the index for the next wall-containing quadrant
				getLineOfVars(&vec, "i,i", &wallGridIndex, &numWalls);
			}
		}
		vec.close();
	}
	free(quadrantsAtThisLevel);

	// Set up the Physics Engine
//	PHY_Reset();
//	PHY_EnableGravity(true);
//	PHY_EnableDamping(true);
	// Add walls
	int numWalls;
	getLineOfVars(&lvl, "i", &numWalls);
	for (int i = 0; i < numWalls; i++) {
		float wallDist, wallXDir, wallYDir, wallZDir;
		getLineOfVars(&lvl, "f,f,f,f", &wallDist, &wallXDir, &wallYDir, &wallZDir);
//		PHY_AddWall(i, -wallDist, wallXDir, -wallZDir, wallYDir);
	}

	Entity *e = NULL;
	Coord c; int v; float f;

	getLineOfVars(&lvl, "i", &w->numEntities);
	w->entities = (Entity *)malloc(sizeof(Entity) * (w->numEntities));
	int entIndex = -1;

	// Read lines until we have as many objects as specified above
	#define IFEQUAL(str, arr)	if (str.compare(0, string(arr).length(), string(arr)) == 0)
	while (entIndex < w->numEntities && getline(lvl, s)) {
		if (s.compare(0, 2, string("//")) != 0) {			// Not a comment? Yay data!
			const char *cArr = s.c_str();
			cArr = &cArr[s.find(':') + 1];
			IFEQUAL(s, "NewObject") {
				e = &(w->entities[++entIndex]);
				createObjectFromTemplate(&e->obj, &objects[indexOfString(cArr, objectNames, numObjects)]);
				e->rotation.enabled = false;
				e->exists = true;
				e->Eid = 0; //PHY_NOPHYSICS;
				e->EnumSprings = 0;
				e->Edensity = 1.0f;
			} else IFEQUAL(s, "Physics")	{
				if (cArr[0] == 'Y' || cArr[0] == 'y') {
//					e->PHY_id = PHY_GetNumberRigidBodies();
//					PHY_AddBody(e->PHY_id, e->PHY_density, 2.0f, 2.0f, 2.0f);
				}
			} else IFEQUAL(s, "PHY_Density") {
				// Must set density before scaling!
//				if (e->PHY_id != PHY_NOPHYSICS) {
					getLineOfVars(NULL, "*f", cArr, &e->Edensity);
//					PHY_SetBody(e->PHY_id, e->PHY_density, 2.0f, 2.0f, 2.0f);
//				} else { MSGBOX("You must enable physics to give this object a density."); }
			} else IFEQUAL(s, "PHY_AddSpring") {
//				if (e->PHY_id != PHY_NOPHYSICS) {
//					e->springs[e->PHY_numSprings] = PHY_GetNumberOfSprings();
					getLineOfVars(NULL, "*i,C", cArr, &v, &c);
//					V3D v3d = V3D(c);
//					PHY_AddSpring(e->springs[e->PHY_numSprings], e->PHY_id, (unsigned int)v, v3d); // V3D(c));
//					e->PHY_numSprings++;
//				} else { MSGBOX("You must enable physics for this object to add a spring to it."); }
			} else IFEQUAL(s, "Texture")	{
				getLineOfVars(NULL, "*i", cArr, &v);
				e->textureId = w->textures[v - 1];
			} else IFEQUAL(s, "Translate")	{
				getLineOfVars(NULL, "*C", cArr, &c);
//				if (e->PHY_id != PHY_NOPHYSICS) {
//					PHY_SetBodyPosn(e->PHY_id, V3D(c).x, V3D(c).y, V3D(c).z);
//				} else {
                translateObject(&e->obj, c);
//				}
			} else IFEQUAL(s, "Rotate")		{
				getLineOfVars(NULL, "*C", cArr, &c);
//				if (e->PHY_id != PHY_NOPHYSICS) {
//					PHY_SetBodyDirn(e->PHY_id, RAD(c.x), RAD(-c.z), RAD(c.y));
//				} else {
                rotateObject(&e->obj, RAD(c.x), RAD(c.y), RAD(c.z));
//				}
			} else IFEQUAL(s, "Scale") {	// Must scale before translating (if using physics)!
				getLineOfVars(NULL, "*f", cArr, &f);
//				if (e->PHY_id != PHY_NOPHYSICS) {
//					PHY_SetBody(e->PHY_id, e->PHY_density, 2.0f * f, 2.0f * f, 2.0f * f);
//				}
				scaleObject(&e->obj, f);	// Need to scale the vertices regardless of whether or not physics is used
			} else IFEQUAL(s, "Skew") {		// Must skew before translating (if using physics)!
				float f2, f3;
				getLineOfVars(NULL, "*f,f,f", cArr, &f, &f2, &f3);
//				if (e->PHY_id != PHY_NOPHYSICS) {
//					PHY_SetBody(e->PHY_id, e->PHY_density, 2.0f * f, 2.0f * f2, 2.0f * f3);
//				}
				skewObject(&e->obj, f, f2, f3);
			} else IFEQUAL(s, "AnimRotate") {
				getLineOfVars(NULL, "*f,f,f", cArr, &e->rotation.angles.z, &e->rotation.angles.y, &e->rotation.angles.x);
				e->rotation.angles *= ONERAD;
				calculateRotationInits(&e->rotation);
				e->rotation.enabled = true;
			} else {
				MSGBOX("Invalid command.");
				return false;
			}
		}
	}

//	PHY_Run(); PHY_Run(); PHY_Run();

	#undef IFEQUAL
	lvl.close();
	return true;
}
#undef	GETVALUE

void freeWorld(World *w) {
	for (int i = 0; i < w->numEntities; i++) { destroyObject(&w->entities[i].obj, false); }
	free(w->entities);
	glDeleteTextures(w->numTextures, w->textures);
	free(w->textures);
	for (int i = 0; i < w->numTerrains; i++) { freeTerrain(&w->terrains[i]); }
	free(w->terrains);
	free(w->quadrants);
//	PHY_Reset();
}

int indexOfString(const char *str, char **list, int n) {
	for (int i = 0; i < n; i++) {
		int c = 0;
		while (str[c] != '\0') {
			if (str[c] != list[i][c]) { break; }
			c++;
		}
		if (str[c] == '\0' && list[i][c] == '\0') { return i; }
	}
	return -1;
}


void createObjectFromTemplate(Object *newObj, Object *obj) {
	newObj->fc_n = obj->fc_n;
	newObj->fc_vt_n = obj->fc_vt_n;
	newObj->fc_edge_n = obj->fc_edge_n;
	newObj->fc_vtI = obj->fc_vtI;
	newObj->fc_edgeI = obj->fc_edgeI;
	newObj->vt_n = obj->vt_n;
	newObj->nv_n = obj->nv_n;
	newObj->nv = (Coord *) malloc(sizeof(Coord) * obj->nv_n);
	memcpy(newObj->nv, obj->nv, sizeof(Coord) * obj->nv_n);
	newObj->sn = (Coord *) malloc(sizeof(Coord) * obj->fc_n);
	memcpy(newObj->sn, obj->sn, sizeof(Coord) * obj->fc_n);
	newObj->tx = obj->tx;
	newObj->cp = obj->cp;
}


/*********************************************************************
 *	destroyObject - Frees all of an object's acquired memory
 *		- Deallocates any memory allocated when the object was created
 *		- If an object is not a template, it will point to memory that
 *		belongs to one... therefore, do not free those resources.
 *	Parameters: Object pointer, isTemplate boolean
 *	Return:		None.
 *	Effect:		Deallocates all memory belonging to the object.
 */
void destroyObject(Object *obj, bool isTemplate) {
	if (isTemplate) {
		for (int i = 0; i < obj->fc_n; i++) {
			free(obj->fc_vtI[i]);
			free(obj->fc_edgeI[i]);
			//free(obj->fc_clr[i]);
		}
		free(obj->fc_vtI);
		free(obj->fc_edgeI);
		//free(obj->fc_clr);
		free(obj->fc_vt_n);
		free(obj->fc_edge_n);
		free(obj->tx);
	}
	free(obj->nv);
	free(obj->sn);
}


//Texture *
bool loadTexture(const char *filename, GLuint *texture) {
    Texture t;
	// Test if the file exists
//	if (FILE *file = fopen(filename, "r")) {	// File exists?
//		fclose(file);
//	} else {
//		MSGBOX("Texture file does not exist.");
//		return false;
//	}
//	RGBImage textureImage;
	// Load the texture image
//	AUX_RGBImageRec *textureImage;
//	if (!(textureImage = auxDIBImageLoad(filename))) {
//		MSGBOX("Texture image could not be read.");
//		return false;
//	}

    // -----------------------------------------------------------
    FILE *file = fopen(filename, "rb");
    if (!file) {
        MSGBOX("Texture image could not be read.");
        return false;
    }
	BitmapHeader header;
	fread(&header, sizeof(BitmapHeader), 1, file);
	if(header.bfType != 0x4D42) {
		MSGBOX("Texture file not a bitmap!");
		fclose(file); return false;
	}
	t.width = (int)header.biWidth;
	t.height = abs((int)header.biHeight);
	if (t.width <= 0 || t.height <= 0) {
		MSGBOX("Invalid texture image size");
		fclose(file); return false;
	}
	t.bytes = ((int)header.biBitCount) >> 3;	// number of bytes... because 2^3 == 8
	t.dlen = t.width * t.height * t.bytes;
    t.data = (unsigned char *)malloc(sizeof(unsigned char) * t.dlen);
	fseek(file, header.bfOffBits, SEEK_SET);
//	t->nv_n = (texheight * texwidth) << 1;
//	t->nv = (Coord *)malloc(sizeof(Coord) * (t->nv_n + extraSpace * 2));
//	unsigned char *line = (unsigned char *)malloc(sizeof(unsigned char) * texwidth * bytes);
	// Set the smaller axis to have points between 0 and 1
//	GLfloat size = max(1.0f / (GLfloat)(texwidth - 1), 1.0f / (GLfloat)(texheight - 1));
	// Read each character and define a vertex
	unsigned char *datap = t.data;
	int linelen = t.width * t.bytes;
	char upsidedown = ((int)header.biHeight < 0) ? 1 : 0;
	if (upsidedown) { datap = &t.data[t.dlen - linelen]; }
    int h = -1;
    while (++h < t.height) {
		fread(datap, 1, linelen, file);
		if (upsidedown) { datap -= linelen; }
        else { datap += linelen; }
	}
	GLenum format = (t.bytes == 3) ? GL_RGB : (t.bytes == 4) ? GL_RGBA : GL_RED;
	glGenTextures(1, texture);
	glBindTexture(GL_TEXTURE_2D, *texture);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, t.width, t.height, 0, format, GL_UNSIGNED_BYTE, t.data);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	free(t.data);
//    return t;
	return true;
}

