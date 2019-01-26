
//#include <windows.h>
//#include <gl\gl.h>
//#include <gl\glu.h>
//#include "glaux.h"
//#include <math.h>
//#include <stdarg.h>


//#include <fstream>
//#include <string>
//using namespace std;

//#include "GenLib.h"
//#include "physlib.h"



enum ReadState { Identity, Face_Vertices, Face_Edges, Vertices, TexCoords };
inline ReadState operator++(ReadState &rs, int) { return rs = (ReadState)(rs + 1); }
inline ReadState operator--(ReadState &rs, int) { return rs = (ReadState)(rs - 1); }


#pragma pack(push)
#pragma pack(1)
typedef struct {
	WORD bfType;
	DWORD bfSize;
	WORD bfReserved1;
	WORD bfReserved2;
	DWORD bfOffBits;
	DWORD biSize;
	LONG biWidth;
	LONG biHeight;
	WORD biPlanes;
	WORD biBitCount;
	DWORD biCompression;
	DWORD biSizeImage;
	LONG biXPelsPerMeter;
	LONG biYPelsPerMeter;
	DWORD biClrUsed;
	DWORD biClrImportant;
} BitmapHeader;
#pragma pack(pop)

typedef struct {
  int width, height;
  int bytes;
  int dlen;
  unsigned char *data;
} Texture;

/*   *******************************************************************
 *	getNextNumber - Gets the next number in a delimited string
 *		- Given a pointer to a pointer to the first character in an
 *		array of characters, this function gets the number existing
 *		at the start of the array and then moves the pointer to the
 *		character existing immediately after the next instance of
 *		the given delimiter. Consequently, each time this function
 *		is called, it will retrieve 'the next number'. If there are
 *		no more delimiters, the pointer will become NULL.
 *	Parameters: Pointer to a pointer to &arr[0].
 *	Return:		The retrieved number.
 *	Effect:		Moves the pointer to indexOf(delim) + 1.
 */
GLfloat getNextNumber(const char **arr, char delim);

char *getLineOfVars(ifstream *f, const char *format, ...);

/*   *******************************************************************
 *	LoadObjects - Loads a VJS file and creates template objects
 *		- Reads the contents of a VJS file and uses the information to
 *		create a series of object templates that can be used by the
 *		LVL (World) file loaded next.
 *	Parameters: Filename of the VJS file.
 *	Return:		True.
 *	Effect:		Populates the objects[] array.
 */
bool loadObjects(const char *filename, Object **objects, int *numObjects, char ***objectNames);

void freeObjects(Object *objects, int numObjects, char **objectNames);

/*   *******************************************************************
 *	LoadTerrain - Loads a bitmap heightmap
 *		- Loads a bitmap file and uses its data to populate the given
 *		terrain structure
 *	Parameters: Filename of bitmap file, pointer to terrain structure
 *	Return:		True if successful, false otherwise.
 *	Effect:		Populates the given terrain structure.
 */
bool loadTerrain(const char *filename, Terrain *t, Coord quadrantSize, int extraSpace);

void freeTerrain(Terrain *t);


/*   *******************************************************************
 *	LoadWorld - Creates objects from templates
 *		- Creates objects by duplicating templates (from the VJS) and
 *		puts them into their world space location
 *		- Defines properties such as rotation for each object
 *	Parameters: Filename (currently unused).
 *	Return:		True.
 *	Effect:		Populates the entities[] array.
 */
bool loadWorld(const char *filename, World *w, Object *objects, int numObjects, char **objectNames);

void freeWorld(World *w);

/* *******************************************************************
 *	indexOfString - Returns the list index corresponding to the given string
 *		- Given a list of strings and a string to search for, this function
 *		traverses the list and, if the string is found, returns the index of
 *		the matching string in the list.
 *	Parameters: String to search for, string list, number of items in list
 *	Return:		Returns the index of str in list, or -1 if not found
 *	Effect:		None.
 */
int indexOfString(const char *str, char **list, int n);


/* *******************************************************************
 *	createObjectFromTemplate - Duplicates a given template
 *		- Given an existing object obj, and an empty object newObj,
 *		this method analyses the existing object, allocates memory
 *		space of the same size and copies the contents across
 *		- This new object is independent of the template's vertices
 *		and surface normals, as these are instance-specific attributes.
 *	Parameters: Empty object pointer, pointer to existing object.
 *	Return:		None.
 *	Effect:		Populates the newObj object.
 */
void createObjectFromTemplate(Object *newObj, Object *obj);


/* *******************************************************************
 *	destroyObject - Frees all of an object's acquired memory
 *		- Deallocates any memory allocated when the object was created
 *		- If an object is not a template, it will point to memory that
 *		belongs to one... therefore, do not free those resources.
 *	Parameters: Object pointer, isTemplate boolean
 *	Return:		None.
 *	Effect:		Deallocates all memory belonging to the object.
 */
void destroyObject(Object *obj, bool isTemplate);


/*   *******************************************************************
 *	loadTexture - Puts a given texture into GL memory
 *		- Given a filename, this function opens the corresponding file
 *		as a bitmap, puts its data into GL memory, and puts the
 *		GL texture index into the given variable.
 *	Parameters: Texture filename, index variable
 *	Return:		True if loaded correctly, false otherwise.
 *	Effect:		None.
 */
bool loadTexture(const char *filename, GLuint *texture);

