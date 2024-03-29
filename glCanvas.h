// ====================================================================
// GLCanvas class by Rob Jagnow.
// ====================================================================

#ifndef _GL_CANVAS_H_
#define _GL_CANVAS_H_

#include <cstdlib>
#include <cassert>

class ArgParser;
class Camera;
class Cloth;
class Fluid;

// ====================================================================
// NOTE:  All the methods and variables of this class are static
// ====================================================================

class GLCanvas
{

public:

    // Set up the canvas and enter the rendering loop
    // Note that this function will not return but can be
    // terminated by calling 'exit(0)'
    static void initialize(ArgParser *_args);
    static void Render();

private:

    static void InitLight();
    static void Load();

    // various static variables
    static ArgParser *args;
    static Camera *camera;
    static Cloth *cloth;
    static Fluid *fluid;
    static int display_list_index;

    // state of the mouse cursor
    static int mouseButton;
    static int mouseX;
    static int mouseY;
    static bool controlPressed;

    // Callback functions for mouse and keyboard events
    static void display(void);
    static void reshape(int w, int h);
    static void mouse(int button, int state, int x, int y);
    static void motion(int x, int y);
    static void keyboard(unsigned char key, int x, int y);
    static void idle();
};

// ====================================================================

int HandleGLError();

// ====================================================================

#endif
