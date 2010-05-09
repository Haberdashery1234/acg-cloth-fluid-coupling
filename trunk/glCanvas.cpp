#include "glCanvas.h"

// Included files for OpenGL Rendering
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "glCanvas.h"
#include "argparser.h"
#include "camera.h"
#include "cloth.h"
#include "fluid.h"
#include "matrix.h"

// ========================================================
// static variables of GLCanvas class

// State of the mouse cursor
int GLCanvas::mouseButton;
int GLCanvas::mouseX;
int GLCanvas::mouseY;
int GLCanvas::display_list_index;
ArgParser* GLCanvas::args;
Camera* GLCanvas::camera;
Cloth* GLCanvas::cloth;
Fluid* GLCanvas::fluid;
bool GLCanvas::controlPressed;

void GLCanvas::Load()
{
    delete cloth;
    cloth = NULL;
    delete fluid;
    fluid = NULL;
    if (args->fluid_file != "")
        fluid = new Fluid(args);
    if (args->cloth_file != "")
    {
        if (fluid)
          cloth = new Cloth(args, fluid);
        else
          cloth = new Cloth(args);
    }
}

// ========================================================
// Initialize all appropriate OpenGL variables, set
// callback functions, and start the main event loop.
// This function will not return but can be terminated
// by calling 'exit(0)'
// ========================================================

void GLCanvas::initialize(ArgParser *_args)
{
    args = _args;
    cloth = NULL;
    fluid = NULL;
    Load();
    Vec3f camera_position = Vec3f(0,0,4);
    Vec3f point_of_interest = Vec3f(0,0,0);
    Vec3f up = Vec3f(0,1,0);
    camera = new PerspectiveCamera(camera_position, point_of_interest, up, 20 * M_PI/180.0);

    // setup glut stuff
    glutInitWindowSize(args->width, args->height);
    glutInitWindowPosition(100,100);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
    glutCreateWindow("OpenGL Viewer");

    // basic rendering
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
    GLfloat ambient[] = { 0.2, 0.2, 0.2, 1.0 };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glCullFace(GL_BACK);
    glDisable(GL_CULL_FACE);

    display_list_index = glGenLists(1);

    // Initialize callback functions
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);

    Render();

    // Enter the main rendering loop
    glutMainLoop();
}

// ========================================================

void GLCanvas::InitLight()
{
    // Set the last component of the position to 0 to indicate
    // a directional light source

    GLfloat position[4] = { 30,30,100, 1};
    GLfloat diffuse[4] = { 0.75,0.75,0.75,1};
    GLfloat specular[4] = { 0,0,0,1};
    GLfloat ambient[4] = { 0.2, 0.2, 0.2, 1.0 };

    GLfloat zero[4] = {0,0,0,0};
    glLightfv(GL_LIGHT1, GL_POSITION, position);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
    glLightfv(GL_LIGHT1, GL_AMBIENT, zero);
    glEnable(GL_LIGHT1);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_COLOR_MATERIAL);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

    GLfloat spec_mat[4] = {1,1,1,1};
    float glexponent = 30;
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &glexponent);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, spec_mat);

    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    float back_color[] = { 0.5,0.5,1.0,1};
    glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, back_color);
    glEnable(GL_LIGHT1);
}


void GLCanvas::display(void)
{
    // Clear the display buffer, set it to the background color
    glClearColor(1,1,1,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set the camera parameters
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    InitLight(); // light will be a headlamp!
    camera->glPlaceCamera();

    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);

    glCallList(display_list_index);
    HandleGLError();

    // Swap the back buffer with the front buffer to display
    // the scene
    glutSwapBuffers();
}

// ========================================================
// Callback function for window resize
// ========================================================

void GLCanvas::reshape(int w, int h)
{
    // Set the OpenGL viewport to fill the entire window
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);

    // Set the camera parameters to reflect the changes
    camera->glInit(w, h);

    args->width = w;
    args->height = h;
}

// ========================================================
// Callback function for mouse click or release
// ========================================================

void GLCanvas::mouse(int button, int state, int x, int y)
{
    // Save the current state of the mouse.  This will be
    // used by the 'motion' function
    mouseButton = button;
    mouseX = x;
    mouseY = y;
    controlPressed = glutGetModifiers() & GLUT_ACTIVE_CTRL;
}

// ========================================================
// Callback function for mouse drag
// ========================================================

void GLCanvas::motion(int x, int y)
{
    // Left button = rotation
    // (rotate camera around the up and horizontal vectors)
    if (mouseButton == GLUT_LEFT_BUTTON)
    {
        camera->rotateCamera(0.005*(mouseX-x), 0.005*(mouseY-y));
        mouseX = x;
        mouseY = y;
    }
    // Middle button = translation
    // (move camera perpendicular to the direction vector)
    else if (mouseButton == GLUT_MIDDLE_BUTTON)
    {
        camera->truckCamera((mouseX-x)*0.5, (y-mouseY)*0.5);
        mouseX = x;
        mouseY = y;
    }
    // Right button = dolly or zoom
    // (move camera along the direction vector)
    else if (mouseButton == GLUT_RIGHT_BUTTON)
    {
        if (controlPressed)
        {
            camera->zoomCamera(mouseY-y);
        }
        else
        {
            camera->dollyCamera(mouseY-y);
        }
        mouseX = x;
        mouseY = y;
    }

    // Redraw the scene with the new camera parameters
    glutPostRedisplay();
}

// ========================================================
// Callback function for keyboard events
// ========================================================

void GLCanvas::keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 'a':
    case 'A':
        // toggle continuous animation
        args->animate = !args->animate;
        if (args->animate)
            printf ("animation started, press 'A' to stop\n");
        else
            printf ("animation stopped, press 'A' to start\n");
        break;
    case ' ':
        // a single step of animation
        if (cloth) cloth->Animate();
        if (fluid) fluid->Animate();
        Render();
        break;
    case 'm':
    case 'M':
        args->particles = !args->particles;
        Render();
        break;
    case 'v':
    case 'V':
        args->velocity = !args->velocity;
        Render();
        break;
    case 'f':
    case 'F':
        args->force = !args->force;
        Render();
        break;
    case 'e':
    case 'E':
        args->edge_velocity = (args->edge_velocity+1)%4;
        Render();
        break;
    case 'd':
    case 'D':
        args->dense_velocity = (args->dense_velocity+1)%4;
        Render();
        break;
    case 's':
    case 'S':
        args->surface = !args->surface;
        Render();
        break;
    case 'w':
    case 'W':
        args->wireframe = !args->wireframe;
        Render();
        break;
    case 'b':
    case 'B':
        args->bounding_box = !args->bounding_box;
        Render();
        break;
    case 'c':
    case 'C':
        args->cubes = !args->cubes;
        Render();
        break;
    case 'p':
    case 'P':
        args->pressure = !args->pressure;
        Render();
        break;
    case 'r':
    case 'R':
        // reset system
        Load();
        Render();
        break;
    case '+':
    case '=':
        args->isosurface += 0.1;
        if (args->isosurface >= 0.999) args->isosurface = 0.99;
        printf ("isosurface %f\n", args->isosurface);
        Render();
        break;
    case '-':
    case '_':
        args->isosurface -= 0.1;
        if (args->isosurface <= 0.001) args->isosurface = 0.01;
        printf ("isosurface %f\n", args->isosurface);
        Render();
        break;
    case 'q':
    case 'Q':
        delete cloth;
        cloth = NULL;
        delete fluid;
        fluid = NULL;
        delete args;
        args = NULL;
        delete camera;
        camera = NULL;
        printf ("program exiting\n");
        exit(0);
        break;
    case 'u':
    case 'U':
        args->updraft = !args->updraft;
        break;        
    default:
        printf("UNKNOWN KEYBOARD INPUT  '%c'\n", key);
    }
}


void GLCanvas::idle()
{
    if (args->animate)
    {
        // do 10 steps of animation before rendering
        for (int i = 0; i < 10; i++)
        {
            if (cloth) cloth->Animate();
            if (fluid) fluid->Animate();
        }
        Render();
    }
}



void GLCanvas::Render()
{
    glNewList(display_list_index, GL_COMPILE_AND_EXECUTE);
    // =========================================================
    // put your GL drawing calls inside the display list for efficiency

    glMatrixMode(GL_MODELVIEW);

    BoundingBox box;
    if (fluid != NULL)
    {
        box = fluid->getBoundingBox();
        if (cloth != NULL)
        {
            box.Extend(cloth->getBoundingBox());
        }
    }
    else
    {
        assert (cloth != NULL);
        box = cloth->getBoundingBox();
    }

    // center the volume in the window
    Matrix m;
    m.SetToIdentity();
    Vec3f center = box.getCenter();
    m *= Matrix::MakeScale(1/float(box.maxDim()));
    m *= Matrix::MakeTranslation(-center);
    glMultMatrixf(m.glGet());

    if (cloth) cloth->Paint();
    if (fluid) fluid->Paint();
    // =========================================================
    glEndList();
    glutPostRedisplay();
}

// ========================================================
// ========================================================

int HandleGLError()
{
    GLenum error;
    int i = 0;
    while ((error = glGetError()) != GL_NO_ERROR)
    {
        printf ("GL ERROR(%d):  %s\n", i, gluErrorString(error));
        i++;
    }
    if (i == 0) return 1;
    return 0;
}
