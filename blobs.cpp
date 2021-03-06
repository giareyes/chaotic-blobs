#include "SETTINGS.h"
#include <cmath>
#include <iostream>
#include <vector>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "TRIANGLE_MESH.h"
#include "STVK.h"
#include "WALL.h"

using namespace std;

// the resolution of the OpenGL window -- independent of the field resolution
int xScreenRes = 800;
int yScreenRes = 800;

// Text for the title bar of the window
string windowLabel("Chaotic Blobs");

// mouse tracking variables
int xMouse         = -1;
int yMouse         = -1;
int mouseButton    = -1;
int mouseState     = -1;
int mouseModifiers = -1;

// animate the current runEverytime()?
bool animate = false;
bool singleStep = false;
// float dt = 1.0/15360.0;
float dt = 1.0/100.0;

// the current viewer eye position
VEC3 eyeCenter(0.9, 0.6, 1);

// current zoom level into the field
float zoom = 2.0;

//Real poissonsRatio = 0.0;
// old:
Real poissonsRatio = 0.4;
Real youngsModulus = 1.0;

// testing:
// Real poissonsRatio = 0.4;
// Real youngsModulus = 0.8;

TRIANGLE_MESH triangleMesh(poissonsRatio, youngsModulus);
TRIANGLE_MESH blob2(poissonsRatio, youngsModulus);
VEC2 bodyForce;

enum SCENE { STRETCH, SQUASH, LSHEAR, RSHEAR, SINGLE, MOTION};
SCENE scene = SINGLE;
int sceneNum = 1;

int meshFlag = 0;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawTriangle(const TRIANGLE& triangle, const VEC4& color)
{
  glColor4f(color[0], color[1], color[2], color[3]);
  glBegin(GL_TRIANGLES);
    for (int x = 0; x < 3; x++)
      glVertex2f(triangle.vertex(x)[0], triangle.vertex(x)[1]);
  glEnd();

  if(meshFlag == 1)
  {
    glColor4f(0,0,0,1);
    glBegin(GL_LINE_LOOP);
      for (int x = 0; x < 3; x++)
        glVertex2f(triangle.vertex(x)[0], triangle.vertex(x)[1]);
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawMesh(const VEC4& color1, const VEC4& color2)
{
  // draw the triangles
  const vector<TRIANGLE>& triangles = triangleMesh.triangles();
  const vector<TRIANGLE>& triangles2 = blob2.triangles();
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    drawTriangle(triangles[x], color1);
    drawTriangle(triangles2[x], color2);
  }


  // get vertices
  vector<VEC2>& vertices = triangleMesh.vertices();
  const vector<int>& constrainedVertices = triangleMesh.constrainedVertices();
  vector<VEC2>& vertices2 = blob2.vertices();
  const vector<int>& constrainedVertices2 = blob2.constrainedVertices();

  //draw eyes
  glPointSize(10.0);
  glColor4f(0,0,0,1);
  glBegin(GL_POINTS);

  glVertex2f(vertices[9][0],
             vertices[9][1]);
  glVertex2f(vertices2[9][0],
             vertices2[9][1]);
  glVertex2f(vertices[11][0],
             vertices[11][1]);
  glVertex2f(vertices2[11][0],
             vertices2[11][1]);

  glEnd();

  // draw the constrained vertices if mesh flag is on
  if(meshFlag == 1)
  {
    glPointSize(5.0);
    glColor4f(1,0,0,1);
    glBegin(GL_POINTS);
    for (unsigned int x = 0; x < constrainedVertices.size(); x++)
    {
      glVertex2f(vertices[constrainedVertices[x]][0],
                 vertices[constrainedVertices[x]][1]);
      glVertex2f(vertices2[constrainedVertices2[x]][0],
                 vertices2[constrainedVertices2[x]][1]);
    }
  }
  glEnd();
}

void drawWalls()
{
  vector<WALL>& walls = triangleMesh.walls();
  // draw the walls
  glColor4f(101.0/255,106.0/255,110.0/255,1);
  for (unsigned int x = 0; x < walls.size(); x++)
  {
    glPushMatrix();
      // translate to the point
      glTranslatef(walls[x].point()[0], walls[x].point()[1], 0);

      // apply a rotation
      float angle = asin(walls[x].normal()[0]) / (2 * M_PI) * 360.0;
      glRotatef(-angle, 0, 0, 1);

      // make it a plane at 0,0
      glTranslatef(0, -0.5, 0);
      glScalef(50,1,1);
      glutSolidCube(1.0);
    glPopMatrix();
  }
  glFlush();
}
///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  float halfZoom = zoom * 0.5;

  glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);

  // set the matrix mode back to modelview
  glMatrixMode(GL_MODELVIEW);

  // set the lookat transform
  glLoadIdentity();
  gluLookAt(eyeCenter[0], eyeCenter[1], 1,  // eye
            eyeCenter[0], eyeCenter[1], 0,  // center
            0, 1, 0);   // up

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawMesh(VEC4(213.0 / 255.0, 246.0 / 255, 247.0 / 255.0,1.0), VEC4(157.0 / 255.0, 230.0 / 255, 156.0 / 255.0,1.0));
  drawWalls();

  glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////
// Map the keyboard keys to something here
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'q':
      exit(0);
      break;

    case 'a':
      animate = !animate;
      break;

    case ' ':
      animate = true;
      singleStep = true;
      break;

    case 'v':
      cout << " eye: " << eyeCenter << endl;
      cout << " zoom: " << zoom << endl;
      break;

    case 's':
      break;

    default:
      break;
  }
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  int modifiers = glutGetModifiers();
  mouseButton = button;
  mouseState = state;
  mouseModifiers = modifiers;

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && modifiers & GLUT_ACTIVE_SHIFT)
  {
    // make sure nothing else is called
    return;
  }

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is clicked and moving
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  if (mouseButton == GLUT_LEFT_BUTTON &&
      mouseState == GLUT_DOWN &&
      mouseModifiers & GLUT_ACTIVE_SHIFT)
  {
    // make sure nothing else is called
    return;
  }

  printf("in mouse motion\n");

  vector<VEC2>& vertices = triangleMesh.vertices();
  vector<int>& unconstrainedVertices = triangleMesh.unconstrainedVertices();
  vector<VEC2>& vertices2 = blob2.vertices();
  vector<int>& unconstrainedVertices2 = blob2.unconstrainedVertices();

  float xDiff = x - xMouse;
  float yDiff = y - yMouse;
  float speed = 0.01;

  if (mouseButton == GLUT_LEFT_BUTTON)
  {
    // TO DO: right now x and y are in viewspace coordiantes. must convert!!!

    // printf("x is: %f\n", (double)x);
    // printf("y is: %f\n", (double)y);
    // eyeCenter[0] -= xDiff * speed;
    // eyeCenter[1] += yDiff * speed;
    double screenX = ((double) x / 400.0) - 0.09;
    double screenY = ((double) (800 - y) / 400.0) - 0.35;
    double mouseX = ((double) xMouse / 400.0) - 0.09;
    double mouseY = ((double) (800 - yMouse) / 400.0) - 0.35;
    // printf("xMouse is: %f\n", screenX);
    // printf("yMouse is: %f\n", screenY);

    for (int i = 0; i < unconstrainedVertices.size(); i++)
    {
      // printf("unconstrained vertex at: (%f %f)\n", vertices[unconstrainedVertices[i]][0], vertices[unconstrainedVertices[i]][1]);
      if (screenX <= vertices[unconstrainedVertices[i]][0] + 0.04
        && screenX >= vertices[unconstrainedVertices[i]][0] - 0.04
        && screenY <= vertices[unconstrainedVertices[i]][1] + 0.03
        && screenY >= vertices[unconstrainedVertices[i]][1] - 0.03)
      {
        // printf("x on vertex\n");
        vertices[unconstrainedVertices[i]][0] = screenX + 0.1;
        vertices[unconstrainedVertices[i]][1] = screenY + 0.1;
      }

      if (screenX <= vertices2[unconstrainedVertices2[i]][0] + 0.04
        && screenX >= vertices2[unconstrainedVertices2[i]][0] - 0.04
        && screenY <= vertices2[unconstrainedVertices2[i]][1] + 0.03
        && screenY >= vertices2[unconstrainedVertices2[i]][1] - 0.03)
      {
        vertices2[unconstrainedVertices2[i]][0] = screenX - 0.1;
        vertices2[unconstrainedVertices2[i]][1] = screenY - 0.1;
      }
    }
  }
  // if (mouseButton == GLUT_RIGHT_BUTTON)
  //   zoom -= yDiff * speed;

  xMouse = x;
  yMouse = y;
}

///////////////////////////////////////////////////////////////////////
// Do something if the mouse is not clicked and moving
///////////////////////////////////////////////////////////////////////
void glutPassiveMouseMotion(int x, int y)
{
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  if (animate)
  {
    static int frame = 0;
    switch (scene) {
      case STRETCH:
        triangleMesh.stretch2(0.01);
        blob2.stretch2(0.01);
        break;
      case RSHEAR:
        triangleMesh.stepShearTest(0.01);
        blob2.stepShearTest(0.01);
        break;
      case LSHEAR:
        triangleMesh.stepShearTest(-0.01);
        blob2.stepShearTest(-0.01);
        break;
      case SQUASH:
        triangleMesh.stretch2(-0.01);
        blob2.stretch2(-0.01);
        break;
      case SINGLE:
        triangleMesh.addBodyForce(bodyForce);
        blob2.addBodyForce(bodyForce);
        break;
      case MOTION:
        break;
    }

    if (scene == MOTION)
    {
      if (frame == 0 && sceneNum == 0) // if we are on the first frame, insert a force so they jump
      {
        triangleMesh.addBodyForce(VEC2(10.0, 200.0));
        blob2.addBodyForce(VEC2(-10.0, 200.0));
      }
      for( int i = 0; i < 15; i ++)
      {
        if (sceneNum == 0)
        {
          blob2.addBodyForce(bodyForce);
          triangleMesh.addBodyForce(bodyForce);
        }
        triangleMesh.stepMotion(dt, bodyForce);
        blob2.stepMotion(dt, bodyForce);
      }
    }
    else
    {
      printf("else\n");
      triangleMesh.stepQuasistatic();
      blob2.stepQuasistatic();
    }
    frame++;

    if (singleStep)
    {
      animate = false;
      singleStep = false;
    }
  }
  glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int glutWindow()
{
  glutInitDisplayMode(GLUT_DOUBLE| GLUT_RGBA);
  glutInitWindowSize(xScreenRes, yScreenRes);
  glutInitWindowPosition(10, 10);
  glutCreateWindow(windowLabel.c_str());

  // set the viewport resolution (w x h)
  glViewport(0, 0, (GLsizei) xScreenRes, (GLsizei) yScreenRes);

  // set the background color to gray
  //glClearColor(0.1, 0.1, 0.1, 0);
  glClearColor(1,1,1,1);

  // register all the callbacks
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);

  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  //glEnable(GL_MULTISAMPLE);
  glLineWidth(1.0);
  glEnable(GL_LINE_SMOOTH);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  // enter the infinite GL loop
  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

string toUpper(const string& input)
{
  string copy(input);
  for (unsigned int x = 0; x < input.length(); x++)
    copy[x] = std::toupper(copy[x]);
  return copy;
}

///////////////////////////////////////////////////////////////////////
// process the command line
///////////////////////////////////////////////////////////////////////
void readCommandLine(int argc, char** argv)
{
  if (argc > 1)
  {
    string sceneType(argv[1]);
    sceneType = toUpper(sceneType);

    if (sceneType.compare("LSHEAR") == 0)
      scene = LSHEAR;
    else if (sceneType.compare("RSHEAR") == 0 || sceneType.compare("SHEAR") == 0 )
      scene = RSHEAR;
    else if (sceneType.compare("SQUASH") == 0)
      scene = SQUASH;
    else if (sceneType.compare("STRETCH") == 0)
      scene = STRETCH;
    else if (sceneType.compare("SINGLE") == 0)
      scene = SINGLE;
    else
    {
      scene = MOTION;
      sceneNum = 0;
    }

    if (argc > 2)
    {
      for(int x = 2; x < argc; x++)
      {
        if (argv[x][1] == 'm')
          meshFlag = 1;
        else if (argv[x][1] == 'b')
          sceneNum = 1;
      }
    }
  }
  else
  {
    scene = MOTION;
    sceneNum = 0;
  }

  // build the scene
  triangleMesh.buildBlob(1.15, sceneNum);
  blob2.buildBlob(0.25, sceneNum);
  bodyForce[0] = 0;
  bodyForce[1] = -0.3;

  triangleMesh.addWall(WALL(VEC2(1,0), VEC2(-0.09,0)));
  triangleMesh.addWall(WALL(VEC2(-1,0), VEC2(1.89,0)));
  triangleMesh.addWall(WALL(VEC2(0,1), VEC2(0,-0.35)));
  blob2.addWall(WALL(VEC2(1,0), VEC2(-0.09,0)));
  blob2.addWall(WALL(VEC2(-1,0), VEC2(1.89,0)));
  blob2.addWall(WALL(VEC2(0,1), VEC2(0,-0.35)));

}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  cout << " Usage: " << argv[0] << " <which scene> " << endl;
  cout << "\t Valid values: " << endl;
  cout << "\t\t <which test>: SINGLE, LSHEAR, RSHEAR, SQUASH, STRETCH, MOTION" << endl;

  readCommandLine(argc, argv);

  // initialize GLUT and GL
  glutInit(&argc, argv);

  // open the GL window
  glutWindow();
  return 1;
}
