#include "boundingbox.h"

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

#include "cloth.h"
#include "argparser.h"
#include <fstream>
#include "vectors.h"
#include <iostream>

using namespace std;

// ================================================================================
// some helper drawing functions

void insertNormal(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3, const Vec3f &p4);
void DrawSpring(const Vec3f &a_o, const Vec3f &b_o, const Vec3f &a, const Vec3f &b, float correction);
void CheckCollision(ClothParticle &p1, ClothParticle &p2);

// ================================================================================
// ================================================================================

Cloth::Cloth(ArgParser *_args)
{
  args =_args;
  fluid = NULL;

  // open the file
  ifstream istr(args->cloth_file.c_str());
  assert (istr != NULL);
  string token;

  // read in the simulation parameters
  istr >> token >> k_structural;
  assert (token == "k_structural");  // (units == N/m)  (N = kg*m/s^2)
  istr >> token >> k_shear;
  assert (token == "k_shear");
  istr >> token >> k_bend;
  assert (token == "k_bend");
  istr >> token >> damping;
  assert (token == "damping");
  // istr >> token >> normal;
  // assert (token == "normal");
  // NOTE: correction factor == .1, means springs shouldn't stretch more than 10%
  //       correction factor == 100, means don't do any correction
  istr >> token >> provot_structural_correction;
  assert (token == "provot_structural_correction");
  istr >> token >> provot_shear_correction;
  assert (token == "provot_shear_correction");

  // the cloth dimensions
  istr >> token >> nx >> ny; // (units == meters)
  assert (token == "m");
  assert (nx >= 2 && ny >= 2);

  // the corners of the cloth
  Vec3f a,b,c,d;
  istr >> token >> a;
  assert (token == "p");
  istr >> token >> b;
  assert (token == "p");
  istr >> token >> c;
  assert (token == "p");
  istr >> token >> d;
  assert (token == "p");

  // fabric weight  (units == kg/m^2)
  // denim ~300 g/m^2
  // silk ~70 g/m^2
  float fabric_weight;
  istr >> token >> fabric_weight;
  assert (token == "fabric_weight");
  float area = AreaOfTriangle(a,b,c) + AreaOfTriangle(a,c,d);

  // create the particles
  particles = new ClothParticle[nx*ny];
  float mass = area*fabric_weight / float(nx*ny);
  for (int i = 0; i < nx; i++)
  {
    float x = i/float(nx-1);
    Vec3f ab = (1-x)*a + x*b;
    Vec3f dc = (1-x)*d + x*c;
    for (int j = 0; j < ny; j++)
    {
      float y = j/float(ny-1);
      ClothParticle &p = getParticle(i,j);
      Vec3f abdc = (1-y)*ab + y*dc;
      p.setOriginalPosition(abdc);
      p.setPosition(abdc);
      p.setVelocity(Vec3f(0,0,0));
      p.setMass(mass);
      p.setFixed(false);
    }
  }
  
  

  // the fixed particles
  while (istr >> token)
  {
    assert (token == "f");
    int i,j;
    float x,y,z;
    istr >> i >> j >> x >> y >> z;
    ClothParticle &p = getParticle(i,j);
    p.setPosition(Vec3f(x,y,z));
    p.setFixed(true);
  }

  computeBoundingBox();
  
  collision_boundary = box.maxDim()/(2*max2(nx,ny));
  collision_boundary *= 1.3;
}

Cloth::Cloth(ArgParser *_args, Fluid *_fluid)
{
  args =_args;
  fluid = _fluid;

  // open the file
  ifstream istr(args->cloth_file.c_str());
  assert (istr != NULL);
  string token;

  // read in the simulation parameters
  istr >> token >> k_structural;
  assert (token == "k_structural");  // (units == N/m)  (N = kg*m/s^2)
  istr >> token >> k_shear;
  assert (token == "k_shear");
  istr >> token >> k_bend;
  assert (token == "k_bend");
  istr >> token >> damping;
  assert (token == "damping");
  // istr >> token >> normal;
  // assert (token == "normal");
  // NOTE: correction factor == .1, means springs shouldn't stretch more than 10%
  //       correction factor == 100, means don't do any correction
  istr >> token >> provot_structural_correction;
  assert (token == "provot_structural_correction");
  istr >> token >> provot_shear_correction;
  assert (token == "provot_shear_correction");

  // the cloth dimensions
  istr >> token >> nx >> ny; // (units == meters)
  assert (token == "m");
  assert (nx >= 2 && ny >= 2);

  // the corners of the cloth
  Vec3f a,b,c,d;
  istr >> token >> a;
  assert (token == "p");
  istr >> token >> b;
  assert (token == "p");
  istr >> token >> c;
  assert (token == "p");
  istr >> token >> d;
  assert (token == "p");

  // fabric weight  (units == kg/m^2)
  // denim ~300 g/m^2
  // silk ~70 g/m^2
  float fabric_weight;
  istr >> token >> fabric_weight;
  assert (token == "fabric_weight");
  float area = AreaOfTriangle(a,b,c) + AreaOfTriangle(a,c,d);

  // create the particles
  particles = new ClothParticle[nx*ny];
  float mass = area*fabric_weight / float(nx*ny);
  for (int i = 0; i < nx; i++)
  {
    float x = i/float(nx-1);
    Vec3f ab = (1-x)*a + x*b;
    Vec3f dc = (1-x)*d + x*c;
    for (int j = 0; j < ny; j++)
    {
      float y = j/float(ny-1);
      ClothParticle &p = getParticle(i,j);
      Vec3f abdc = (1-y)*ab + y*dc;
      p.setOriginalPosition(abdc);
      p.setPosition(abdc);
      p.setVelocity(Vec3f(0,0,0));
      p.setMass(mass);
      p.setFixed(false);
    }
  }
  
  

  // the fixed particles
  while (istr >> token)
  {
    assert (token == "f");
    int i,j;
    float x,y,z;
    istr >> i >> j >> x >> y >> z;
    ClothParticle &p = getParticle(i,j);
    p.setPosition(Vec3f(x,y,z));
    p.setFixed(true);
  }

  computeBoundingBox();
  
  // assign cloth particles to the proper cells in the voxel grid
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      ClothParticle &p = getParticle(i,j);
      Vec3f pos = p.getPosition();
      Cell *cell = fluid->getCell(int(pos.x()/fluid->getDX()),int(pos.y()/fluid->getDY()),int(pos.z()/fluid->getDZ()));
      cell->addClothParticle(&p);
    }
  }
  
  collision_boundary = box.maxDim()/(2*max2(nx,ny));
  collision_boundary *= 1.3;
}

// ================================================================================

void Cloth::computeBoundingBox()
{
    box = BoundingBox(getParticle(0,0).getPosition());
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            box.Extend(getParticle(i,j).getPosition());
        }
    }
}

// ================================================================================
// ================================================================================
// ================================================================================

void DrawForce(const Vec3f &p, const Vec3f &f);

// ================================================================================
// Calculate "Force" without spring constant
// ================================================================================
Vec3f CalculateForce(const ClothParticle &p1, const ClothParticle &p2, float kF)
{
    Vec3f p1Orig = p1.getOriginalPosition();
    Vec3f p2Orig = p2.getOriginalPosition();
    Vec3f p1Curr = p1.getPosition();
    Vec3f p2Curr = p2.getPosition();

    Vec3f currDiff = p1Curr - p2Curr;
    Vec3f origDiff = p1Orig - p2Orig;
    Vec3f currDiffDir = Vec3f(currDiff.x()/currDiff.Length(), currDiff.y()/currDiff.Length(), currDiff.z()/currDiff.Length());
    Vec3f force = (currDiff.Length() - origDiff.Length()) * currDiffDir;

    return (force*(-kF));
}

Vec3f GetUpdraft()
{
    return Vec3f(rand()%30-15, 15+(rand()%10-5), rand()%20-10);
}

void Cloth::Paint() const
{

    // =====================================================================================
    // render the bounding box
    // =====================================================================================
    if (args->bounding_box)
        box.Paint();

    // =====================================================================================
    // render the particles
    // =====================================================================================

    if (args->particles)
    {
        glDisable(GL_LIGHTING);
        glPointSize(3);
        glBegin(GL_POINTS);
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                const ClothParticle &p = getParticle(i,j);
                const Vec3f &pos = p.getPosition();
                bool fixed = p.isFixed();
                if (fixed)
                    glColor3f(0,1,0);
                else
                    glColor3f(0,0,0);
                glVertex3f(pos.x(),pos.y(),pos.z());
            }
        }
        glEnd();
        glEnable(GL_LIGHTING);
    }

    // =====================================================================================
    // render the velocity at each particle
    // =====================================================================================

    if (args->velocity)
    {
        glDisable(GL_LIGHTING);
        glLineWidth(2);
        glBegin(GL_LINES);
        glColor3f(1,0,0);
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                const ClothParticle &p = getParticle(i,j);
                const Vec3f &pos = p.getPosition();
                const Vec3f &vel = p.getVelocity();
                Vec3f v = pos + vel;
                glVertex3f(pos.x(),pos.y(),pos.z());
                glVertex3f(v.x(),v.y(),v.z());
            }
        }
        glEnd();
        glEnable(GL_LIGHTING);
    }

    // =====================================================================================
    // render the force at each particle
    // =====================================================================================

    if (args->force)
    {
        glDisable(GL_LIGHTING);
        glLineWidth(2);
        glBegin(GL_LINES);
        glColor3f(0,0,1);
        Vec3f sumForce;
        Vec3f updraftF;
        if (args->updraft)
        {
            updraftF = GetUpdraft();
        }
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                // ASSIGNMENT:
                // visualize the forces acting on each particle
                const ClothParticle &p = getParticle(i,j);
                if (args->updraft)
                {
                    sumForce = updraftF*p.getMass();
                }
                else
                {
                    sumForce = Vec3f(0,0,0);
                }

                sumForce += args->gravity*p.getMass();
                // calculate Structural forces
                if (i > 0)
                {
                    const ClothParticle &p1 = getParticle(i-1,j);
                    sumForce += CalculateForce(p, p1,k_structural);
                }
                if (j > 0)
                {
                    const ClothParticle &p2 = getParticle(i,j-1);
                    sumForce += CalculateForce(p, p2,k_structural);
                }
                if (i < nx-1)
                {
                    const ClothParticle &p3 = getParticle(i+1,j);
                    sumForce += CalculateForce(p, p3,k_structural);
                }
                if (j < ny-1)
                {
                    const ClothParticle &p4 = getParticle(i,j+1);
                    sumForce += CalculateForce(p, p4,k_structural);
                }
                // calculate Flexion forces
                if (i > 1)
                {
                    const ClothParticle &p5 = getParticle(i-2,j);
                    sumForce += CalculateForce(p, p5,k_bend);
                }
                if (j > 1)
                {
                    const ClothParticle &p6 = getParticle(i,j-2);
                    sumForce += CalculateForce(p, p6,k_bend);
                }
                if (i < nx-2)
                {
                    const ClothParticle &p7 = getParticle(i+2,j);
                    sumForce += CalculateForce(p, p7,k_bend);
                }
                if (j < ny-2)
                {
                    const ClothParticle &p8 = getParticle(i,j+2);
                    sumForce += CalculateForce(p, p8,k_bend);
                }
                // calculate Shear forces
                if (i > 0 && j > 0)
                {
                    const ClothParticle &p9 = getParticle(i-1,j-1);
                    sumForce += CalculateForce(p, p9,k_shear);
                }
                if (i > 0 && j < ny-1)
                {
                    const ClothParticle &p10 = getParticle(i-1,j+1);
                    sumForce += CalculateForce(p, p10,k_shear);
                }
                if (i < nx-1 && j > 0)
                {
                    const ClothParticle &p11 = getParticle(i+1,j-1);
                    sumForce += CalculateForce(p, p11,k_shear);
                }
                if (i < nx-1 && j < ny-1)
                {
                    const ClothParticle &p12 = getParticle(i+1,j+1);
                    sumForce += CalculateForce(p, p12,k_shear);
                }
                DrawForce(p.getPosition(), sumForce);
            }
        }
        glEnd();
        glEnable(GL_LIGHTING);
    }

    // =====================================================================================
    // render the cloth surface
    // =====================================================================================

    if (args->surface)
    {
        glColor3f(1,1,1);
        glBegin(GL_QUADS);
        for (int i = 0; i < nx-1; i++)
        {
            for (int j = 0; j < ny-1; j++)
            {
                const Vec3f &a = getParticle(i,j).getPosition();
                const Vec3f &b = getParticle(i,j+1).getPosition();
                const Vec3f &c = getParticle(i+1,j+1).getPosition();
                const Vec3f &d = getParticle(i+1,j).getPosition();
                insertNormal(a,b,c,d);
                glVertex3f(a.x(),a.y(),a.z());
                glVertex3f(b.x(),b.y(),b.z());
                glVertex3f(c.x(),c.y(),c.z());
                glVertex3f(d.x(),d.y(),d.z());
            }
        }
        glEnd();
    }

    // =====================================================================================
    // render the wireframe cloth (the structural and shear springs)
    // =====================================================================================

    if (args->wireframe)
    {
        glDisable(GL_LIGHTING);
        glLineWidth(1);
        glBegin(GL_LINES);
        for (int i = 0; i < nx-1; i++)
        {
            for (int j = 0; j < ny-1; j++)
            {
                const Vec3f &a_o = getParticle(i,j).getOriginalPosition();
                const Vec3f &b_o = getParticle(i,j+1).getOriginalPosition();
                const Vec3f &c_o = getParticle(i+1,j+1).getOriginalPosition();
                const Vec3f &d_o = getParticle(i+1,j).getOriginalPosition();
                const Vec3f &a = getParticle(i,j).getPosition();
                const Vec3f &b = getParticle(i,j+1).getPosition();
                const Vec3f &c = getParticle(i+1,j+1).getPosition();
                const Vec3f &d = getParticle(i+1,j).getPosition();
                DrawSpring(a_o,b_o,a,b, provot_structural_correction);
                DrawSpring(a_o,c_o,a,c, provot_shear_correction);
                DrawSpring(a_o,d_o,a,d, provot_structural_correction);
                DrawSpring(b_o,c_o,b,c, provot_structural_correction);
                DrawSpring(b_o,d_o,b,d, provot_shear_correction);
                DrawSpring(c_o,d_o,c,d, provot_structural_correction);
            }
        }
        glEnd();
        glEnable(GL_LIGHTING);
    }
}

// ================================================================================
// ================================================================================

void Cloth::Animate()
{
    // ASSIGNMENT:
    // move the vertices of the mesh

    Vec3f sumForce;
    glDisable(GL_LIGHTING);
    glPointSize(3);
    glBegin(GL_POINTS);

    Vec3f updraftF;
    if (args->updraft)
    {
        updraftF = GetUpdraft();
    }
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            ClothParticle &p1 = getParticle(i,j);
            if (args->updraft)
            {
                sumForce = updraftF*p1.getMass();
            }
            else
            {
                sumForce = Vec3f(0,0,0);
            }
            if (!p1.isFixed())
            {
                sumForce += args->gravity*p1.getMass();
                // calculate Structural forces
                if (i > 0)
                {
                    const ClothParticle &p2 = getParticle(i-1,j);
                    sumForce += CalculateForce(p1, p2, k_structural);
                }
                if (j > 0)
                {
                    const ClothParticle &p3 = getParticle(i,j-1);
                    sumForce += CalculateForce(p1, p3, k_structural);
                }
                if (i < nx-1)
                {
                    const ClothParticle &p4 = getParticle(i+1,j);
                    sumForce += CalculateForce(p1, p4, k_structural);
                }
                if (j < ny-1)
                {
                    const ClothParticle &p5 = getParticle(i,j+1);
                    sumForce += CalculateForce(p1, p5, k_structural);
                }
                // calculate Flexion forces
                if (i > 1)
                {
                    const ClothParticle &p6 = getParticle(i-2,j);
                    sumForce += CalculateForce(p1, p6, k_bend);
                }
                if (j > 1)
                {
                    const ClothParticle &p7 = getParticle(i,j-2);
                    sumForce += CalculateForce(p1, p7, k_bend);
                }
                if (i < nx-2)
                {
                    const ClothParticle &p8 = getParticle(i+2,j);
                    sumForce += CalculateForce(p1, p8, k_bend);
                }
                if (j < ny-2)
                {
                    const ClothParticle &p9 = getParticle(i,j+2);
                    sumForce += CalculateForce(p1, p9, k_bend);
                }
                // calculate Shear forces
                if (i > 0 && j > 0)
                {
                    const ClothParticle &p10 = getParticle(i-1,j-1);
                    sumForce += CalculateForce(p1, p10, k_shear);
                }
                if (i > 0 && j < ny-1)
                {
                    const ClothParticle &p11 = getParticle(i-1,j+1);
                    sumForce += CalculateForce(p1, p11, k_shear);
                }
                if (i < nx-1 && j > 0)
                {
                    const ClothParticle &p12 = getParticle(i+1,j-1);
                    sumForce += CalculateForce(p1, p12, k_shear);
                }
                if (i < nx-1 && j < ny-1)
                {
                    const ClothParticle &p13 = getParticle(i+1,j+1);
                    sumForce += CalculateForce(p1, p13, k_shear);
                }
                sumForce += (-damping * p1.getVelocity());
                p1.setAcceleration(Vec3f(sumForce.x()/p1.getMass(), sumForce.y()/p1.getMass(), sumForce.z()/p1.getMass()));
                p1.setVelocity((p1.getAcceleration()*args->timestep)+p1.getVelocity());
                p1.setPosition(p1.getPosition()+(p1.getVelocity()*args->timestep));
            }

            if (i > 0)
            {
                ClothParticle &p2 = getParticle(i-1, j);
                makeCorrection(p1, p2, provot_structural_correction);
            }
            if (j > 0)
            {
                ClothParticle &p2 = getParticle(i, j-1);
                makeCorrection(p1, p2, provot_structural_correction);
            }
            if (i > 0 && j > 0)
            {
                ClothParticle &p2 = getParticle(i-1, j-1);
                makeCorrection(p1, p2, provot_shear_correction);
            }
            if (i > 0 && j < ny-1)
            {
                ClothParticle &p2 = getParticle(i-1, j+1);
                makeCorrection(p1, p2, provot_shear_correction);
            }
        }
    }
    if (fluid) AdjustCells();
    if (!strcmp(args->collision, "both") || !strcmp(args->collision, "cloth"))CheckCollision();
    if (fluid) AdjustCells();
    glEnd();
    glEnable(GL_LIGHTING);
}

// puts each cloth particle into its proper cell if it has been moved out of its old one.
void Cloth::AdjustCells()
{
  for (int i = 0; i < fluid->getNX(); i++)
  {
    for (int j = 0; j < fluid->getNY(); j++)
    {
      for (int k = 0; k < fluid->getNZ(); k++)
      {
        Cell *currentcell = fluid->getCell(i,j,k);
        std::vector<ClothParticle*> pvec = currentcell->getClothParticles();
        for (unsigned int q = 0; q < pvec.size(); q++)
        {
          ClothParticle *p = pvec[q];
          Vec3f pos = p->getPosition();
          Cell *cell = fluid->getCell(int(pos.x()/fluid->getDX()),int(pos.y()/fluid->getDY()),int(pos.z()/fluid->getDZ()));
          if (cell != currentcell)
          {
            currentcell->removeClothParticle(p);
            cell->addClothParticle(p);
          }
        }
      }
    }
  }
}

void Cloth::makeCorrection(ClothParticle &p1, ClothParticle &p2, float corr)
{
    Vec3f len_diff, len_diff1, len_diff2;
    Vec3f ab_o;
    Vec3f ab;
    float length_o;     // the original length
    float length;       // the current length


    if (!p1.isFixed() && !p2.isFixed())
    {
        ab_o = p1.getOriginalPosition()-p2.getOriginalPosition();
        ab = p1.getPosition()-p2.getPosition();
        length_o = ab_o.Length();
        length = ab.Length();
        ab.Normalize();
        if (length > (1+corr) * length_o)
        {
            len_diff = (length-((1+corr)*length_o))*ab;
            len_diff1 = Vec3f(len_diff.x()/2, len_diff.y()/2, len_diff.z()/2);
            len_diff2 = -1*(len_diff1);
            p2.setPosition(p2.getPosition()+len_diff1);
            p1.setPosition(p1.getPosition()+len_diff2);
        }
        else if (length < (1-corr) * length_o)
        {
            len_diff = (((1-corr)*length_o) - length)*ab;
            len_diff1 = Vec3f(len_diff.x()/2, len_diff.y()/2, len_diff.z()/2);
            len_diff2 = -1*(len_diff1);
            p2.setPosition(p2.getPosition()-len_diff1);
            p1.setPosition(p1.getPosition()-len_diff2);
        }
    }
    else if (!p1.isFixed() && p2.isFixed())
    {
        ab_o = p1.getOriginalPosition()-p2.getOriginalPosition();
        ab = p1.getPosition()-p2.getPosition();
        length_o = ab_o.Length();
        length = ab.Length();
        ab.Normalize();
        if (length > (1+corr) * length_o)
        {
            len_diff = (length-((1+corr)*length_o))*ab;
            p1.setPosition(p1.getPosition()-len_diff);
        }
        else if (length < (1-corr) * length_o)
        {
            len_diff = (((1-corr)*length_o) - length)*ab;
            p1.setPosition(p1.getPosition()+len_diff);
        }
    }
    else if (p1.isFixed() && !p2.isFixed())
    {
        ab_o = p1.getOriginalPosition()-p2.getOriginalPosition();
        ab = p1.getPosition()-p2.getPosition();
        length_o = ab_o.Length();
        length = ab.Length();
        ab.Normalize();
        if (length > (1+corr) * length_o)
        {
            len_diff = (length-((1+corr)*length_o))*ab;
            p2.setPosition(p2.getPosition()+len_diff);
        }
        else if (length < (1-corr) * length_o)
        {
            len_diff = (((1-corr)*length_o) - length)*ab;
            p2.setPosition(p2.getPosition()-len_diff);
        }
    }
}

void Cloth::CheckCollision()
{
  for (int i = 0; i < fluid->getNX(); i++)
  {
    for (int j = 0; j < fluid->getNY(); j++)
    {
      for (int k = 0; k < fluid->getNZ(); k++)
      {
        std::vector<Cell*> cells;
        if (i > 0)
        {
          if (j > 0)
          {
            if (k > 0)  cells.push_back(fluid->getCell(i-1,j-1,k-1));          
            cells.push_back(fluid->getCell(i-1,j-1,k));
            if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i-1,j-1,k+1));
          }
          
          if (k > 0)  cells.push_back(fluid->getCell(i-1,j,k-1));
          cells.push_back(fluid->getCell(i-1,j,k));
          if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i-1,j,k+1));
          
          if (j < fluid->getNY()-1)
          {
            if (k > 0)  cells.push_back(fluid->getCell(i-1,j+1,k-1));
            cells.push_back(fluid->getCell(i-1,j+1,k));
            if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i-1,j+1,k+1));
          }
        }
        else if (i < fluid->getNX()-1)
        {
          if (j > 0)
          {
            if (k > 0)  cells.push_back(fluid->getCell(i+1,j-1,k-1));
            cells.push_back(fluid->getCell(i+1,j-1,k));
            if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i+1,j-1,k+1));
          }
          
          if (k > 0)  cells.push_back(fluid->getCell(i+1,j,k-1));
          cells.push_back(fluid->getCell(i+1,j,k));
          if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i+1,j,k+1));
          
          if (j < fluid->getNY()-1)
          { 
            if (k > 0)  cells.push_back(fluid->getCell(i+1,j+1,k-1));
            cells.push_back(fluid->getCell(i+1,j+1,k));
            if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i+1,j+1,k+1));
          }
        }
        else
        {
          if (j > 0)
          {
            if (k > 0)  cells.push_back(fluid->getCell(i,j-1,k-1)); 
            cells.push_back(fluid->getCell(i,j-1,k));
            if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i,j-1,k+1));
          }
          
          if (j < fluid->getNY()-1)
          {
            if (k > 0)  cells.push_back(fluid->getCell(i,j+1,k-1));
            cells.push_back(fluid->getCell(i,j+1,k));
            if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i,j+1,k+1)); 
          }
          
          if (k > 0)  cells.push_back(fluid->getCell(i,j,k-1));
          if (k < fluid->getNZ()-1) cells.push_back(fluid->getCell(i,j,k+1));
        }
        Cell *currentcell = fluid->getCell(i,j,k);
        if (!strcmp(args->collision, "cloth"))
        {
          std::vector<ClothParticle*> pvec = currentcell->getClothParticles();
          for (unsigned int q = 0; q < pvec.size(); q++)
          {
            ClothParticle *p = pvec[q];
            for (unsigned int c = 0; c < cells.size(); c++)
            {
              std::vector<ClothParticle*> pvec2 = cells[c]->getClothParticles();
              for (unsigned int cp = 0; cp < pvec2.size(); cp++)
              {
                ClothParticle *p1 = pvec2[cp];
                if (abs((p1->getPosition() - p->getPosition()).Length()) <= collision_boundary)
                {
                  Vec3f temp = p1->getVelocity();
                  p1->setVelocity(p->getVelocity());
                  p->setVelocity(temp);
                }
              }
            }
            for (unsigned int c1 = 0; c1 < pvec.size(); c1++)
            {
              ClothParticle *p1 = pvec[c1];
              if (q == c1)
                continue;
              if (abs((p1->getPosition() - p->getPosition()).Length()) <= collision_boundary)
              {
                Vec3f p1v = p1->getVelocity()*p1->getMass();
                Vec3f pv = p->getVelocity()*p->getMass();
                Vec3f sum = Vec3f(p1v.x()+pv.x(),p1v.y()+pv.y(),p1v.z()+pv.z());
                p1->setVelocity(sum);
                p->setVelocity(sum);
              }
            }
          }
          cells.clear();
        }
        else if (!strcmp(args->collision, "both"))
        {
          cells.push_back(currentcell);
          std::vector<FluidParticle*> pvec = currentcell->getParticles();
          for (unsigned int q = 0; q < pvec.size(); q++)
          {
            FluidParticle *p = pvec[q];
            for (unsigned int c = 0; c < cells.size(); c++)
            {
              std::vector<ClothParticle*> pvec2 = cells[c]->getClothParticles();
              for (unsigned int cp = 0; cp < pvec2.size(); cp++)
              {
                ClothParticle *p1 = pvec2[cp];
                if (abs((p1->getPosition() - p->getPosition()).Length()) <= collision_boundary)
                {
                  double rate = ((currentcell->getParticles().size()-1)/currentcell->getParticles().size());
                  // std::cout << "fluid-cloth collision\n";
                  Vec3f p1v = p1->getVelocity()*p1->getMass();
                  Vec3f pv = p->getVelocity()*fluid->getMass();
                  Vec3f sum = Vec3f(p1v.x()+pv.x(),p1v.y()+pv.y(),0);
                  p1->setVelocity(sum);
                  p->setVelocity(sum);
                  currentcell->set_u_plus(currentcell->get_u_plus()*rate);
                  currentcell->set_v_plus(currentcell->get_v_plus()*rate);
                  currentcell->set_w_plus(currentcell->get_w_plus()*rate);
                }
              }
            }
          }
          cells.clear();
        }
      }
    }
  }
}

// ================================================================================
// ================================================================================
// some helper drawing functions
// ================================================================================

void DrawSpring(const Vec3f &a_o, const Vec3f &b_o, const Vec3f &a, const Vec3f &b, float correction)
{
    Vec3f ab_o = b_o-a_o;
    Vec3f ab = b-a;
    float length_o = ab_o.Length(); // the original length
    float length = ab.Length();     // the current length
    if (length > (1+correction) * length_o)
    {
        // draw the spring in cyan if it's over-stretched
        glColor3f(0,1,1);
    }
    else if (length < (1-correction) * length_o)
    {
        // draw the spring in magenta if it's under-stretched
        glColor3f(1,0,1);
    }
    else
    {
        // otherwise draw it black
        glColor3f(0,0,0);
    }
    glVertex3f(a.x(),a.y(),a.z());
    glVertex3f(b.x(),b.y(),b.z());
}

void insertNormal(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3, const Vec3f &p4)
{
    // compute the normal of a non-planar quadrilateral (not well-defined)
    Vec3f v12 = p2 - p1;
    Vec3f v23 = p3 - p2;
    Vec3f v34 = p4 - p3;
    Vec3f v41 = p1 - p4;

    // compute normals at each corner and average them
    Vec3f normal1, normal2, normal3, normal4;
    Vec3f::Cross3(normal1,v41,v12);
    Vec3f::Cross3(normal2,v12,v23);
    Vec3f::Cross3(normal3,v23,v34);
    Vec3f::Cross3(normal4,v34,v41);
    Vec3f normal = normal1+normal2+normal3+normal4;

    normal.Normalize();
    glNormal3f(normal.x(), normal.y(), normal.z());
}

void DrawForce(const Vec3f &p, const Vec3f &f)
{
    Vec3f tmp = p+f;
    glVertex3f(p.x(),p.y(),p.z());
    glVertex3f(tmp.x(),tmp.y(),tmp.z());
}



// ================================================================================
// ================================================================================

