#ifndef _FLUID_H_
#define _FLUID_H_

#include <cassert>
#include <cstdio>
#include <string>

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

#include "argparser.h"
#include "boundingbox.h"
#include "vectors.h"
#include "cell.h"

class ArgParser;
class MarchingCubes;

// ========================================================================
// ========================================================================

class Fluid {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Fluid(ArgParser *_args);
  ~Fluid() { delete [] cells; }
  void Load();

  // ===============================
  // ANIMATION & RENDERING FUNCTIONS
  void Animate();
  void Paint() const; 
  BoundingBox getBoundingBox() const {
    return BoundingBox(Vec3f(0,0,0),Vec3f(nx*dx,ny*dy,nz*dz)); }
    
  // ==============
  // CELL ACCESSORS
  int Index(int i, int j, int k) const {
    assert (i >= -1 && i <= nx);
    assert (j >= -1 && j <= ny);
    assert (k >= -1 && k <= nz);
    return (i+1)*(ny+2)*(nz+2) + (j+1)*(nz+2) + (k+1);
  }
  Cell* getCell(int i, int j, int k) const { return &cells[Index(i,j,k)]; }
    
  // ===============
  // OTHER ACCESSORS 
  float getNX() { return nx; }
  float getNY() { return ny; }
  float getNZ() { return nz; }
  float getDX() { return dx; }
  float getDY() { return dy; }
  float getDZ() { return dz; }

private:

  // =================
  // ANIMATION HELPERS
  void AddNewParticles();
  void RemoveParticles();
  void ComputeNewVelocities();
  void SetBoundaryVelocities();
  void CopyVelocities();
  float AdjustForIncompressibility();
  void UpdatePressures();
  void MoveParticles();
  void ReassignParticles();
  void SetEmptySurfaceFull();

  // =====================
  // NAVIER-STOKES HELPERS
  Vec3f getInterpolatedVelocity(const Vec3f &pos) const;
  float getPressure(int i, int j, int k) const { return getCell(i,j,k)->getPressure(); }
  // velocity accessors
  float get_u_plus(int i, int j, int k) const { return getCell(i,j,k)->get_u_plus(); }
  float get_v_plus(int i, int j, int k) const { return getCell(i,j,k)->get_v_plus(); }  
  float get_w_plus(int i, int j, int k) const { return getCell(i,j,k)->get_w_plus(); }  
  float get_new_u_plus(int i, int j, int k) const { return getCell(i,j,k)->get_new_u_plus(); }  
  float get_new_v_plus(int i, int j, int k) const { return getCell(i,j,k)->get_new_v_plus(); }  
  float get_new_w_plus(int i, int j, int k) const { return getCell(i,j,k)->get_new_w_plus(); }  
  float get_u_avg(int i, int j, int k) const { return 0.5*(get_u_plus(i-1,j,k)+get_u_plus(i,j,k)); }
  float get_v_avg(int i, int j, int k) const { return 0.5*(get_v_plus(i,j-1,k)+get_v_plus(i,j,k)); }
  float get_w_avg(int i, int j, int k) const { return 0.5*(get_w_plus(i,j,k-1)+get_w_plus(i,j,k)); }
  float get_uv_plus(int i, int j, int k) const { 
    return 0.5*(get_u_plus(i,j,k) + get_u_plus(i,j+1,k)) * 0.5*(get_v_plus(i,j,k) + get_v_plus(i+1,j,k)); }
  float get_uw_plus(int i, int j, int k) const { 
    return 0.5*(get_u_plus(i,j,k) + get_u_plus(i,j,k+1)) * 0.5*(get_w_plus(i,j,k) + get_w_plus(i+1,j,k)); }
  float get_vw_plus(int i, int j, int k) const { 
    return 0.5*(get_v_plus(i,j,k) + get_v_plus(i,j,k+1)) * 0.5*(get_w_plus(i,j,k) + get_w_plus(i,j+1,k)); }
  // velocity modifiers
  void set_new_u_plus(int i, int j, int k, float f) { getCell(i,j,k)->set_new_u_plus(f); }
  void set_new_v_plus(int i, int j, int k, float f) { getCell(i,j,k)->set_new_v_plus(f); }
  void set_new_w_plus(int i, int j, int k, float f) { getCell(i,j,k)->set_new_w_plus(f); }
  void adjust_new_u_plus(int i, int j, int k, float f) { getCell(i,j,k)->adjust_new_u_plus(f); }
  void adjust_new_v_plus(int i, int j, int k, float f) { getCell(i,j,k)->adjust_new_v_plus(f); }
  void adjust_new_w_plus(int i, int j, int k, float f) { getCell(i,j,k)->adjust_new_w_plus(f); }

  // ========================================
  // RENDERING SURFACE (using Marching Cubes)
  float interpolateIsovalue(const Vec3f &c) const;
  float getIsovalue(int i, int j, int k) const;

  // ============
  // LOAD HELPERS
  bool inShape(Vec3f &pos, const std::string &shape);
  void GenerateParticles(const std::string &shape, const std::string &placement);
  
private:

  // don't use this constructor
  Fluid() { assert(0); }
  
  // ==============
  // REPRESENTATION
  ArgParser *args;

  // fluid parameters
  int nx,ny,nz;    // number of grid cells in each dimension
  float dx,dy,dz;  // dimensions of each grid cell
  Cell *cells;     // NOTE: padded with extra cells on each side

  // simulation parameters
  bool xy_free_slip;
  bool yz_free_slip;
  bool zx_free_slip;
  bool compressible;
  float viscosity;
  float density; // average # of particles initialized in each "Full" cell
  
  bool source;
  std::string sourcetype;
  float sourcedensity;
  float sourcevelocity;
  
  bool sink;
  std::string sinktype;

  MarchingCubes *marchingCubes;  // to display an isosurface 
};

#endif
