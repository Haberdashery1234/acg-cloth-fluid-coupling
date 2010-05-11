#include "fluid.h"
#include "argparser.h"
#include "boundingbox.h"
#include "vectors.h"
#include "matrix.h"
#include "marching_cubes.h"
#include <fstream>
#include <iomanip>

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

#define BETA_0 1.7
#define EPSILON 0.0001

void CubePaint();
void ClearCubePaint();

// ==============================================================
// ==============================================================
// CONSTRUCTOR
// ==============================================================
// ==============================================================

Fluid::Fluid(ArgParser *_args) {
  args = _args;
  Load();
  marchingCubes = new MarchingCubes(nx+1,ny+1,nz+1,dx,dy,dz);
  SetEmptySurfaceFull();
}

// ==============================================================

void Fluid::Load() {    

  // open the file
  assert (args->fluid_file != "");
  std::ifstream istr(args->fluid_file.c_str());
  assert (istr != NULL);
  std::string token, token2, token3;

  // load in the grid size & dimensions
  istr >> token >> nx >> ny >> nz;  assert (token=="grid");
  assert (nx > 0 && ny > 0 && nz > 0);
  istr >> token >> dx >> dy >> dz; assert (token=="cell_dimensions");
  cells = new Cell[(nx+2)*(ny+2)*(nz+2)];

  // simulation parameters
  istr >> token >> token2;  assert (token=="flow");
  if (token2 == "compressible") compressible = true;
  else { assert (token2 == "incompressible"); compressible = false; }
  istr >> token >> token2;  assert (token=="xy_boundary");
  if (token2 == "free_slip") xy_free_slip = true;
  else { assert  (token2 == "no_slip"); xy_free_slip = false; }
  istr >> token >> token2;  assert (token=="yz_boundary");
  if (token2 == "free_slip") yz_free_slip = true;
  else { assert  (token2 == "no_slip"); yz_free_slip = false; }
  istr >> token >> token2;  assert (token=="zx_boundary");
  if (token2 == "free_slip") zx_free_slip = true;
  else { assert  (token2 == "no_slip"); zx_free_slip = false; }
  istr >> token >> viscosity;  assert (token=="viscosity");
  float gravity;
  istr >> token >> gravity;  assert (token=="gravity");
  args->gravity = Vec3f(0,-9.8,0) * gravity;
  
  // initialize marker particles 
  istr >> token >> token2 >> token3;  assert (token=="initial_particles");
  istr >> token >> density;  assert (token=="density");
  GenerateParticles(token2,token3);

  // initialize velocities
  istr >> token >> token2;  assert (token=="initial_velocity");
  if (token2 == "zero") {
    // default is zero
  } else {
    assert (token2 == "random");
    int i,j,k;
    float max_dim = max3(dx,dy,dz);
    for (i = -1; i <= nx; i++) {
      for (j = -1; j <= ny; j++) {
        for (k = -1; k <= nz; k++) {
          getCell(i,j,k)->set_u_plus((2*(rand() / double(RAND_MAX))-1)*max_dim);
          getCell(i,j,k)->set_v_plus((2*(rand() / double(RAND_MAX))-1)*max_dim);
          getCell(i,j,k)->set_w_plus((2*(rand() / double(RAND_MAX))-1)*max_dim);
        }
      }
    }
  }
  // read in sources and sinks
  istr >> token >> token2; assert (token =="source");
  if (token2 == "none") {
    source = false;
  } else if (token2 == "spout" || token2 == "rain" || token2 == "fountain" || token2 == "pour") {
    source = true;
    istr >> sourcedensity;
    istr >> sourcevelocity;
    sourcetype = token2;
  } else {
    assert(0);
  }
  istr >> token >> token2; assert (token =="sink");
  if (token2 == "none") {
    sink = false;
  } else if (token2 == "end" || token2 == "bottom") {
    sink = true;
    sinktype = token2;
  } else {
    assert(0);
  }
  // read in custom velocities
  while(istr >> token) {
    int i,j,k;
    float velocity;
    assert (token == "u" || token == "v" || token == "w");
    istr >> i >> j >> k >> velocity;
    assert(i >= 0 && i < nx);
    assert(j >= 0 && j < ny);
    assert(k >= 0 && k < nz);
    if      (token == "u") getCell(i,j,k)->set_u_plus(velocity);
    else if (token == "v") getCell(i,j,k)->set_v_plus(velocity);
    else if (token == "w") getCell(i,j,k)->set_w_plus(velocity);
    else assert(0);
  }
  SetBoundaryVelocities();
}

// ==============================================================

bool Fluid::inShape(Vec3f &pos, const std::string &shape) {
  // return true if this point is inside the "shape"
  // defined procedurally (using an implicit surface)
  if (shape == "everywhere") {
    return true;
  } else if (shape == "left") {
    // a blob of particles on the lower left (for the dam)
    return (pos.x() < 0.2*nx*dx && pos.y() < 0.5*ny*dy);
  } else if (shape == "drop") {
    // a shallow pool of particles on the bottom
    float h = ny*dy/6.0;
    if (pos.y() < 2*h) return true;
    // and a sphere of particles above
    Vec3f center = Vec3f(nx*dx*0.5, 5*h,nz*dz*0.5);
    float length = (center-pos).Length();
    if (length < 0.8*h) return true;
    return false;
  } else if (shape == "droponly") {
    // and a sphere of particles above
    float h = ny*dy/6.0;
    Vec3f center = Vec3f(nx*dx*0.5, 5*h,nz*dz*0.5);
    float length = (center-pos).Length();
    if (length < 0.8*h) return true;
    return false;
  } else if (shape == "pool") {
    // a shallow pool of particles on the bottom
    float h = ny*dy/6.0;
    if (pos.y() < 2*h) return true;
    return false;
  } else if (shape == "none") {
    return false;
  }
  else {
    std::cout << "unknown shape: " << shape << std::endl;
    exit(0);
  }
}

// ==============================================================

void Fluid::GenerateParticles(const std::string &shape, const std::string &placement) {
  // create a set of points according to the "placement" token,
  // then check whether they are inside of the "shape"
  if (placement == "uniform") {
    int dens = (int)pow(density,0.334);
    std::cout << "density " << dens << " " << density << std::endl;
    assert (dens*dens*dens == density);
    // the uniform grid spacing
    float spacing = 1/float(dens);
    for (float x = 0.5*spacing*dx; x < nx*dx; x += spacing*dx) {
      for (float y = 0.5*spacing*dy; y < ny*dy; y += spacing*dy) {
        for (float z = 0.5*spacing*dz; z < nz*dz; z += spacing*dz) {
          Vec3f pos = Vec3f(x,y,z);
          if (inShape(pos,shape)) {
            Cell *cell = getCell(int(x/dx),int(y/dy),int(z/dz));
            FluidParticle *p = new FluidParticle();
            p->setPosition(pos);
            cell->addParticle(p);
          }
        }
      }
    }
  } else {
    assert (placement == "random");
    // note: we don't necessarily have the same number of particles in each cell
    for (int n = 0; n < nx*ny*nz*density; n++) {
      Vec3f pos = Vec3f(rand()/double(RAND_MAX)*nx*dx,
                        rand()/double(RAND_MAX)*ny*dy,
                        rand()/double(RAND_MAX)*nz*dz);
      if (inShape(pos,shape)) {      
        Cell *cell = getCell(int(pos.x()/dx),int(pos.y()/dy),int(pos.z()/dz));
        FluidParticle *p = new FluidParticle();
        p->setPosition(pos);
        cell->addParticle(p);
      }
    }
  }
}

// ==============================================================
// ==============================================================
// ANIMATION
// ==============================================================
// ==============================================================

void Fluid::Animate() {

  // the animation manager:  this is what gets done each timestep!
  if (source)
    AddNewParticles();
  
  if (sink)
    RemoveParticles();

  ComputeNewVelocities();
  SetBoundaryVelocities();
  
  // compressible / incompressible flow
  if (compressible == false) {
    for (int iters = 0; iters < 20; iters++) {
      float max_divergence = AdjustForIncompressibility();
      SetBoundaryVelocities();
      if (max_divergence < EPSILON) break;
    }
  }

  UpdatePressures();
  CopyVelocities();

  // advanced the particles through the fluid
  MoveParticles();
  ReassignParticles();
  SetEmptySurfaceFull();
}

// ==============================================================

void Fluid::AddNewParticles() {
  if (sourcetype == "rain")
  {
    for (int n = 0; n < nx*sourcedensity; n++) {
      Vec3f pos = Vec3f(rand()/double(RAND_MAX)*nx*dx,
                        ((rand()%(RAND_MAX/4))/double(RAND_MAX)+0.75)*ny*dy,
                        rand()/double(RAND_MAX)*nz*dz);
      Cell *cell = getCell(int(pos.x()/dx),int(pos.y()/dy),int(pos.z()/dz));
      FluidParticle *p = new FluidParticle();
      p->setPosition(pos);
      cell->addParticle(p);
    }
  }
  else if (sourcetype == "spout")
  {
    for (int n = 0; n < nx*ny*nz*sourcedensity; n++) {
      Vec3f pos = Vec3f(rand()/double(RAND_MAX)*2*dx,
                        ((rand()%(RAND_MAX/10))/double(RAND_MAX)+0.9)*ny*dy,
                        rand()/double(RAND_MAX)*dz);
      Cell *cell = getCell(int(pos.x()/dx),int(pos.y()/dy),int(pos.z()/dz));
      FluidParticle *p = new FluidParticle();
      p->setPosition(pos);
      cell->addParticle(p);
    }
  }
  else if (sourcetype == "pour")
  {
    for (int n = 0; n < nx*ny*nz*sourcedensity; n++) {
      Vec3f pos = Vec3f(rand()/double(RAND_MAX)*2*dx+(dx*nx*0.5),
                        ((rand()%(RAND_MAX/10))/double(RAND_MAX)+0.9)*ny*dy,
                        rand()/double(RAND_MAX)*dz);
      Cell *cell = getCell(int(pos.x()/dx),int(pos.y()/dy),int(pos.z()/dz));
      FluidParticle *p = new FluidParticle();
      p->setPosition(pos);
      cell->addParticle(p);
    }
  }
  else if (sourcetype == "fountain")
  {
    for (int n = 0; n < nx*ny*nz*sourcedensity; n++) {
      Vec3f pos = Vec3f(rand()/double(RAND_MAX)*2*dx+(dx*nx*0.5),
                        ((rand()%(RAND_MAX/10))/double(RAND_MAX))*ny*dy,
                        rand()/double(RAND_MAX)*dz);
      Cell *cell = getCell(int(pos.x()/dx),int(pos.y()/dy),int(pos.z()/dz));
      FluidParticle *p = new FluidParticle();
      p->setPosition(pos);
      cell->addParticle(p);
    }
  }
  else
    assert(0);
}

// ==============================================================

void Fluid::RemoveParticles() {
  if (sinktype == "end")
  {
    for (int j = 0; j < ny; j++)
    {
      Cell *cell = getCell(nx-1, j, 0);
      std::vector<FluidParticle*> p = cell->getParticles();
      for (int i = 0; i < int(p.size()); i++)
      {
        cell->removeParticle(p[i]);
      }
    }
  }
  else if (sinktype == "bottom")
  {
    for (int j = 0; j < nx; j++)
    {
      Cell *cell = getCell(j, 0, 0);
      std::vector<FluidParticle*> p = cell->getParticles();
      for (int i = 0; i < int(p.size()); i++)
      {
        cell->removeParticle(p[i]);
      }
    }
  }
  else
    assert(0);
}

// ==============================================================

void Fluid::ComputeNewVelocities() {
  float dt = args->timestep;
  int i,j,k;

  // using the formulas from Foster & Metaxas
  
  for (i = 0; i < nx-1; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        float new_u_plus =
          get_u_plus(i,j,k) +            
          dt * ((1/dx) * (square(get_u_avg(i,j,k)) - square(get_u_avg(i+1,j,k))) +
                (1/dy) * (get_uv_plus(i,j-1,k) - get_uv_plus(i,j,k)) + 
                (1/dz) * (get_uw_plus(i,j,k-1) - get_uw_plus(i,j,k)) +
                args->gravity.x() +
                (1/dx) * (getPressure(i,j,k)-getPressure(i+1,j,k)) +
                (viscosity/square(dx)) * (get_u_plus(i+1,j  ,k  ) - 2*get_u_plus(i,j,k) + get_u_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_u_plus(i  ,j+1,k  ) - 2*get_u_plus(i,j,k) + get_u_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_u_plus(i  ,j  ,k+1) - 2*get_u_plus(i,j,k) + get_u_plus(i  ,j  ,k-1)) );
        cell->set_new_u_plus(new_u_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny-1; j++) {
      for (k = 0; k < nz; k++) {	
        Cell *cell = getCell(i,j,k);
        float new_v_plus =
          get_v_plus(i,j,k) +
          dt * ((1/dx) * (get_uv_plus(i-1,j,k) - get_uv_plus(i,j,k)) +
                (1/dy) * (square(get_v_avg(i,j,k)) - square(get_v_avg(i,j+1,k))) +
                (1/dz) * (get_vw_plus(i,j,k-1) - get_vw_plus(i,j,k)) +
                args->gravity.y() +
                (1/dy) * (getPressure(i,j,k)-getPressure(i,j+1,k)) +
                (viscosity/square(dx)) * (get_v_plus(i+1,j  ,k  ) - 2*get_v_plus(i,j,k) + get_v_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_v_plus(i  ,j+1,k  ) - 2*get_v_plus(i,j,k) + get_v_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_v_plus(i  ,j  ,k+1) - 2*get_v_plus(i,j,k) + get_v_plus(i  ,j  ,k-1)) );
        cell->set_new_v_plus(new_v_plus);
      }
    }
  }

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz-1; k++) {
        Cell *cell = getCell(i,j,k);
        float new_w_plus =
          get_w_plus(i,j,k) +
          dt * ((1/dx) * (get_uw_plus(i-1,j,k) - get_uw_plus(i,j,k)) +
                (1/dy) * (get_vw_plus(i,j-1,k) - get_vw_plus(i,j,k)) +
                (1/dz) * (square(get_w_avg(i,j,k)) - square(get_w_avg(i,j,k+1))) +
                args->gravity.z() +
                (1/dz) * (getPressure(i,j,k)-getPressure(i,j,k+1)) +
                (viscosity/square(dx)) * (get_w_plus(i+1,j  ,k  ) - 2*get_w_plus(i,j,k) + get_w_plus(i-1,j  ,k  )) +
                (viscosity/square(dy)) * (get_w_plus(i  ,j+1,k  ) - 2*get_w_plus(i,j,k) + get_w_plus(i  ,j-1,k  )) +
                (viscosity/square(dz)) * (get_w_plus(i  ,j  ,k+1) - 2*get_w_plus(i,j,k) + get_w_plus(i  ,j  ,k-1)) );
        cell->set_new_w_plus(new_w_plus);
      }
    }
  }
  
  if (source)
  {
    if (sourcetype == "spout")
    {
      i = int((dy*ny*9)/10);
      while (i < ny)
      {
        Cell *cell = getCell(0,i,0);
        cell->adjust_new_u_plus(sourcevelocity);
        i++;
      }
    }
    else if (sourcetype == "fountain")
    {
      i = int(nx / 2);
      i--;
      Cell *cell = getCell(i,0,0);
      cell->adjust_new_v_plus(sourcevelocity);
      i++;
      cell = getCell(i,0,0);
      cell->adjust_new_v_plus(sourcevelocity);
      i++;
      cell = getCell(i,0,0);
      cell->adjust_new_v_plus(sourcevelocity);
    }
  }
}

// ==============================================================

void Fluid::SetBoundaryVelocities() {

  // zero out flow perpendicular to the boundaries (no sources or sinks)
  for (int j = -1; j <= ny; j++) {
    for (int k = -1; k <= nz; k++) {
      getCell(-1  ,j,k)->set_u_plus(0);
      getCell(nx-1,j,k)->set_u_plus(0);
      getCell(nx  ,j,k)->set_u_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int k = -1; k <= nz; k++) {
      getCell(i,-1  ,k)->set_v_plus(0);
      getCell(i,ny-1,k)->set_v_plus(0);
      getCell(i,ny  ,k)->set_v_plus(0);
    }
  }
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i,j,-1  )->set_w_plus(0);
      getCell(i,j,nz-1)->set_w_plus(0);
      getCell(i,j,nz  )->set_w_plus(0);
    }
  }

  // free slip or no slip boundaries (friction with boundary)
  float xy_sign = (xy_free_slip) ? 1 : -1;
  float yz_sign = (yz_free_slip) ? 1 : -1;
  float zx_sign = (zx_free_slip) ? 1 : -1;
  for (int i = 0; i < nx; i++) {
    for (int j = -1; j <= ny; j++) {
      getCell(i,j,-1)->set_u_plus(xy_sign*getCell(i,j,0)->get_u_plus());
      getCell(i,j,nz)->set_u_plus(xy_sign*getCell(i,j,nz-1)->get_u_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(i,-1,k)->set_u_plus(zx_sign*getCell(i,0,k)->get_u_plus());
      getCell(i,ny,k)->set_u_plus(zx_sign*getCell(i,ny-1,k)->get_u_plus());
    }
  }
  for (int j = 0; j < ny; j++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i,j,-1)->set_v_plus(xy_sign*getCell(i,j,0)->get_v_plus());
      getCell(i,j,nz)->set_v_plus(xy_sign*getCell(i,j,nz-1)->get_v_plus());
    }
    for (int k = -1; k <= nz; k++) {
      getCell(-1,j,k)->set_v_plus(yz_sign*getCell(0,j,k)->get_v_plus());
      getCell(nx,j,k)->set_v_plus(yz_sign*getCell(nx-1,j,k)->get_v_plus());
    }
  }
  for (int k = 0; k < nz; k++) {
    for (int i = -1; i <= nx; i++) {
      getCell(i,-1,k)->set_w_plus(zx_sign*getCell(i,0,k)->get_w_plus());
      getCell(i,ny,k)->set_w_plus(zx_sign*getCell(i,ny-1,k)->get_w_plus());
    }
    for (int j = -1; j <= ny; j++) {
      getCell(-1,j,k)->set_w_plus(yz_sign*getCell(0,j,k)->get_w_plus());
      getCell(nx,j,k)->set_w_plus(yz_sign*getCell(nx-1,j,k)->get_w_plus());
    }
  }
}

// ==============================================================

void Fluid::CopyVelocities() {
  float dt = args->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *c = getCell(i,j,k);
        c->copyVelocity();
        if (fabs(c->get_u_plus()) > 0.5*dx/dt ||
            fabs(c->get_v_plus()) > 0.5*dy/dt ||
            fabs(c->get_w_plus()) > 0.5*dz/dt) {
          // velocity has exceeded reasonable threshhold
          std::cout << "velocity has exceeded reasonable threshhold, stopping animation" << std::endl;
          args->animate=false;
        }
      }
    }
  }
}

// ==============================================================

float Fluid::AdjustForIncompressibility() {


  // ASSIGNMENT: 
  // adjust the face velocities of each cell so that each cell has
  // zero divergence
  
  float maxd = 0;
  
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        float divergence = 
          - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
        (1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
        (1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );
        if (fabs(divergence) > maxd)
          maxd = fabs(divergence);
      }
    }
  }
  
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        if (getCell(i,j,k)->getStatus() == CELL_EMPTY)
        {
          set_new_u_plus(i,j,k,0);
          set_new_u_plus(i-1,j,k,0);
          set_new_v_plus(i,j,k,0);
          set_new_v_plus(i,j-1,k,0);
          getCell(i,j,k)->setPressure(0);
        }
      }
    }
  }

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nz; k++)
      {
        if (getCell(i,j,k)->getStatus() == CELL_FULL)
        {
          float divergence = 
            - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
          (1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
          (1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );
          float dt = args->timestep;
          float beta = BETA_0/((2*dt) * (1/square(dx) + 1/square(dy) + 1/square(dz)));
          float dp = beta*divergence;
          set_new_u_plus(i,j,k,get_new_u_plus(i,j,k) + (dt/dx) * dp);
          set_new_u_plus(i-1,j,k,get_new_u_plus(i-1,j,k) - (dt/dx) * dp);
          set_new_v_plus(i,j,k,get_new_v_plus(i,j,k) + (dt/dy) * dp);
          set_new_v_plus(i,j-1,k,get_new_v_plus(i,j-1,k) - (dt/dy) * dp);
          set_new_w_plus(i,j,k,get_new_w_plus(i,j,k) + (dt/dz) * dp);
          set_new_w_plus(i,j,k-1,get_new_w_plus(i,j,k-1) - (dt/dz) * dp);
          getCell(i,j,k)->setPressure(0);
        }
        else if (getCell(i,j,k)->getStatus() == CELL_SURFACE)
        {
          if (getCell(i+1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() != CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() != CELL_EMPTY)
          {
            set_new_u_plus(i,j,k,get_new_u_plus(i-1,j,k) - (dx/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)));
          }
          else if (getCell(i+1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() != CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() != CELL_EMPTY)
          {
            set_new_u_plus(i-1,j,k,get_new_u_plus(i,j,k) + (dx/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)));
          }
          else if (getCell(i+1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() == CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() != CELL_EMPTY)
          {
            set_new_v_plus(i,j,k,get_new_v_plus(i,j-1,k) - (dy/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)));
          }
          else if (getCell(i+1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() != CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() == CELL_EMPTY)
          {
            set_new_v_plus(i,j-1,k,get_new_v_plus(i,j,k) + (dy/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)));
          }
          else if (getCell(i+1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() == CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() != CELL_EMPTY)
          {
            set_new_u_plus(i,j,k,get_new_u_plus(i-1,j,k));
            set_new_v_plus(i,j,k,get_new_v_plus(i,j-1,k));
          }
          else if (getCell(i+1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() == CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() != CELL_EMPTY)
          {
            set_new_u_plus(i-1,j,k,get_new_u_plus(i,j,k));
            set_new_v_plus(i,j,k,get_new_v_plus(i,j-1,k));
          }
          else if (getCell(i+1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() != CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() == CELL_EMPTY)
          {
            set_new_u_plus(i-1,j,k,get_new_u_plus(i,j,k));
            set_new_v_plus(i,j-1,k,get_new_v_plus(i,j,k));
          }
          else if (getCell(i+1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() != CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() == CELL_EMPTY)
          {
            set_new_u_plus(i,j,k,get_new_u_plus(i-1,j,k));
            set_new_v_plus(i,j-1,k,get_new_v_plus(i,j,k));
          }
          else if (getCell(i+1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() == CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() == CELL_EMPTY)
          {
            float splitVelocity = (get_new_u_plus(i,j,k) + get_new_u_plus(i-1,j,k))/2;
            set_new_v_plus(i,j,k,splitVelocity);
            set_new_v_plus(i,j-1,k,splitVelocity);
          }
          else if (getCell(i+1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() != CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() != CELL_EMPTY)
          {
            float splitVelocity = (get_new_v_plus(i,j,k) + get_new_v_plus(i,j-1,k))/2;
            set_new_u_plus(i,j,k,splitVelocity);
            set_new_u_plus(i-1,j,k,splitVelocity);
          }
          else if (getCell(i+1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() == CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() == CELL_EMPTY)
          {
            set_new_u_plus(i,j,k,get_new_u_plus(i-1,j,k));
            //set_new_v_plus(i,j,k,0);
            //set_new_v_plus(i,j-1,k,0);
          }
          else if (getCell(i+1,j,k)->getStatus() != CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() == CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() == CELL_EMPTY)
          {
            set_new_u_plus(i-1,j,k,get_new_u_plus(i,j,k));
            //set_new_v_plus(i,j,k,0);
            //set_new_v_plus(i,j-1,k,0);
          }
          else if (getCell(i+1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() == CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() != CELL_EMPTY)
          {
            set_new_v_plus(i,j,k,get_new_v_plus(i,j-1,k));
            //set_new_u_plus(i,j,k,0);
            //set_new_u_plus(i-1,j,k,0);
          }
          else if (getCell(i+1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() != CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() == CELL_EMPTY)
          {
            set_new_v_plus(i,j-1,k,get_new_v_plus(i,j,k));
            //set_new_u_plus(i,j,k,0);
            //set_new_u_plus(i-1,j,k,0);
          }
          else if (getCell(i+1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j+1,k)->getStatus() == CELL_EMPTY &&
              getCell(i-1,j,k)->getStatus() == CELL_EMPTY &&
              getCell(i,j-1,k)->getStatus() == CELL_EMPTY)
          {
            /*
            set_new_v_plus(i,j,k,0);
            set_new_v_plus(i,j-1,k,0);
            set_new_u_plus(i,j,k,0);
            set_new_u_plus(i-1,j,k,0);
            */
            //adjust_new_v_plus(i,j-1,k,args->gravity.y());
          }
          getCell(i,j,k)->setPressure(0);
        }
      }
    }
  }
  
  // should the maximum divergence in any cell before adjustment
  return maxd;
}

// ==============================================================

void Fluid::UpdatePressures() {
  for (int i = -1; i <= nx; i++) {
    for (int j = -1; j <= ny; j++) {
      for (int k = -1; k <= nz; k++) {
	Cell *c = getCell(i,j,k);
	if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz) {
	  // compute divergence and increment/decrement pressure
	  float pressure = c->getPressure();
	  float divergence = 
	    - ( (1/dx) * (get_new_u_plus(i,j,k) - get_new_u_plus(i-1,j,k)) +
		(1/dy) * (get_new_v_plus(i,j,k) - get_new_v_plus(i,j-1,k)) +
		(1/dz) * (get_new_w_plus(i,j,k) - get_new_w_plus(i,j,k-1)) );
	  float dt = args->timestep;
	  float beta = BETA_0/((2*dt) * (1/square(dx) + 1/square(dy) + 1/square(dz)));
	  float dp = beta*divergence;
	  c->setPressure(pressure + dp);
	} else {
	  // zero out boundary cells (just in case)
	  c->setPressure(0);
	}

	// also, zero out empty cells (from Foster 2001 paper)
	if (c->getStatus() == CELL_EMPTY) {
	  c->setPressure(0);
	}

      }
    }
  }
}

// ==============================================================

void Fluid::MoveParticles() {
  float dt = args->timestep;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
	std::vector<FluidParticle*> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          Vec3f pos = p->getPosition();
          Vec3f vel = getInterpolatedVelocity(pos);
          Vec3f pos2 = pos + vel*dt;
#if 0
          // euler integration
          p->setPosition(pos2);
#else
          // trapezoid integration
          Vec3f vel2 = getInterpolatedVelocity(pos2);
          Vec3f pos3 = pos + 0.5*(vel+vel2)*dt;
          p->setPosition(pos3);
#endif
        }
      }
    }
  }
}

// ==============================================================

void Fluid::ReassignParticles() {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
	std::vector<FluidParticle*> &particles = cell->getParticles();
        for (unsigned int iter = 0; iter < particles.size(); iter++) {
          FluidParticle *p = particles[iter];
          Vec3f pos = p->getPosition();
          int i2 = (int)min2((nx-1),max2(0,floor(pos.x()/dx)));
          int j2 = (int)min2((ny-1),max2(0,floor(pos.y()/dy)));
          int k2 = (int)min2((nz-1),max2(0,floor(pos.z()/dz)));
          // if the particle has crossed one of the cell faces 
          // assign it to the new cell
          if (i != i2 || j != j2 || k != k2) {
            cell->removeParticle(p);
            getCell(i2,j2,k2)->addParticle(p);
          } 
        }
      }
    }
  }
}

// ==============================================================

void Fluid::SetEmptySurfaceFull() {
  int i,j,k;
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        if (cell->numParticles() == 0)
	  cell->setStatus(CELL_EMPTY);
        else 
          cell->setStatus(CELL_FULL);
      }
    }
  }

  // pick out the boundary cells
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
        Cell *cell = getCell(i,j,k);
        if (cell->getStatus() == CELL_FULL &&
            (getCell(i-1,j,k)->getStatus() == CELL_EMPTY ||
             getCell(i+1,j,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j-1,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j+1,k)->getStatus() == CELL_EMPTY ||
             getCell(i,j,k-1)->getStatus() == CELL_EMPTY ||
             getCell(i,j,k+1)->getStatus() == CELL_EMPTY)) {
          cell->setStatus(CELL_SURFACE);
        }
      }
    }
  }
}

// ==============================================================

Vec3f Fluid::getInterpolatedVelocity(const Vec3f &pos) const {

  // Get the cell the given particle is in
  int i = int(floor(pos.x()/dx)); if (i < 0) i = 0; if (i >= nx) i = nx-1;
  int j = int(floor(pos.y()/dy)); if (j < 0) j = 0; if (j >= ny) j = ny-1;
  int k = int(floor(pos.z()/dz)); if (k < 0) k = 0; if (k >= nz) k = nz-1;
  
  // Get the relative position of the particle in the cell
  float a = (pos.x() - dx*i)/dx;
  float b = (pos.y() - dy*j)/dy;
  float c = (pos.z() - dz*k)/dz;
  
  // Stores the 8 velocities affecting the particle
  float v1 = get_u_plus(i,j,k);
  float v2 = get_u_plus(i-1,j,k);
  float v3 = 0;
  float v4 = 0;
  float v5 = 0;
  float v6 = 0;
  float v7 = 0;
  float v8 = 0;
  
  // Stores the weights associated with each of the 8 velocities
  float w1 = 0;
  float w2 = 0;
  float w3 = 0;
  float w4 = 0;
  float w5 = 0;
  float w6 = 0;
  float w7 = 0;
  float w8 = 0;
  
  // Find the velocities and their respective weights in the x direction
  
  // This is done by essentially splitting the cell into quadrants, and which quadrant
  // the particle is in dictates which velocities are used and the equations to calculate
  // their weights.
  if (b > 0.5)
  {
    v3 = get_u_plus(i,j+1,k);
    v4 = get_u_plus(i-1,j+1,k);
    if (c > 0.5)
    {
      v5 = get_u_plus(i,j,k+1);
      v6 = get_u_plus(i-1,j,k+1);
      v7 = get_u_plus(i,j+1,k+1);
      v8 = get_u_plus(i-1,j+1,k+1);
      w1 = a * b * c;
      w2 = (1-a) * b * c;
      w3 = a * (b-0.5) * c;
      w4 = (1-a) * (b-0.5) * c;
      w5 = a * b * (c-0.5);
      w6 = (1-a) * b * (c-0.5);
      w7 = a * (b-0.5) * (c-0.5);
      w8 = (1-a) * (b-0.5) * (c-0.5);
    }
    else
    {
      v5 = get_u_plus(i,j,k-1);
      v6 = get_u_plus(i-1,j,k-1);
      v7 = get_u_plus(i,j+1,k-1);
      v8 = get_u_plus(i-1,j+1,k-1);
      w1 = a * b * (c+0.5);
      w2 = (1-a) * b * (c+0.5);
      w3 = a * (b-0.5) * (c+0.5);
      w4 = (1-a) * (b-0.5) * (c+0.5);
      w5 = a * b * c;
      w6 = (1-a) * b * c;
      w7 = a * (b-0.5) * c;
      w8 = (1-a) * (b-0.5) * c;
    }
  }
  else
  {
    v3 = get_u_plus(i,j-1,k);
    v4 = get_u_plus(i-1,j-1,k);
    if (c > 0.5)
    {
      v5 = get_u_plus(i,j,k+1);
      v6 = get_u_plus(i-1,j,k+1);
      v7 = get_u_plus(i,j-1,k+1);
      v8 = get_u_plus(i-1,j-1,k+1);
      w1 = a * (b+0.5) * c;
      w2 = (1-a) * (b+0.5) * c;
      w3 = a * b * c;
      w4 = (1-a) * b * c;
      w5 = a * (b+0.5) * (c-0.5);
      w6 = (1-a) * (b+0.5) * (c-0.5);
      w7 = a * b * (c-0.5);
      w8 = (1-a) * b * (c-0.5);
    }
    else
    {
      v5 = get_u_plus(i,j,k-1);
      v6 = get_u_plus(i-1,j,k-1);
      v7 = get_u_plus(i,j-1,k-1);
      v8 = get_u_plus(i-1,j-1,k-1);
      w1 = a * (b+0.5) * (c+0.5);
      w2 = (1-a) * (b+0.5) * (c+0.5);
      w3 = a * b * (c+0.5);
      w4 = (1-a) * b * (c+0.5);
      w5 = a * (b+0.5) * c;
      w6 = (1-a) * (b+0.5) * c;
      w7 = a * b * c;
      w8 = (1-a) * b * c;
    }
  }
  
  // Calculate the weighted average in the x direction
  float xavg = (w1*v1+w2*v2+w3*v3+w4*v4+w5*v5+w6*v6+w7*v7+w8*v8)
    /(w1+w2+w3+w4+w5+w6+w7+w8);
  
  // Reset all velocity and weight values for y calculation
  v1 = get_v_plus(i,j,k);
  v2 = get_v_plus(i,j-1,k);
  v3 = 0;
  v4 = 0;
  v5 = 0;
  v6 = 0;
  v7 = 0;
  v8 = 0;
  
  w1 = 0;
  w2 = 0;
  w3 = 0;
  w4 = 0;
  w5 = 0;
  w6 = 0;
  w7 = 0;
  w8 = 0;
  
  // Find the velocities and their respective weights in the y direction
  if (a > 0.5)
  {
    v3 = get_v_plus(i+1,j,k);
    v4 = get_v_plus(i+1,j-1,k);
    if (c > 0.5)
    {
      v5 = get_v_plus(i,j,k+1);
      v6 = get_v_plus(i,j-1,k+1);
      v7 = get_v_plus(i+1,j,k+1);
      v8 = get_v_plus(i+1,j-1,k+1);
      w1 = a * b * c;
      w2 = a * (1-b) * c;
      w3 = (a-0.5) * b * c;
      w4 = (a-0.5) * (1-b) * c;
      w5 = a * b * (c-0.5);
      w6 = a * (1-b) * (c-0.5);
      w7 = (a-0.5) * b * (c-0.5);
      w8 = (a-0.5) * (1-b) * (c-0.5);
    }
    else
    {
      v5 = get_v_plus(i,j,k-1);
      v6 = get_v_plus(i,j-1,k-1);
      v7 = get_v_plus(i+1,j,k-1);
      v8 = get_v_plus(i+1,j-1,k-1);
      w1 = a * b * (c+0.5);
      w2 = a * (1-b) * (c+0.5);
      w3 = (a-0.5) * b * (c+0.5);
      w4 = (a-0.5) * (1-b) * (c+0.5);
      w5 = a * b * c;
      w6 = a * (1-b) * c;
      w7 = (a-0.5) * b * c;
      w8 = (a-0.5) * (1-b) * c;
    }
  }
  else
  {
    v3 = get_v_plus(i-1,j,k);
    v4 = get_v_plus(i-1,j-1,k);
    if (c > 0.5)
    {
      v5 = get_v_plus(i,j,k+1);
      v6 = get_v_plus(i,j-1,k+1);
      v7 = get_v_plus(i-1,j,k+1);
      v8 = get_v_plus(i-1,j-1,k+1);
      w1 = (a+0.5) * b * c;
      w2 = (a+0.5) * (1-b) * c;
      w3 = a * b * c;
      w4 = a * (1-b) * c;
      w5 = (a+0.5) * b * (c-0.5);
      w6 = (a+0.5) * (1-b) * (c-0.5);
      w7 = a * b * (c-0.5);
      w8 = a * (1-b) * (c-0.5);
    }
    else
    {
      v5 = get_v_plus(i,j,k-1);
      v6 = get_v_plus(i,j-1,k-1);
      v7 = get_v_plus(i-1,j,k-1);
      v8 = get_v_plus(i-1,j-1,k-1);
      w1 = (a+0.5) * b * (c+0.5);
      w2 = (a+0.5) * (1-b) * (c+0.5);
      w3 = a * b * (c+0.5);
      w4 = a * (1-b) * (c+0.5);
      w5 = (a+0.5) * b * c;
      w6 = (a+0.5) * (1-b) * c;
      w7 = a * b * c;
      w8 = a * (1-b) * c;
    }
  }
  
  // Calculate the weighted average in the y direction
  float yavg = (w1*v1+w2*v2+w3*v3+w4*v4+w5*v5+w6*v6+w7*v7+w8*v8)
    /(w1+w2+w3+w4+w5+w6+w7+w8);
  
  // Reset all velocity and weight values for z calculation
  v1 = get_w_plus(i,j,k);
  v2 = get_w_plus(i,j,k-1);
  v3 = 0;
  v4 = 0;
  v5 = 0;
  v6 = 0;
  v7 = 0;
  v8 = 0;
  
  w1 = 0;
  w2 = 0;
  w3 = 0;
  w4 = 0;
  w5 = 0;
  w6 = 0;
  w7 = 0;
  w8 = 0;
  
    // Find the velocities and their respective weights in the z direction
  if (a > 0.5)
  {
    v3 = get_w_plus(i+1,j,k);
    v4 = get_w_plus(i+1,j,k-1);
    if (b > 0.5)
    {
      v5 = get_w_plus(i,j+1,k);
      v6 = get_w_plus(i,j+1,k-1);
      v7 = get_w_plus(i+1,j+1,k);
      v8 = get_w_plus(i+1,j+1,k-1);
      w1 = a * b * c;
      w2 = a * b * (1-c);
      w3 = (a-0.5) * b * c;
      w4 = (a-0.5) * b * (1-c);
      w5 = a * (b-0.5) * c;
      w6 = a * (b-0.5) * (1-c);
      w7 = (a-0.5) * (b-0.5) * c;
      w8 = (a-0.5) * (b-0.5) * (1-c);
    }
    else
    {
      v5 = get_w_plus(i,j-1,k);
      v6 = get_w_plus(i,j-1,k-1);
      v7 = get_w_plus(i+1,j-1,k);
      v8 = get_w_plus(i+1,j-1,k-1);
      w1 = a * (b+0.5) * c;
      w2 = a * (b+0.5) * (1-c);
      w3 = (a-0.5) * (b+0.5) * c;
      w4 = (a-0.5) * (b+0.5) * (1-c);
      w5 = a * b * c;
      w6 = a * b * (1-c);
      w7 = (a-0.5) * b * c;
      w8 = (a-0.5) * b * (1-c);
    }
  }
  else
  {
    v3 = get_w_plus(i-1,j,k);
    v4 = get_w_plus(i-1,j,k-1);
    if (b > 0.5)
    {
      v5 = get_w_plus(i,j+1,k);
      v6 = get_w_plus(i,j+1,k-1);
      v7 = get_w_plus(i-1,j+1,k);
      v8 = get_w_plus(i-1,j+1,k-1);
      w1 = (a+0.5) * b * c;
      w2 = (a+0.5) * b * (1-c);
      w3 = a * b * c;
      w4 = a * b * (1-c);
      w5 = (a+0.5) * (b-0.5) * c;
      w6 = (a+0.5) * (b-0.5) * (1-c);
      w7 = a * (b-0.5) * c;
      w8 = a * (b-0.5) * (1-c);
    }
    else
    {
      v5 = get_w_plus(i,j-1,k);
      v6 = get_w_plus(i,j-1,k-1);
      v7 = get_w_plus(i-1,j-1,k);
      v8 = get_w_plus(i-1,j-1,k-1);
      w1 = (a+0.5) * (b+0.5) * c;
      w2 = (a+0.5) * (b+0.5) * (1-c);
      w3 = a * (b+0.5) * c;
      w4 = a * (b+0.5) * (1-c);
      w5 = (a+0.5) * b * c;
      w6 = (a+0.5) * b * (1-c);
      w7 = a * b * c;
      w8 = a * b * (1-c);
    }
  }

  // Calculate the weighted average in the z direction  
  float zavg = (w1*v1+w2*v2+w3*v3+w4*v4+w5*v5+w6*v6+w7*v7+w8*v8)
    /(w1+w2+w3+w4+w5+w6+w7+w8);
  
  return Vec3f(xavg,yavg,zavg);
}


// ==============================================================
// ==============================================================
// RENDERING
// ==============================================================
// ==============================================================

void Fluid::Paint() const {

  // paint the boundaries of the volume
  if (args->bounding_box) {
    BoundingBox box = getBoundingBox();
    box.Paint();
  }

  // =====================================================================================
  // render the particles
  // =====================================================================================
  if (args->particles) 
  {
    glDisable(GL_LIGHTING);
    glColor3f(0,0,0);
    glPointSize(3);
    glBegin(GL_POINTS);
    glColor3f(0,0,.7);
    for (int x = 0; x < nx; x++) 
    {
      for (int y = 0; y < ny; y++) 
      {
        for (int z = 0; z < nz; z++) 
        {
          Cell *cell = getCell(x,y,z);
          std::vector<FluidParticle*> &particles = cell->getParticles();
          for (unsigned int iter = 0; iter < particles.size(); iter++) 
          {
            FluidParticle *p = particles[iter];
            Vec3f v = p->getPosition();
            glVertex3f(v.x(),v.y(),v.z());
          }
        } 
      }
    }
    glEnd();
    glEnable(GL_LIGHTING);
  } 

  // =====================================================================================
  // render stubby triangles to visualize the u, v, and w velocities between cell faces
  // =====================================================================================
  if (args->edge_velocity > 0) {
    glDisable(GL_LIGHTING);
    glLineWidth(3);
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          float u = get_u_plus(i,j,k);
          float v = get_v_plus(i,j,k);
          float w = get_w_plus(i,j,k);
          float x = i*dx;
          float y = j*dy;
          float z = k*dz;
          if (args->edge_velocity == 1) {
            if (fabs(u) > 0) {
              glColor3f(1,0,0);
              glBegin(GL_TRIANGLE_FAN);
              glVertex3f(x+dx+u,y+0.5*dy,z+0.5*dz);
              glVertex3f(x+dx,y+0.55*dy,z+0.55*dz);
              glVertex3f(x+dx,y+0.55*dy,z+0.45*dz);
              glVertex3f(x+dx,y+0.45*dy,z+0.45*dz);
              glVertex3f(x+dx,y+0.45*dy,z+0.55*dz);
              glVertex3f(x+dx,y+0.55*dy,z+0.55*dz);
              glEnd();
            }
          }
          if (args->edge_velocity == 2) {
            if (fabs(v) > 0) {
              glColor3f(0,1,0);
              glBegin(GL_TRIANGLE_FAN);
              glVertex3f(x+0.5*dx,y+dy+v,z+0.5*dz);
              glVertex3f(x+0.55*dx,y+dy,z+0.55*dz);
              glVertex3f(x+0.55*dx,y+dy,z+0.45*dz);
              glVertex3f(x+0.45*dx,y+dy,z+0.45*dz);
              glVertex3f(x+0.45*dx,y+dy,z+0.55*dz);
              glVertex3f(x+0.55*dx,y+dy,z+0.55*dz);
              glEnd();
            }
          }
          if (args->edge_velocity == 3) {
            if (fabs(w) > 0) {
              glColor3f(0,0,1);
              glBegin(GL_TRIANGLE_FAN);
              glVertex3f(x+0.5*dx,y+0.5*dy,z+dz+w);
              glVertex3f(x+0.55*dx,y+0.55*dy,z+dz);
              glVertex3f(x+0.55*dx,y+0.45*dy,z+dz);
              glVertex3f(x+0.45*dx,y+0.45*dy,z+dz);
              glVertex3f(x+0.45*dx,y+0.55*dy,z+dz);
              glVertex3f(x+0.55*dx,y+0.55*dy,z+dz);
              glEnd();
            }
          }
        }
      }
    }

    glEnable(GL_LIGHTING);
  }

  // =====================================================================================
  // visualize the average cell velocity
  // =====================================================================================
  if (args->velocity) {
    glDisable(GL_LIGHTING);
    glLineWidth(3);
    glBegin(GL_LINES);
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          glColor3f(1,0,0);
          glVertex3f((i+0.5)*dx,(j+0.5)*dy,(k+0.5)*dz);
          glColor3f(1,1,1);
          glVertex3f((i+0.5)*dx+get_u_avg(i,j,k),
                     (j+0.5)*dy+get_v_avg(i,j,k),
                     (k+0.5)*dz+get_w_avg(i,j,k));
        } 
      }
    }
    glEnd();
    glEnable(GL_LIGHTING);
  } 

  // =====================================================================================
  // visualize the average cell velocity densely
  // =====================================================================================
  if (args->dense_velocity > 0) {
    glDisable(GL_LIGHTING);
    glLineWidth(3);
    glBegin(GL_LINES);
    if (args->dense_velocity == 1) {
      float z = nz*dz / 2.0;
      for (float x = 0; x <= (nx+0.01)*dx; x+=0.25*dx) {
        for (float y = 0; y <= (ny+0.01)*dy; y+=0.25*dy) {
          glColor3f(1,0,0);
          glVertex3f(x,y,z);
          Vec3f vel = getInterpolatedVelocity(Vec3f(x,y,z));
          glColor3f(1,1,1);
          glVertex3f(x+vel.x(),y+vel.y(),z+vel.z());
        }
      } 
    }
    if (args->dense_velocity == 2) {
      float y = ny*dy / 2.0;
      for (float x = 0; x <= (nx+0.01)*dx; x+=0.25*dx) {
        for (float z = 0; z <= (nz+0.01)*dz; z+=0.25*dz) {
          glColor3f(1,0,0);
          glVertex3f(x,y,z);
          Vec3f vel = getInterpolatedVelocity(Vec3f(x,y,z));
          glColor3f(1,1,1);
          glVertex3f(x+vel.x(),y+vel.y(),z+vel.z());
        }
      } 
    }
    if (args->dense_velocity == 3) {
      float x = nx*dx / 2.0;
      for (float y = 0; y <= (ny+0.01)*dy; y+=0.25*dy) {
        for (float z = 0; z <= (nz+0.01)*dz; z+=0.25*dz) {
          glColor3f(1,0,0);
          glVertex3f(x,y,z);
          Vec3f vel = getInterpolatedVelocity(Vec3f(x,y,z));
          glColor3f(1,1,1);
          glVertex3f(x+vel.x(),y+vel.y(),z+vel.z());
        }
      } 
    }
    glEnd();
    glEnable(GL_LIGHTING);
  } 

  // =====================================================================================
  // render a marching cubes representation of the surface
  //    note: make sure you set the Fluid::getIsovalue() to what you want
  // =====================================================================================
  if (args->surface) {
    for (int i = 0; i <= nx; i++) {
      for (int j = 0; j <= ny; j++) {
	for (int k = 0; k <= nz; k++) {
	  marchingCubes->set(i,j,k,interpolateIsovalue(Vec3f((i-0.5),(j-0.5),(k-0.5))));
	} 
      }
    }
    marchingCubes->Paint(args->isosurface);
  } 

  // =====================================================================================
  // render the MAC cells (FULL, SURFACE, or EMPTY)
  // =====================================================================================
  if (args->cubes) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          glPushMatrix();
          Matrix m;
          m.SetToIdentity();
          m *= Matrix::MakeTranslation(Vec3f((i+0.1)*dx,(j+0.1)*dy,(k+0.1)*dz));
          m *= Matrix::MakeScale(0.8);
          m *= Matrix::MakeScale(Vec3f(dx,dy,dz));
          glMultMatrixf(m.glGet());
          Cell *cell = getCell(i,j,k);
          if (cell->getStatus() == CELL_BOUNDARY) {
            glColor3f(0,0,0);
            CubePaint(); 
          } else if (cell->getStatus() == CELL_FULL) {
            glColor3f(1,0,0);
            CubePaint(); 
          } else if (cell->getStatus() == CELL_SURFACE) {
            glColor3f(0,0,1);
            CubePaint(); 
          }
          glPopMatrix();
        } 
      }
    }
  }
  
  if (args->voxels) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          glPushMatrix();
          Matrix m;
          m.SetToIdentity();
          m *= Matrix::MakeTranslation(Vec3f((i+0.1)*dx,(j+0.1)*dy,(k+0.1)*dz));
          m *= Matrix::MakeScale(0.8);
          m *= Matrix::MakeScale(Vec3f(dx,dy,dz));
          glMultMatrixf(m.glGet());
          Cell *cell = getCell(i,j,k);
          glColor3f(0,0,0);
          ClearCubePaint(); 
          glPopMatrix();
        } 
      }
    }
  }
  
  // =====================================================================================
  // render a visualization of the pressure
  // =====================================================================================
  if (args->pressure) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            glPushMatrix();
            Matrix m;
            m.SetToIdentity();
            m *= Matrix::MakeTranslation(Vec3f((i+0.1)*dx,(j+0.1)*dy,(k+0.1)*dz));
            m *= Matrix::MakeScale(0.8);
            m *= Matrix::MakeScale(Vec3f(dx,dy,dz));
            glMultMatrixf(m.glGet());
            float p = getCell(i,j,k)->getPressure();
            p *= 0.1;
            if (p > 1) p = 1;
            if (p < -1) p = -1;
            assert(p >= -1 && p <= 1);
            if (p < 0) {
              glColor3f(1+p,1+p,1);
            } else {
              glColor3f(1,1-p,1-p);
            }
            CubePaint(); 
            glPopMatrix();
        } 
      }
    }
  }
}

// ==============================================================

float Fluid::getIsovalue(int i, int j, int k) const {
  i = max2(0,(min2(i,nx-1)));
  j = max2(0,(min2(j,ny-1)));
  k = max2(0,(min2(k,nz-1)));
  Cell *c = getCell(i,j,k);
  if (c->getStatus() == CELL_EMPTY) return 0;
  // note: this is technically not a correct thing to do
  //       the number of particles is not an indication of it's "fullness"
  if (c->getStatus() == CELL_SURFACE) return 0.5 + c->numParticles()/float(density);
  if (c->getStatus() == CELL_FULL) return 2;
  assert(0);
  return 0;
}

// ==============================================================

float Fluid::interpolateIsovalue(const Vec3f &v) const {

  float x = v.x();
  float y = v.y();
  float z = v.z();

  // get the values at the corners
  float a = getIsovalue(int(floor(x)),int(floor(y)),int(floor(z)));
  float b = getIsovalue(int(floor(x)),int(floor(y)),int( ceil(z)));
  float c = getIsovalue(int(floor(x)),int( ceil(y)),int(floor(z)));
  float d = getIsovalue(int(floor(x)),int( ceil(y)),int( ceil(z)));
  float e = getIsovalue(int( ceil(x)),int(floor(y)),int(floor(z)));
  float f = getIsovalue(int( ceil(x)),int(floor(y)),int( ceil(z)));
  float g = getIsovalue(int( ceil(x)),int( ceil(y)),int(floor(z)));
  float h = getIsovalue(int( ceil(x)),int( ceil(y)),int( ceil(z)));

  float x_frac = x - (floor(x));
  float y_frac = y - (floor(y));
  float z_frac = z - (floor(z));

  assert (x_frac >= 0 && x_frac <= 1);
  assert (y_frac >= 0 && y_frac <= 1);
  assert (z_frac >= 0 && z_frac <= 1);
  
  float answer = triInterpolate(x_frac,y_frac,z_frac,a,b,c,d,e,f,g,h);
  return answer;
}

// ==============================================================

void CubePaint() {

  glBegin (GL_QUADS);

  // front 
  glNormal3f(0,0,1.0);
  glVertex3f(0.0,0.0,1.0);
  glVertex3f(1.0,0.0,1.0);
  glVertex3f(1.0,1.0,1.0);
  glVertex3f(0.0,1.0,1.0);

  // top 
  glNormal3f(0,1.0,0);
  glVertex3f(0.0,1.0,0.0);
  glVertex3f(0.0,1.0,1.0);
  glVertex3f(1.0,1.0,1.0);
  glVertex3f(1.0,1.0,0.0);

  // right
  glNormal3f(1.0,0,0);
  glVertex3f(1.0,0.0,0.0);
  glVertex3f(1.0,1.0,0.0);
  glVertex3f(1.0,1.0,1.0);
  glVertex3f(1.0,0.0,1.0);

  // back 
  glNormal3f(0,0,-1.0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,1.0,0.0);
  glVertex3f(1.0,1.0,0.0);
  glVertex3f(1.0,0.0,0.0);

  // bottom
  glNormal3f(0,-1.0,0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(1.0,0.0,0.0);
  glVertex3f(1.0,0.0,1.0);
  glVertex3f(0.0,0.0,1.0);

  // left 
  glNormal3f(-1.0,0,0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,0.0,1.0);
  glVertex3f(0.0,1.0,1.0);
  glVertex3f(0.0,1.0,0.0);

  glEnd();
}

// ==============================================================

void ClearCubePaint() {

  glBegin (GL_LINES);

  // front 
  glNormal3f(0,0,1.0);
  glVertex3f(0.0,0.0,1.0);
  glVertex3f(1.0,0.0,1.0);
  glVertex3f(1.0,1.0,1.0);
  glVertex3f(0.0,1.0,1.0);

  // top 
  glNormal3f(0,1.0,0);
  glVertex3f(0.0,1.0,0.0);
  glVertex3f(0.0,1.0,1.0);
  glVertex3f(1.0,1.0,1.0);
  glVertex3f(1.0,1.0,0.0);

  // right
  glNormal3f(1.0,0,0);
  glVertex3f(1.0,0.0,0.0);
  glVertex3f(1.0,1.0,0.0);
  glVertex3f(1.0,1.0,1.0);
  glVertex3f(1.0,0.0,1.0);

  // back 
  glNormal3f(0,0,-1.0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,1.0,0.0);
  glVertex3f(1.0,1.0,0.0);
  glVertex3f(1.0,0.0,0.0);

  // bottom
  glNormal3f(0,-1.0,0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(1.0,0.0,0.0);
  glVertex3f(1.0,0.0,1.0);
  glVertex3f(0.0,0.0,1.0);

  // left 
  glNormal3f(-1.0,0,0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(0.0,0.0,1.0);
  glVertex3f(0.0,1.0,1.0);
  glVertex3f(0.0,1.0,0.0);

  glEnd();
}
