#ifndef _MARCHING_CUBES_
#define _MARCHING_CUBES_

#include "vectors.h"

// ==================================================================================
// A helper class for the Marching Cubes algorithm

class GridValue
{
public:
    GridValue() {}
    GridValue(const Vec3f &p, const Vec3f &n, float v) : position(p), normal(n), value(v) { }
    Vec3f position;
    Vec3f normal;
    float value;
};

// ==================================================================================
class MarchingCubes
{

public:
    MarchingCubes(int _nx, int _ny, int _nz, float _dx, float _dy, float _dz) :
        nx(_nx), ny(_ny), nz(_nz), dx(_dx), dy(_dy), dz(_dz)
    {
        values = new float[nx*ny*nz];
    }
    ~MarchingCubes()
    {
        delete [] values;
    }

    // position
    float get(int x, int y, int z) const
    {
        assert (x >= 0 && x < nx);
        assert (y >= 0 && y < ny);
        assert (z >= 0 && z < nz);
        return values[z*nx*ny + y*nx + x];
    }
    // normal
    Vec3f getNormal(int x, int y, int z) const
    {
        assert (x >= 0 && x < nx);
        assert (y >= 0 && y < ny);
        assert (z >= 0 && z < nz);
        float dx;
        float dy;
        float dz;
        if (x == 0)         dx = get(x+1,y,z)-get(x  ,y,z);
        else if (x == nx-1) dx = get(x  ,y,z)-get(x-1,y,z);
        else                dx = get(x+1,y,z)-get(x-1,y,z);
        if (y == 0)         dy = get(x,y+1,z)-get(x,y  ,z);
        else if (y == ny-1) dy = get(x,y  ,z)-get(x,y-1,z);
        else                dy = get(x,y+1,z)-get(x,y-1,z);
        if (z == 0)         dz = get(x,y,z+1)-get(x,y,z  );
        else if (z == nz-1) dz = get(x,y,z  )-get(x,y,z-1);
        else                dz = get(x,y,z+1)-get(x,y,z-1);
        Vec3f norm(-dx,-dy,-dz);
        norm.Normalize();
        return norm;
    }
    void set(int x, int y, int z, float v) const
    {
        assert (x >= 0 && x < nx);
        assert (y >= 0 && y < ny);
        assert (z >= 0 && z < nz);
        values[z*nx*ny + y*nx + x] = v;
    }

    // =============
    // THE DRAW CODE
    void Paint(float isosurface) const;

private:

    // private helper functions
    void PaintTetra(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, float isosurface) const;
    void PaintTetraHelper4(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, float isosurface) const;
    void PaintTetraHelper3(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, float isosurface) const;
    void PaintTetraHelper2(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, float isosurface) const;
    void PaintTetraHelper1(const GridValue &a, const GridValue &b, const GridValue &c, const GridValue &d, float isosurface) const;
    void drawIfBoundary(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) const;

    // ==============
    // REPRESENTATION
    int nx, ny, nz;
    float dx, dy, dz;
    float* values;
};

// ==================================================================================

void insertNormal(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3);
void drawTriangle(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3);
void drawTriangleWithNormals(const Vec3f &n1, const Vec3f &p1,
                             const Vec3f &n2, const Vec3f &p2,
                             const Vec3f &n3, const Vec3f &p3);

#endif

