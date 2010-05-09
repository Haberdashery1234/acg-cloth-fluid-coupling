#ifndef _CLOTH_H_
#define _CLOTH_H_

#include <stdio.h>
#include <vector>
#include "argparser.h"
#include "boundingbox.h"
#include "fluid.h"
#include "clothparticle.h"

// =====================================================================================
// Cloth System
// =====================================================================================

class Cloth
{
public:
    Cloth(ArgParser *args);
    Cloth(ArgParser *args, Fluid *fluid);
    ~Cloth()
    {
        delete [] particles;
    }

    // ACCESSORS
    const BoundingBox& getBoundingBox() const
    {
        return box;
    }

    // PAINTING & ANIMATING
    void Paint() const;
    void Animate();
    void CheckCollision();

private:

    // PRIVATE ACCESSORS
    const ClothParticle& getParticle(int i, int j) const
    {
        assert (i >= 0 && i < nx && j >= 0 && j < ny);
        return particles[i + j*nx];
    }
    ClothParticle& getParticle(int i, int j)
    {
        assert (i >= 0 && i < nx && j >= 0 && j < ny);
        return particles[i + j*nx];
    }

    // HELPER FUNCTION
    void computeBoundingBox();
    void makeCorrection(ClothParticle &p1, ClothParticle &p1, float corr);

    // REPRESENTATION
    ArgParser *args;
    // grid data structure
    int nx, ny;
    ClothParticle *particles;
    BoundingBox box;
    // simulation parameters
    float damping;
    // spring constants
    float k_structural;
    float k_shear;
    float k_bend;
    // correction thresholds
    float provot_structural_correction;
    float provot_shear_correction;
    float collision_boundary;
    
    Fluid *fluid;
};

#endif
