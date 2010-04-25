#ifndef _CLOTH_H_
#define _CLOTH_H_

#include <stdio.h>
#include <vector>
#include "argparser.h"
#include "edge.h"
#include "boundingbox.h"

// =====================================================================================
// Cloth Particles
// =====================================================================================

class ClothParticle
{
public:
    // ACCESSORS
    const Vec3f& getOriginalPosition() const
    {
        return original_position;
    }
    const Vec3f& getPosition() const
    {
        return position;
    }
    const Vec3f& getVelocity() const
    {
        return velocity;
    }
    const Vec3f& getAcceleration() const
    {
        return acceleration;
    }
    Vec3f getForce() const
    {
        return mass*acceleration;
    }
    float getMass() const
    {
        return mass;
    }
    bool isFixed() const
    {
        return fixed;
    }
    std::vector<Edge*> getEdges()
    {
      return edges;
    }
    // MODIFIERS
    void setOriginalPosition(const Vec3f &p)
    {
        original_position = p;
    }
    void setPosition(const Vec3f &p)
    {
        position = p;
    }
    void setVelocity(const Vec3f &v)
    {
        velocity = v;
    }
    void setAcceleration(const Vec3f &a)
    {
        acceleration = a;
    }
    void setMass(float m)
    {
        mass = m;
    }
    void setFixed(bool b)
    {
        fixed = b;
    }
    void addEdge(Edge *e)
    {
      assert(e!=NULL);
      edges.push_back(e);
    }
    
private:
    // REPRESENTATION
    Vec3f original_position;
    Vec3f position;
    Vec3f velocity;
    Vec3f acceleration;
    std::vector<Edge*> edges;
    float mass;
    bool fixed;
};

// =====================================================================================
// Cloth System
// =====================================================================================

class Cloth
{
public:
    Cloth(ArgParser *args);
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
    void addEdgeToCloth(Edge *e)
    {
      assert(e!=NULL);
      cloth_edges.push_back(e);
      printf("x:%f y:%f z:%f", e->getStartVertex()->x(), e->getStartVertex()->y(), e->getStartVertex()->z());
      printf("     x:%f y:%f z:%f\n", e->getEndVertex()->x(), e->getEndVertex()->y(), e->getEndVertex()->z());
    }
    
    void printEdges()
    {
      Edge *e;
      for (unsigned int i = 0; i < cloth_edges.size(); i++)
      {
        e = cloth_edges[i];
        std::cout << "EDGE: " << e << " SV: " << e->getStartVertex() << " EV: " << e->getEndVertex() << "\n";
        printf("Edge %d/%d: ", i, cloth_edges.size());
        printf("(S)x:%f y:%f z:%f", e->getStartVertex()->x(), e->getStartVertex()->y(), e->getStartVertex()->z());
        printf("     (E)x:%f y:%f z:%f\n\n", e->getEndVertex()->x(), e->getEndVertex()->y(), e->getEndVertex()->z());
      }
    }

private:

    // PRIVATE ACCESSORS
    const ClothParticle& getParticle(int i, int j) const
    {
        assert (i >= 0 && i < nx && j >= 0 && j < ny);
        return particles[i + j*nx];
    }
    ClothParticle& getParticle(int i, int j)
    {
        //assert (i >= 0 && i < nx && j >= 0 && j < ny);
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
    std::vector<Edge*> cloth_edges;
    // simulation parameters
    float damping;
    // spring constants
    float k_structural;
    float k_shear;
    float k_bend;
    // correction thresholds
    float provot_structural_correction;
    float provot_shear_correction;
};

#endif
