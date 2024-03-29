#ifndef _CELL_H_
#define _CELL_H_

#include <cassert>
#include <vector>
#include "clothparticle.h"

// ==============================================================================

class FluidParticle
{
public:
    // accessor
    Vec3f getPosition() const
    {
      return position;
    }
    Vec3f getVelocity() const
    {
      return Vec3f(u_vel, v_vel, w_vel);;
    }
    // modifer
    void setPosition(Vec3f p)
    {
      position = p;
    }
    void set_u_vel(float u)
    {
      u_vel = u;
    }
    void set_v_vel(float v)
    {
      v_vel = v;
    }
    void set_w_vel(float w)
    {
      w_vel = w;
    }
    void setVelocity(Vec3f v)
    {
      u_vel = v.x();
      v_vel = v.y();
      w_vel = v.z();
    }
private:
    // representation
    Vec3f position;
    Vec3f velocity;
    float u_vel, v_vel, w_vel;
};

// ==============================================================================
enum CELL_STATUS { CELL_BOUNDARY, CELL_EMPTY, CELL_SURFACE, CELL_FULL };

class Cell
{
public:

    // ========================
    // CONSTRUCTOR & DESTRUCTOR
    Cell()
    {
        pressure = 0;
        status = CELL_BOUNDARY;
        u_plus = 0;
        v_plus = 0;
        w_plus = 0;
        new_u_plus = 0;
        new_v_plus = 0;
        new_w_plus = 0;
    }
    ~Cell()
    {
        for (unsigned int i = 0; i < particles.size(); i++)
        {
            delete particles[i];
        }
    }

    // =========
    // ACCESSORS
    float getPressure() const
    {
        return pressure;
    }
    enum CELL_STATUS getStatus() const
    {
        return status;
    }
    float get_u_plus() const
    {
        return u_plus;
    }
    float get_v_plus() const
    {
        return v_plus;
    }
    float get_w_plus() const
    {
        return w_plus;
    }
    float get_new_u_plus() const
    {
        return new_u_plus;
    }
    float get_new_v_plus() const
    {
        return new_v_plus;
    }
    float get_new_w_plus() const
    {
        return new_w_plus;
    }
    int numParticles() const
    {
        return particles.size();
    }
    std::vector<FluidParticle*>& getParticles()
    {
        return particles;
    }
    std::vector<ClothParticle*>& getClothParticles()
    {
        return cparticles;
    }

    // =========
    // MODIFIERS
    void setPressure(float p)
    {
        pressure = p;
    }
    void setStatus(enum CELL_STATUS s)
    {
        status = s;
    }
    void set_u_plus(float f)
    {
        u_plus = new_u_plus = f;
    }
    void set_v_plus(float f)
    {
        v_plus = new_v_plus = f;
    }
    void set_w_plus(float f)
    {
        w_plus = new_w_plus = f;
    }
    void set_new_u_plus(float f)
    {
        new_u_plus = f;
    }
    void set_new_v_plus(float f)
    {
        new_v_plus = f;
    }
    void set_new_w_plus(float f)
    {
        new_w_plus = f;
    }
    void adjust_new_u_plus(float f)
    {
        new_u_plus += f;
    }
    void adjust_new_v_plus(float f)
    {
        new_v_plus += f;
    }
    void adjust_new_w_plus(float f)
    {
        new_w_plus += f;
    }
    void copyVelocity()
    {
        u_plus = new_u_plus;
        new_u_plus = 0;
        v_plus = new_v_plus;
        new_v_plus = 0;
        w_plus = new_w_plus;
        new_w_plus = 0;
    }
    void addParticle(FluidParticle *p)
    {
        assert(p != NULL);
        particles.push_back(p);
    }
    void removeParticle(FluidParticle *p)
    {
        assert(p != NULL);
        for (std::vector<FluidParticle*>::iterator i = particles.begin(); i != particles.end(); i++)
        {
            if (*i == p)
            {
                particles.erase(i);
                return;
            }
        }
        assert (0);
    }
    void addClothParticle(ClothParticle *p)
    {
        assert(p != NULL);
        cparticles.push_back(p);
    }
    void removeClothParticle(ClothParticle *p)
    {
        assert(p != NULL);
        for (std::vector<ClothParticle*>::iterator i = cparticles.begin(); i != cparticles.end(); i++)
        {
            if (*i == p)
            {
                cparticles.erase(i);
                return;
            }
        }
        assert (0);
    }

private:

    // ==============
    // REPRESENTATION
    float pressure;
    enum CELL_STATUS status;

    // velocities at the center of each face (flowing in the positive direction)
    float u_plus,v_plus,w_plus;
    float new_u_plus,new_v_plus,new_w_plus;

    std::vector<FluidParticle*> particles;
    std::vector<ClothParticle*> cparticles;
};

#endif
