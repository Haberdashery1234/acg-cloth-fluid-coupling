#ifndef _CLOTH_PARTICLE_H_
#define _CLOTH_PARTICLE_H_

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

private:
    // REPRESENTATION
    Vec3f original_position;
    Vec3f position;
    Vec3f velocity;
    Vec3f acceleration;
    float mass;
    bool fixed;
};

#endif
