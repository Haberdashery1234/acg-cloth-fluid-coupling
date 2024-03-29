#ifndef _MATRIX_H_
#define _MATRIX_H_

//
// originally implemented by Justin Legakis
//

#include <cmath>
#include <cassert>

#include "vectors.h"

// ====================================================================
// ====================================================================

class Matrix
{

public:

    // CONSTRUCTORS & DESTRUCTOR
    Matrix()
    {
        Clear();
    }
    Matrix(const Matrix& m);
    Matrix& operator=(const Matrix& m);
    ~Matrix() {}

    // ACCESSORS
    float Get(int row, int col) const
    {
        assert (row >= 0 && row < 4);
        assert (col >= 0 && col < 4);
        return data[row+col*4];
    }
    const float* glGet(void) const
    {
        return data;
    }

    // MODIFIERS
    void Set(int row, int col, float v)
    {
        assert (row >= 0 && row < 4);
        assert (col >= 0 && col < 4);
        data[row+col*4] = v;
    }
    void SetToIdentity();
    void Clear();

    void Transpose(Matrix &m) const;
    void Transpose()
    {
        Transpose(*this);
    }

    int Inverse(Matrix &m, float epsilon = 1e-08) const;
    int Inverse(float epsilon = 1e-08)
    {
        return Inverse(*this,epsilon);
    }

    // OVERLOADED OPERATORS
    friend Matrix operator+(const Matrix &m1, const Matrix &m2);
    friend Matrix operator-(const Matrix &m1, const Matrix &m2);
    friend Matrix operator*(const Matrix &m1, const Matrix &m2);
    friend Matrix operator*(const Matrix &m1, float f);
    friend Matrix operator*(float f, const Matrix &m)
    {
        return m * f;
    }
    Matrix& operator+=(const Matrix& m)
    {
        *this = *this + m;
        return *this;
    }
    Matrix& operator-=(const Matrix& m)
    {
        *this = *this - m;
        return *this;
    }
    Matrix& operator*=(const float f)
    {
        *this = *this * f;
        return *this;
    }
    Matrix& operator*=(const Matrix& m)
    {
        *this = *this * m;
        return *this;
    }

    // TRANSFORMATIONS
    static Matrix MakeTranslation(const Vec3f &v);
    static Matrix MakeScale(const Vec3f &v);
    static Matrix MakeScale(float s)
    {
        return MakeScale(Vec3f(s,s,s));
    }
    static Matrix MakeXRotation(float theta);
    static Matrix MakeYRotation(float theta);
    static Matrix MakeZRotation(float theta);
    static Matrix MakeAxisRotation(const Vec3f &v, float theta);

    // Use to transform a point with a matrix
    // that may include translation
    void Transform(Vec4f &v) const;
    void Transform(Vec3f &v) const
    {
        Vec4f v2 = Vec4f(v.x(),v.y(),v.z(),1);
        Transform(v2);
        v.Set(v2.x(),v2.y(),v2.z());
    }
    void Transform(Vec2f &v) const
    {
        Vec4f v2 = Vec4f(v.x(),v.y(),1,1);
        Transform(v2);
        v.Set(v2.x(),v2.y());
    }

    // Use to transform the direction of the ray
    // (ignores any translation)
    void TransformDirection(Vec3f &v) const
    {
        Vec4f v2 = Vec4f(v.x(),v.y(),v.z(),0);
        Transform(v2);
        v.Set(v2.x(),v2.y(),v2.z());
    }

    // INPUT / OUTPUT
    friend std::ostream& operator<<(std::ostream &ostr, const Matrix &m);
    friend std::istream& operator>>(std::istream &istr, Matrix &m);

private:

    // REPRESENTATION
    // column-major order
    float	data[16];

};

// ====================================================================
// ====================================================================

#endif
