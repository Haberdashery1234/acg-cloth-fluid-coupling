#ifndef _UTILS_H
#define _UTILS_H

#include <cassert>

#define square(x) ((x)*(x))

#define max2(x,y) (((x)>=(y))?(x):(y))
#define min2(x,y) (((x)>=(y))?(y):(x))
#define max3(x,y,z) (((x)>=(y)&&(x)>=(z))?(x):(((y)>=(x)&&(y)>=(z))?(y):(z)))
#define min3(x,y,z) (((x)<=(y)&&(x)<=(z))?(x):(((y)<=(x)&&(y)<=(z))?(y):(z)))
#define max4(x,y,z,w) (max2(max2((x),(y)),max2((z),(w))))
#define min4(x,y,z,w) (min2(min2((x),(y)),min2((z),(w))))
#define mid3(x,y,z) ( (((y)<(x))&&((x)<(z))) ? (x) : \
		      ( (((z)<(x))&&((x)<(y))) ? (x) : \
		        ( (((x)<(y))&&((y)<(z))) ? (y) : \
		          ( (((z)<(y))&&((y)<(x))) ? (y) : (z) ))))


inline float triInterpolate(float x_frac, float y_frac, float z_frac,
                            float a, float b, float c, float d, float e, float f, float g, float h)
{

    assert (x_frac >= 0 && x_frac <= 1);
    assert (y_frac >= 0 && y_frac <= 1);
    assert (z_frac >= 0 && z_frac <= 1);

    // trilinear interpolation
    float ab = (1-z_frac)*a + z_frac*b;
    float cd = (1-z_frac)*c + z_frac*d;
    float ef = (1-z_frac)*e + z_frac*f;
    float gh = (1-z_frac)*g + z_frac*h;
    float abcd = (1-y_frac)*ab + y_frac*cd;
    float efgh = (1-y_frac)*ef + y_frac*gh;
    float abcdefgh = (1-x_frac)*abcd + x_frac*efgh;
    return abcdefgh;
}


inline double AreaOfTriangle(double a, double b, double c)
{
    // Area of Triangle =  (using Heron's Formula)
    //  sqrt[s*(s-a)*(s-b)*(s-c)]
    //    where s = (a+b+c)/2
    // also... Area of Triangle = 0.5 * x * c
    double s = (a+b+c) / (double)2;
    double tmp = s*(s-a)*(s-b)*(s-c);
    if (tmp < 0) return 0;
    double answer = sqrt(tmp);
    return answer;
}

inline double AreaOfTriangle(const Vec3f &a, const Vec3f &b, const Vec3f &c)
{
    return AreaOfTriangle((a-b).Length(), (b-c).Length(), (c-a).Length());
}

#endif

