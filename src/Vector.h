
#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>

const float PI = 3.14159265358979f;

class Vector  
{
public:
    Vector();
    Vector(float x, float y, float z);
    Vector(float[]);
    Vector(int[]);
    ~Vector();

    float x;
    float y;
    float z;

    /* Operators for vector addition. Note that these
       create a new Vector.
    */
    Vector operator+(const Vector& rhs) const;
    Vector operator-(const Vector& rhs) const;

    /* Operators for scalar multiplication. Note that
       these create a new Vector
    */
    friend Vector  operator*    (const Vector& rhs, float fScale);
    friend Vector  operator*    (float fScale, const Vector& rhs);
    Vector operator/(const float factor) const;

    /* 
        Operators for vector arithmetic updates.
    */
    Vector&  operator*=   (float fScale);
    Vector&  operator/=   (float fScale);
    Vector&  operator+=   (const Vector& rhs);
    Vector&  operator-=   (const Vector& rhs);

    /* 
        Operators for vector comparisons.
    */
    bool operator==   (const Vector& rhs);
    bool operator!=   (const Vector& rhs);

    /* Array operator allows access to x, y and z
       members through indices VECTOR_X, VECTOR_Y,
       and VECTOR_Z, respectively (0, 1, and 2).
       Only these three values should be used. If
       another index is used, the Z-value is returned.
       (RHS version)
    */
    float operator[](const int index) const;

    /* Array operator allows access to x, y and z
       members through indices VECTOR_X, VECTOR_Y,
       and VECTOR_Z, respectively (0, 1, and 2).
       Only these three values should be used. If
       another index is used, the Z-value is returned.
       (LHS version)
    */
    float& operator[](const int index);

    /* Operation for vector addition. This operation
       changes the original Vector. If the original is
       not needed after the operation, then this is faster
       than creating a new Vector object.
    */
    void translateBy(const Vector& rhs);

    /* Operation for scalar multiplication. This operation
       changes the original Vector. If the original is
       not needed after the operation, then this is faster
       than creating a new Vector object.
    */
    void scaleBy(const float factor);

    /* Normalizes the vector, if the vector is not the
       zero vector. If it is the zero vector, then it is
       unchanged.
    */
    float normalize();

    /* Zeros out the vector.
    */
    void zero();
    
    /* Return the squared length of the vector. Slightly
       faster, and we may only need to compare relative lengths.
    */
    float getSquaredLength() const;
    
    /* Return the magnitude of the vector
     */
    float getMagnitude() const;

    /* Return the distance between this vector and the
     specified right-hand side vector.
     */
    float distance(const Vector& rhs) const;
    
    /* Return the dot product of this vector with the
       specified right-hand side vector.
    */
    float dot(const Vector& rhs) const;

    /* Return the cross product of this vector with the
       specified right-hand side vector.
    */
    Vector cross(const Vector& rhs) const;

    /* Put this vector in the specified float array. The
       array must be at least 3 elements long. Each component
       is cast to a float, so loss of precision is likely.
    */
    void toArray(float array[]) const;

    /* Set this vector to the contents of the first three
       elements of the specified float array. The
       array must be at least 3 elements long.
    */
    void fromArray(float array[]);

    /* Rotate this vector about an axis. The rotation is specified
       in degrees.
    */
    void rotateX(const float degrees);
    void rotateY(const float degrees);
    void rotateZ(const float degrees);

    /* Rotate this vector about an axis. The rotation is specified
       in radians. This is marginally faster than the rotations specified in degrees, if
       you already know the radians.
    */
    void radianRotateX(const float radians);
    void radianRotateY(const float radians);
    void radianRotateZ(const float radians);

    /* Rotate this vector about an arbitrary axis. The rotation
       is specified in degrees.
    */
    void rotateAxis(const Vector& axis, const float degrees);

    /* Return a linearly interpolated vector between this
       vector and an endpoint. t should vary between 0 and 1.
    */
    Vector interpolate1(const Vector& endPoint, const float t) const;

    /* Return a qudratic Bezier interpolated vector with the
       three controls points this vector, midControl, and endControl.
       t should vary between 0 and 1.
    */
    Vector interpolate2(const Vector& midControl, const Vector& endControl, const float t) const;

    /* Return a cubic Bezier interpolated vector with the four
       controls points this vector, leftControl, rightControl,
       and endControl. t should vary between 0 and 1.
    */
    Vector interpolate3(const Vector& leftControl, const Vector& rightControl, const Vector& endControl, const float t) const;

};

#endif //end of VECTOR_H
