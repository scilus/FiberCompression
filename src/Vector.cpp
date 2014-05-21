
#include <assert.h>

#include "Vector.h"

Vector::Vector() :
x(0),
y(0),
z(0)
{
}

Vector::Vector(float newX, float newY, float newZ) :
x(newX),
y(newY),
z(newZ)
{
}


Vector::Vector(float array[])
{
    x = array[0];
    y = array[1];
    z = array[2];
}

Vector::Vector(int array[])
{
    x = (float)array[0];
    y = (float)array[1];
    z = (float)array[2];
}

Vector::~Vector()
{

}

Vector& Vector::operator*=(float fScale)
{
    x *= fScale;
    y *= fScale;
    z *= fScale;

    return *this;
}

Vector& Vector::operator/=(float fScale)
{
    x /= fScale;
    y /= fScale;
    z /= fScale;

    return *this;
}

Vector& Vector::operator+=(const Vector& rhs)
{
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;

    return *this;
}

Vector& Vector::operator-=(const Vector& rhs)
{
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;

    return *this;
}

bool Vector::operator==(const Vector& rhs)
{
    return (x == rhs.x && y == rhs.y && z == rhs.z);
}

bool Vector::operator!=(const Vector& rhs)
{
    return !(*this == rhs);
}

Vector Vector::operator+(const Vector& rhs) const
{
    return Vector(x + rhs.x, y + rhs.y, z + rhs.z);
}

Vector Vector::operator-(const Vector& rhs) const
{
    return Vector(x - rhs.x, y - rhs.y, z - rhs.z);
}

Vector operator*(const Vector& rhs, float fScale)
{
    Vector result;

    result.x = rhs.x*fScale;
    result.y = rhs.y*fScale;
    result.z = rhs.z*fScale;

    return result;
}

Vector operator*(float fScale, const Vector& rhs)
{
    Vector result;

    result.x = rhs.x*fScale;
    result.y = rhs.y*fScale;
    result.z = rhs.z*fScale;

    return result;
}

Vector Vector::operator/(const float factor) const
{
    return Vector(x / factor, y / factor, z / factor);
}

float Vector::operator[](const int index) const
{
    assert(index >= 0 && index <= 2);
    
    switch (index)
    {
        case 0:
            return x;
            break;
        case 1:
            return y;
            break;
        case 2:
            return z;
            break;
        default :
            return z;
            break;
    }
}

float& Vector::operator[](const int index)
{
    assert(index >= 0 && index <= 2);
    
    switch (index)
    {
        case 0:
            return x;
            break;
        case 1:
            return y;
            break;
        case 2:
            return z;
            break;
        default :
            return z;
            break;
    }
}

void Vector::translateBy(const Vector& rhs)
{
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
}

void Vector::scaleBy(const float factor)
{
    x *= factor;
    y *= factor;
    z *= factor;
}

float Vector::normalize()
{
    float length = sqrt(x * x + y * y + z * z);
    if (length > 0.000001)
    {
        x /= length;
        y /= length;
        z /= length;
    }
    return length;
}

void Vector::zero()
{
    x = y = z = 0.0;
}

float Vector::getSquaredLength() const
{
    return (x * x + y * y + z * z);
}

float Vector::getMagnitude() const
{
    return sqrt(x * x + y * y + z * z);
}

float Vector::distance(const Vector& rhs) const
{
    return sqrt(pow(x - rhs.x, 2) + pow(y - rhs.y, 2) + pow(z - rhs.z, 2));
}

float Vector::dot(const Vector& rhs) const
{
    return (x * rhs.x + y * rhs.y + z * rhs.z);
}

Vector Vector::cross(const Vector& rhs) const
{
    return Vector((y * rhs.z) - (z * rhs.y),
                  (z * rhs.x) - (x * rhs.z),
                  (x * rhs.y) - (y * rhs.x));
}

void Vector::toArray(float array[]) const
{
    array[0] = x;
    array[1] = y;
    array[2] = z;
}

void Vector::fromArray(float array[])
{
    x = array[0];
    y = array[1];
    z = array[2];
}

void Vector::rotateX(const float degrees) {
    float radians = degrees *  PI/ 180.0f;
    float cosAngle = cos(radians);
    float sinAngle = sin(radians);
    float origY = y;
    y =    y * cosAngle - z * sinAngle;
    z = origY * sinAngle + z * cosAngle;
}

void Vector::rotateY(const float degrees) {
    float radians = degrees * PI / 180.0f;
    float cosAngle = cos(radians);
    float sinAngle = sin(radians);
    float origX = x;
    x =    x * cosAngle + z * sinAngle;
    z = z * cosAngle - origX * sinAngle;
}

void Vector::rotateZ(const float degrees) {
    float radians = degrees * PI / 180.0f;
    float cosAngle = cos(radians);
    float sinAngle = sin(radians);
    float origX = x;
    x =    x * cosAngle - y * sinAngle;
    y = origX * sinAngle + y * cosAngle;
}

void Vector::radianRotateX(const float radians) {
    float cosAngle = cos(radians);
    float sinAngle = sin(radians);
    float origY = y;
    y =    y * cosAngle - z * sinAngle;
    z = origY * sinAngle + z * cosAngle;
}

void Vector::radianRotateY(const float radians) {
    float cosAngle = cos(radians);
    float sinAngle = sin(radians);
    float origX = x;
    x =    x * cosAngle + z * sinAngle;
    z = z * cosAngle - origX * sinAngle;
}

void Vector::radianRotateZ(const float radians) {
    float cosAngle = cos(radians);
    float sinAngle = sin(radians);
    float origX = x;
    x =    x * cosAngle - y * sinAngle;
    y = origX * sinAngle + y * cosAngle;
}

void Vector::rotateAxis(const Vector& axis, float degrees)
{
    Vector rotator;
    Vector rx;
    Vector ry;
    Vector rz;

    // Normalize the axis and store in rz for later
    // potential use.
    rz = axis / axis.getMagnitude();

    // If this vector is parallel to the axis (Now in rz),
    // don't bother rotating it. Check with Cauchy-Schwartz
    // u dot v == length of u * length of v iff u is linearly
    // dependant on v. Length of rz is one, because it is
    // normalized.
    if (!(fabs(this->dot(rz)) == this->getMagnitude()))
    {
        // If we're not already rotating around the Z, transform to Z
        if (axis.x == 0 && axis.y == 0) {
            // In this case, the axis is along Z already, so we
            // wont bother rotating our axis of rotation to Z.
            // Since we're checking directly for 0, this is going
            // to happen extremely rarely. However, if it does
            // happen, some of the below math falls apart. So,
            // we check.
            rotator = *this;
        }
        else
        {
            // Build the rotation matrix
            // rz was assigned already while normalizing the axis.
            rx = this->cross(axis);
            rx.normalize();
            ry = rz.cross(rx);

            // Move this vector such that the axis would be in Z
            rotator = Vector(rx.dot(*this), ry.dot(*this), rz.dot(*this));
        }

        // Rotate this vector around Z
        rotator.rotateZ(degrees);

        if (axis.x == 0 && axis.y == 0)
        {
            *this = rotator;
        }
        else
        {
            // Move back so axis is in original location.
            this->x = rotator.x * rx.x +
                      rotator.y * ry.x +
                      rotator.z * rz.x;
            this->y = rotator.x * rx.y +
                      rotator.y * ry.y +
                      rotator.z * rz.y;
            this->z = rotator.x * rx.z +
                      rotator.y * ry.z +
                      rotator.z * rz.z;
        }
    }
}

Vector Vector::interpolate1(const Vector& endPoint, const float t) const
{
    return Vector((1.0f - t) * x + t * endPoint.x,
                  (1.0f - t) * y + t * endPoint.y,
                  (1.0f - t) * z + t * endPoint.z);
}

Vector Vector::interpolate2(const Vector& midControl, const Vector& endControl, const float t) const
{
    Vector left = this->interpolate1(midControl, t);
    Vector right = midControl.interpolate1(endControl, t);
    return left.interpolate1(right, t);
}

Vector Vector::interpolate3(const Vector& leftControl, const Vector& rightControl, const Vector& endControl, const float t) const
{
    Vector begin = this->interpolate1(leftControl, t);
    Vector mid = leftControl.interpolate1(rightControl, t);
    Vector end = rightControl.interpolate1(endControl, t);
    return begin.interpolate2(mid, end, t);
}
