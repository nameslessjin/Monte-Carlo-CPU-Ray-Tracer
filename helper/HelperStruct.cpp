#include "HelperStruct.hpp"

Color::Color()
{
    r = 0.0f;
    g = 0.0f;
    b = 0.0f;
}

Color::Color(float _r, float _g, float _b)
{
    r = _r;
    g = _g;
    b = _b;
}

Color::Color(glm::vec3 c)
{
    r = c.r;
    g = c.g;
    b = c.b;
}

void Color::clamp()
{
    if (r > 1.0)
        r = 1.0;
    if (r < 0.0)
        r = 0.0;
    if (g > 1.0)
        g = 1.0;
    if (g < 0.0)
        g = 0.0;
    if (b > 1.0)
        b = 1.0;
    if (b < 0.0)
        b = 0.0;
}

void Color::print()
{
    std::cout << "r: " << r << " b: " << b << " g: " << g << '\n';
}

Color &Color::operator+=(Color const &c)
{
    r += c.r;
    g += c.g;
    b += c.b;
    return *this;
}

Color Color::operator+(Color const &c)
{
    return Color(r + c.r, g + c.g, b + c.b);
}

Color Color::operator*(glm::vec3 vec3)
{
    return Color(r * vec3.r, g * vec3.g, b * vec3.b);
}

Color Color::operator/=(float f)
{
    r /= f;
    g /= f;
    b /= f;
    return *this;
}

Ray::Ray()
{
    dir = glm::vec3(0, 0, 0);
    pos = glm::vec3(0, 0, 0);
}

Ray::Ray(glm::vec3 _dir, glm::vec3 _pos)
{
    dir = _dir;
    pos = _pos;
}

glm::vec3 vec3(double *v3)
{
  return glm::vec3(v3[0], v3[1], v3[2]);
}

glm::vec3 vec3(float a, float b, float c)
{
  return glm::vec3(a, b, c);
}

void clamp(float &f)
{
  if (f > 1.0f)
    f = 1.0f;
  else if (f < 0.0f)
    f = 0.0f;
}

bool AABB::intersect(const Ray &ray, float t_min, float t_max) {

    

    return false;
}