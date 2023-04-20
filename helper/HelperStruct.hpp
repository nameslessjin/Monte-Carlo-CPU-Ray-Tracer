#ifndef HELPER_STRUCT_HPP
#define HELPER_STRUCT_HPP

#include <glm/glm.hpp>
#include <iostream>

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct GLM_Vertex
{
  glm::vec3 position;
  glm::vec3 kd;
  glm::vec3 ks;
  glm::vec3 n;
  float shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};


struct Color
{
  float r, g, b;

  Color();
  Color(float _r, float _g, float _b);
  Color(glm::vec3 c);

  void clamp();
  void print();

  Color &operator+=(Color const &c);
  Color operator+(Color const &c);
  Color operator*(glm::vec3 vec3);
  Color operator/=(float f);

};

struct Ray
{
  glm::vec3 dir;
  glm::vec3 pos;

  Ray();
  Ray(glm::vec3 _dir, glm::vec3 _pos);
};

glm::vec3 vec3(double *v3);

glm::vec3 vec3(float a, float b, float c);

void clamp(float &f);

#endif