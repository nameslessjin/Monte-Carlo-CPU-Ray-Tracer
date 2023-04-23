#ifndef HELPER_STRUCT_HPP
#define HELPER_STRUCT_HPP

#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <stdio.h>

struct Vertex
{
  double position[3];
  double normal[3];

  double color_diffuse[3];
  double roughness;
  double metallic;

  double color_specular[3]; // not used in monte carlo
  double shininess; // not used in monte carlo
};

struct GLM_Vertex
{
  glm::vec3 position;
  glm::vec3 kd;
  glm::vec3 n;
  float roughness;
  float metallic;

  glm::vec3 ks; // not used in monte carlo
  float shininess; // not used in monte carlo
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double radius;

  double color_diffuse[3];
  double roughness;
  double metallic;

  double color_specular[3]; // not used in monte carlo
  double shininess; // not used in monte carlo
};

struct Light
{
  double position[3];
  double normal[3];
  double color[3];
  double p[4][3];
};

struct MonteCarlo {
  glm::vec3 p, w_i, w_o, albedo, n, F0;
  Light light;
  float metallic, roughness;
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
  Color operator/(Color const &c);
};

struct Ray
{
  glm::vec3 dir;
  glm::vec3 pos;

  Ray();
  Ray(glm::vec3 _dir, glm::vec3 _pos);

  void print();
};

class AABB
{
public:
  glm::vec3 min, max;
  int sphere_i;
  int triangle_i;
  AABB *left, *right; // Binary Tree like HBV
  AABB() : min(glm::vec3()), max(glm::vec3()), left(nullptr), right(nullptr), sphere_i(-1), triangle_i(-1) {}
  AABB(const glm::vec3 &min, const glm::vec3 &max) : min(min), max(max), left(nullptr), right(nullptr), sphere_i(-1), triangle_i(-1) {}
  AABB(const glm::vec3 &min, const glm::vec3 &max, int sphere_i, int triangle_i) : min(min), max(max), left(nullptr), right(nullptr), sphere_i(sphere_i), triangle_i(triangle_i) {}
  AABB(const glm::vec3 &min, const glm::vec3 &max, int sphere_i, int triangle_i, AABB *left, AABB *right) : min(min), max(max), sphere_i(sphere_i), triangle_i(triangle_i), left(left), right(right) {}

  bool intersect(const Ray &ray);
  float surfaceArea();
  void print();
  bool operator>(const AABB *(&aabb));
};

glm::vec3 mergedMin(const AABB &aabb1, const AABB &aabb2);
glm::vec3 mergedMax(const AABB &aabb1, const AABB &aabb2);
float mergedSurfaceArea(const AABB &aabb1, const AABB &aabb2);
float sahCost(AABB &aabb1, AABB &aabb2);
std::vector<AABB> buildHVB(std::vector<AABB> &aabbs);
std::pair<int, int> findBestPair(std::vector<AABB> &aabbs);

glm::vec3 vec3(double *v3);

glm::vec3 vec3(float a, float b, float c);

void clamp(float &f);

#endif