#ifndef HW3_HPP
#define HW3_HPP

#include <vector>
#include <array>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifdef WIN32
#include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
#include <GL/gl.h>
#include <GL/glut.h>
#elif defined(__APPLE__)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif


#include <random>
#include <glm/gtx/string_cast.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mutex>
#include <atomic>
#include <queue>
#include "../helper/ThreadPool.hpp"
#include "../helper/HelperStruct.hpp"

#if defined(WIN32) || defined(_WIN32)
#  define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char *filename = NULL;

// different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

// the field of view of the camera
#define fov 60.0

#define MAX_REFLECT 0
#define ANTI_ALIASING_SAMPLE 16 // has to be a square number
#define LIGHT_SAMPLES 100

#define ASERT(cond)                                                      \
  do {                                                                   \
    if ((cond) == false) {                                               \
      std::cerr << #cond << " failed at line " << __LINE__ << std::endl; \
      exit(1);                                                           \
    }                                                                    \
  } while (0)

using minAABBHeap = std::priority_queue<AABB *, std::vector<AABB *>, std::greater<AABB *>>;
std::mutex mtx;
unsigned char buffer[HEIGHT][WIDTH][3];
float img[HEIGHT][WIDTH][3];

std::vector<Triangle> triangles;
std::vector<Sphere> spheres;
std::vector<Light> lights;
double F0[3];
double ambient_light[3];

std::vector<AABB> all_aabbs;
AABB *HBV;

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
float e = 1e-4;
float sigma = 1e-5;
float e5 = 1e-5;

std::vector<glm::vec3> randomPointsInQuadrilateral(const int light_samples);
float calcLightArea(const Light &light);
glm::vec3 calculateF(const MonteCarlo &mc);
float calculateD(const MonteCarlo &mc, const glm::vec3 &m);
float findAngleRad(const glm::vec3 &u, const glm::vec3 &v);
float positiveChar(float t);
float calculateG1(const MonteCarlo &mc, const glm::vec3 &v, const glm::vec3 &m);
glm::vec3 calculateFD(const MonteCarlo &mc);
glm::vec3 calculateFS(const MonteCarlo &mc);
glm::vec3 calculateBRDF(const MonteCarlo &mc);
Color calculateMonteCarlo(const GLM_Vertex &v, const Light &light);
glm::vec3 shadowRayColor(const Ray &shadow_ray, const Light &light);
glm::vec3 maxPointTriangle(Triangle &t);
glm::vec3 minPointTriangle(Triangle &t);
glm::vec3 maxPointSphere(Sphere &s);
glm::vec3 minPointSphere(Sphere &s);
AABB *pickCloserAABB(AABB *aabb1, AABB *aabb2, const Ray &ray);
AABB *checkIntersectionWithAABB(AABB *aabb, const Ray &ray);
void contructHVB();
void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);
GLM_Vertex calcGLMVertex(Triangle &t, glm::vec3 &intersection);
GLM_Vertex calcGLMVertex(Sphere &s, glm::vec3 &intersection);
Color phong_shading(Sphere &s, Light &l, glm::vec3 intersection);
Color phong_shading(Triangle &t, Light &l, glm::vec3 &intersection);
Color calc_phong_shading(glm::vec3 l_dir, glm::vec3 l_color, GLM_Vertex vertex);
bool check_single_sphere_intersection(const Ray &r, const Sphere &s, float &t, glm::vec3 &intersection);
bool check_single_triangle_intersection(const Ray &r, const Triangle &tri, float &t, glm::vec3 &intersection);
float calc_triangle_area_xy(glm::vec3 a, glm::vec3 b, glm::vec3 c);
float calc_triangle_area_xz(glm::vec3 a, glm::vec3 b, glm::vec3 c);
float calc_triangle_area(glm::vec3 a, glm::vec3 b, glm::vec3 c);
bool point_triangle_xy(glm::vec3 c0, glm::vec3 c1, glm::vec3 c2, glm::vec3 c);
bool point_triangle_xz(glm::vec3 c0, glm::vec3 c1, glm::vec3 c2, glm::vec3 c);
glm::vec3 doubleToVec(const double *d);
bool check_block(const Ray &r, const Sphere &s, const Light &l);
bool check_block(const Ray &r, const Triangle &t, const Light &l);
Color check_intersection(Ray &ray, int time);
glm::vec3 calc_reflect_dir(int sphere_i, int triangle_i, glm::vec3 dir, glm::vec3 intersection);
Color calc_shadow_ray(int sphere_i, int triangle_i, glm::vec3 &intersection);
Color calc_ray_color(int sphere_i, int triangle_i, glm::vec3 intersection, Ray &ray, int time);
void generate_ray(Ray &ray, int x, int y);
void generate_ray_antialiasing(std::vector<Ray> &rays, int x, int y, int num_samples);
void stratified_sampling(std::vector<Ray> &rays, int x, int y, int num_samples);
void draw_pixel(int x, int y, std::atomic<int> &finished);


#endif