#ifndef HW3_HPP
#define HW3_HPP

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
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <mutex>
#include <atomic>
#include <queue>
#include "../helper/ThreadPool.hpp"
#include "../helper/HelperStruct.hpp"

#ifdef WIN32
#define strcasecmp _stricmp
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

#define MAX_REFLECT 3
#define ANTI_ALIASING_SAMPLE 16

using minAABBHeap = std::priority_queue<AABB *, std::vector<AABB *>, std::greater<AABB *>>;
std::mutex mtx;
unsigned char buffer[HEIGHT][WIDTH][3];
float img[HEIGHT][WIDTH][3];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
std::vector<AABB> all_aabbs;
AABB *HBV;
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
float e = 1e-4;
float sigma = 1e-5;

float calculateD(const MonteCarlo &mc, const glm::vec3 &m);
float findAngleRad(const glm::vec3 &u, const glm::vec3 &v);
float positiveChar(float t);
float calculateG1(const MonteCarlo &mc, glm::vec3 &v, glm::vec3 &m);
float calculateFD(const MonteCarlo &mc);
float calculateFS(const MonteCarlo &mc);
glm::vec3 calculateBRDF(const MonteCarlo &mc);
Color calculateMonteCarlo();
glm::vec3 shadowRayColor(Ray &shadow_ray, Light &light);
glm::vec2 findPixelTopLeftCorner(int x, int y);
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
GLM_Vertex calc_barycentric_interpolation(Triangle &t, glm::vec3 &intersection);
Color phong_shading(Sphere &s, Light &l, glm::vec3 intersection);
Color phong_shading(Triangle &t, Light &l, glm::vec3 &intersection);
Color calc_phong_shading(glm::vec3 l_dir, glm::vec3 l_color, GLM_Vertex vertex);
bool check_single_sphere_intersection(const Ray &r, Sphere &s, float &t, glm::vec3 &intersection);
bool check_single_triangle_intersection(const Ray &r, Triangle &tri, float &t, glm::vec3 &intersection);
float calc_triangle_area_xy(glm::vec3 a, glm::vec3 b, glm::vec3 c);
float calc_triangle_area_xz(glm::vec3 a, glm::vec3 b, glm::vec3 c);
float calc_triangle_area(glm::vec3 a, glm::vec3 b, glm::vec3 c);
bool point_triangle_xy(glm::vec3 c0, glm::vec3 c1, glm::vec3 c2, glm::vec3 c);
bool point_triangle_xz(glm::vec3 c0, glm::vec3 c1, glm::vec3 c2, glm::vec3 c);
bool check_block(Ray &r, Sphere &s, Light &l);
bool check_block(Ray &r, Triangle &t, Light &l);
Color check_intersection(Color &c, Ray &ray, int time);
glm::vec3 calc_reflect_dir(int sphere_i, int triangle_i, glm::vec3 dir, glm::vec3 intersection);
void calc_shadow_ray(Color &color, int sphere_i, int triangle_i, glm::vec3 &intersection);
void calc_ray_color(Color &c, int sphere_i, int triangle_i, glm::vec3 intersection, Ray &ray, int time);
void generate_ray(Ray &ray, int x, int y);
void generate_ray_antialiasing(std::vector<Ray> &rays, int x, int y, int num_samples);
void draw_pixel(int x, int y, std::atomic<int> &finished);

#endif