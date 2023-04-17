/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Your name here>
 * *************************
 */

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

#include <thread>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

#define MAX_REFLECT 7

unsigned char buffer[HEIGHT][WIDTH][3];

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

  Color()
  {
    r = 0.0f;
    g = 0.0f;
    b = 0.0f;
  }
  Color(float _r, float _g, float _b)
  {
    r = _r;
    g = _g;
    b = _b;
  }
  Color(glm::vec3 c)
  {
    r = c.r;
    g = c.g;
    b = c.b;
  }

  void clamp()
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

  Color &operator+=(Color const &c)
  {
    r += c.r;
    g += c.g;
    b += c.b;
    return *this;
  }

  Color &operator+(Color const &c)
  {
    r += c.r;
    g += c.g;
    b += c.b;
    return *this;
  }

  Color operator*(glm::vec3 vec3)
  {

    Color new_c = Color(r * vec3.r, g * vec3.g, b * vec3.b);
    return new_c;
  }

  Color operator/=(float f)
  {
    r /= f;
    g /= f;
    b /= f;
    return *this;
  }

  void print()
  {
    std::cout << "r: " << r << " b: " << b << " g: " << g << '\n';
  }
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
float e = 1e-13;

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

struct Ray
{
  glm::vec3 dir;
  glm::vec3 pos;

  Ray()
  {
    dir = glm::vec3(0, 0, 0);
    pos = glm::vec3(0, 0, 0);
  }
  Ray(glm::vec3 _dir, glm::vec3 _pos)
  {
    dir = _dir;
    pos = _pos;
  }

};

void clamp(float &f);
glm::vec3 vec3(double *v3);
glm::vec3 vec3(float a, float b, float c);
GLM_Vertex calc_barycentric_interpolation(Triangle &t, glm::vec3 &intersection);
Color phong_shading(Sphere &s, Light &l, glm::vec3 intersection);
Color phong_shading(Triangle &t, Light &l, glm::vec3 &intersection);
Color calc_phong_shading(glm::vec3 l_dir, glm::vec3 l_color, GLM_Vertex vertex);
bool check_single_sphere_intersection(Ray &r, Sphere &s, float &t, glm::vec3 &intersection);
bool check_single_triangle_intersection(Ray &r, Triangle &tri, float &t, glm::vec3 &intersection);
float calc_triangle_area_xy(glm::vec3 a, glm::vec3 b, glm::vec3 c);
float calc_triangle_area_xz(glm::vec3 a, glm::vec3 b, glm::vec3 c);
float calc_triangle_area(glm::vec3 a, glm::vec3 b, glm::vec3 c);
bool point_triangle_xy(glm::vec3 c0, glm::vec3 c1, glm::vec3 c2, glm::vec3 c);
bool point_triangle_xz(glm::vec3 c0, glm::vec3 c1, glm::vec3 c2, glm::vec3 c);
bool check_block(Ray &r, Sphere &s, Light &l, int i, int j);
bool check_block(Ray &r, Triangle &t, Light &l, int i, int j);
Color check_intersection(Color &c, Ray &ray, int s_i, int t_i, int time);
glm::vec3 calc_reflect_dir(int sphere_i, int triangle_i, glm::vec3 dir, glm::vec3 intersection);
void calc_shadow_ray(Color &color, int sphere_i, int triangle_i, glm::vec3 intersection);
void calc_ray_color(Color &c, int sphere_i, int triangle_i, glm::vec3 intersection, Ray &ray, int time);
void generate_ray(Ray &ray, int x, int y);
void generate_ray_antialiasing(std::vector<Ray> &rays, int x, int y);

void clamp(float &f)
{
  if (f > 1.0f)
    f = 1.0f;
  else if (f < 0.0f)
    f = 0.0f;
}

glm::vec3 vec3(double *v3)
{
  return glm::vec3(v3[0], v3[1], v3[2]);
}

glm::vec3 vec3(float a, float b, float c)
{
  return glm::vec3(a, b, c);
}

void generate_ray(Ray &ray, int x, int y)
{
  float aspect_ratio = WIDTH * 1.0f / HEIGHT;
  float fov_rad = glm::radians(fov);
  float half_fov_tan = std::tan(fov_rad * 0.5f);

  // convert pixel coordinates to normalized device coordinates (NDC)
  // generate a ray to go through the center of a pixel
  float ndcX = (2.0f * (x + 0.5f) / (WIDTH * 1.0f) - 1.0f) * aspect_ratio * half_fov_tan;
  float ndcY = (1.0f - 2.0f * (y + 0.5f) / (HEIGHT * 1.0f)) * half_fov_tan;

  glm::vec3 dirs(ndcX, -ndcY, -1.0f);
  ray.dir = glm::normalize(dirs);
}

void generate_ray_antialiasing(std::vector<Ray> &rays, int x, int y) {

  float aspect_ratio = WIDTH * 1.0f / HEIGHT;
  float fov_rad = glm::radians(fov);
  float half_fov_tan = std::tan(fov_rad * 0.5f);

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      float ndcX = (2.0f * (x + 0.125f + (i * 0.25f)) / (WIDTH * 1.0f) - 1.0f) * aspect_ratio * half_fov_tan;
      float ndcY = (1.0f - 2.0f * (y + 0.125f + (j * 0.25f)) / (HEIGHT * 1.0f)) * half_fov_tan;
      glm::vec3 dirs(ndcX, -ndcY, -1.0f);
      rays[i * 4 + j].dir = glm::normalize(dirs);
    }
  }

}

float calc_triangle_area_xy(glm::vec3 a, glm::vec3 b, glm::vec3 c)
{
  return 0.5f * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

float calc_triangle_area_xz(glm::vec3 a, glm::vec3 b, glm::vec3 c)
{
  return 0.5f * ((b.x - a.x) * (c.z - a.z) - (c.x - a.x) * (b.z - a.z));
}

float calc_triangle_area(glm::vec3 a, glm::vec3 b, glm::vec3 c)
{
  return 0.5f * glm::length(glm::cross(b - a, c - a));
}

bool point_triangle_xy(glm::vec3 c0, glm::vec3 c1, glm::vec3 c2, glm::vec3 c)
{
  float c_c1_c2 = calc_triangle_area_xy(c, c1, c2);
  float c0_c_c2 = calc_triangle_area_xy(c0, c, c2);
  float c0_c1_c = calc_triangle_area_xy(c0, c1, c);
  float c0_c1_c2 = calc_triangle_area_xy(c0, c1, c2);
  float alpha = c_c1_c2 / c0_c1_c2, beta = c0_c_c2 / c0_c1_c2, gamma = 1.0f - alpha - beta;
  return alpha >= 0 && beta >= 0 && gamma >= 0;
}

bool point_triangle_xz(glm::vec3 c0, glm::vec3 c1, glm::vec3 c2, glm::vec3 c)
{
  float c_c1_c2 = calc_triangle_area_xz(c, c1, c2);
  float c0_c_c2 = calc_triangle_area_xz(c0, c, c2);
  float c0_c1_c = calc_triangle_area_xz(c0, c1, c);
  float c0_c1_c2 = calc_triangle_area_xz(c0, c1, c2);
  float alpha = c_c1_c2 / c0_c1_c2, beta = c0_c_c2 / c0_c1_c2, gamma = 1.0f - alpha - beta;
  return alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f;
}

bool check_single_triangle_intersection(Ray &r, Triangle &tri, float &t, glm::vec3 &intersection)
{

  glm::vec3 c0 = vec3(tri.v[0].position);
  glm::vec3 c1 = vec3(tri.v[1].position);
  glm::vec3 c2 = vec3(tri.v[2].position);

  // check if there is an intersection between the plane and ray
  glm::vec3 plane_n = glm::normalize(glm::cross(c1 - c0, c2 - c1));
  float n_dot_d = glm::dot(plane_n, r.dir);

  // light direction is parallel to the plane
  if (abs(n_dot_d) < e)
    return false;

  // calculate the intersection between ray and plane
  float d = -glm::dot(plane_n, c0);
  t = -(glm::dot(plane_n, r.pos) + d) / n_dot_d;
  if (t <= 0.0f)
    return false;
  intersection = r.pos + r.dir * t;

  // check if the intersection is inside the polygon with 2d projection against x-y plane
  glm::vec3 c = intersection;
  bool in_xy = point_triangle_xy(c0, c1, c2, c);
  bool in_xz = point_triangle_xz(c0, c1, c2, c);
  return in_xy || in_xz;
}

bool check_single_sphere_intersection(Ray &r, Sphere &s, float &t, glm::vec3 &intersection)
{

  float x0_xc = r.pos.x - s.position[0];
  float y0_yc = r.pos.y - s.position[1];
  float z0_zc = r.pos.z - s.position[2];
  float a = 1.0f;
  float b = 2 * (r.dir[0] * x0_xc + r.dir[1] * y0_yc + r.dir[2] * z0_zc);
  float c = pow(x0_xc, 2) + pow(y0_yc, 2) + pow(z0_zc, 2) - pow(s.radius, 2);
  float bc = pow(b, 2) - 4 * c;

  if (bc < 0)
    return false;

  float t0 = (-b + sqrt(bc)) / (2 * a);
  float t1 = (-b - sqrt(bc)) / (2 * a);

  if (t0 >= 0.0f && t1 >= 0.0f)
  {
    t = std::min(t0, t1);
  }
  else if (t0 < 0.0f && t1 < 0.0f)
  {
    return false;
  }
  else
  {
    t = std::max(t0, t1);
  }

  intersection = r.pos + r.dir * t;

  return true;
}

Color calc_phong_shading(glm::vec3 l_dir, glm::vec3 l_color, GLM_Vertex vertex)
{

  glm::vec3 n = vertex.n;
  glm::vec3 kd = vertex.kd;
  glm::vec3 ks = vertex.ks;
  glm::vec3 intersection = vertex.position;
  float shininess = vertex.shininess;

  // diffuse
  float l_dot_n = glm::dot(l_dir, n);
  clamp(l_dot_n);

  // specular
  glm::vec3 r = glm::normalize(2.0f * l_dot_n * n - l_dir);
  glm::vec3 v = glm::normalize(-intersection); // eye_dir 0 - intersection;
  float r_dot_v = glm::dot(r, v);
  clamp(r_dot_v);
  float specular = pow(r_dot_v, shininess);

  glm::vec3 c = l_color * (kd * l_dot_n + ks * specular);

  return Color(c);
}

Color phong_shading(Sphere &s, Light &l, glm::vec3 intersection)
{

  GLM_Vertex v;
  glm::vec3 l_pos = vec3(l.position);
  glm::vec3 l_color = vec3(l.color);
  glm::vec3 l_dir = glm::normalize(l_pos - intersection);
  glm::vec3 s_pos = vec3(s.position);
  v.kd = vec3(s.color_diffuse);
  v.ks = vec3(s.color_specular);
  v.n = glm::normalize(intersection - s_pos);
  v.shininess = s.shininess;
  v.position = intersection;

  return calc_phong_shading(l_dir, l_color, v);
}

Color phong_shading(Triangle &t, Light &l, glm::vec3 &intersection)
{

  glm::vec3 l_pos = vec3(l.position);
  glm::vec3 l_color = vec3(l.color);
  glm::vec3 l_dir = glm::normalize(l_pos - intersection);

  // barycentric interpolation
  GLM_Vertex v = calc_barycentric_interpolation(t, intersection);

  return calc_phong_shading(l_dir, l_color, v);
}

bool check_block(Ray &r, Sphere &s, Light &l, int i, int j)
{

  glm::vec3 shadow_ray_intersection;
  float tmp_t = -1.0f;
  glm::vec3 l_pos = vec3(l.position);

  // if the shadow ray hits a sphere that is not the current sphere and there is an intersection
  if (check_single_sphere_intersection(r, s, tmp_t, shadow_ray_intersection) && i != j)
  {
    float intersection_to_light = glm::length(l_pos - r.pos);
    float intersection_to_shadow_intersect = glm::length(shadow_ray_intersection - r.pos);

    // if the light is beyond the object with shadow ray intersection, the object is blocking the light
    if (intersection_to_light > intersection_to_shadow_intersect)
    {
      return true;
    }
  }

  return false;
}

bool check_block(Ray &r, Triangle &t, Light &l, int i, int j)
{

  glm::vec3 shadow_ray_intersection;
  float tmp_t = -1.0f;
  glm::vec3 l_pos = vec3(l.position);

  // if the shadow ray hits a sphere that is not the current sphere and there is an intersection
  if (check_single_triangle_intersection(r, t, tmp_t, shadow_ray_intersection) && i != j)
  {
    float intersection_to_light = glm::length(l_pos - r.pos);
    float intersection_to_shadow_intersect = glm::length(shadow_ray_intersection - r.pos);

    // if the light is beyond the object with shadow ray intersection, the object is blocking the light
    if (intersection_to_light > intersection_to_shadow_intersect)
    {
      return true;
    }
  }

  return false;
}

GLM_Vertex calc_barycentric_interpolation(Triangle &t, glm::vec3 &intersection)
{

  GLM_Vertex vertex;

  Vertex v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
  glm::vec3 c0 = vec3(v0.position), c1 = vec3(v1.position), c2 = vec3(v2.position), c = intersection;

  // barycentric interpolation
  float c_c1_c2 = calc_triangle_area(c, c1, c2);
  float c0_c_c2 = calc_triangle_area(c0, c, c2);
  float c0_c1_c = calc_triangle_area(c0, c1, c);
  float c0_c1_c2 = calc_triangle_area(c0, c1, c2);

  float alpha = c_c1_c2 / c0_c1_c2;
  float beta = c0_c_c2 / c0_c1_c2;
  float gamma = 1.0f - alpha - beta;

  vertex.n = glm::normalize(alpha * vec3(v0.normal) + beta * vec3(v1.normal) + gamma * vec3(v2.normal));
  vertex.kd = alpha * vec3(v0.color_diffuse) + beta * vec3(v1.color_diffuse) + gamma * vec3(v2.color_diffuse);
  vertex.ks = alpha * vec3(v0.color_specular) + beta * vec3(v1.color_specular) + gamma * vec3(v2.color_specular);
  vertex.shininess = alpha * v0.shininess + beta * v1.shininess + gamma * v2.shininess;
  vertex.position = intersection;

  return vertex;
}

void calc_shadow_ray(Color &color, int sphere_i, int triangle_i, glm::vec3 intersection)
{

  color = Color();

  // check out each light source
  for (int i = 0; i < num_lights; ++i)
  {

    // create shadow ray
    Light &light = lights[i];
    glm::vec3 l_pos = vec3(light.position);
    glm::vec3 shadow_ray_dir = glm::normalize(l_pos - intersection);
    Ray shadow_ray(shadow_ray_dir, intersection);
    bool blocked = false;

    // check to see if shadow_ray is blocked by any spheres
    for (int j = 0; j < num_spheres && !blocked; ++j)
    {
      if (check_block(shadow_ray, spheres[j], light, sphere_i, j))
      {
        blocked = true;
      }
    }

    // check to see if shadow ray is blocked by any triangles
    for (int j = 0; j < num_triangles && !blocked; ++j)
    {
      if (check_block(shadow_ray, triangles[j], light, triangle_i, j))
      {
        blocked = true;
      }
    }

    // if the ray is not block then calculate phong shading
    if (!blocked)
    {
      if (triangle_i == -1)
        color += phong_shading(spheres[sphere_i], light, intersection);
      else
        color += phong_shading(triangles[triangle_i], light, intersection);
    }
  }
}

void calc_ray_color(Color &c, int sphere_i, int triangle_i, glm::vec3 intersection, Ray &ray, int time)
{
  calc_shadow_ray(c, sphere_i, triangle_i, intersection);

  if (time > 0)
  {
    glm::vec3 r = calc_reflect_dir(sphere_i, triangle_i, ray.dir, intersection);
    Ray reflect_ray(r, intersection);
    Color color;
    Color reflect_color = check_intersection(color, reflect_ray, sphere_i, triangle_i, time - 1);
    glm::vec3 ks;

    if (sphere_i != -1)
    {
      Sphere &s = spheres[sphere_i];
      ks = vec3(s.color_specular);
    }
    else
    {
      Triangle &t = triangles[triangle_i];
      GLM_Vertex v = calc_barycentric_interpolation(t, intersection);
      ks = v.ks;
    }

    c = c * (1.0f - ks) + reflect_color * ks;
  }
}

Color check_intersection(Color &c, Ray &ray, int s_i, int t_i, int time)
{
  float sphere_t = 0.0f, triangle_t = 0.0f;
  int sphere_i = -1, triangle_i = -1;
  Color &color = c;
  glm::vec3 s_intersection, t_intersection;

  // find the closest sphere intersect with ray
  for (int i = 0; i < num_spheres; ++i)
  {
    float tmp_t = -1.0f;

    if (i == s_i)
      continue;
    if (!check_single_sphere_intersection(ray, spheres[i], tmp_t, s_intersection))
      continue;

    // calculate the intersection in the smallest t
    if (tmp_t < sphere_t || sphere_i == -1)
    {
      sphere_t = tmp_t;
      sphere_i = i;
    }
  }

  // find the cloest triangle intersection with ray
  for (int i = 0; i < num_triangles; ++i)
  {
    float tmp_t = -1.0f;

    if (i == t_i)
      continue;
    if (!check_single_triangle_intersection(ray, triangles[i], tmp_t, t_intersection))
      continue;

    // calculate the intersection in the smallest t
    if (tmp_t < triangle_t || triangle_i == -1)
    {
      triangle_t = tmp_t;
      triangle_i = i;
    }
  }

  s_intersection = ray.pos + ray.dir * sphere_t;
  t_intersection = ray.pos + ray.dir * triangle_t;

  if (triangle_i == -1 && sphere_i != -1)
  {
    calc_ray_color(color, sphere_i, -1, s_intersection, ray, time);
  }
  else if (triangle_i != -1 && sphere_i == -1)
  {
    calc_ray_color(color, -1, triangle_i, t_intersection, ray, time);
  }
  else if (triangle_i != -1 && sphere_i != -1)
  {
    if (sphere_t < triangle_t)
    {
      calc_ray_color(color, sphere_i, -1, s_intersection, ray, time);
    }
    else
    {
      calc_ray_color(color, -1, triangle_i, t_intersection, ray, time);
    }
  }

  return color;
}

glm::vec3 calc_reflect_dir(int sphere_i, int triangle_i, glm::vec3 dir, glm::vec3 intersection)
{

  glm::vec3 n;

  if (sphere_i != -1)
  {
    glm::vec3 s_pos = vec3(spheres[sphere_i].position);
    n = glm::normalize(intersection - s_pos);
  }
  else
  {
    Triangle &t = triangles[triangle_i];

    GLM_Vertex v = calc_barycentric_interpolation(t, intersection);
    n = v.n;
  }

  glm::vec3 l_dir = -dir;
  float l_dot_n = glm::dot(l_dir, n);
  clamp(l_dot_n);
  glm::vec3 r = glm::normalize(2.0f * l_dot_n * n - l_dir);

  return r;
}

Color tracing(int x, int y)
{
  // single ray per pixel
  // Color color(1.0f, 1.0f, 1.0f);

  // Ray ray;
  // generate_ray(ray, x, y);
  // check_intersection(color, ray, -1, -1, MAX_REFLECT);

  // anti-aliasing with 16 rays per pixel
  Color color(0.0f, 0.0f, 0.0f);
  std::vector<Ray> rays;
  for (int i = 0; i < 16; ++i) {
    Ray ray;
    rays.push_back(ray);
  }
  generate_ray_antialiasing(rays, x, y);

  for (int i = 0; i < 16; ++i) {
    Color c(1.0f, 1.0f, 1.0f);
    check_intersection(c, rays[i], -1, -1, MAX_REFLECT);
    color += c;
  }

  color /= 16.0f;

  return color;
}

void draw_pixel(int x, int y) {
  Color color = tracing(x, y);
  color += Color(0.1f, 0.1f, 0.1f); // add ambient light
  color.clamp();
  plot_pixel(x, y, color.r * 255, color.g * 255, color.b * 255);
}

// MODIFY THIS FUNCTION
void draw_scene()
{

  for (unsigned int x = 0; x < WIDTH; x++)
  {
    glPointSize(2.0);
    glBegin(GL_POINTS);

    for (unsigned int y = 0; y < HEIGHT; y++)
    {
      draw_pixel(x, y);
    }

    glEnd();
    glFlush();
  }



  printf("Done!\n");
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x, y, r, g, b);
  if (mode == MODE_JPEG)
    plot_pixel_jpeg(x, y, r, g, b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if (strcasecmp(expected, found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE *file, const char *check, double p[3])
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check(check, str);
  fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
  printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check("rad:", str);
  fscanf(file, "%lf", r);
  printf("rad: %f\n", *r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file, "%s", s);
  parse_check("shi:", s);
  fscanf(file, "%lf", shi);
  printf("shi: %f\n", *shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv, "r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file, "%i", &number_of_objects);

  printf("number of objects: %i\n", number_of_objects);

  parse_doubles(file, "amb:", ambient_light);

  for (int i = 0; i < number_of_objects; i++)
  {
    fscanf(file, "%s\n", type);
    printf("%s\n", type);
    if (strcasecmp(type, "triangle") == 0)
    {
      printf("found triangle\n");
      for (int j = 0; j < 3; j++)
      {
        parse_doubles(file, "pos:", t.v[j].position);
        parse_doubles(file, "nor:", t.v[j].normal);
        parse_doubles(file, "dif:", t.v[j].color_diffuse);
        parse_doubles(file, "spe:", t.v[j].color_specular);
        parse_shi(file, &t.v[j].shininess);
      }

      if (num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if (strcasecmp(type, "sphere") == 0)
    {
      printf("found sphere\n");

      parse_doubles(file, "pos:", s.position);
      parse_rad(file, &s.radius);
      parse_doubles(file, "dif:", s.color_diffuse);
      parse_doubles(file, "spe:", s.color_specular);
      parse_shi(file, &s.shininess);

      if (num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if (strcasecmp(type, "light") == 0)
    {
      printf("found light\n");
      parse_doubles(file, "pos:", l.position);
      parse_doubles(file, "col:", l.color);

      if (num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n", type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  // hack to make it only draw once
  static int once = 0;
  if (!once)
  {
    draw_scene();
    if (mode == MODE_JPEG)
      save_jpg();
  }
  once = 1;
}

int main(int argc, char **argv)
{
  if ((argc < 2) || (argc > 3))
  {
    printf("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if (argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if (argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc, argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(WIDTH, HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
#ifdef __APPLE__
  // This is needed on recent Mac OS X versions to correctly display the window.
  glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
#endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
