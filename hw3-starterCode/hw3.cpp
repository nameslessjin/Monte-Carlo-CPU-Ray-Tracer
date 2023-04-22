/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Your name here>
 * *************************
 */

#include "hw3.hpp"

glm::vec3 maxPointSphere(Sphere &s)
{
  glm::vec3 center = vec3(s.position);
  float r = s.radius;
  glm::vec3 addition = vec3(r, r, r);
  glm::vec3 error = vec3(sigma, sigma, sigma);
  return center + addition + error;
}

glm::vec3 minPointSphere(Sphere &s)
{
  glm::vec3 center = vec3(s.position);
  float r = s.radius;
  glm::vec3 subtraction = vec3(-r, -r, -r);
  glm::vec3 error = vec3(sigma, sigma, sigma);
  return center + subtraction - error;
}

glm::vec3 maxPointTriangle(Triangle &t)
{
  Vertex v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
  glm::vec3 pos(std::max({v0.position[0], v1.position[0], v2.position[0]}),
                std::max({v0.position[1], v1.position[1], v2.position[1]}),
                std::max({v0.position[2], v1.position[2], v2.position[2]}));
  glm::vec3 error = vec3(sigma, sigma, sigma);

  return pos + sigma;
}

glm::vec3 minPointTriangle(Triangle &t)
{
  Vertex v0 = t.v[0], v1 = t.v[1], v2 = t.v[2];
  glm::vec3 pos(
      std::min({v0.position[0], v1.position[0], v2.position[0]}),
      std::min({v0.position[1], v1.position[1], v2.position[1]}),
      std::min({v0.position[2], v1.position[2], v2.position[2]}));
  glm::vec3 error = vec3(sigma, sigma, sigma);
  return pos - error;
}

int countHBV(AABB *hbv, int level)
{

  int count = 0;

  if (!hbv)
    return count;

  std::cout << "parent: , level: " << level << '\n';
  hbv->print();
  if (hbv->left)
  {
    std::cout << "left children: \n";
    hbv->left->print();
  }
  if (hbv->right)
  {
    std::cout << "right children: \n";
    hbv->right->print();
  }

  ++count;
  count += countHBV(hbv->left, level + 1);
  count += countHBV(hbv->right, level + 1);

  return count;
}

void contructHVB()
{
  for (int i = 0; i < num_spheres; ++i)
  {
    Sphere &s = spheres[i];
    glm::vec3 sphere_min = minPointSphere(s);
    glm::vec3 sphere_max = maxPointSphere(s);
    AABB aabb(sphere_min, sphere_max, i, -1);
    all_aabbs.push_back(aabb);
  }

  for (int i = 0; i < num_triangles; ++i)
  {
    Triangle &t = triangles[i];
    glm::vec3 t_min = minPointTriangle(t);
    glm::vec3 t_max = maxPointTriangle(t);
    AABB aabb(t_min, t_max, -1, i);
    all_aabbs.push_back(aabb);
  }

  // for (int i = 0; i < all_aabbs.size(); ++i) {
  //   all_aabbs[i].print();
  // }

  std::cout << "num all_aabbs: " << all_aabbs.size() << '\n';
  std::vector<AABB> hbv = buildHVB(all_aabbs);

  if (hbv.size() == 1)
  {
    HBV = new AABB(hbv[0]);
    // int count = countHBV(HBV, 0);
    // std::cout << "succ " << count << '\n';
  }
  // else if (hbv.size() == 0) {
  //   std::cout << "empty\n";
  // } else {
  //   std::cout << "suck\n";
  // }
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

void generate_ray_antialiasing(std::vector<Ray> &rays, int x, int y, int num_samples)
{

  float aspect_ratio = WIDTH * 1.0f / HEIGHT;
  float fov_rad = glm::radians(fov);
  float half_fov_tan = std::tan(fov_rad * 0.5f);
  int rt_num_samples = sqrt(num_samples);

  float increment = 1.0f / rt_num_samples;
  float start = increment / 2;

  for (int i = 0; i < rt_num_samples; ++i)
  {
    for (int j = 0; j < rt_num_samples; ++j)
    {
      float ndcX = (2.0f * (x + start + (i * increment)) / (WIDTH * 1.0f) - 1.0f) * aspect_ratio * half_fov_tan;
      float ndcY = (1.0f - 2.0f * (y + start + (j * increment)) / (HEIGHT * 1.0f)) * half_fov_tan;
      glm::vec3 dirs(ndcX, -ndcY, -1.0f);
      rays[i * rt_num_samples + j].dir = glm::normalize(dirs);
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

  if (abs(glm::length(intersection) - glm::length(r.pos)) < e)
    return false;

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
  float bc = pow(b, 2) - 4 * a * c;

  if (bc < 0.0f)
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

  if (abs(glm::length(intersection) - glm::length(r.pos)) < e)
    return false;

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

bool check_block(Ray &r, Sphere &s, Light &l)
{

  glm::vec3 shadow_ray_intersection;
  float tmp_t = -1.0f;
  glm::vec3 l_pos = vec3(l.position);

  // if the shadow ray hits a sphere that is not the current sphere and there is an intersection
  if (check_single_sphere_intersection(r, s, tmp_t, shadow_ray_intersection))
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

bool check_block(Ray &r, Triangle &t, Light &l)
{

  glm::vec3 shadow_ray_intersection;
  float tmp_t = -1.0f;
  glm::vec3 l_pos = vec3(l.position);

  // if the shadow ray hits a sphere that is not the current sphere and there is an intersection
  if (check_single_triangle_intersection(r, t, tmp_t, shadow_ray_intersection))
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

void calc_shadow_ray(Color &color, int sphere_i, int triangle_i, glm::vec3 &intersection)
{

  color = Color();

  // check out each light source
  for (int i = 0; i < num_lights; ++i)
  {

    // create area light for softshadow
    Light &light = lights[i];

    int x_count = 10, y_count = 1, z_count = 10;
    float multiplier = 1.0f * x_count * y_count * z_count;
    for (float x = -1.0f; x < 1.0f; x += 0.2)
    {
      for (float y = 0.0f; y < 1.0f; y += 10.0f)
      {
        for (float z = -1.0f; z < 1.0f; z += 0.2)
        {
          Light new_light;
          new_light.position[0] = light.position[0] + x;
          new_light.position[1] = light.position[1] + y;
          new_light.position[2] = light.position[2] + z;
          new_light.color[0] = light.color[0] / multiplier;
          new_light.color[1] = light.color[1] / multiplier;
          new_light.color[2] = light.color[2] / multiplier;

          // create shadow ray
          glm::vec3 l_pos = vec3(new_light.position);
          glm::vec3 shadow_ray_dir = glm::normalize(l_pos - intersection);
          Ray shadow_ray(shadow_ray_dir, intersection);
          bool blocked = false;

          // check to see if shadow_ray is blocked by any spheres
          for (int j = 0; j < num_spheres; ++j)
          {
            if (blocked) break;
            if (check_block(shadow_ray, spheres[j], new_light))
            {
              blocked = true;
            }
          }

          // check to see if shadow ray is blocked by any triangles
          for (int j = 0; j < num_triangles; ++j)
          {
            if (blocked) break;
            if (check_block(shadow_ray, triangles[j], new_light))
            {
              blocked = true;
            }
          }

          // if the ray is not block then calculate phong shading
          if (!blocked)
          {
            if (triangle_i == -1)
              color += phong_shading(spheres[sphere_i], new_light, intersection);
            else
              color += phong_shading(triangles[triangle_i], new_light, intersection);
          }
        }
      }
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
    Color reflect_color = check_intersection(color, reflect_ray, time - 1);
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

bool isAABB1Closer(AABB *aabb1, AABB *aabb2, Ray &ray) {


  int sphere_t_1 = aabb1->sphere_i;
  int triangle_t_1 = aabb1->triangle_i;
  bool aabb1_intersect = false;
  float intersection_t_1;
  float distance1;
  glm::vec3 interaction1;

  if (sphere_t_1 != -1) {
    aabb1_intersect = check_single_sphere_intersection(ray, spheres[sphere_t_1], intersection_t_1, interaction1);
  }
  if (triangle_t_1 != -1) {
    aabb1_intersect = check_single_triangle_intersection(ray, triangles[triangle_t_1], intersection_t_1, interaction1);
  }

  int sphere_t_2 = aabb2->sphere_i;
  int triangle_t_2 = aabb2->triangle_i;
  bool aabb2_intersect = false;
  float intersection_t_2;
  float distance2;
  glm::vec3 interaction2;

  if (sphere_t_2 != -1) {
    aabb2_intersect = check_single_sphere_intersection(ray, spheres[sphere_t_2], intersection_t_2, interaction2);
  }
  if (triangle_t_2 != -1) {
    aabb2_intersect = check_single_triangle_intersection(ray, triangles[triangle_t_2], intersection_t_2, interaction2);
  }

  // if ray intersects with both objects, pick the one takes less time
  if (aabb1_intersect && aabb2_intersect) {
    return intersection_t_1 < intersection_t_2;
  }

  // if ray doesn't not intersect with object 1 return false
  if (!aabb1_intersect && aabb2_intersect) return false;

  // if ray doesn't not intersect with object 2 return false
  if (aabb1_intersect && !aabb2_intersect) return true;

  return true;

}

bool checkIntersectionWithAABB(AABB *aabb, AABB *(&intersected_aabb), Ray &ray)
{

  // if aabb is null or ray doesn't intersection with it return false
  if (!aabb || !aabb->intersect(ray))
    return false;

  // if ray intersecting with a left return true
  if (!aabb->left && !aabb->right)
  {
    if (!intersected_aabb || (intersected_aabb && !isAABB1Closer(intersected_aabb, aabb, ray))) {
      intersected_aabb = aabb;
    }
    return true;
  }

  bool left_intersection = false, right_intersection = false;

  if (aabb->left)
  {
    left_intersection = checkIntersectionWithAABB(aabb->left, intersected_aabb, ray);
  }
  if (aabb->right)
  {
    right_intersection = checkIntersectionWithAABB(aabb->right, intersected_aabb, ray);
  }

  return left_intersection || right_intersection;
}

Color check_intersection(Color &c, Ray &ray, int time)
{
  float sphere_t = 0.0f, triangle_t = 0.0f;
  int sphere_i = -1, triangle_i = -1;
  Color &color = c;
  glm::vec3 s_intersection, t_intersection;
  float tmp_t = -1.0f;

  AABB *intersected_aabb;
  bool has_intersection = checkIntersectionWithAABB(HBV, intersected_aabb, ray);

  if (!has_intersection) return color;

  sphere_i = intersected_aabb->sphere_i;
  triangle_i = intersected_aabb->triangle_i;

  if (sphere_i != -1 && check_single_sphere_intersection(ray, spheres[sphere_i], tmp_t, s_intersection)) {
      sphere_t = tmp_t;
  }

  if (triangle_i != -1 && check_single_triangle_intersection(ray, triangles[triangle_i], tmp_t, t_intersection)) {
    triangle_t = tmp_t;
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
  int num_samples = ANTI_ALIASING_SAMPLE;
  std::vector<Ray> rays;

  for (int i = 0; i < num_samples; ++i)
  {
    Ray ray;
    rays.push_back(ray);
  }
  generate_ray_antialiasing(rays, x, y, num_samples);

  for (int i = 0; i < num_samples; ++i)
  {
    Color c(1.0f, 1.0f, 1.0f);
    check_intersection(c, rays[i], MAX_REFLECT);
    color += c;
  }

  color /= num_samples * 1.0f;

  return color;
}

void draw_pixel(int x, int y, std::atomic<int> &finished)
{
  Color color = tracing(x, y);
  color += Color(0.1f, 0.1f, 0.1f); // add ambient light
  color.clamp();

  img[y][x][0] = color.r;
  img[y][x][1] = color.g;
  img[y][x][2] = color.b;

  ++finished;
  print_progress(finished, WIDTH * HEIGHT, mtx);
}

void fill_image_plane()
{
  int num_threads = std::thread::hardware_concurrency();
  std::atomic<int> finished(0);
  int cols_thread = WIDTH / num_threads;
  std::thread threads[num_threads];

  print_progress(finished, WIDTH * HEIGHT, mtx);

  ThreadPool thread_pool(num_threads);

  // HBV->print();
  for (unsigned int x = 0; x < WIDTH; ++x)
  {
    for (unsigned int y = 0; y < HEIGHT; ++y)
    {
      thread_pool.enqueue(draw_pixel, x, y, std::ref(finished));
    }
  }

  thread_pool.wait();

  // print final progress
  print_progress(finished, WIDTH * HEIGHT, mtx);
  std::cout << '\n';
}

// MODIFY THIS FUNCTION
void draw_scene()
{
  contructHVB();
  fill_image_plane();

  glPointSize(2.0);
  glBegin(GL_POINTS);

  for (unsigned int x = 0; x < WIDTH; x++)
  {

    for (unsigned int y = 0; y < HEIGHT; y++)
    {
      plot_pixel(x, y, img[y][x][0] * 255, img[y][x][1] * 255, img[y][x][2] * 255);
    }
  }

  glEnd();
  glFlush();

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
