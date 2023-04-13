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

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

#define DEG_TO_RAD(degrees) (degrees * M_PI / 180.0);

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
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

struct Color {
  float r, g, b;

  Color() {r = 0; g = 0; b = 0;}
  Color(float _r, float _g, float _b) {r = _r, g =_g; b = _b;}

  Color &clamp() {
    if (r > 1.0) r = 1.0;
    if (r < 0.0) r = 0.0;
    if (g > 1.0) g = 1.0;
    if (g < 0.0) g = 0.0;
    if (b > 1.0) b = 1.0;
    if (b < 0.0) b = 0.0;

    return *this;
  }

  Color &operator+(Color &c) {return Color(c.r + r, c.g + g, c.b + b).clamp();}
  Color &operator+=(Color const &c) {return Color(c.r + r, c.g + g, c.b + b).clamp();}
  void print() {
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

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);



struct Ray {
  glm::vec3 dir;
  glm::vec3 pos;

  Ray() {dir = glm::vec3(0,0,0);}
  Ray(glm::vec3 _dir, glm::vec3 _pos) { dir = _dir; pos = _pos;}

  void generate_ray(int x, int y) {
    float aspect_ratio = WIDTH * 1.0f / HEIGHT;
    float fov_rad = DEG_TO_RAD(fov);
    float half_fov_tan = std::tan(fov_rad * 0.5f);

    // convert pixel coordinates to normalized device coordinates (NDC)
    float ndcX = (2.0f * (x + 0.5f) / (WIDTH * 1.0f) - 1.0f) * aspect_ratio * half_fov_tan;
    float ndcY = (1.0f - 2.0f * (y + 0.5f) / (HEIGHT * 1.0f)) * half_fov_tan;

    glm::vec3 dirs (ndcX, ndcY, -1.0f);
    dir = glm::normalize(dirs);

  }

  bool check_single_sphere_intersection(Sphere &s, float &t, glm::vec3 &intersection) {

      float x0_xc = pos.x - s.position[0];
      float y0_yc = pos.y - s.position[1];
      float z0_zc = pos.z - s.position[2];
      float a = 1.0f;
      float b = 2 * (dir[0] * x0_xc + y0_yc + z0_zc);
      float c = pow(x0_xc, 2) + pow(y0_yc, 2) + pow(z0_zc, 2) - pow(s.radius, 2);

      if (pow(b, 2) - 4 * c < 0) return false;

      float t0 = (-b + sqrt(pow(b, 2) - 4 * c)) / 2;
      float t1 = (-b - sqrt(pow(b, 2) - 4 * c)) / 2;

      if (t0 >= 0 && t1 >= 0) {
        t = std::min(t0, t1);
      } else if (t0 < 0 && t1 < 0) {
        return false;
      } else {
        t = std::max(t0, t1);
      }

      intersection = pos + dir * t;


      return true;
  }

};


void clamp(float &f);
Color phong_shading(Sphere &s, Light &l, glm::vec3 intersection);
Color check_spheres_intersection(Color &c, Ray &ray);

void clamp(float &f) {
  if (f > 1.0) f = 1.0;
  else if (f < 0.0f) f = 0.0;
}

Color phong_shading(Sphere &s, Light &l, glm::vec3 intersection) {

    // diffuse
    glm::vec3 l_pos(l.position[0], l.position[1], l.position[2]);
    glm::vec3 s_pos(s.position[0], s.position[1], s.position[2]);
    glm::vec3 l_color(l.color[0], l.color[1], l.color[2]);
    glm::vec3 kd(s.color_diffuse[0], s.color_diffuse[1], s.color_diffuse[2]);
    glm::vec3 ks(s.color_specular[0], s.color_specular[1], s.color_specular[2]);
    glm::vec3 l_dir = glm::normalize(l_pos - intersection);
    glm::vec3 n = glm::normalize(intersection - s_pos);
    float l_dot_n = glm::dot(l_dir, n);
    clamp(l_dot_n);

    // specular
    glm::vec3 r = glm::normalize(2 * l_dot_n * n - l_dir);
    glm::vec3 v = glm::normalize(-intersection); // eye_dir 0 - intersection;
    float r_dot_v = glm::dot(r, v);
    float specular = pow(r_dot_v, s.shininess);

    glm::vec3 c = l_color * (kd * l_dot_n + ks * specular);

    return Color(c.r, r.g, r.b);      
}

Color check_spheres_intersection(Color &c, Ray &ray) {

  float t = 0.0f;
  int intersection_index = -1;
  Color &color = c;
  glm::vec3 intersection;

  // find the closest sphere intersect with ray
  for (int i = 0; i < num_spheres; ++i) {

    float tmp_t = -1.0f;

    if (!ray.check_single_sphere_intersection(spheres[i], tmp_t, intersection)) continue;

    // calculate the intersection in the smallest t
    if (tmp_t < t || intersection_index == -1) {
      t = tmp_t;
      intersection_index = i;
      intersection = ray.pos + ray.dir * t;
    }
  }

  if (intersection_index != -1) {

    // check out each light source
    for (int i = 0; i < num_lights; ++i) {

      // create shadow ray
      Light &light = lights[i];
      glm::vec3 l_pos(light.position[0], light.position[1], light.position[2]);
      glm::vec3 shadow_ray_dir = l_pos - intersection;
      Ray shadow_ray(shadow_ray_dir, intersection);
      bool in_shadow = false;

      // check to see if shadow_ray is blocked by any spheres
      for (int j = 0; j < num_spheres; ++j) {
        glm::vec3 shadow_ray_intersection;
        float tmp_t = -1.0f;
        
        // if the shadow ray hits a sphere that is not the current sphere and there is an intersection
        if (shadow_ray.check_single_sphere_intersection(spheres[j], tmp_t, shadow_ray_intersection) && intersection_index != j) {
          float intersection_to_light = glm::length(l_pos - intersection);
          float intersection_to_shadow_intersect = glm::length(shadow_ray_intersection - intersection);

          // if the light is beyond the object with shadow ray intersection, the object is blocking the light
          if (intersection_to_light - intersection_to_shadow_intersect) {
            in_shadow = true;
            break;
          }
        }
      }

      if (!in_shadow) {
        color += phong_shading(spheres[intersection_index], lights[i], intersection);
      }

    }
    
  }

  return color;
}

Color tracing(int x, int y) {

  Ray ray;
  Color color;

  ray.generate_ray(x, y);
  check_spheres_intersection(color, ray);

  return color;

}


//MODIFY THIS FUNCTION
void draw_scene()
{

  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);

    for(unsigned int y=0; y<HEIGHT; y++)
    {
      Color color = tracing(x, y);
      color.print();
      plot_pixel(x, y, color.r * 255, color.g * 255, color.b * 255);
    }


    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
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
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
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
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
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

