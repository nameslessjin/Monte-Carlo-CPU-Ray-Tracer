# CSCI 420 Programming Assignment 3: Ray Tracing - Monte Carlo Sampling

Subject             : CSCI420 - Computer Graphics <br>
Assignment 3        : Ray Tracing <br>
Author              : Jinsen Wu                   <br>
USC ID              : 9338292958                  <br>

## Introduction
In this assignment, we create a ray tracer that utilize backward ray casting from camera.  We will shoot rays to the image panel on each pixels and trace the rays to their intersection with the cloest objects.  After we find the intersections, we will cast shadow ray to determine the color through phong shading and obsticle checking.  If the shadow ray is blocked by another object then the color will be black which becomes shadow else we use phong shading to find the correct color.

## Core Credit Features
1. Send out rays from the camera at (0, 0, 0).  Use backwards ray tracing on each pixel.    -   Y
2. Write the intersection code.                                                             -   Y
    - Triangle Intersection                                                                 -   Y
    - Sphere Intersection                                                                   -   Y
3. Implement the illumination equations                                                     -   Y
    - Triangle Phong Shading                                                                -   Y
    - Sphere Phong Shading                                                                  -   Y
4. Shadow rays                                                                              -   Y
5. Create still images showing off your ray tracer                                          -   Y

## Extra Credits Features
1. Recursive Reflection                                                                     -   Y
2. Good antialiasing                                                                        -   Y
3. Soft shadows                                                                             -   Y
4. Hierachical Bounding Volume with AABB                                                    -   Y
5. Monte-Carlo Sampling                                                                     -   Y
6. Multi-threading                                                                          -   Y

## Note
This program is for monte carlo sampling.  For the core phong shading, please check out the other program in the submission named "USC-CSCI420-HW3-core-credit".  I will submit a JPEG folders to contain still images from both programs.  The example of monte carlo sampling is also provided in JPEG folder under hw3-starterCode.
