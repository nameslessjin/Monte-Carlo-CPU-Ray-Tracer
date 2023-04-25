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
3. Implement the illumination equations                                                     -   Y
4. Create still images showing off your ray tracer                                          -   Y

## Extra Credits Features
