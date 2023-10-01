Monte-Carlo CPU Ray Tracer



https://github.com/nameslessjin/Monte-Carlo-CPU-Ray-Tracer/assets/42518630/34708312-e49b-45d0-b424-d45dcbdd3dfb


## Introduction
In this project, we created a ray tracer that utilizes backward ray casting from the camera.  We shoot rays to the image panel on each pixel and trace the rays to their intersection with the closest objects.  After we find the intersections, we will cast a shadow ray to determine the color through phong shading and obstacle  checking.  If the shadow ray is blocked by another object then the color will be black which becomes shadow else we use Phong shading to find the correct color.

## Extra Features
1. Recursive Reflection                                                                  
2. Good antialiasing with stratified_sampling                                            
3. Soft shadows with random sampling for each area light                               
4. Bounding Volume Hierarchy with AABB                                               
5. Monte-Carlo Sampling                                                                
6. Multi-threading                                                               
