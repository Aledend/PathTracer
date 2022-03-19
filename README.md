# PathTracer

August 2021

/* Disclaimer
This project is an adaptation of Peter Shirley's Ray Tracing in One Weekend: https://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf. The document provided comes with a clear picture of what a Ray Tracer is, but lacks in completeness, as not everything is provided code wise. I highly recommend this document for learning and finding/filling in the missing parts yourself, but for a start to finish code walk through, I would recommend https://raytracing.github.io/books/RayTracingInOneWeekend.html.
*/

This is my attempt at learning compute shaders and using them to accelerate the rendering of a path-tracer.

I was quickly fascinated in the beautiful images this path tracer could provide. But there was one thing wrong with it. As complexity increased, it grew dreadfully slow. My own twist on this project is the converting of everything to using compute shaders. Offloading all that cpu work load to the gpu was challenging, but provided amazing results. My frame-time went from minutes per frame to frames per second. It even provides beautiful results at steady framerate allowing for interactive camera movement.

Example of rendering a 1200x800 window, running 60 passes per pixel with a bounce depth of 50:
![pathtrace](https://user-images.githubusercontent.com/16718975/130252463-94a20c81-f02e-47c2-ad6e-fed54764e80a.gif)



However, there is still much that can be done. The compute shader was thrown together as a proof of concept, there are very likely many better ways to provide uniform random behaviour, and a lot can be done to relieve the computing that has to be done. I would start by building a collision tree and relieve the gpu from passing through as many spheres as possible. Most performance costs originate in the light-bouncing however, so time will tell whether a collision tree would pay off in such a compact area, or if it only adds complexity. I would also look into sorting the spheres by material, and see if the compute shader can be relieved from the severe branching that comes with handling different materials.

I had great fun playing around with this project, and will for sure revisit it in the future!

![pathtrace image](https://user-images.githubusercontent.com/16718975/130251813-46e7e57b-5a5e-48f7-8b43-649a01afdfe7.png)
