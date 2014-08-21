RayTracing
==========

Simple Ray Tracer

![alt text](sample_output.png "Sample output")

Compile It
----------

make

Run It
------

./RayTrace ./input/simple/input12.txt


About
-----

A Phong reflection model is used for local illumination.  Shadows are present, although simple.  If you would like to run the program on an entire directory:

ls ./input/simple/*.txt | xargs -I {} ./RayTrace {}
