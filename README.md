# Volume Path Tracer (Final)

#### Group Member: Ruan Hang, Yijin Guo, Yuhao Wang

### **Summary**

In this project, we will build a volumetric path tracer that can handle scattering and absorption inside participating media inside lajolla. We will split the development into 5 steps and build 5 volumetric path tracers. Each has more features than the previous ones.After complementing the renderer, we introduce the atmospheric medium to simulate the natural phenomenon.Then we render some scene files for the simulation of sunrise and generate a video() for it.

### Build

Use CMake to build. 

```
mkdir build
cd build
cmake ..
cmake --build .
```

It requires compilers that support C++17 (gcc version >= 8, clang  version >= 7, Apple Clang version >= 11.0, MSVC version >=  19.14).

### Run

To render a single file, 

```
bash run.sh -i <input_file.xml> [-o file_name]
# example: bash run.sh -i ./scenes/volpath_test/volpath_test1.xml -o test1 => get a rendered result "text1.exr"
```

To generate a video,

```
bash animation.sh <scene_folder>
```

This will output a video `output.mp4`. 

### **Resources List**

[Volumetric Path Tracing](https://cseweb.ucsd.edu/~tzli/cse272/wi2023/homework2.pdf) - Jan Novák,liyan Georgiev,Johannes Hanika,Jaroslav Křiváne, Wojciec

[Monte Carlo methods for physically based volume rendering](https://cs.dartmouth.edu/~wjarosz/publications/novak18monte-sig.html)

[Display of The Earth Taking into Account Atmospheric Scattering](http://nishitalab.org/user/nis/cdrom/sig93_nis.pdf)

[Precomputed Atmospheric Scattering](https://hal.inria.fr/inria-00288758v1/document)