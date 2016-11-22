Mac OS X   :  
Before you can use a dynamic library as a dependent library, 
the library and its header files must be installed on your computer. 
The standard locations for header files are ~/include, /usr/local/include and /usr/include. 
The standard locations for dynamic libraries are ~/lib, /usr/local/lib, and /usr/lib.
% g++ main.cpp -o magnoom -O3 -Wall -fno-strict-aliasing -I ./include/ -L ./lib -lAntTweakBar -framework GLUT -pthread -framework OpenGL -Wno-deprecated-declarations
% echo "" | gcc -xc - -v -E
% sudo cp ~/..../AntTweakBar/lib/libAntTweakBar.dylib /usr/local/lib
% sudo cp ~/..../AntTweakBar/include/AntTweakBar.h /usr/local/include/
% sudo cp -r ~/..../AntTweakBar /usr/local/
% g++ main.cpp -o magnoom -O3 -Wall -fno-strict-aliasing -lAntTweakBar -framework GLUT -pthread -framework OpenGL -Wno-deprecated-declarations
Note, in OS X: https://lukecyca.com/2008/glui-235-framework-for-mac-os-x.html
DYLD_LIBRARY_PATH=DYLD_LIBRARY_PATH=/usr/local/opt/gcc/lib/gcc/6/    

Raspberry Pi3:
pi@raspberrypi: sudo apt-get install libgl1-mesa-dev libgles2-mesa-dev libglew-dev:armhf libglewmx-dev:armhf libglib2.0-dev libglu1-mesa-dev
pi@raspberrypi: g++ main.cpp -o magnoom -pthread -O3 -Wall -fno-strict-aliasing -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL

Ubuntu   :  
$ g++ main.cpp -o magnoom -pthread -O3 -Wall -fno-strict-aliasing -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL 
$ sudo export PATH="/usr/bin:$PATH" {/usr/local/bin/ld: this linker was not configured to use sysroots}
CentOS   : 
$ g++ JS2v7.cpp -o JS2v7 -pthread -O3 -Wall -fno-strict-aliasing -I ./include -L ./lib -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL
     LD_LIBRARY_PATH=~/Desktop/JSpinx4/lib/
     LD_LIBRARY_PATH=~/Desktop/JSpinx4/include/