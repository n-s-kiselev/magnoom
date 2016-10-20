#Magnoom

Advanced software for Atomistic Spin Dynamics

--------------------

Main features:
  * platform independent code writen in C/C++
  * minimal dependencies on external libraries
  * GUI entirely based on OpenGL 
  * real time control of the parameters

##Dependencies:
  * gcc >=4.8.1 (possibly can be compiled also with earlier versions of gcc)
  * OpenGL >= 3.0 On Ubuntu you may check your version of OpenGL with `glxinfo | grep "OpenGL version"`
  * [AntTweakBar] (http://anttweakbar.sourceforge.net/)

##Building from source:
  * install a C compiler
    - Ubuntu: `sudo apt-get install gcc`
    - MacOSX: [link] (https://developer.apple.com/xcode/download/)
    - Windows: [link] (http://sourceforge.net/projects/mingw-w64/)
  * install AntTweakBar library and its header files 
  * copile with 
    - Ubuntu:`g++ main.cpp -o magnoom -pthread -O3 -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL`
    - MacOSX: `g++ main.cpp -o magnoom -pthread -O3 -lAntTweakBar -framework GLUT  -framework OpenGL`
    - Windows:`not known`
    
##Contributing

Contributions are welcome!!! 

To contribute to the code development:
  - __fork__ our repository on github
  - send a __pull request__ when your corrections are done and tested.
