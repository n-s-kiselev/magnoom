#Magnoom

Advanced software for Atomistic Spin Dynamics allowing in real-time visualization and control

###Developers:
  1. [Nikolai S. Kiselev](http://www.fz-juelich.de/SharedDocs/Personen/PGI/PGI-1/EN/Kiselev_N.html?nn=758466)
  2. [Jonathan Chico](http://www.fz-juelich.de/SharedDocs/Personen/PGI/PGI-1/EN/Chico_J.html?nn=758466)

###Main features:
  * platform independent code writen in C/C++
  * minimal dependencies on external libraries
  * GUI entirely based on OpenGL 
  * real time control of the parameters

###Dependencies:
  * gcc >=4.8.1 (possibly can be compiled also with earlier versions of gcc)
  * OpenGL >= 2.1 (On Ubuntu and MacOS check your version of OpenGL with `glxinfo | grep "OpenGL version"`
  * [AntTweakBar](http://anttweakbar.sourceforge.net/)

###Building from source:
  * install a C compiler
    - Ubuntu: `sudo apt-get install gcc`
    - MacOSX: [link] (https://developer.apple.com/xcode/download/)
    - Windows: [link] (http://sourceforge.net/projects/mingw-w64/)
  * install [AntTweakBar](http://anttweakbar.sourceforge.net/) library and its header files such that gcc can find them
  * compile with 
    - Ubuntu:`g++ main.cpp -o magnoom -pthread -O3 -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL`
    - MacOSX: `g++ main.cpp -o magnoom -pthread -O3 -lAntTweakBar -framework GLUT  -framework OpenGL -Wno-deprecated-declarations`
    - Windows:`not known`
    
###Contributing:

Contributions are welcome!!! 

To contribute to the code development:
  - [__fork__](https://help.github.com/articles/fork-a-repo/) our repository on github
  - send a [__pull request__](https://help.github.com/articles/about-pull-requests/) when your corrections are done and tested.
