![Alt text](https://github.com/n-s-kiselev/magnoom/blob/master/MagnoomWiki/TitleImage.png)

# Magnoom

Advanced software for Atomistic Spin Dynamics (see also [Wikipage](https://github.com/n-s-kiselev/magnoom/wiki) )

### Developer:
  [Nikolai S. Kiselev](http://www.fz-juelich.de/SharedDocs/Personen/PGI/PGI-1/EN/Kiselev_N.html?nn=758466) 

### Contributors:
  1. [Filipp N. Rybakov](http://www.hopfion.com) (Implementation for Windows OS and alternative threads semaphores)
  2. [Andriy S. Savchenko](https://www.researchgate.net/profile/A_Savchenko) (Implementation of a variety of initial states)

### Main features:
  * platform independent code writen in C/C++
  * efficiently CPU parallelized, multithreading code
  * minimal dependencies on external libraries
  * real time control of the parameters
  * GUI entirely based on OpenGL 
  * advanced visualisation and post processing tools:
    - slicing 
    - filters
  * reading vector field from [OVF files](http://math.nist.gov/oommf/doc/userguide11b2/userguide/vectorfieldformat.html) ([OOMMF](http://math.nist.gov/oommf/) and [MuMax3](http://mumax.github.io/))  
  


### Dependencies:
  * gcc >=4.8.1 (possibly can be compiled also with earlier versions of gcc)
  * OpenGL >= 2.1 (On Ubuntu and MacOS check your version of OpenGL with `glxinfo | grep "OpenGL version"`
  * [AntTweakBar](http://anttweakbar.sourceforge.net/) by Philippe Decaudin wich is cool C/C++ library allowing graphical user interface in OpenGL applications 

### Building from source:
  * install a C compiler
    - Ubuntu: `sudo apt-get install gcc`
    - MacOSX: [link] (https://developer.apple.com/xcode/download/)
    - Windows: [Microsoft Visual Studio Community 12.0] (https://www.visualstudio.com/vs/community/)
  * install [AntTweakBar](http://anttweakbar.sourceforge.net/) library and its header files such that gcc (or Windows VS) can find them, for detail see [likn](http://anttweakbar.sourceforge.net/doc/tools:anttweakbar:download#contact)
  * compile with 
    - Ubuntu:`g++ main.cpp -o magnoom -pthread -O3 -lAntTweakBar -lpthread  -lglut -lGLU -lGLEW -lGL`
    - MacOSX: `g++ main.cpp -o magnoom -pthread -O3 -lAntTweakBar -framework GLUT  -framework OpenGL -Wno-deprecated-declarations`
    - Windows: via System Preferences add the path to the `vcvars64.bat`, for example `"C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\amd64\vcvars64.bat"`, then compile with: `cl main.cpp -o magnoom /O3`
    
### Contributing:

Contributions are welcome!!! 

To contribute to the code development:
  - [__fork__](https://help.github.com/articles/fork-a-repo/) our repository on github
  - send a [__pull request__](https://help.github.com/articles/about-pull-requests/) when your corrections are done and tested.
