Separator Library
========================

## Background

The backened solver [GLPK](https://www.gnu.org/software/glpk/) can be installed by following the instructions located [here](https://en.wikibooks.org/wiki/GLPK/Linux_OS#Install).

GLPK Reference Manual (and API documentation) is available [here](http://www.chiark.greenend.org.uk/doc/glpk-doc/glpk.pdf).

## Sample Usage (ROS / catkin)

In CMakeLists.txt add the library name to `find_package()`.

```cmake
find_package(catkin REQUIRED COMPONENTS separator)
```

Example: see test.cpp

## Credits
Part of the code is based on the ACL motoralloc library done by Brett Lopez