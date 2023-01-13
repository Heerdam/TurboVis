# MdVis
## Overview
The TurboVis project and the Turbodorn library is the result of my [Bachleor Thesis]() and found it's final form in the  [Semester Thesis]() as part of my Master's program. <br/>

## Building and Usage
This project contains a library, an example implementation and tests. Following a short outline on how to use any of those three targets.
### The Libary
The Turbodorn library is a header-only library that can be very easily added to your project
by leveraging the FetchContent functionality of CMake. <br/>
Simply add following Code to your CMake file. The library will fetch all needed dependencies with the exception of the [hdf5 library](https://www.hdfgroup.org/solutions/hdf5/) which needs to be installed on the system for Turbodorn to compile.
```
include(FetchContent)

FetchContent_Declare(turbodorn
    GIT_REPOSITORY      https://github.com/Heerdam/TurboVis.git
	GIT_TAG				v0.2
    GIT_SHALLOW         TRUE
)
FetchContent_MakeAvailable(turbodorn)

target_link_libraries(TARGET turbodorn)
```

### The Example Implementation
The usage of the example implementation is documented [here]().
```
git clone https://github.com/Heerdam/TurboVis.git
mkdir build
cd build
cmake ..
make tv
```

### The Tests
```
git clone https://github.com/Heerdam/TurboVis.git
mkdir build
cd build
cmake ..
make turbodorn_test
```