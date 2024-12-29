<br />
<div align="center">
    <a href="https://github.com/alexandremnd/AtomicOrbitals">
        <img src="img/illustration.png" alt="Logo" width="180" height="180">
    </a>
    <h3 align="center">Atomic orbitals</h3>
    <p align="center">
        An Hartree-Fock approach to find atomic ground state and electronic density using Slater type orbitals (STO) and Gaussian type orbitals (GTO).
    </p>
</div>

# Getting started
## Prerequisites
1. A C++ compiler (this project is tested on g++ and clang++) supporting C++20 standard
2. [CMake](https://cmake.org/) (version at least 3.12, checks using ```cmake --version```)
3. [Boost](https://www.boost.org/)
4. (Optional) [Doxygen](https://www.doxygen.nl)
5. (Optional) [GoogleTest](https://github.com/google/googletest)

You may install Doxygen and GoogleTest with your package manager (apt/dnf/yum/pacman/... for Linux, chocolatey for Windows, Homebrew for MacOS).

### Ubuntu
To install dependencies on your Ubuntu-bases distribution:
```bash
    sudo apt install cmake libgtest-dev doxygen
```

### Fedora
To install dependencies on your Fedora-based distribution:
```bash
    sudo dnf install cmake gtest-devel doxygen g++
```

## Cloning and building the project

1. Clone the repository in the directory of your choice:
```bash
    git clone https://github.com/alexandremnd/AtomicOrbitals.git
    cd AtomicOrbitals
```
2. Build the project without tests:
```bash
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DENABLE_DOCS=OFF -DENABLE_TESTS=OFF
    cmake --build build --config release
    ./build/AtomicOrbitals
```

*Note:* You can switch ```ENABLE_DOCS``` and ```ENABLE_TESTS``` if you want to build the documentation and test cases.
When setting ```ENABLE_TESTS```, ```CMAKE_BUILD_TYPE``` is set to Debug for every future build. If you want to rebuild in Release mode, do not forget to add ```-DCMAKE_BUILD_TYPE=Release``` and ```-DENABLE_TESTS=OFF``` to reset cache file to next build !

1. Compiled project is located in ```build/AtomicOrbitals``` and (optionnaly) generated documentation in ```build/docs/html/index.html```

## Generating and running test cases
To generate test cases, [GoogleTest](https://github.com/google/googletest) is required.

1. Build the project with ```-DENABLE_TESTS=ON``` (build type will be set to Debug automatically):
```bash
    cmake -B build -DENABLE_TESTS=ON
    cmake --build build --config debug
```

2. Execute test with:
```bash
    cmake --build build --target test --output-on-failure
```

## Generate documentation

This project uses Doxygen to generate documentation. Make sure you have [Doxygen](https://www.doxygen.nl) installed.

1. ```cd``` into the root folder.
2. Run doxygen to create the documentation files:
   ```bash
   doxygen
   ```
3. After the build, the generated documentation will be available in ```build/docs/html/index.html```