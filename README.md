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
2. CMake (version at least 3.12)
3. Boost
4. (Optional) [Doxygen](https://www.doxygen.nl)
5. (Optional) [GoogleTest](https://github.com/google/googletest)

You may install Doxygen and GoogleTest with your package manager (apt/dnf/yum/pacman/... for Linux, chocolatey for Windows, Homebrew for MacOS).

### Ubuntu
To install dependencies in your Ubuntu-bases distribution:
```bash
    sudo apt install cmake libgtest-dev doxygen
```

### Fedora
To install dependencies in your Fedora-based distribution:
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
    mkdir build && cd build
    cmake ..
    make
```

3. Compiled project is in ```build/AtomicOrbitals```.

## Generate test cases
To generate test cases, Google Test is required.
1. Build the project with ```-DENABLE_TESTS=ON```:
```bash
    mkdir build
    cd build
    cmake .. -DENABLE_TESTS=ON
    make
```

2. Execute test with either:
```bash
    ctest
```
or
```bash
    make test
```

## Generate documentation

This project uses Doxygen to generate documentation. Make sure you have [Doxygen](https://www.doxygen.nl) installed.

1. ```cd``` into the root folder.
2. Run doxygen to create the documentation files:
   ```bash
   doxygen
   ```
3. After the build, the generated documentation will be available in ```build/docs/html/index.html```