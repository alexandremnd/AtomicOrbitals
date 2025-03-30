<br />
<div align="center">
    <a href="https://github.com/alexandremnd/AtomicOrbitals">
        <img src="img/illustration.png" alt="Logo" width="180" height="180">
    </a>
    <h3 align="center">Atomic orbitals</h3>
    <p align="center">
        An Hartree-Fock approach to find atomic ground state and electronic density using Gaussian type orbitals (GTO).
    </p>
</div>

# Getting started
## Prerequisites
1. A C++ compiler (this project is tested on g++ and clang++) supporting C++20 standard
2. [CMake](https://cmake.org/) (version at least 3.12, checks using ```cmake --version```)
3. [Boost](https://www.boost.org/)
4. [Python](https://www.python.org/) with necessary packages in [requirements](requirements.txt)
5. (Optional) [OpenMP](https://www.openmp.org/) : highly recommanded for large basis set.
6. (Optional) [Doxygen](https://www.doxygen.nl)
7. (Optional) [GoogleTest](https://github.com/google/googletest)

We assume that you already installed a C++ compiler supporting at least the C++20 standard such as [g++](https://gcc.gnu.org/)/[clang](https://clang.llvm.org/).

### Ubuntu (tested on Ubuntu 24.04.01 LTS)
To install dependencies on your Ubuntu-based distribution:
```bash
    sudo apt install cmake python3-dev libboost-all-dev libgtest-dev doxygen
```

### Fedora (tested on Fedora 40 and 41)
To install dependencies on your Fedora-based distribution:
```bash
    sudo dnf install cmake python3-devel boost-devel gtest-devel doxygen
```

### Arch (tested on rolling release of February)
To install dependencies on your Arch-based distribution:
```bash
    sudo pacman -Syu cmake base-devel python3 boost doxygen gtest
```

## Cloning and building the project

1. Clone the repository in the directory of your choice:
```bash
    git clone https://github.com/alexandremnd/AtomicOrbitals.git
    cd AtomicOrbitals
```

3. Download requirements for executing python scripts from (using a virtual environnement is recommended):
```bash
    pip install -r requirements.txt
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
    cmake --build build --target test
```

# Authors
- [Alexandre Menard](https://github.com/alexandremnd/)
- [Ewan Bataille](https://github.com/EwanBat)
- Jérémy Atané

# Contributing
Feel free to contribute for this project. Currently, we are still working on laying solid foundations for this project, so many design choice may change rapidly.

# References
- [Binding energies](data/experimental_binding_energies.csv) from the [NIST](https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html).
- Gaussian basis set from the [Basis Set Exchange](https://www.basissetexchange.org/) website.

# License
This project is licensed under the GNU GPLv3 - see [LICENCE](LICENSE) file for details
