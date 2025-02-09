#pragma once

typedef unsigned long size_t;

class System {
  public:
    virtual double overlap(size_t i, size_t j) const = 0;
    virtual double kinetic(size_t i, size_t j) const = 0;
    virtual double electron_nucleus(size_t i, size_t j) const = 0;
    virtual double electron_electron(size_t i, size_t j, size_t k,
                                     size_t l) const = 0;
    virtual double nuclear_energy() const = 0;

    virtual size_t size() const = 0;
};