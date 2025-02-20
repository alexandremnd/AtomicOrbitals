#pragma once

#include <iostream>
#include <vector>

template <typename T> class Tensor3D {
  public:
    Tensor3D(){};

    Tensor3D(int size_1index, int size_2index, int size_3index)
        : m_size_1index(size_1index), m_size_2index(size_2index),
          m_size_3index(size_3index) {
        m_values.resize(size_1index * size_2index * size_3index);
    }

    Tensor3D(Tensor3D const &other)
        : m_size_1index(other.m_size_1index),
          m_size_2index(other.m_size_2index),
          m_size_3index(other.m_size_3index), m_values(other.m_values) {}

    Tensor3D(Tensor3D &&other)
        : m_size_1index(other.m_size_1index),
          m_size_2index(other.m_size_2index),
          m_size_3index(other.m_size_3index),
          m_values(std::move(other.m_values)) {}

    T operator()(int i, int j, int k) const {
        return m_values.at(m_size_3index * (j + m_size_2index * i) + k);
    };

    T &operator()(int i, int j, int k) {
        return m_values.at(m_size_3index * (j + m_size_2index * i) + k);
    };

    void dim() {
        std::cout << "Dimension : " << m_size_1index << "x" << m_size_2index
                  << "x" << m_size_3index << std::endl;
    }

    // Setters & Getters
    void set_size_1index(int size_1index) { m_size_1index = size_1index; }
    void set_size_2index(int size_2index) { m_size_2index = size_2index; }
    void set_size_3index(int size_3index) { m_size_3index = size_3index; }
    void set_values(std::vector<T> values) { m_values = values; }

    int size_1index() const { return m_size_1index; }
    int size_2index() const { return m_size_2index; }
    int size_3index() const { return m_size_3index; }

  private:
    int m_size_1index;
    int m_size_2index;
    int m_size_3index;
    std::vector<T> m_values;
};