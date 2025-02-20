#pragma once

#include <iostream>
#include <vector>

template <typename T> class Tensor4D {
  public:
    Tensor4D(){};

    Tensor4D(int size_1index, int size_2index, int size_3index, int size_4index)
        : m_size_1index(size_1index), m_size_2index(size_2index),
          m_size_3index(size_3index), m_size_4index(size_4index) {
        m_values.resize(size_1index * size_2index * size_3index * size_4index);
    }

    Tensor4D(Tensor4D const &other)
        : m_size_1index(other.m_size_1index),
          m_size_2index(other.m_size_2index),
          m_size_3index(other.m_size_3index),
          m_size_4index(other.m_size_4index), m_values(other.m_values) {}

    Tensor4D(Tensor4D &&other)
        : m_size_1index(other.m_size_1index),
          m_size_2index(other.m_size_2index),
          m_size_3index(other.m_size_3index),
          m_size_4index(other.m_size_4index),
          m_values(std::move(other.m_values)) {}

    T operator()(int i, int j, int k, int l) const {
        return m_values.at(
            m_size_4index * (k + m_size_3index * (j + m_size_2index * i)) + l);
    };

    T &operator()(int i, int j, int k, int l) {
        return m_values.at(
            m_size_4index * (k + m_size_3index * (j + m_size_2index * i)) + l);
    };

    void dim() {
        std::cout << "Dimension : " << m_size_1index << "x" << m_size_2index
                  << "x" << m_size_3index << "x" << m_size_4index << std::endl;
    }

    // Setters & Getters
    void set_size_1index(int size_1index) { m_size_1index = size_1index; }
    void set_size_2index(int size_2index) { m_size_2index = size_2index; }
    void set_size_3index(int size_3index) { m_size_3index = size_3index; }
    void set_size_4index(int size_4index) { m_size_4index = size_4index; }
    void set_values(std::vector<T> values) { m_values = values; }

    int size_1index() const { return m_size_1index; }
    int size_2index() const { return m_size_2index; }
    int size_3index() const { return m_size_3index; }
    int size_4index() const { return m_size_4index; }

  private:
    int m_size_1index;
    int m_size_2index;
    int m_size_3index;
    int m_size_4index;
    std::vector<T> m_values;
};