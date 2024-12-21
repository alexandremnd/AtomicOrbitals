#include <stdexcept>
#include <vector>
#include <iostream>

template <typename T>
class Yoshimine {
public:
    /**
     * @brief Construct a new Yoshimine container using Yoshimine sort
     *
     * @param size The size of the Yoshimine container
     */
    Yoshimine(int size) {
        m_yoshimine.resize(size);
    }

    inline T operator()(int a, int b, int c, int d) const {
        return this->operator()(a, b, c, d);
    }

    inline T& operator()(int a, int b, int c, int d) {
        int abcd = index(a, b, c, d);

        if (abcd >= m_yoshimine.size()) {
            throw std::out_of_range("Yoshimine index out of range");
        }

        return m_yoshimine.at(abcd);
    }

    void print_content() {
        for (int i = 0; i < m_yoshimine.size(); i++) {
            std::cout << i << ": " << m_yoshimine.at(i) << std::endl;
        }
    }

    int size() const { return m_yoshimine.size(); }

private:
    inline int index(int a, int b, int c, int d) const {
        int ab, cd, abcd;

        if (a > b) {
            ab = (a*(a+1))/2 + b;
        } else {
            ab = (b*(b+1))/2 + a;
        }

        if (c > d) {
            cd = (c*(c+1))/2 + d;
        } else {
            cd = (d*(d+1))/2 + c;
        }

        if (ab > cd) {
            abcd = (ab*(ab+1))/2 + cd;
        } else {
            abcd = (cd*(cd+1))/2 + ab;
        }

        return abcd;
    }

    std::vector<T> m_yoshimine;
};