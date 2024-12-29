#pragma once

#include "BasisSet/gaussian_contracted.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Eigen/Dense"
#include <stdexcept>
#include <vector>

template <typename T>
class Atom {
    public:
        /**
         * Throws an invalid_argument exception if Z < 1.
         *
         * @param Z Number of protons in the nucleus.
         * @param position Position of the nucleus.
         */
        Atom(int Z, Eigen::Vector3d position) : m_Z(Z), m_position(position) {
            if (Z < 1) {
                throw std::invalid_argument("Atom: The atomic number must be greater than 0.");
            }
        };

        /**
         * @brief Builds an atom at the origin.
         *
         * Throws an invalid_argument exception if Z < 1.
         *
         * @param Z Number of protons in the nucleus.
         */
        Atom(int Z) : Atom(Z, Eigen::Vector3d{0., 0., 0.}) {};

        /**
         * Throws an invalid_argument exception if Z < 1.
         *
         * @param Z Number of protons in the nucleus.
         */
        Atom(int Z, float x, float y, float z) : Atom(Z, Eigen::Vector3d{x, y, z}) {};

        Atom(Atom& atom) : m_Z(atom.Z()), m_position(atom.m_position) {};
        Atom(Atom&& atom) : m_Z(atom.Z()), m_position(atom.m_position), m_orbital(std::move(atom.m_orbital)) {};

        /**
         * @param orbital The orbital to add to the atom.
         */
        void add_orbital(const T& orbital) {
            m_orbital.push_back(orbital);
        }

        /**
         * @param n The principal quantum number. (0 < n)
         * @param l The azimuthal quantum number. (0 <= l < n)
         * @param m The magnetic quantum number. (-l <= m <= l)
         * @param alpha The radial exponential decay rate. (alpha > 0)
         */
        void add_slater_orbital(const int n, const int l, const int m, const double alpha) requires std::is_same_v<T, SlaterPrimitive> {
            m_orbital.push_back(SlaterPrimitive(n, l, m, alpha));
        }

        /**
         * @brief Adds a contracted gaussian orbital (s type) to the atom.
         *
         * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
         *
         * Throws an invalid_argument exception if (weight.size() != decay.size())
         *
         * @param weight Coefficient of the primitive gaussians in the linear combination.
         * @param decay Exponential decay rate of the primitive gaussians.
         */
        void add_gaussian_orbital_stype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<T, ContractedGaussian> {
            if (weight.size() != decay.size()) {
                throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
            }

            ContractedGaussian cg{};
            for (size_t i = 0; i < weight.size(); i++) {
                cg.add_primitive(weight[i], decay[i], 0, 0, 0);
            }

            m_orbital.push_back(cg);
        }

        /**
         * @brief Add a contracted gaussian orbital (p type) to the atom.
         *
         * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
         *
         * Throws an invalid_argument exception if (weight.size() != decay.size())
         * @param weight Coefficient of the primitive gaussians in the linear combination.
         * @param decay Exponential decay rate of the primitive gaussians.
         */
        void add_gaussian_orbital_ptype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<T, ContractedGaussian> {
            if (weight.size() != decay.size()) {
                throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
            }

            ContractedGaussian cg_x{};
            ContractedGaussian cg_y{};
            ContractedGaussian cg_z{};

            for (size_t i = 0; i < weight.size(); i++) {
                cg_x.add_primitive(weight[i], decay[i], 1, 0, 0);
                cg_y.add_primitive(weight[i], decay[i], 0, 1, 0);
                cg_z.add_primitive(weight[i], decay[i], 0, 0, 1);
            }

            m_orbital.push_back(cg_x);
            m_orbital.push_back(cg_y);
            m_orbital.push_back(cg_z);
        }

        /**
         * @brief Add a contracted gaussian orbital (d type) to the atom.
         *
         * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
         *
         * Throws an invalid_argument exception if (weight.size() != decay.size())
         * @param weight Coefficient of the primitive gaussians in the linear combination.
         * @param decay Exponential decay rate of the primitive gaussians.
         */
        void add_gaussian_orbital_dtype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<T, ContractedGaussian> {
            if (weight.size() != decay.size()) {
                throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
            }

            ContractedGaussian cg_xx{};
            ContractedGaussian cg_yy{};
            ContractedGaussian cg_zz{};
            ContractedGaussian cg_xy{};
            ContractedGaussian cg_xz{};
            ContractedGaussian cg_yz{};

            for (size_t i = 0; i < weight.size(); i++) {
                cg_xx.add_primitive(weight[i], decay[i], 2, 0, 0);
                cg_yy.add_primitive(weight[i], decay[i], 0, 2, 0);
                cg_zz.add_primitive(weight[i], decay[i], 0, 0, 2);
                cg_xy.add_primitive(weight[i], decay[i], 1, 1, 0);
                cg_xz.add_primitive(weight[i], decay[i], 1, 0, 1);
                cg_yz.add_primitive(weight[i], decay[i], 0, 1, 1);
            }

            m_orbital.push_back(cg_xx);
            m_orbital.push_back(cg_yy);
            m_orbital.push_back(cg_zz);
            m_orbital.push_back(cg_xy);
            m_orbital.push_back(cg_xz);
            m_orbital.push_back(cg_yz);
        }

        /**
         * @brief Add a contracted gaussian orbital (f type ie l = 3) to the atom.
         *
         * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
         *
         * Throws an invalid_argument exception if (weight.size() != decay.size())
         * @param weight Coefficient of the primitive gaussians in the linear combination.
         * @param decay Exponential decay rate of the primitive gaussians.
         */
        void add_gaussian_orbital_ftype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<T, ContractedGaussian> {
            if (weight.size() != decay.size()) {
                throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
            }

            ContractedGaussian cg_xxx{};
            ContractedGaussian cg_yyy{};
            ContractedGaussian cg_zzz{};
            ContractedGaussian cg_xxy{};
            ContractedGaussian cg_xxz{};
            ContractedGaussian cg_xyy{};
            ContractedGaussian cg_yyz{};
            ContractedGaussian cg_xzz{};
            ContractedGaussian cg_yzz{};
            ContractedGaussian cg_xyz{};

            for (size_t i = 0; i < weight.size(); i++) {
                cg_xxx.add_primitive(weight[i], decay[i], 3, 0, 0);
                cg_yyy.add_primitive(weight[i], decay[i], 0, 3, 0);
                cg_zzz.add_primitive(weight[i], decay[i], 0, 0, 3);
                cg_xxy.add_primitive(weight[i], decay[i], 2, 1, 0);
                cg_xxz.add_primitive(weight[i], decay[i], 2, 0, 1);
                cg_xyy.add_primitive(weight[i], decay[i], 1, 2, 0);
                cg_yyz.add_primitive(weight[i], decay[i], 0, 2, 1);
                cg_xzz.add_primitive(weight[i], decay[i], 1, 0, 2);
                cg_yzz.add_primitive(weight[i], decay[i], 0, 1, 2);
                cg_xyz.add_primitive(weight[i], decay[i], 1, 1, 1);
            }

            m_orbital.push_back(cg_xxx);
            m_orbital.push_back(cg_yyy);
            m_orbital.push_back(cg_zzz);
            m_orbital.push_back(cg_xxy);
            m_orbital.push_back(cg_xxz);
            m_orbital.push_back(cg_xyy);
            m_orbital.push_back(cg_yyz);
            m_orbital.push_back(cg_xzz);
            m_orbital.push_back(cg_yzz);
            m_orbital.push_back(cg_xyz);
        }

        inline int Z() const { return m_Z; }
        inline Eigen::Vector3d position() const { return m_position; }

        inline const int n_orbitals() const { return m_orbital.size(); }
        inline const std::vector<T>& get_orbitals() const { return m_orbital; }
        inline const T& get_orbital(size_t i) const { return m_orbital[i]; }
        inline const T& operator()(size_t i) const { return m_orbital[i]; }

    private:
        int m_Z;
        Eigen::Vector3d m_position;
        std::vector<T> m_orbital;
};