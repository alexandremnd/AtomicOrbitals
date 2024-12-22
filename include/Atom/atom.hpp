#pragma once

#include "BasisSet/gaussian_contracted.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Eigen/Dense"
#include <stdexcept>
#include <vector>

template <typename T>
class Atom {
    public:
        Atom(int Z) : m_Z(Z), m_position(0., 0., 0.) {};
        Atom(int Z, Eigen::Vector3f position) : m_Z(Z), m_position(position) {};
        Atom(int Z, float x, float y, float z) : m_Z(Z), m_position(x, y, z) {};
        Atom(Atom& atom) : m_Z(atom.Z()), m_position(atom.m_position) {};

        /**
         * @brief Add a basis orbital of any type to the atom.
         *
         * @param orbital The orbital to add to the atom.
         */
        void add_orbital(const T& orbital) {
            m_orbital.push_back(orbital);
        }

        /**
         * @brief Add a slater orbital to the atom.
         *
         * @param n The principal quantum number.
         * @param l The azimuthal quantum number.
         * @param m The magnetic quantum number.
         * @param alpha The radial exponential decay rate.
         */
        void add_slater_orbital(const int n, const int l, const int m, const double alpha) requires std::is_same_v<T, SlaterPrimitive> {
            m_orbital.push_back(SlaterPrimitive(n, l, m, alpha));
        }

        /**
         * @brief Add a contracted gaussian orbital (s type) to the atom.
         * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
         * @param weight Coefficient of the primitive gaussians in the linear combination.
         * @param decay Exponential decay rate of the primitive gaussians.
         */
        void add_gaussian_orbital_stype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<T, ContractedGaussian> {
            if (weight.size() != decay.size()) {
                throw std::invalid_argument("[Atom] The weight and decay vectors must have the same size.");
            }

            ContractedGaussian cg{};
            for (size_t i = 0; i < weight.size(); i++) {
                cg.add_primitive(weight[i], decay[i], 0, 0, 0);
            }

            m_orbital.push_back(cg);
        }

        /**
         * @brief Add a contracted gaussian orbital (p type) to the atom.
         * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
         * @param weight Coefficient of the primitive gaussians in the linear combination.
         * @param decay Exponential decay rate of the primitive gaussians.
         */
        void add_gaussian_orbital_ptype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<T, ContractedGaussian> {
            if (weight.size() != decay.size()) {
                throw std::invalid_argument("[Atom] The weight and decay vectors must have the same size.");
            }

            ContractedGaussian cg_x{};
            ContractedGaussian cg_y{};
            ContractedGaussian cg_z{};

            for (size_t i = 0; i < weight.size(); i++) {
                cg_x.add_primitive(weight[i], decay[i], 1, 0, 0);
                cg_y.add_primitive(weight[i], decay[i], 0, 1, 0);
                cg_z.add_primitive(weight[i], decay[i], 0, 0, 1);
            }

            m_orbital.insert(m_orbital.end(), {cg_x, cg_y, cg_z});
        }

        /**
         * @brief Add a contracted gaussian orbital (d type) to the atom.
         * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
         * @param weight Coefficient of the primitive gaussians in the linear combination.
         * @param decay Exponential decay rate of the primitive gaussians.
         */
        void add_gaussian_orbital_dtype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<T, ContractedGaussian> {
            if (weight.size() != decay.size()) {
                throw std::invalid_argument("[Atom] The weight and decay vectors must have the same size.");
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

            m_orbital.insert(m_orbital.end(), {cg_xx, cg_yy, cg_zz, cg_xy, cg_xz, cg_yz});
        }

        /**
         * @brief Add a contracted gaussian orbital (f type) to the atom.
         * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
         * @param weight Coefficient of the primitive gaussians in the linear combination.
         * @param decay Exponential decay rate of the primitive gaussians.
         */
        void add_gaussian_orbital_ftype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<T, ContractedGaussian> {
            if (weight.size() != decay.size()) {
                throw std::invalid_argument("[Atom] The weight and decay vectors must have the same size.");
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

            m_orbital.insert(m_orbital.end(), {cg_xxx, cg_yyy, cg_zzz, cg_xxy, cg_xxz, cg_xyy, cg_yyz, cg_xzz, cg_yzz, cg_xyz});
        }

        inline int Z() const { return m_Z; }
        inline Eigen::Vector3f position() const { return m_position; }

    private:
        int m_Z;
        Eigen::Vector3f m_position;
        std::vector<T> m_orbital;
};