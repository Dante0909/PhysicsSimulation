#pragma once

/**
 * @file Solvers.hpp
 *
 * @brief Implémentation de plusieurs algorihtmes de solveurs pour un système
 *        d'équations linéaires
 *
 * Nom:
 * Code permanent :
 * Email :
 *
 */

#include "Math3D.h"

namespace gti320
{
    // Identification des solveurs
    enum eSolverType { kNone, kGaussSeidel, kColorGaussSeidel, kCholesky };

    // Paramètres de convergences pour les algorithmes itératifs
    static const float eps = 1e-7f;
    static const float tau = 1e-5f;

    /**
     * Résout Ax = b avec la méthode Gauss-Seidel
     */
    static void gaussSeidel(const Matrix<float, Dynamic, Dynamic>& A,
        const Vector<float, Dynamic>& b,
        Vector<float, Dynamic>& x, int k_max)
    {
        // TODO 
        //
        // Implémenter la méthode de Gauss-Seidel
        
        int const n = A.rows();
        int k = 0;
        Vector<float, Dynamic> x_prec;
        Vector<float, Dynamic> r;
        auto const bnorm = b.norm();
        x.resize(n);
        do {
            x_prec = x;
            for (int i = 0; i < n; ++i) {
                x(i) = b(i);
                for (int j = 0; j <= i - 1; ++j) {
                    x(i) -= A(i, j) * x(j);
                }
                for (int j = i + 1; j < n; ++j) {
                    x(i) -= A(i, j) * x(j);
                }
                x(i) /= A(i, i);
            }
            k++;
            r = A * x - b;

        } while (k < k_max && ((x-x_prec).norm() / x.norm()) > tau && ((r.norm() / bnorm) > eps));
       
    }

    /**

     * Résout Ax = b avec la méthode Gauss-Seidel (coloration de graphe)
     */
    static void gaussSeidelColor(const Matrix<float, Dynamic, Dynamic>& A, const Vector<float, Dynamic>& b, Vector<float, Dynamic>& x, const Partitions& P, const int maxIter)
    {
        // TODO 
        //
        // Implémenter la méthode de Gauss-Seidel avec coloration de graphe.
        // Les partitions avec l'index de chaque particule sont stockées dans la table des tables, P.
        gaussSeidel(A, b, x, maxIter);
        std::cout << P.size();
        

    }

    /**
     * Résout Ax = b avec la méthode de Cholesky
     */
    static void cholesky(const Matrix<float, Dynamic, Dynamic>& A,
        const Vector<float, Dynamic>& b,
        Vector<float, Dynamic>& x)
    {
        // TODO 
        //
        // Calculer la matrice L de la factorisation de Cholesky

        int const rows = A.rows();
        Matrix<float, Dynamic, Dynamic> L(rows, rows);
        
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j <= i; ++j) {
                float s = 0;

                if (i == j)
                {
                    for (int k = 0; k < j; ++k)
                        s += L(j, k) * L(j, k);
                    L(j,j) = std::sqrt(A(j,j) - s);
                }
                else 
                {
                    for (int k = 0; k < j; k++)
                        s += L(i,k) * L(j,k);
                    L(i,j) = (A(i,j) - s) / L(j,j);
                }
            }
        }


        // TODO
        //
        // Résoudre Ly = b
        
        Vector<float, Dynamic> y(rows);
        for (int i = 0; i < rows; ++i) {
            y(i) = b(i);
            for (int j = 0; j < i; ++j) {
                y(i) -= L(i, j) * y(j);
                
            }
            y(i) /= L(i, i);
        }
        // TODO
        //
        // Résoudre L^t x = y
        // 
        // Remarque : ne pas caculer la transposer de L, c'est inutilement
        // coûteux.
        x.resize(rows);
        for (int i = rows - 1; i >= 0; --i) {
            x(i) = y(i);
            for (int j = i + 1; j < rows; ++j) {
                x(i) -= L(j, i) * x(j);
            }
            x(i) /= L(i, i);
        }
    }

}
