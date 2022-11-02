/**
 * @file main.cpp
 *
 * @brief Point d'entrée pour exécuter les tests de la librairie d'algèbre linéaire.
 *
 * Nom: Mathias Glorieux
 * Code permanent : AR53050
 * Email : mathias.glorieux.1@ens.etsmtl.ca
 *
 */

#include <gtest/gtest.h>

int main(int argc, char** argv)
{
    // Execute tous les tests unitaires.
    //
    // Les tests sont écrits dans les fichiers :
    //   tests/Tests1a.cpp
    //   tests/Tests1b.cpp
    //   tests/TestsSupplementaire1a.cpp
    //   tests/TestsSupplementaire1b.cpp
    //
    ::testing::InitGoogleTest(&argc, argv);
    const int ret = RUN_ALL_TESTS();
    return ret;
}
