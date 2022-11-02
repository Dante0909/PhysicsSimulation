#include "ParticleSystem.h"

using namespace gti320;

/**
 * Calcule des forces qui affectent chacune des particules.
 *
 * Étant donné une particule p, la force est stockée dans le vecteur p.f
 * Les forces prisent en compte sont : la gravité et la force des ressorts.
 */
void ParticleSystem::computeForces()
{
    // TODO 
    //
    // Calcul de la force gravitationnelle sur chacune des particules

    static const Vector2f g(0, -9.81);
    for (Particle& p : m_particles)
    {
        p.f = p.m * g;
    }

    // TODO
    //
    // Pour chaque ressort, ajouter la force exercée à chacune des exptrémités sur
    // les particules appropriées. Pour un ressort s, les deux particules
    // auxquelles le ressort est attaché sont m_particles[s.index0] et
    // m_particles[s.index1]. On rappelle que les deux forces sont de même
    // magnitude mais dans des directions opposées.
    for (const Spring& s : m_springs)
    {
        Vector2f ab = m_particles[s.index1].x - m_particles[s.index0].x;
        float const l = ab.norm();
        
        Vector2f force = s.k * (1.f - (s.l0 / l)) * ab;
        
        m_particles[s.index0].f = m_particles[s.index0].f + force;
        m_particles[s.index1].f = m_particles[s.index1].f + (-1.0f * force);
    }
}

/**
 * Assemble les données du système dans les vecteurs trois vecteurs d'état _pos,
 * _vel et _force.
 */
void ParticleSystem::pack(Vector<float, Dynamic>& _pos,
    Vector<float, Dynamic>& _vel,
    Vector<float, Dynamic>& _force)
{
    // TODO 
    //
    // Copier les données des particules dans les vecteurs. Attention, si on a
    // changé de modèle, il est possible que les vecteurs n'aient pas la bonne
    // taille. Rappel : la taille de ces vecteurs doit être 2 fois le nombre de
    // particules.
    int size = m_particles.size();
    _pos.resize(size * 2);
    _vel.resize(size * 2);
    _force.resize(size * 2);

    for (int i = 0; i < size; ++i) {
        _pos(2 * i) = m_particles[i].x(0);
        _pos(2 * i + 1) = m_particles[i].x(1);

        _vel(2 * i) = m_particles[i].v(0);
        _vel(2 * i + 1) = m_particles[i].v(1);
        
        _force(2 * i) = m_particles[i].f(0);
        _force(2 * i + 1) = m_particles[i].f(1);
    }
}

/**
 * Copie les données des vecteurs d'états dans le particules du système.
 */
void ParticleSystem::unpack(const Vector<float, Dynamic>& _pos,
    const Vector<float, Dynamic>& _vel)
{
    // TODO 
    //
    // Mise à jour des propriétés de chacune des particules à partir des valeurs
    // contenues dans le vecteur d'état.
    int size = m_particles.size();
    assert(_pos.size() == 2 * size);

    for (int i = 0; i < size; ++i) {
        m_particles[i].x(0) = _pos(2 * i);
        m_particles[i].x(1) = _pos(2 * i + 1);

        m_particles[i].v(0) = _vel(2 * i);
        m_particles[i].v(1) = _vel(2 * i + 1);
    }
}



/**
 * Construction de la matirce de masses.
 */
void ParticleSystem::buildMassMatrix(Matrix<float, Dynamic, Dynamic>& M)
{
    const int numParticles = m_particles.size();
    const int dim = 2 * numParticles;
    M.resize(dim, dim);
    M.setZero();

    // TODO 
    //
    // Inscrire la masse de chacune des particules dans la matrice de masses M
    //
    for (int i = 0; i < numParticles; ++i)
    {
        float m = m_particles[i].fixed ? 1e9f : m_particles[i].m;
 
        M(2 * i, 2 * i) = m;
        M(2 * i + 1, 2 * i + 1) = m;
    }
}


/**
 * Construction de la matrice de rigidité.
 */
void ParticleSystem::buildDfDx(Matrix<float, Dynamic, Dynamic>& dfdx)
{
    const int numParticles = m_particles.size();
    const int numSprings = m_springs.size();
    const int dim = 2 * numParticles;
    dfdx.resize(dim, dim);
    dfdx.setZero();

    // Pour chaque ressort...
    for (const Spring& spring : m_springs)
    {
        auto ab = m_particles[spring.index1].x - m_particles[spring.index0].x;
        auto distance = ab.norm();
        Matrix<float, 2, 1> upperDyadique;
        upperDyadique(0, 0) = ab(0);
        upperDyadique(1, 0) = ab(1);
        
        auto produitDyadique = (1.f / (distance * distance)) * (upperDyadique * upperDyadique.transpose<float, 1, 2, ColumnStorage>());

        auto const alpha = spring.k * (1 - spring.l0 / distance);
        Matrix<float, 2, 2> alphaMatrix;
        alphaMatrix(0, 0) = alpha;
        alphaMatrix(1, 1) = alpha;

        auto const dfixj = alphaMatrix + (spring.k * (spring.l0 / distance) * produitDyadique);
     

        int const df1 = spring.index0 * 2;
        int const df2 = df1 + 1;
        int const df3 = spring.index1 * 2;
        int const df4 = df3 + 1;

        //dfdx.block(df1, df1, 2, 2) += dfixi;

        dfdx(df1, df1) -= dfixj(0, 0);
        dfdx(df1, df2) -= dfixj(0, 1);
        dfdx(df2, df1) -= dfixj(1, 0);
        dfdx(df2, df2) -= dfixj(1, 1);

        //dfdx.block(df3, df3, 2, 2) += dfixi;

        dfdx(df3, df3) -= dfixj(0, 0);
        dfdx(df3, df4) -= dfixj(0, 1);
        dfdx(df4, df3) -= dfixj(1, 0);
        dfdx(df4, df4) -= dfixj(1, 1);

        //dfdx.block(df1, df3, 2, 2) += dfixj;

        dfdx(df1, df3) += dfixj(0, 0);
        dfdx(df1, df4) += dfixj(0, 1);
        dfdx(df2, df3) += dfixj(1, 0);
        dfdx(df2, df4) += dfixj(1, 1);

        //dfdx.block(df3, df1, 2, 2) += dfixj;

        dfdx(df3, df1) += dfixj(0, 0);
        dfdx(df4, df1) += dfixj(0, 1);
        dfdx(df3, df2) += dfixj(1, 0);
        dfdx(df4, df2) += dfixj(1, 1);
        // TODO
        //
        // Calculer le coefficients alpha et le produit dyadique tel que décrit au cours.
        // Utiliser les indices spring.index0 et spring.index1 pour calculer les coordonnées des endroits affectés.
        //
        // Astuce: créer une matrice de taille fixe 2 par 2 puis utiliser la classe SubMatrix pour accumuler 
        // les modifications sur la diagonale (2 endroits) et pour mettre à jour les blocs non diagonale (2 endroits).


    }
}
