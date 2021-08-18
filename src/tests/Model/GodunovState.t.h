#include <cxxtest/TestSuite.h>
#include <cassert>
#include <cmath>

#include <Equations/datastructures.hpp>
#include <Model/common.hpp>
#include <Equations/Setup.h>

#include "values.h"

namespace seissol {
  namespace unit_test {
    class GodunovStateTestSuite;
  }
}

class seissol::unit_test::GodunovStateTestSuite : public CxxTest::TestSuite
{
  public:
    const real epsilon = std::numeric_limits<real>::epsilon();

    void testGodunovState()
    {
      real localData[seissol::tensor::QgodLocal::size()];
      real neighborData[seissol::tensor::QgodLocal::size()];
      init::QgodLocal::view::type QgodLocal = init::QgodLocal::view::create(localData);
      init::QgodNeighbor::view::type QgodNeighbor = init::QgodNeighbor::view::create(neighborData); 
      QgodLocal.setZero();
      QgodNeighbor.setZero();

      //test homogeneous material
#ifdef USE_ANISOTROPIC
      seissol::model::AnisotropicMaterial local(materialVal_1, 22);
      seissol::model::AnisotropicMaterial neighbor(materialVal_1, 22);
#elif defined USE_POROELASTIC
      seissol::model::PoroElasticMaterial local(materialVal_1, 10);
      seissol::model::PoroElasticMaterial neighbor(materialVal_1, 10);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
      seissol::model::ViscoElasticMaterial local(materialVal_1, 3 + NUMBER_OF_RELAXATION_MECHANISMS*4);
      seissol::model::ViscoElasticMaterial neighbor(materialVal_1, 3 + NUMBER_OF_RELAXATION_MECHANISMS*4);
#else
      seissol::model::ElasticMaterial local(materialVal_1, 3);
      seissol::model::ElasticMaterial neighbor(materialVal_1, 3);
#endif

      seissol::model::getTransposedGodunovState( local,
          neighbor,
          FaceType::regular,
          QgodLocal,
          QgodNeighbor );
      test_matrix(QgodLocal, solution_homogeneous_local);
      test_matrix(QgodNeighbor, solution_homogeneous_neighbor);

      //test free surface
      seissol::model::getTransposedGodunovState( local,
          neighbor,
          FaceType::freeSurface,
          QgodLocal,
          QgodNeighbor );
      test_matrix(QgodLocal, solution_boundary);
      test_nan(QgodNeighbor);



      //test heterogeneous material
#ifdef USE_ANISOTROPIC
      neighbor = seissol::model::AnisotropicMaterial(materialVal_2, 22);
#elif defined USE_POROELASTIC
      neighbor = seissol::model::PoroElasticMaterial(materialVal_2, 10);
#elif defined USE_VISCOELASTIC || defined USE_VISCOELASTIC2
      neighbor = seissol::model::ViscoElasticMaterial(materialVal_2, 3 + NUMBER_OF_RELAXATION_MECHANISMS * 4);
#else
      neighbor = seissol::model::ElasticMaterial(materialVal_2, 3);
#endif

      seissol::model::getTransposedGodunovState( local,
          neighbor,
          FaceType::regular,
          QgodLocal,
          QgodNeighbor );
      test_matrix(QgodLocal, solution_heterogeneous_local);
      test_matrix(QgodNeighbor, solution_heterogeneous_neighbor);

    }

    template<typename T>
    void test_matrix(init::QgodLocal::view::type QgodLocal, T solution) {
      for (unsigned int i = 0; i < solution[0].size(); i++) {
        for (unsigned int j = 0; j < solution.size(); j++) {
          test_relative(QgodLocal(i,j), solution[j][i]);
        }
      }
    }

    template<typename T>
    void test_local(init::QgodLocal::view::type QgodLocal, T solution) {
      for (unsigned int i = 0; i < solution[0].size(); i++) {
        for (unsigned int j = 0; j < solution.size(); j++) {
          if (i == j) {
            test_relative(QgodLocal(i,j), 1-solution[j][i]);
          }
          else {
            test_relative(QgodLocal(i,j), -solution[j][i]);
          }
        }
      }
    }

    template<typename T>
    void test_neighbor(init::QgodLocal::view::type QgodNeighbor, T solution) {
      for (unsigned int i = 0; i < solution[0].size(); i++) {
        for (unsigned int j = 0; j < solution.size(); j++) {
          test_relative(QgodNeighbor(i,j), solution[j][i]);
        }
      }
    }

    void test_nan(init::QgodLocal::view::type QgodNeighbor) {
      for (unsigned int i = 0; i < 9; i++) {
        for (unsigned int j = 0; j < 9; j++) {
          //CXXTEST TS_ASSERT_IS_NAN
          TS_ASSERT(std::isnan(QgodNeighbor(i,j)));
        }
      }
    }

    void test_relative(real a, real b) {
      if (std::abs(a) < epsilon) {
        TS_ASSERT_DELTA(0, b, epsilon);
      } else if (std::abs(b) < epsilon) {
        TS_ASSERT_DELTA(0, a, epsilon);
      } else {
        TS_ASSERT_DELTA(a, b, std::abs(a)*1e2*epsilon);
        TS_ASSERT_DELTA(a, b, std::abs(b)*1e2*epsilon);
      }
    }

};
