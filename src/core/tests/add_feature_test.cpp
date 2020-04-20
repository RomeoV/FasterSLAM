#include "add_feature.h"  // import file to test

#include <cmath>

#include "linalg.h"
#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
  "add feature"_test = [] {
    given(
        "A particle with three known features and vehicle position "
        "{1,1,acos(4/5)}) "  // the angle is so that we can have a triangle with 3,4,5 side lengths for the landmark
        "and two features at {1,3} and {5,4} in local coordinates") = [] {
      /*  y
       *  |  F2 (1,3)
       *  |  |            F1 (5,4)
       *  |  |      5____/|
       *  |  |  ____/     | 3
       *  |  | /__________|
       *  |  R(1,1)  4           // Robot R points straight at F1
       *  ------------------------x
       */
      Particle p;
      initParticle(&p, 10);
      Vector3d pos = {1., 1., acos(4. / 5.)};
      copy(pos, 3, p.xv);
      p.Nfa = 3;

      Vector2d landmarks[2] = {
          {5, 0}, {2, M_PI / 2 - acos(4. / 5.)}};  // local robot frame

      Matrix2d R = {pow(0.1, 2), 0,
                    0, pow(1.0 * M_PI / 180, 2)};  // this is from the yglee config

      when("I add the new features to the particle") = [&]() mutable {
        add_feature(&p, landmarks, 2, R);

        then("I expect there to be feature means at (1,3) and (5,4)") = [&] {
          "F1"_test = [&] {
            expect(that % fabs(p.xf[2 * 3 + 0] - 5) < 1e-14) << "F1 - distance";
            expect(that % fabs(p.xf[2 * 3 + 1] - 4) < 1e-14) << "F1 - angle";
          };
          "F2"_test = [&] {
            expect(that % fabs(p.xf[2 * 4 + 0] - 1) < 1e-14) << "F2 - distance";
            expect(that % fabs(p.xf[2 * 4 + 1] - 3) < 1e-14) << "F2 - angle";
          };
        };

        then(
            "I expect there to be two more features and the same number of max "
            "features") = [&] {
          "Nfa"_test = [&] {
            expect(that % p.Nfa == 3 + 2) << "number of actual features";
            expect(that % p.Nf == 10) << "number of max features";
          };
        };
      };
      delParticleMembers(&p);
    };
  };
}
