#ifndef _PARTICLE_H_
#define _PARTICLE_H_

namespace Aperture {

struct Particle {
  Vec3<double> x = {0.0, 0.0, 0.0};
  Vec3<double> u = {0.0, 0.0, 0.0};
  double u0 = 0.0;
  double e_over_m = -1.0;
};  // ----- end of struct Particle -----


}

#endif  // _PARTICLE_H_
