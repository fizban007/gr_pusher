#ifndef _PARTICLE_H_
#define _PARTICLE_H_

namespace Aperture {

template <typename Data>
struct Particle {
  typedef Particle<Data> self_type;

  Vec3<Data> x = {0.0, 0.0, 0.0};
  Vec3<Data> u = {0.0, 0.0, 0.0};
  bool is_photon = false;
  // double u0 = 0.0;
  // double e_over_m = -1.0;

  friend std::ostream& operator<<(std::ostream& os, self_type& ptc) {
    os << "x = ( " << ptc.x[0].x() << ", " << ptc.x[1].x() << ", "
       << ptc.x[2].x() << " ), u = ( " << ptc.u[0].x() << ", " << ptc.u[1].x()
       << ", " << ptc.u[2].x() << " )";
    return os;
  }

};  // ----- end of struct Particle -----

}  // namespace Aperture

#endif  // _PARTICLE_H_
