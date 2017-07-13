#include "metrics.h"
#include "metrics/metric_base.h"

namespace Aperture {

namespace metric {

metric::metric_cartesian g_metric_cartesian;
metric::metric_spherical g_metric_spherical;
metric::metric_log_spherical g_metric_log_spherical;
metric::metric_schwarzschild g_metric_schwarzschild;
metric::metric_kerr_schild g_metric_kerr_schild;
metric::metric_boyer_lindquist g_metric_boyer_lindquist;
metric::metric_log_boyer_lindquist g_metric_log_boyer_lindquist;

}

namespace detail {

HOST_DEVICE void
general_sph_pos_to_cart(Vec3<Scalar>& pos) {
  double x1n = pos.x, x2n = pos.y, x3n = pos.z;
  // These become x, y, z
  pos.x = x1n * sin(x2n) * cos(x3n);
  pos.y = x1n * sin(x2n) * sin(x3n);
  pos.z = x1n * cos(x2n);
}

HOST_DEVICE void
general_sph_pos_from_cart(Vec3<Scalar>& pos) {
  double x1n = pos.x, x2n = pos.y, x3n = pos.z;
  // These become r, theta, phi
  pos.x = sqrt(x1n * x1n + x2n * x2n + x3n * x3n);
  pos.y = acos(x3n / pos.x);  // acos(z / r)
  // TODO: Check correctness of this statement!
  pos.z = atan(x2n / x1n) + (x1n < 0) * sgn(x2n) * CONST_PI;  // atan(y / x)
}

HOST_DEVICE void
general_sph_vec_to_cart(Vec3<Scalar>& v, const Vec3<Scalar>& pos) {
  double v1n = v.x, v2n = v.y, v3n = v.z;
  // Now these become vx, vy, vz
  double sy = sin(pos.y), cy = cos(pos.y);
  double sz = sin(pos.z), cz = cos(pos.z);
  v.x = v1n * sy * cz + v2n * cy * cz -
        v3n * sz;
  v.y = v1n * sy * sz + v2n * cy * sz +
        v3n * cz;
  v.z = v1n * cy - v2n * sy;
}

HOST_DEVICE void
general_sph_vec_from_cart(Vec3<Scalar>& v, const Vec3<Scalar>& pos) {
  double v1n = v.x, v2n = v.y, v3n = v.z;
  double sy = sin(pos.y), cy = cos(pos.y);
  double sz = sin(pos.z), cz = cos(pos.z);
  // These become vr, vtheta, vphi
  v.x = v1n * sy * cz + v2n * sy * sz +
        v3n * cy;
  v.y = v1n * cy * cz + v2n * cy * sz -
        v3n * sy;
  v.z = -v1n * sz + v2n * cz;
}

HOST_DEVICE void
log_sph_pos_to_cart(Vec3<Scalar>& pos) {
  double x1n = pos.x, x2n = pos.y, x3n = pos.z;
  // These become x, y, z
  pos.x = exp(x1n) * sin(x2n) * cos(x3n);
  pos.y = exp(x1n) * sin(x2n) * sin(x3n);
  pos.z = exp(x1n) * cos(x2n);
}

HOST_DEVICE void
log_sph_pos_from_cart(Vec3<Scalar>& pos) {
  double x1n = pos.x, x2n = pos.y, x3n = pos.z;
  // These become r, theta, phi
  pos.x = sqrt(x1n * x1n + x2n * x2n + x3n * x3n);
  pos.y = acos(x3n / pos.x);  // acos(z / r)
  // TODO: Check correctness of this statement!
  pos.z = atan(x2n / x1n) + (x1n < 0) * CONST_PI;  // atan(y / x)
  pos.x = log(pos.x);
}

#if defined(__AVX2__) && (defined(__ICC) || defined(__INTEL_COMPILER))
void
log_sph_pos_to_cart(__m256d* pos_x, __m256d* pos_y, __m256d* pos_z) {
  // __m256d x = *pos_x, y = *pos_y, z = *pos_z;
  __m256d r = _mm256_exp_pd(*pos_x);
  __m256d x2n = *pos_y, x3n = *pos_z;
  __m256d sin_theta = _mm256_sin_pd(*pos_y);
  *pos_x = _mm256_mul_pd(r, _mm256_mul_pd(sin_theta, _mm256_cos_pd(x3n)));
  *pos_y = _mm256_mul_pd(r, _mm256_mul_pd(sin_theta, _mm256_sin_pd(x3n)));
  *pos_z = _mm256_mul_pd(r, _mm256_cos_pd(x2n));
}

void
log_sph_pos_from_cart(__m256d* pos_x, __m256d* pos_y, __m256d* pos_z) {
  __m256d x = *pos_x, y = *pos_y, z = *pos_z;
  *pos_x = _mm256_sqrt_pd(_mm256_fmadd_pd(x, x, _mm256_fmadd_pd(y, y, _mm256_mul_pd(z, z))));
  *pos_y = _mm256_acos_pd(_mm256_div_pd(z, *pos_x));
  __m256d zero = _mm256_set1_pd(0.0);
  *pos_z = _mm256_add_pd(_mm256_atan_pd(_mm256_div_pd(y, x)),
                         _mm256_blendv_pd(zero, _mm256_set1_pd(CONST_PI), _mm256_cmp_pd(x, zero, _CMP_LT_OQ)));
  *pos_x = _mm256_log_pd(*pos_x);
}

void
general_sph_vec_to_cart(__m256d* v1, __m256d* v2, __m256d* v3,
                        const __m256d* pos_x, const __m256d* pos_y, const __m256d* pos_z) {
  __m256d sin_theta = _mm256_sin_pd(*pos_y);
  __m256d cos_theta = _mm256_cos_pd(*pos_y);
  __m256d sin_phi = _mm256_sin_pd(*pos_z);
  __m256d cos_phi = _mm256_cos_pd(*pos_z);
  __m256d v1n = *v1, v2n = *v2, v3n = *v3;

  *v1 = _mm256_fmadd_pd(v1n, _mm256_mul_pd(sin_theta, cos_phi),
                        _mm256_sub_pd(_mm256_mul_pd(v2n, _mm256_mul_pd(cos_theta, cos_phi)),
                                      _mm256_mul_pd(v3n, sin_phi)));
  *v2 = _mm256_fmadd_pd(v1n, _mm256_mul_pd(sin_theta, sin_phi),
                        _mm256_fmadd_pd(v3n, cos_phi, _mm256_mul_pd(v2n, _mm256_mul_pd(cos_theta, sin_phi))));
  *v3 = _mm256_sub_pd(_mm256_mul_pd(v1n, cos_theta), _mm256_mul_pd(v2n, sin_theta));
}

void
general_sph_vec_from_cart(__m256d* v1, __m256d* v2, __m256d* v3,
                          const __m256d* pos_x, const __m256d* pos_y, const __m256d* pos_z) {
  __m256d sin_theta = _mm256_sin_pd(*pos_y);
  __m256d cos_theta = _mm256_cos_pd(*pos_y);
  __m256d sin_phi = _mm256_sin_pd(*pos_z);
  __m256d cos_phi = _mm256_cos_pd(*pos_z);
  __m256d v1n = *v1, v2n = *v2, v3n = *v3;

  *v2 = _mm256_fmadd_pd(v1n, _mm256_mul_pd(cos_theta, cos_phi),
                        _mm256_sub_pd(_mm256_mul_pd(v2n, _mm256_mul_pd(cos_theta, sin_phi)),
                                      _mm256_mul_pd(v3n, sin_theta)));
  *v1 = _mm256_fmadd_pd(v1n, _mm256_mul_pd(sin_theta, cos_phi),
                        _mm256_fmadd_pd(v3n, cos_theta, _mm256_mul_pd(v2n, _mm256_mul_pd(sin_theta, sin_phi))));
  *v3 = _mm256_sub_pd(_mm256_mul_pd(v2n, cos_phi), _mm256_mul_pd(v1n, sin_phi));
}

#endif

}

MetricType parse_metric(const std::string& str) {
  if (metric::g_metric_cartesian.parse(str)) return MetricType::Cartesian;
  if (metric::g_metric_spherical.parse(str)) return MetricType::spherical;
  if (metric::g_metric_log_spherical.parse(str)) return MetricType::log_spherical;
  if (metric::g_metric_schwarzschild.parse(str)) return MetricType::Schwarzschild;
  if (metric::g_metric_kerr_schild.parse(str)) return MetricType::Kerr_Schild;
  if (metric::g_metric_boyer_lindquist.parse(str)) return MetricType::Boyer_Lindquist;
  if (metric::g_metric_log_boyer_lindquist.parse(str)) return MetricType::log_Boyer_Lindquist;
  return MetricType::Cartesian;
}



}
