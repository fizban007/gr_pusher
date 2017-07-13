#include "utils/memory.h"
#include "catch.hpp"

using namespace Aperture;

TEST_CASE("Trying nonvoid pointer") {
  void* p = aligned_malloc(100 * sizeof(int), 32u);
  int* p_int = reinterpret_cast<int*>(p);
  p_int[30] = 30;
  // EXPECT_EQ(&p, &p_int);
  aligned_free(reinterpret_cast<void*>(p_int));
  p = nullptr;
}
