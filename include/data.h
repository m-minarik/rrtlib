#ifndef __DATA__
#define __DATA__

#include <string>

#include "RAPID.H"
#include "types.h"

#include "matrix.h" // from libcip

bool off_to_rapid(const std::string &off, RAPID_model *model, BoundingBox &bb,
                  const Matrix *R_init = nullptr, const Matrix *t_init = nullptr);

#endif // __DATA__
