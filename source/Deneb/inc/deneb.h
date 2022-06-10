#pragma once

#include "deneb_basis.h"
#include "deneb_config_macro.h"
#include "deneb_contour.h"
#include "deneb_data.h"
#include "deneb_element.h"
#include "deneb_equation.h"
#include "deneb_header.h"
#include "deneb_limiter.h" 
#include "deneb_pressurefix.h"
#include "deneb_artificial_viscosity.h" 
#include "deneb_point.h"
#include "deneb_saveload.h"
#include "deneb_timescheme.h"

namespace deneb {
void DENEB_Initialize();
void DENEB_Finalize();
}  // namespace deneb