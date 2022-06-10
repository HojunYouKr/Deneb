#include "deneb.h"

#include "avocado.h"

namespace deneb {
void DENEB_Initialize() {
  MASTER_MESSAGE(avocado::GetTitle("Deneb Initialization"));
  std::string message = DENEB_VERSION;
  MASTER_MESSAGE(message + "\n");

  auto& config = AVOCADO_CONFIG;
  std::string resource_dir = config->GetConfigValue(RESOURCE_DIR);
  avocado::ToAbsolutePath(resource_dir);
  DENEB_ELEMENT_INITIALIZE(resource_dir);
  DENEB_POINT_INITIALIZE(resource_dir);
  DENEB_SAVELOAD_INITIALIZE();

  DENEB_DATA_INITIALIZE();
  DENEB_EQUATION_INITIALIZE(config->GetConfigValue(EQUATION));
  DENEB_TIMESCHEME_INITIALIZE(config->GetConfigValue(TIMESCHEME));

  const int cell_post_order =
      std::stoi(config->GetConfigValue(CELL_POST_ORDER));
  const int face_post_order =
      std::stoi(config->GetConfigValue(FACE_POST_ORDER));
  DENEB_CONTOUR_INITIALIZE(cell_post_order, face_post_order);
  DENEB_CONTOUR->SetStrandID(0);

  DENEB_LIMITER_INITIALIZE(config->GetConfigValue(LIMITER)); 
  DENEB_PRESSUREFIX_INITIALIZE(config->GetConfigValue(PRESSUREFIX)); 
  DENEB_ARTIFICIAL_VISCOSITY_INITIALIZE(
      config->GetConfigValue(ARTIFICIAL_VISCOSITY));
}
void DENEB_Finalize() {
  DENEB_ARTIFICIAL_VISCOSITY_FINALIZE(); 
  DENEB_PRESSUREFIX_FINALIZE();          
  DENEB_LIMITER_FINALIZE();              
  DENEB_CONTOUR_FINALIZE();
  DENEB_TIMESCHEME_FINALIZE();
  DENEB_EQUATION_FINALIZE();
  DENEB_DATA_FINALIZE();

  DENEB_SAVELOAD_FINALIZE();
  DENEB_POINT_FINALIZE();
  DENEB_ELEMENT_FINALIZE();
}
}  // namespace deneb