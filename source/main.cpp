
#include "avocado.h"
#include "deneb.h"

#include <petsc.h>

int main(int argc, char* argv[]) {
  PetscInitialize(&argc, &argv, (char*)0, (char*)0);

  const std::string codename =
      DENEB_VERSION + " (Last update: " + DENEB_LAST_UPDATE + ")";
  std::vector<std::string> add_argv;
  avocado::AVOCADO_Initialize(argc, argv, codename.c_str(), add_argv);
  deneb::DENEB_Initialize();

  // build data (BE AWARE the dependency)
  DENEB_DATA->BuildData();
  DENEB_EQUATION->BuildData();
  DENEB_TIMESCHEME->BuildData();
  DENEB_CONTOUR->BuildData();
  DENEB_LIMITER->BuildData();             
  DENEB_PRESSUREFIX->BuildData();         
  DENEB_ARTIFICIAL_VISCOSITY->BuildData();
  AVOCADO_MEMORY_PROFILER->PrintMemoryUsage();

  // time marching
  DENEB_TIMESCHEME->Marching();

  SYNCRO();
  MASTER_MESSAGE("Program normal exit.\n");
  deneb::DENEB_Finalize();
  avocado::AVOCADO_Finalize();
  PetscFinalize();
  return 0;
}
