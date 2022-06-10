#include "deneb_pressurefix.h"

#include "avocado.h"
#include "deneb_data.h"
#include "deneb_equation.h"

namespace deneb {
std::shared_ptr<PressureFix> DENEB_PRESSUREFIX_NAME = nullptr;
std::shared_ptr<PressureFix> PressureFix::GetPressureFix(
    const std::string& name) {
  if (!name.compare("Off"))
    return std::make_shared<PressureFixOff>();
  else if (!name.compare("On"))
    return std::make_shared<PressureFixOn>();
  ERROR_MESSAGE("Wrong pressure fix option (On/Off):" + name + "\n");
  return nullptr;
}

PressureFixOff::PressureFixOff() {
  MASTER_MESSAGE(avocado::GetTitle("PressureFixOff"));
}
void PressureFixOff::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("PressureFixOff::BuildData()"));
}

PressureFixOn::PressureFixOn() : iteration_(10) {
  MASTER_MESSAGE(avocado::GetTitle("PressureFixOn"));
}
void PressureFixOn::BuildData(void) {
  MASTER_MESSAGE(avocado::GetTitle("PressureFixOn::BuildData()"));

  int max_num_points = 0;
  int max_num_cell_points = 0;
  int max_num_face_points = 0;
  int max_num_bdry_points = 0;
  const std::vector<int>& num_cell_points = DENEB_DATA->GetNumCellPoints();
  const std::vector<int>& num_face_points = DENEB_DATA->GetNumFacePoints();
  const std::vector<int>& num_bdry_points = DENEB_DATA->GetNumBdryPoints();
  if (num_cell_points.size() != 0)
    max_num_cell_points =
        *std::max_element(num_cell_points.begin(), num_cell_points.end());
  if (num_face_points.size() != 0)
    max_num_face_points =
        *std::max_element(num_face_points.begin(), num_face_points.end());
  if (num_bdry_points.size() != 0)
    max_num_bdry_points =
        *std::max_element(num_bdry_points.begin(), num_bdry_points.end());
  max_num_face_points = std::max(max_num_bdry_points, max_num_face_points);
  max_num_points = std::max(max_num_cell_points, max_num_face_points);

  const int& num_states = DENEB_EQUATION->GetNumStates();
  solution_.resize(num_states * max_num_points);
}
void PressureFixOn::Execute(double* solution) {
  static const int& num_states = DENEB_EQUATION->GetNumStates();
  static const int& num_bases = DENEB_DATA->GetNumBases();
  static const int sb = num_states * num_bases;

  // inner face sweep
  static const int& num_inner_faces = DENEB_DATA->GetNumInnerFaces();
  static const auto& num_face_points = DENEB_DATA->GetNumFacePoints();
  static const auto& face_owner_cell = DENEB_DATA->GetFaceOwnerCell();
  static const auto& face_neighbor_cell = DENEB_DATA->GetFaceNeighborCell();
  static const auto& face_owner_basis_value =
      DENEB_DATA->GetFaceOwnerBasisValue();
  static const auto& face_neighbor_basis_value =
      DENEB_DATA->GetFaceNeighborBasisValue();
  for (int iface = 0; iface < num_inner_faces; iface++) {
    const int& num_points = num_face_points[iface];
    const int& owner_cell = face_owner_cell[iface];
    const int& neighbor_cell = face_neighbor_cell[iface];

    int i = 0;
    for (; i < iteration_; i++) {
      avocado::Kernel0::f4(&solution[owner_cell * sb],
                           &face_owner_basis_value[iface][0], &solution_[0],
                           num_states, num_bases, num_points, 1.0, 0.0);

      bool flag = false;
      for (int ipoint = 0; ipoint < num_points; ipoint++) {
        const std::vector<double>& pressure_fix_values =
            DENEB_EQUATION->ComputePressureFixValues(
                &solution_[ipoint * num_states]);
        for (int ivalue = 0; ivalue < pressure_fix_values.size(); ivalue++) {
          if (pressure_fix_values[ivalue] <= 0.0) {
            flag = true;
            break;
          }
        }
        if (flag) {
          break;
        }
      }
      if (flag) {
        for (int istate = 0; istate < num_states; istate++)
          cblas_dscal(num_bases - 1, 0.5,
                      &solution[owner_cell * sb + istate * num_bases + 1], 1);
      } else
        break;
    }
    if (i == iteration_) {
      for (int istate = 0; istate < num_states; istate++)
        memset(&solution[owner_cell * sb + istate * num_bases + 1], 0,
               (num_bases - 1) * sizeof(double));
    }

    i = 0;
    for (; i < iteration_; i++) {
      avocado::Kernel0::f4(&solution[neighbor_cell * sb],
                           &face_neighbor_basis_value[iface][0], &solution_[0],
                           num_states, num_bases, num_points, 1.0, 0.0);

      bool flag = false;
      for (int ipoint = 0; ipoint < num_points; ipoint++) {
        const std::vector<double>& pressure_fix_values =
            DENEB_EQUATION->ComputePressureFixValues(
                &solution_[ipoint * num_states]);
        for (int ivalue = 0; ivalue < pressure_fix_values.size(); ivalue++) {
          if (pressure_fix_values[ivalue] <= 0.0) {
            flag = true;
            break;
          }
        }
        if (flag) {
          break;
        }
      }
      if (flag) {
        for (int istate = 0; istate < num_states; istate++)
          cblas_dscal(num_bases - 1, 0.5,
                      &solution[neighbor_cell * sb + istate * num_bases + 1],
                      1);
      } else
        break;
    }
    if (i == iteration_) {
      for (int istate = 0; istate < num_states; istate++)
        memset(&solution[neighbor_cell * sb + istate * num_bases + 1], 0,
               (num_bases - 1) * sizeof(double));
    }
  }

  // bdry sweep
  static const int& num_bdries = DENEB_DATA->GetNumBdries();
  static const auto& num_bdry_points = DENEB_DATA->GetNumBdryPoints();
  static const auto& bdry_owner_cell = DENEB_DATA->GetBdryOwnerCell();
  static const auto& bdry_owner_basis_value =
      DENEB_DATA->GetBdryOwnerBasisValue();
  static const auto& bdry_owner_coefficients =
      DENEB_DATA->GetBdryOwnerCoefficients();
  for (int ibdry = 0; ibdry < num_bdries; ibdry++) {
    const int& num_points = num_bdry_points[ibdry];
    const int& owner_cell = bdry_owner_cell[ibdry];

    int i = 0;
    for (; i < iteration_; i++) {
      avocado::Kernel0::f4(&solution[owner_cell * sb],
                           &bdry_owner_basis_value[ibdry][0], &solution_[0],
                           num_states, num_bases, num_points, 1.0, 0.0);

      bool flag = false;
      for (int ipoint = 0; ipoint < num_points; ipoint++) {
        const std::vector<double>& pressure_fix_values =
            DENEB_EQUATION->ComputePressureFixValues(
                &solution_[ipoint * num_states]);
        for (int ivalue = 0; ivalue < pressure_fix_values.size(); ivalue++) {
          if (pressure_fix_values[ivalue] <= 0.0) {
            flag = true;
            break;
          }
        }
        if (flag) {
          break;
        }
      }
      if (flag) {
        for (int istate = 0; istate < num_states; istate++)
          cblas_dscal(num_bases - 1, 0.5,
                      &solution[owner_cell * sb + istate * num_bases + 1], 1);
      } else
        break;

      if (i == iteration_) {
        for (int istate = 0; istate < num_states; istate++)
          memset(&solution[owner_cell * sb + istate * num_bases + 1], 0,
                 (num_bases - 1) * sizeof(double));
      }
    }
  }

  // outer face sweep
  static const int& num_faces = DENEB_DATA->GetNumFaces();
  for (int iface = num_inner_faces; iface < num_faces; iface++) {
    const int& num_points = num_face_points[iface];
    const int& owner_cell = face_owner_cell[iface];

    int i = 0;
    for (; i < iteration_; i++) {
      avocado::Kernel0::f4(&solution[owner_cell * sb],
                           &face_owner_basis_value[iface][0], &solution_[0],
                           num_states, num_bases, num_points, 1.0, 0.0);

      bool flag = false;
      for (int ipoint = 0; ipoint < num_points; ipoint++) {
        const std::vector<double>& pressure_fix_values =
            DENEB_EQUATION->ComputePressureFixValues(
                &solution_[ipoint * num_states]);
        for (int ivalue = 0; ivalue < pressure_fix_values.size(); ivalue++) {
          if (pressure_fix_values[ivalue] <= 0.0) {
            flag = true;
            break;
          }
        }
        if (flag) {
          break;
        }
      }
      if (flag) {
        for (int istate = 0; istate < num_states; istate++)
          cblas_dscal(num_bases - 1, 0.5,
                      &solution[owner_cell * sb + istate * num_bases + 1], 1);
      } else
        break;

      if (i == iteration_) {
        for (int istate = 0; istate < num_states; istate++)
          memset(&solution[owner_cell * sb + istate * num_bases + 1], 0,
                 (num_bases - 1) * sizeof(double));
      }
    }
  }

  // cell sweep
  static const int& num_cells = DENEB_DATA->GetNumCells();
  static const auto& num_cell_points = DENEB_DATA->GetNumCellPoints();
  static const auto& cell_basis_value = DENEB_DATA->GetCellBasisValue();
  for (int icell = 0; icell < num_cells; icell++) {
    int i = 0;
    for (; i < iteration_; i++) {
      avocado::Kernel0::f4(&solution[icell * sb], &cell_basis_value[icell][0],
                           &solution_[0], num_states, num_bases,
                           num_cell_points[icell], 1.0, 0.0);

      bool flag = false;
      for (int ipoint = 0; ipoint < num_cell_points[icell]; ipoint++) {
        const std::vector<double>& pressure_fix_values =
            DENEB_EQUATION->ComputePressureFixValues(
                &solution_[ipoint * num_states]);
        for (int ivalue = 0; ivalue < pressure_fix_values.size(); ivalue++) {
          if (pressure_fix_values[ivalue] <= 0.0) {
            flag = true;
            break;
          }
        }
        if (flag) {
          break;
        }
      }
      if (flag) {
        for (int istate = 0; istate < num_states; istate++)
          cblas_dscal(num_bases - 1, 0.5,
                      &solution[icell * sb + istate * num_bases + 1], 1);
      } else
        break;
    }
    if (i == iteration_) {
      for (int istate = 0; istate < num_states; istate++)
        memset(&solution[icell * sb + istate * num_bases + 1], 0,
               (num_bases - 1) * sizeof(double));
    }
  }
}
}  // namespace deneb
