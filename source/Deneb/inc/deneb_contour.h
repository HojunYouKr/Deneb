#pragma once

#include <memory>
#include <string>
#include <vector>

#define DENEB_CONTOUR_NAME contour_global_ptr
#define DENEB_CONTOUR deneb::DENEB_CONTOUR_NAME
#define DENEB_CONTOUR_INITIALIZE(cell_post_order, face_post_order) \
  DENEB_CONTOUR =                                                  \
      std::make_shared<deneb::Contour>(cell_post_order, face_post_order)
#define DENEB_CONTOUR_FINALIZE() DENEB_CONTOUR.reset()

namespace deneb {
class Contour {
 private:
  enum class FileType : int { FULL = 0, GRIDONLY = 1, SOLUTION = 2 };
  enum class ZoneType : int { FELINE = 1, FETRIS = 2, FETETS = 4 };
  struct ContourData {
    int num_total_nodes_ = 0;
    int num_total_cells_ = 0;
    std::vector<int> num_nodes_;
    std::vector<int> connectivity_;
    std::vector<std::vector<float>> coords_;
    std::vector<std::vector<double>> basis_;
    std::vector<std::vector<double>> basis_grad_;

    int num_target_cells_ = 0;
    std::vector<int> target_cells_;
    std::vector<std::vector<double>> normal_;
  };

 private:
  int strandid_;
  int cell_post_order_;
  int face_post_order_;

  int dimension_;
  std::vector<std::string> cell_variables_;
  std::vector<std::string> face_variables_;

  ContourData cell_;

  int num_faces_;
  std::vector<int> face_tags_;
  std::vector<std::string> face_names_;
  std::vector<ContourData> faces_;

  std::vector<double> sol_;
  std::vector<double> sol_grad_;
  std::vector<double> sol_post_;

 public:
  Contour(const int cell_post_order, const int face_post_order);
  ~Contour(){};

  inline void SetStrandID(const int strandid) { strandid_ = strandid; };
  inline int GetStrandID(void) const { return strandid_; };

  void BuildData(void);
  void CellGrid(const std::string& filename);
  void CellSolution(const std::string& filename, const double* solution,
                    const double time);
  void FaceGrid(const std::string& filename);
  void FaceSolution(const std::string& filename, const double* solution,
                    const double time);

 private:
  void BuildCellData(void);
  void BuildFaceData(void);

  void WriteTecplotVersion(std::vector<int>& buffer) const;
  void WriteIntData(std::vector<int>& buffer, const int data) const;
  void WriteReal32Data(std::vector<int>& buffer, const float data) const;
  void WriteReal64Data(std::vector<int>& buffer, const double data) const;
  void WriteStringData(std::vector<int>& buffer, const std::string& data) const;
};
extern std::shared_ptr<Contour> DENEB_CONTOUR_NAME;
}  // namespace deneb