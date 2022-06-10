#include "deneb_timescheme.h"

#include <algorithm>

#include "avocado.h"
#include "deneb_config_macro.h"
#include "deneb_contour.h"
#include "deneb_equation.h"
#include "deneb_saveload.h"
#include "deneb_timescheme_impeuler.h"
#include "deneb_timescheme_rosenbrock.h"
#include "deneb_timescheme_ssprk.h"
#include "deneb_timescheme_steady.h"

namespace deneb {
bool Interrupt::CheckIterationFinish(const int iteration) {
  if (control_ == Control::ITER)
    if (finish_ <= iteration) return true;
  return false;
}
bool Interrupt::CheckTimeFinish(const unsigned __int64& time,
                                unsigned __int64& timestep) {
  if (control_ == Control::TIME)
    if (finish_ <= (time + timestep)) {
      timestep = finish_ - time;
      return true;
    }
  return false;
}
bool Interrupt::CheckIterationEvent(const int iteration) {
  if (control_ == Control::ITER)
    if (begin_ <= iteration &&
        (finish_ == static_cast<unsigned __int64>(-1) || iteration <= finish_))
      if ((iteration - begin_) % interval_ == 0) return true;
  return false;
}
bool Interrupt::CheckTimeEvent(const unsigned __int64& time,
                               unsigned __int64& timestep) {
  if (control_ == Control::TIME) {
    if (step_ == -1) {
      step_ = 0;
      while (!(time < begin_ + step_ * interval_)) step_++;
    }

    if (begin_ <= time &&
        (finish_ == static_cast<unsigned __int64>(-1) || time <= finish_)) {
      if (begin_ + step_ * interval_ <= time + timestep) {
        timestep = begin_ + step_ * interval_ - time;
        step_++;
        return true;
      }
    }
  }
  return false;
}

std::shared_ptr<Timescheme> DENEB_TIMESCHEME_NAME = nullptr;
std::shared_ptr<Timescheme> Timescheme::GetTimescheme(const std::string& name) {
  if (!name.compare("SSPRK11") || !name.compare("ExplicitEuler"))
    return std::make_shared<TimeschemeSSPRK11>();
  else if (!name.compare("SSPRK33") || !name.compare("TVDRK"))
    return std::make_shared<TimeschemeSSPRK33>();
  else if (!name.compare("SSPRK54") || !name.compare("SSPRK"))
    return std::make_shared<TimeschemeSSPRK54>();
  else if (!name.compare("ImplicitEuler"))
    return std::make_shared<TimeschemeImpEuler>();
  else if (!name.compare("ROS3P") || !name.compare("ROSEN33"))
    return std::make_shared<TimeschemeROS3P>();
  else if (!name.compare("RODAS3") || !name.compare("ROSEN43"))
    return std::make_shared<TimeschemeRODAS3>();
  else if (!name.compare("ROS4") || !name.compare("ROSEN44"))
    return std::make_shared<TimeschemeROS4>();
  else if (!name.compare("RODASP") || !name.compare("ROSEN64"))
    return std::make_shared<TimeschemeRODASP>();
  else if (!name.compare("RODAS5") || !name.compare("ROSEN85"))
    return std::make_shared<TimeschemeRODAS5>();
  else if (!name.compare("ROW6A") || !name.compare("ROSEN66"))
    return std::make_shared<TimeschemeROW6A>();
  else if (!name.compare("Steady"))
    return std::make_shared<TimeschemeSteady>();
  ERROR_MESSAGE("Wrong time scheme (no-exist):" + name + "\n");
  return nullptr;
}

Timescheme::Timescheme(const bool is_implicit) : is_implicit_(is_implicit) {
  auto& config = AVOCADO_CONFIG;
  order_ = std::stoi(config->GetConfigValue(ORDER));
  max_iteration_ = std::stoi(config->GetConfigValue(MAX_ITERATION));
  iteration_ = 0;

  time_resolution_ = std::stoi(config->GetConfigValue(TIME_RESOLUTION));
  time_scale_ = static_cast<unsigned __int64>(std::pow(10, time_resolution_));
  current_time_ = 0;

  const std::string& timestep_control =
      config->GetConfigValue(TIMESTEP_CONTROL);
  timestep_control_value_ =
      std::stod(config->GetConfigValue(TIMESTEP_CONTROL_VALUE));
  if (!timestep_control.compare("CFL"))
    timestep_control_ = TimestepControl::CFL;
  else if (!timestep_control.compare("dt"))
    timestep_control_ = TimestepControl::DT;
  else
    ERROR_MESSAGE("Wrong time step control (no-exist): " + timestep_control +
                  "\n");

  {
    const std::string& text = config->GetConfigValue(STOP_CONDITION);
    Interrupt::Control control;
    unsigned __int64 finish;
    if (!text.compare("Iter")) {
      control = Interrupt::Control::ITER;
      finish = static_cast<__int64>(
          std::stoi(config->GetConfigValue(STOP_CONTROL_VALUE)));
    } else if (!text.compare("Time")) {
      control = Interrupt::Control::TIME;
      finish =
          ConvertTime(std::stod(config->GetConfigValue(STOP_CONTROL_VALUE)));
    } else
      ERROR_MESSAGE("Wrong stop control (no-exist): " + text + "\n");
    stop_.Initialize(control, 0, finish, 0);
  }

  {
    const std::string& text = config->GetConfigValue(POST_CONTROL);
    Interrupt::Control control;
    std::vector<unsigned __int64> values(3);
    if (!text.compare("Iter")) {
      control = Interrupt::Control::ITER;
      for (int i = 0; i < 3; i++)
        values[i] = static_cast<__int64>(
            std::stoi(config->GetConfigValue(POST_CONTROL_VALUE(i))));
    } else if (!text.compare("Time")) {
      control = Interrupt::Control::TIME;
      for (int i = 0; i < 3; i++)
        values[i] = ConvertTime(
            std::stod(config->GetConfigValue(POST_CONTROL_VALUE(i))));
    } else
      ERROR_MESSAGE("Wrong post control (no-exist): " + text + "\n");
    post_.Initialize(control, values[0], values[1], values[2]);
  }

  {
    const std::string& text = config->GetConfigValue(SAVE_CONTROL);
    Interrupt::Control control;
    std::vector<unsigned __int64> values(3);
    if (!text.compare("Iter")) {
      control = Interrupt::Control::ITER;
      for (int i = 0; i < 3; i++)
        values[i] = static_cast<__int64>(
            std::stoi(config->GetConfigValue(SAVE_CONTROL_VALUE(i))));
    } else if (!text.compare("Time")) {
      control = Interrupt::Control::TIME;
      for (int i = 0; i < 3; i++)
        values[i] = ConvertTime(
            std::stod(config->GetConfigValue(SAVE_CONTROL_VALUE(i))));
    } else
      ERROR_MESSAGE("Wrong save control (no-exist): " + text + "\n");
    save_.Initialize(control, values[0], values[1], values[2]);
  }

  const std::string& restart = config->GetConfigValue(RESTART);
  if (!restart.compare("On")) {
    const std::string& restart_filepath = config->GetConfigValue(RESTART_PATH);
    if (avocado::pathattr::FILE != avocado::GetPathAttribute(restart_filepath))
      ERROR_MESSAGE("Restart file path is wrong: " + restart_filepath + "\n");
  }
}
void Timescheme::ComputeLocalTimestep(
    const double* solution, std::vector<double>& local_timestep) const {
  if (timestep_control_ == TimestepControl::DT) {
    for (auto&& timestep : local_timestep) timestep = timestep_control_value_;
    return;
  }

  DENEB_EQUATION->ComputeLocalTimestep(solution, local_timestep);
  avocado::VecScale(static_cast<int>(local_timestep.size()),
                    timestep_control_value_, &local_timestep[0]);
}
double Timescheme::ComputeGlobalTimestep(
    const double* solution, std::vector<double>& local_timestep) const {
  if (timestep_control_ == TimestepControl::DT) return timestep_control_value_;

  DENEB_EQUATION->ComputeLocalTimestep(solution, local_timestep);
  const double min_dt = AVOCADO_MPI->Reduce(
      *std::min_element(local_timestep.begin(), local_timestep.end()),
      avocado::MPI::Op::MIN);

  if (!std::isnormal(min_dt))
    ERROR_MESSAGE("Wrong time step.\n\tdt = " + std::to_string(min_dt) + "\n");
  return min_dt * timestep_control_value_;
}
void Timescheme::InitSolution(void) {
  START_TIMER();
  auto& config = AVOCADO_CONFIG;
  const std::string& restart = config->GetConfigValue(RESTART);
  if (!restart.compare("On")) {
    const std::string& restart_filepath = config->GetConfigValue(RESTART_PATH);
    DENEB_SAVELOAD->Load(restart_filepath);
    solution_ = DENEB_SAVELOAD->GetSolution();
    DENEB_SAVELOAD->ClearSolution();
    const SaveLoad::SaveData& data = DENEB_SAVELOAD->GetSaveData();
    current_time_ = data.time_;
    iteration_ = data.iteration_;
    DENEB_CONTOUR->SetStrandID(data.strandid_);
  } else {
    const double t = GetCurrentTime();
    MASTER_MESSAGE("Initilizing solution at t = " + std::to_string(t) +
                   " ... ");
    solution_.resize(length_);
    DENEB_EQUATION->ComputeInitialSolution(&solution_[0], t);
  }
  SYNCRO();
  MASTER_MESSAGE("Complete. (Time: " + std::to_string(STOP_TIMER()) + "s)\n");
}
}  // namespace deneb