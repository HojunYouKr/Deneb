#include "avocado_memory_profiler.h"

#include "avocado_string.h"
#include "avocado_mpi.h"

#define BUFFER_SIZE 100

namespace avocado {
std::shared_ptr<MemoryProfiler> AVOCADO_MEMORY_PROFILER_NAME = nullptr;
MemoryProfiler::MemoryProfiler()
    : hostname_(GetHostName()),
      coretag_(0),
      memory_used_(),
      hostname_table_(),
      coretag_table_(),
      sort_mapping_(),
      host_ptr_() {
  hostname_table_.push_back(hostname_);
  coretag_table_.push_back(coretag_);
  sort_mapping_.push_back(0);
  host_ptr_.push_back(0);
  host_ptr_.push_back(1);

  const int ndomain = NDOMAIN;
  const int myrank = MYRANK;
  if (ndomain == 1) return;

  const int data_send = static_cast<int>(hostname_.length());
  std::vector<int> num_data_recv(NDOMAIN);
  MPI_Gather(&data_send, 1, MPI_INT, &num_data_recv[0], 1, MPI_INT, MASTER_NODE,
             MPI_COMM_WORLD);

  std::vector<std::vector<char>> recv_names(ndomain);
  if (myrank == MASTER_NODE)
    for (int idomain = 0; idomain < ndomain; idomain++)
      recv_names[idomain].resize(num_data_recv[idomain]);

  std::vector<MPI_Request> recv_req(ndomain);
  std::vector<MPI_Status> recv_sta(ndomain);
  if (myrank == MASTER_NODE) {
    for (int idomain = 0; idomain < ndomain; idomain++)
      MPI_Irecv(&recv_names[idomain][0], num_data_recv[idomain], MPI_CHAR,
                idomain, 78, MPI_COMM_WORLD, &recv_req[idomain]);
  }
  MPI_Request send_req;
  MPI_Status send_sta;
  MPI_Isend(hostname_.c_str(), data_send, MPI_CHAR, MASTER_NODE, 78,
            MPI_COMM_WORLD, &send_req);
  if (myrank == MASTER_NODE) MPI_Waitall(ndomain, &recv_req[0], &recv_sta[0]);
  MPI_Wait(&send_req, &send_sta);

  if (myrank == MASTER_NODE) {
    hostname_table_.resize(ndomain);
    for (int idomain = 0; idomain < ndomain; idomain++) {
      recv_names[idomain].push_back(0);
      hostname_table_[idomain] = std::string(&recv_names[idomain][0]);
    }

    std::vector<std::string> table;
    for (int idomain = 0; idomain < ndomain; idomain++) {
      const int size = table.size();
      if (size == 0)
        table.push_back(hostname_table_[idomain]);
      else
        for (int i = 0; i < size; i++) {
          const int var = table[i].compare(hostname_table_[idomain]);
          if (var < 0 && i == table.size() - 1) {
            table.push_back(hostname_table_[idomain]);
            break;
          }
          if (var > 0) {
            table.insert(table.begin() + i, hostname_table_[idomain]);
            break;
          }
          if (var == 0) break;
        }
    }

    const int name_types = table.size();
    std::vector<int> index(name_types, 0);
    coretag_table_.resize(ndomain);
    for (int idomain = 0; idomain < ndomain; idomain++) {
      int i = 0;
      for (i = 0; i < name_types; i++)
        if (table[i].compare(hostname_table_[idomain]) == 0) break;
      coretag_table_[idomain] = index[i]++;
    }

    sort_mapping_.resize(ndomain);
    int ind = 0;
    host_ptr_.clear();
    host_ptr_.push_back(0);
    for (auto&& tab : table) {
      for (int idomain = 0; idomain < ndomain; idomain++)
        if (tab.compare(hostname_table_[idomain]) == 0)
          sort_mapping_[ind++] = idomain;
      host_ptr_.push_back(ind);
    }
  }
  SYNCRO();
}
void MemoryProfiler::PrintMemoryUsage(void) {
  const double mem_used =
      static_cast<double>(GetCurrentRSS()) / 1024.0 / 1024.0;
  const double mem_size =
      static_cast<double>(GetMemorySize()) / 1024.0 / 1024.0;

  if (NDOMAIN > 1) {
    std::vector<double> mem_used_gather(NDOMAIN);
    std::vector<double> mem_size_gather(NDOMAIN);
    MPI_Gather(&mem_used, 1, MPI_DOUBLE, &mem_used_gather[0], 1, MPI_DOUBLE,
               MASTER_NODE, MPI_COMM_WORLD);
    MPI_Gather(&mem_size, 1, MPI_DOUBLE, &mem_size_gather[0], 1, MPI_DOUBLE,
               MASTER_NODE, MPI_COMM_WORLD);

    MASTER_MESSAGE(avocado::GetTitle("Memory Usage Table (Total)"));
    std::stringstream str;
    if (MYRANK == MASTER_NODE) {
      double total_mem_used = 0.0;
      for (int idomain = 0; idomain < NDOMAIN; idomain++)
        total_mem_used += mem_used_gather[idomain];

      double total_mem_size = 0.0;
      for (int i = 0; i < host_ptr_.size() - 1; i++)
        total_mem_size += mem_size_gather[sort_mapping_[host_ptr_[i]]];

      str << "Total memory size: " << static_cast<int>(total_mem_size) << " Mib"
          << std::endl;
      str << "Total memory usage: " << static_cast<int>(total_mem_used)
          << " Mib (" << std::to_string(total_mem_used / total_mem_size * 100.0)
          << " %)" << std::endl;
    }
    MASTER_MESSAGE(str.str());
    
    MASTER_MESSAGE(avocado::GetTitle("Memory Usage Table (Summary)"));
    str.str("");
    if (MYRANK == MASTER_NODE) {
      for (int i = 0; i < host_ptr_.size() - 1; i++) {
        double mem_sum = 0.0;
        for (int idomain = host_ptr_[i]; idomain < host_ptr_[i + 1]; idomain++)
          mem_sum += mem_used_gather[sort_mapping_[idomain]];
        str << "Host: " << hostname_table_[sort_mapping_[host_ptr_[i]]]
            << "\tMem used: " << static_cast<int>(mem_sum) << " Mib / "
            << static_cast<int>(mem_size_gather[sort_mapping_[host_ptr_[i]]])
            << " Mib";
        str << " ("
            << std::to_string(mem_sum /
                              mem_size_gather[sort_mapping_[host_ptr_[i]]] *
                              100.0)
            << " %)" << std::endl;
      }
    }
    MASTER_MESSAGE(str.str());

    MASTER_MESSAGE(avocado::GetTitle("Memory Usage Table (Details)"));
    str.str("");
    if (MYRANK == MASTER_NODE) {
      for (int i = 0; i < host_ptr_.size() - 1; i++) {
        double mem_sum = 0.0;
        for (int idomain = host_ptr_[i]; idomain < host_ptr_[i + 1]; idomain++)
          mem_sum += mem_used_gather[sort_mapping_[idomain]];
        str << "Host: " << hostname_table_[sort_mapping_[host_ptr_[i]]]
            << "\tMem used: " << static_cast<int>(mem_sum) << " Mib / "
            << static_cast<int>(mem_size_gather[sort_mapping_[host_ptr_[i]]])
            << " Mib";
        str << " ("
            << std::to_string(mem_sum /
                              mem_size_gather[sort_mapping_[host_ptr_[i]]] *
                              100.0)
            << " %)" << std::endl;
        for (int idomain = host_ptr_[i]; idomain < host_ptr_[i + 1];
             idomain++) {
          str << "\tRank: " << sort_mapping_[idomain] << "\tMem used: "
              << static_cast<int>(mem_used_gather[sort_mapping_[idomain]])
              << " Mib";
          str << " ("
              << std::to_string(mem_used_gather[sort_mapping_[idomain]] /
                                mem_size_gather[sort_mapping_[host_ptr_[i]]] *
                                100.0)
              << " %)" << std::endl;
        }
        if (i != host_ptr_.size() - 2) str << std::endl;
      }
    }
    MASTER_MESSAGE(str.str());
    SYNCRO();
  } else {
    MASTER_MESSAGE(avocado::GetTitle("Memory Usage Table (Single Core)"));
    std::stringstream str;
    str << "Total memory size: " << static_cast<int>(mem_size) << " Mib"
        << std::endl;
    str << "Total memory usage: " << static_cast<int>(mem_used) << " Mib ("
        << std::to_string(mem_used / mem_size * 100.0) << " %)" << std::endl;
    MASTER_MESSAGE(str.str());

    MASTER_MESSAGE(avocado::GetTitle("Memory Usage Table (Details)"));
    str.str("");
    str << "Host: " << hostname_ << "\tMem used: " << static_cast<int>(mem_used)
        << " Mib / " << static_cast<int>(mem_size) << " Mib";
    str << " (" << std::to_string(mem_used / mem_size * 100.0) << " %)"
        << std::endl;
    MASTER_MESSAGE(str.str());
  }
}
void MemoryProfiler::MemoryCheck(void) {
  memory_used_[0] = memory_used_[1];
  const double mem_used =
      static_cast<double>(GetCurrentRSS()) / 1024.0 / 1024.0;
  MPI_Allreduce(&mem_used, &memory_used_[1], 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
}

const std::string MemoryProfiler::GetHostName(void) {
  char buf[BUFFER_SIZE];
#ifdef AVOCADO_USE_WINDOWS
  DWORD bufCharCount = BUFFER_SIZE;
  GetComputerName(buf, &bufCharCount);
#else
  gethostname(buf, BUFFER_SIZE);
#endif
  return std::string(buf);
}

size_t MemoryProfiler::GetMemorySize(void) {
#if defined(_WIN32) && (defined(__CYGWIN__) || defined(__CYGWIN32__))
  /* Cygwin under Windows. ------------------------------------ */
  /* New 64-bit MEMORYSTATUSEX isn't available.  Use old 32.bit */
  MEMORYSTATUS status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatus(&status);
  return (size_t)status.dwTotalPhys;

#elif defined(_WIN32)
  /* Windows. ------------------------------------------------- */
  /* Use new 64-bit MEMORYSTATUSEX, not old 32-bit MEMORYSTATUS */
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  return (size_t)status.ullTotalPhys;

#elif defined(__unix__) || defined(__unix) || defined(unix) || \
    (defined(__APPLE__) && defined(__MACH__))
  /* UNIX variants. ------------------------------------------- */
  /* Prefer sysctl() over sysconf() except sysctl() HW_REALMEM and HW_PHYSMEM */

#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))
  int mib[2];
  mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
  mib[1] = HW_MEMSIZE; /* OSX. --------------------- */
#elif defined(HW_PHYSMEM64)
  mib[1] = HW_PHYSMEM64; /* NetBSD, OpenBSD. --------- */
#endif
  int64_t size = 0;    /* 64-bit */
  size_t len = sizeof(size);
  if (sysctl(mib, 2, &size, &len, NULL, 0) == 0) return (size_t)size;
  return 0L; /* Failed? */

#elif defined(_SC_AIX_REALMEM)
  /* AIX. ----------------------------------------------------- */
  return (size_t)sysconf(_SC_AIX_REALMEM) * (size_t)1024L;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
  /* FreeBSD, Linux, OpenBSD, and Solaris. -------------------- */
  return (size_t)sysconf(_SC_PHYS_PAGES) * (size_t)sysconf(_SC_PAGESIZE);

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
  /* Legacy. -------------------------------------------------- */
  return (size_t)sysconf(_SC_PHYS_PAGES) * (size_t)sysconf(_SC_PAGE_SIZE);

#elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))
  /* DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX. -------- */
  int mib[2];
  mib[0] = CTL_HW;
#if defined(HW_REALMEM)
  mib[1] = HW_REALMEM;   /* FreeBSD. ----------------- */
#elif defined(HW_PYSMEM)
  mib[1] = HW_PHYSMEM; /* Others. ------------------ */
#endif
  unsigned int size = 0; /* 32-bit */
  size_t len = sizeof(size);
  if (sysctl(mib, 2, &size, &len, NULL, 0) == 0) return (size_t)size;
  return 0L; /* Failed? */
#endif /* sysctl and sysconf variants */

#else
  return 0L; /* Unknown OS. */
#endif
}

size_t MemoryProfiler::GetCurrentRSS(void) {
#if defined(_WIN32)
  /* Windows -------------------------------------------------- */
  PROCESS_MEMORY_COUNTERS info;
  GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
  return (size_t)info.WorkingSetSize;
#elif defined(__linux__) || defined(__linux) || defined(linux) || \
    defined(__gnu_linux__)
  /* Linux ---------------------------------------------------- */
  long rss = 0L;
  FILE* fp = NULL;
  if ((fp = fopen("/proc/self/statm", "r")) == NULL)
    return (size_t)0L; /* Can't open? */
  if (fscanf(fp, "%*s%ld", &rss) != 1) {
    fclose(fp);
    return (size_t)0L; /* Can't read? */
  }
  fclose(fp);
  return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);

#else
  /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
  return (size_t)0L; /* Unsupported. */
#endif
}
}  // namespace avocado