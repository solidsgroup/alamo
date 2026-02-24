

#ifndef SCIMITARX_NUMERIC_DEBUGPROBE_H_
#define SCIMITARX_NUMERIC_DEBUGPROBE_H_


#pragma once
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <string>

struct DebugProbe {
  static inline int enabled = 1;     // set 0 to disable
  static inline int lev = 0;
  static inline int i = 32, j = 7, k = 8; // choose a cell on the artifact
  static inline int comp = 3;        // component to print
  static inline int stage = -1;      // RK stage
};

AMREX_FORCE_INLINE
bool ProbeCell(int lev, int i, int j, int k) {
  return DebugProbe::enabled && lev==DebugProbe::lev
      && i==DebugProbe::i && j==DebugProbe::j && k==DebugProbe::k;
}

AMREX_FORCE_INLINE
void ProbePrint(const std::string& tag, amrex::Real v) {
  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << tag << " = " << v << "\n";
  }
}

#endif