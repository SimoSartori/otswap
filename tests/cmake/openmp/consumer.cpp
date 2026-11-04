// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

#include <omp.h>

extern "C" int openmp_consumer_threads ()
{
  int threads = 0;
#pragma omp parallel
  {
#pragma omp atomic
    ++threads;
  }
  return threads;
}
