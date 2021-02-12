#pragma once

#include "exdot_serial.h"
#include "thrust/device_vector.h"
#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_CUDA
#include "exdot_cuda.cuh" // accumulate.cuh , config.h, mylibm.cuh
#else
#include "exdot_omp.h" //accumulate.h, mylibm.hpp
#endif

#ifdef MPI_VERSION
#include "mpi_accumulate.h"
#endif //MPI_VERSION
