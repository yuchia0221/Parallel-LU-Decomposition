## Parallel LU Decomposition with OpenMP

#### Table of Contents

-   [Folder Structure](#folder-structure)
-   [Commands](#commands)

---

### Folder Structure

    .
    ├── result                       # Experiment Results
    │   ├── openmp_worker.result     # Result of execution time on each of a different number of threads
    │   ├── openmp_numa.result       # Result of execution time on each of a different number of threads with NUMA
    │   └── data_races_check.result  # Data races check result of Intel's inspector
    ├── lu-omp.cpp                   # Parallelized Version
    ├── lu-omp-numa.cpp              # Parallelized Version with NUMA
    ├── Makefile                     # Recipes for building and running your program
    └── README.md

---

## Commands

Makefile:

> a Makefile that includes recipes for building and running your program. Note: this command will only work on the login node displaying using X11

```bash
make        # builds your code
make runp   # runs a parallel version of your code on W workers
make runs   # runs a serial version of your code on one worker
make check  # runs your parallel code with Intel Thread Checker
make view   # inspect the results of make check with Intel's GUI.
make run-hpc        # creates a HPCToolkit database for performance measurements
make clean          # removes all executable files
make clean-hpc      # removes all HPCToolkit-related files
```

submit.sbatch:

> a script that you can use to launch a batch job that will execute a series of tests on 1..32 threads on a compute node

```bash
sbatch submit.sbatch
```
