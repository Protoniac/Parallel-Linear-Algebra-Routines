# Parallel-Linear-Algebra-Routines
In this repository are some tools for parallel linear algebra routines using both C/C++ and MPI.

## Parallel Block Matrices Multiplication

Multiplication of two matrices A and B by scattering them in multiples block matrices are well suited for parallelism.

Below is explained the different step of this method : 

#### 1) Scatter Matrices A and B into multiple blocks

In each thread will be allocated a corresponding part of *A* and *B* called *A_block* and *B_block*, respectively, following the thread's *rank*. The game of this first step is to find the perfect way to scatter *A* and *B* in blocks, following their row and column dimensions.

#### 2) Compute each block 

In coming ...

#### 3) Gather the computed blocks into the final matrix AB

In coming ...

### Prerequisites
* MPIC++
* MPIEXEC
* Electricity

### Install MPI on Ubuntu / Debian

```shell
sudo apt install mpich
```

### Compilation

Clone the git project with : 
```shell
git clone https://github.com/Protoniac/Parallel-Linear-Algebra-Routines.git
```
Compile and execute :

```shell
mpic++ -c *.cpp -Wall
mpic++ -o main *.o -Wall
mpirun -n NUMBEROFTHREAD ./main
```
The square root of the number of thread given as parameter should be a nutural number.

If the above commands are not working, it might be possible that the symobolic links have not been created properly.

```shell
ln -s /usr/bin/mpic++ mpic++
ln -s /usr/bin/mpirun mpirun
```

Generate NxM matrix :
```shell
cd matrix_samples
python matrix_generator.py N M 1
```
