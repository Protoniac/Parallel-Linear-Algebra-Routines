# Parallel-Linear-Algebra-Routines
In this repository are some tools for parallel linear algebra routines using both C/C++ and MPI.

## Parallel Block Matrices Multiplication

### Prerequisites
* MPIC++
* MPIEXEC
* Electricity

### Install MPI on Ubuntu

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
