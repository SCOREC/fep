## BUILD INSTRUCTIONS

1. allocate a himem or debug node to compile


```
salloc -n 1 -N 1 -t 60 -p himem
```

or

```
salloc -n 1 -N 1 -t 60 -p debug
```

2. get the node name with 'squeue -l'; i.e., 'erp14' and ssh to it

```
mynode=`squeue -l | awk '/erp/ {print $9}'`
ssh $mynode
```

3. setup the environment by running

```
source erp_env_setup.sh.sh
```

5. build all the examples by running

```
MFEM_INSTALL_DIR=$MFEM_ROOT make
```

to clean your build you can use
```
MFEM_INSTALL_DIR=$MFEM_ROOT make clean
```

## RUN INSTRUCTIONS

to run an executables use `sbatch` and the provided `run.sh` script as follows

```
sbatch ./run.sh .executable arguments
```

for example you can have

```
sbatch ./run.sh ./a4_projection --mesh ./data/1x1_square_quad.mesh --order 1
```
