## BUILD INSTRUCTIONS

1. connect to the high memory node

```
ssh erp14
```

2. setup the environment by running

```
source erp_env_setup.sh
```

3. build all the examples by running

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
sbatch ./run.sh ./executable arguments
```

for example you can have

```
sbatch ./run.sh ./a4_projection --mesh ./data/1x1_square_quad.mesh --order 1
```
