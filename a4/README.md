## BUILD INSTRUCTIONS

1. connect - do this each time you start work

```
ssh erp01
```

2. setup the environment by running

```
source erp_env_setup.sh
```

3. build all the examples by running

```
make
```

By this point you should see the executables in the folder, which will be used to run the cases.


Note that if you need to clean your build (i.e., removing the executables) you can use
```
make clean
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
