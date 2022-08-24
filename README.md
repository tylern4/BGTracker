# BGTracker

* build with make (must have nvcc eg module load cuda for nsls2)

## Run a test

```
./test.sh
```

## Using with pybind11

Building with pybind11 requires an extra git submodule to be built. To download with the submodule `git clone --recurse-submodules https://github.com/tylern4/BGTracker.git`. Once downloaded you can follow a similar cmake build setup.

```sh
mkdir build; cd build
cmake .. 
make -j2
```

That will build a python module called `pyGBCuda.cpython-ver-arch-os-etc.so` which needs to be included in your `$PYTHONPATH` for you to be able to run the module. Then to access and run the modules:

```python
import pyGBCuda as pylatt

quad = pylatt.Quad("qh1g2c30a", 0.268, -0.641957314648, NKick=100)
# Runs on GPU with #particles, #nturns
quad.run_qsympass4(1000, 1000)
```