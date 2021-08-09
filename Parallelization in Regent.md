# Parallelization in Regent



## Introduction

To make `task`s run in parallel in Regent, you annotate them with `__demand(__cuda)`. This does most of the work for you, but you do have to ensure that the Cuda code generation is possible in the first place. To do so, there are certain patterns you must follow. We will go over these in this guide.



## For Loops

Every loop has to be in for-each style. They cannot use numerical bounds.

```c++
-- DON'T --
for i = 0, nCells do
	for j = 0, nVertLevels do
		cr[{i, j}].dss = 5
  end
end
```

```c++
-- DO --
var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
for iCell in cell_range do
	cr[iCell].dss = 5
end
```

Note that for loops can be combined as in the example above. Additionally, rewriting the for loop only needs to be done for the outermost loop in the `task`.



## Loop-carried dependency

A loop-carried dependency is if you read and write to the same memory location in multiple iterations of the loop. There is no easy fix for this. Sometimes the function needs to be massively rewritten.

```c++
-- Example --
var cell_range = rect2d { int2d {0, 0}, int2d {nCells - 1, nVertLevels - 1} }
for iCell in cell_range do
  if (iCell.y > 0) then
		cr[iCell].dss = cr[iCell].dss + cr[{iCell.x, iCell.y - 1 }]
  end
end
```

For more information check out slides 57 and onwards [here](https://www3.nd.edu/~zxu2/acms60212-40212/Lec-12-OpenMP.pdf).

### Common patterns

#### Min / Max of an array

```c++
-- ERROR --
var max_dss = 0
for iCell in cell_range do
	max_dss = max(max_dss, cr[iCell].dss)
end
```

```c++
-- FIX --
var max_dss = 0
for iCell in cell_range do
	max_dss max= cr[iCell].dss
end
```



## Stack-allocated arrays

The following will not compile since there is both region access and an assignment to a stack-allocated array in the same for loop.

```c++
__demand(__cuda)
task example2(reg : region(ispace(int1d), double), n : int)
where
  writes (reg)
do
  var range_1d = rect1d {0, n - 1}
  var a : double[n]
 
  --__demand(__openmp)
  for i in range_1d do
    a[int(i)] = 1
    reg[int(i)] = 1
  end
end
```

Stack-allocated arrays don’t get assigned to a GPU, only regions do. This means that…

1. For loops containing only stack-allocated arrays don’t get run in parallel since they are run “host side”.
2. CUDA simply ignores certain ineligible loops instead of converting them, whereas OpenMP would complain. (OpenMP operates on a per-loop basis instead of a per-task basis). That's why sometimes `__demand(__openmp)` will raise an error while `__demand(__cuda)` will not.
3. One cannot write to a stack-allocated array in a loop that runs in parallel.

For more info, read this [issue](https://github.com/StanfordLegion/legion/issues/1124).



## Debugging Tips

As one will quickly be able to tell, the compiler's error messages are not very useful. We go over some common error messages here.  

- `CUDA code generation failed: found a region access outside parallelizable loops`. 

  - This error usually points to the beginning of a for loop. Use the annotation `__demand(__openmp)` directly in front of the for loop, and run your program again. This time the error message will be more informative. 
  - In the case this actually points to a region access outside a for loop, you will have to wrap the block of code with a dummy for loop that only runs once.

- `option __demand(__cuda) is not permitted for non-leaf task`. This error is thrown when you try to annotate a `task` that calls other tasks. The Cuda code generation does not work for such tasks since the overhead for each function is too big. 

- Math statements cause an `inadmissable statement` error. 

  ```c++
  -- Error --
  cmath = terralib.includec("math.h")
  for iCell in cell_range_1d do
    cr[iCell].hx = cmath.pow(etavs, 1.5)
  end
  ```

  ```c++
  -- Fix --
  local pow = regentlib.pow(double)
  for iCell in cell_range_1d do
    cr[iCell].hx = pow(etavs, 1.5)
  end
  ```

- Currently there is an [issue](https://github.com/StanfordLegion/legion/issues/1121) when using a `rect1d` to loop over a region. Explicitly cast the iterator to an int with `int(x)` as a temporary fix.  

 
