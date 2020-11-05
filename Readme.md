# BESOneo2
A fast and efficient BESO topology optimization function for 2D minimum compliance topology optimization subject to a volume constraint.

## About

2D Bi-Directional Evolutionary Topology Optimization (BESO) for computing the optimal minimum compliance structure subject to a given volume constraint. BESO is based on the code presented by [Huang, X., & Xie, Y. M. (2010)](https://doi.org/10.1002/9780470689486), and significant speedups are realized through reduced indexing and efficient matrix construction as detailed in [Ferrari, F., & Sigmund, O. (2020)](https://doi.org/10.1007/s00158-020-02629-w).

### inputs
Parameter | Description
--------|-------
`nx` `ny` | design domain size in the x & y direction
`volfrac` | volume constraint where *0 < `volfrac` <= 1*
`er` | Evolutionary rate
`rmin` | Sensitivity filter radius to solve mesh dependency issues

### features

Specification of non-design regions can be achieved by adding elements to the `pasS` and `pasV` sets in code section C. This will exclude the specified elements from the optimization process and set them permanently as solid or void elements.

#### Example: Non-designable circular region
Add the following to code section C
```matlab
    [cx, cy, cr] = deal(nx/2, ny/2, ny/3);
    [dy, dx] = meshgrid(1:nx, 1:ny);
    f = sqrt((dx-cy).^2+(dy-cx).^2) < cr;
    [pasS, pasV] = deal([],[elNrs(f)]);
```
Calling BESOneo with `[x, obj] = BESOneo2(800,400,0.3,0.02,8)` results in:
![Passive Void example](https://octodex.github.com/images/yaktocat.png)

## Dependencies

BESOneo uses the fast sparse matrix construction routine `fsparse` [Engblom, S., & Lukarski, D. (2016)](https://doi.org/10.1016/j.parco.2016.04.001) which is downloadable from the first authors GitHub page at: https://github.com/stefanengblom/stenglib.
Alternately, the code can be run without dependencies by changing all `fsparse` calls to MATLABs inbuilt `sparse` function.

## Use

`nx` and `ny` specify the number of elements in the x and y direction. `volfrac` can take any positive value up to 1.0 and specifies the volume constraint as a percentage of the total design domain. `er` is the evolutionary ratio and is generally set to 2% (0.02). `rmin` is the filter radius and can take any positive value >= 1
The cantilever example cal be called using `[x, obj] = BESOneo2(400,200,0.3,0.02,3)`
Definitions for a frame reinforcement problem and an L-bracket are included in BESOneo.m and are commented after the main code

## References

Huang, X., & Xie, Y. M. (2010). Evolutionary Topology Optimization of Continuum Structures: Methods and Applications. https://doi.org/10.1002/9780470689486

Ferrari, F., & Sigmund, O. (2020). A new generation 99 line Matlab code for compliance Topology Optimization and its extension to 3D. ArXiv Preprint ArXiv:2005.05436.

Engblom, S., & Lukarski, D. (2016). Fast Matlab compatible sparse assembly on multicore computers. Parallel Computing, 56, 1–17. https://doi.org/10.1016/j.parco.2016.04.001
