# Optimize-ode
Computational optimization of  larger non-linear differential equation system using DifferentialEquations.jl

### Mathematical Formulation of the ODE system for a Stick-slip model using Rate-and-State Friction Velocity

The system models the interaction of <a href="https://www.codecogs.com/eqnedit.php?latex=l&space;\times&space;d" target="_blank"><img src="https://latex.codecogs.com/gif.latex?l&space;\times&space;d" title="l \times d" /></a>  blocks. Each block has 3 degrees of freedom corresponding position **u** , velocity **v** and state **&theta;** . A figure of the physical setup is provided in [Stick-Slip Model](/figures/RSA.pdf).

The governing equations of motions for each block are prescribed as in the following:

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{cases}&space;\dot{u_i}=v_i-v_0\\&space;\dot{v_i}=(-1/m)\Big(k_c(u_{i&plus;1}&plus;u_{j&plus;1}-4u_i&plus;u_{i-1}&plus;u_{j-1})&plus;k_pu_i&plus;\sigma_{n}a\Big[\sinh^{-1}\left(\frac{V}{2v_0}e^{\frac{\theta}{a}}\right)\Big]\Big)\\&space;\dot{\theta_i}=\frac{bv_0}{D_c}\left(e^{\frac{\tau_0-\theta}{b}}-\frac{v}{v_0}\right)\\&space;\end{cases}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{cases}&space;\dot{u_i}=v_i-v_0\\&space;\dot{v_i}=(-1/m)\Big(k_c(u_{i&plus;1}&plus;u_{j&plus;1}-4u_i&plus;u_{i-1}&plus;u_{j-1})&plus;k_pu_i&plus;\sigma_{n}a\Big[\sinh^{-1}\left(\frac{V}{2v_0}e^{\frac{\theta}{a}}\right)\Big]\Big)\\&space;\dot{\theta_i}=\frac{bv_0}{D_c}\left(e^{\frac{\tau_0-\theta}{b}}-\frac{v}{v_0}\right)\\&space;\end{cases}" title="\begin{cases} \dot{u_i}=v_i-v_0\\ \dot{v_i}=(-1/m)\Big(k_c(u_{i+1}+u_{j+1}-4u_i+u_{i-1}+u_{j-1})+k_pu_i+\sigma_{n}a\Big[\sinh^{-1}\left(\frac{V}{2v_0}e^{\frac{\theta}{a}}\right)\Big]\Big)\\ \dot{\theta_i}=\frac{bv_0}{D_c}\left(e^{\frac{\tau_0-\theta}{b}}-\frac{v}{v_0}\right)\\ \end{cases}" /></a>


where m is the block mass, k<sub>p</sub> the block to driving plate spring stiffness, k<sub>c</sub> the inter-block spring stiffness, &nu; the driving plate velocity, &sigma; the applied normal stress, &tau;<sub>0</sub> the static friction force, D<sub>c</sub> the characteristic Length Scale and A and B are friction coefficient of the rate-and-state friction law. For further reading, refer to [Burridge-Knopoff Model](https://pubs.geoscienceworld.org/ssa/bssa/article/57/3/341/116471/model-and-theoretical-seismicity).

N.B. The elastic linear interactions correspond to the terms containing k<sub>c</sub> and k<sub>p</sub> whereas the non-linear part of the ODE is due to the friction law dependent on &nu;, &sigma;, &tau;<sub>0</sub>, D<sub>c</sub>, A and B


### Numerical Implementation

The differential equations are solved using [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl). Solving the differential Equations include 2 steps:

1. Passing the governing equations into a `ODEProblem(f,u0,t,p)`

1. Choosing the Solver Algorithms, tolerances and parameters to solve the system

For the first step, two options have been considered: Either write out the equations for each block and use [ParameterizedFunctions.jl](https://github.com/JuliaDiffEq/ParameterizedFunctions.jl) to pass the equations efficiently (ref. *Model1*) or use a matrix implementation as <a href="https://www.codecogs.com/eqnedit.php?latex=A\times&space;\textbf{x}&space;&plus;&space;f(\textbf{x})&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A\times&space;\textbf{x}&space;&plus;&space;f(\textbf{x})&space;=&space;0" title="A\times \textbf{x} + f(\textbf{x}) = 0" /></a> (ref. *Model2* ). Within each approach, different implementations are suggested to optimize the computations. 

The second step is mainly defined by computational accuracy and has been previously checked. As a conclusion, the `Rodas5()` algorithm has most efficient for larger systems using `reltol=1e-6` and `abstol=1e-6` (and `save_everystep=false`). All functions use in-place allocations.

While the actual problem uses 10x20 blocks, simulates 20000 timesteps and uses `callback` functions to store intermediate results, a benchmark is shown for a 10x10 block system for 100 timesteps. The simulations were run on a local machine.

Model |  CPU Time (avg)| No Alloc | Memory | Setup 
----- | --------- | -------- | ------ | -----
Model 1 | 0.68 sec | 20.24 k | 86.6 MiB | Std
Model 1 | ||| linsolve=LinSolveGMRES()
Model 2 | 74 ms | 54.60 k | 12.75 MiB | Std
Model 2 | 69 ms | 54.99 k | 11.175 MiB | linsolve=LinSolveGMRES()
Model 2 | 74 ms | 52.49 k | 7.63 MiB | Array_Slicing + linsolve
Model 2 | 73.7 ms | 53.73 k | 11.79 MiB | @. + linsolve
Model 2 | 76.3 ms | 51.24 k | 7.67 MiB | @. + Array_Slicing + linsolve
Model 2 | 391.67 ms | 4.7 M | 91.8 MiB | Pre_def
Model 2 | 2.87 s | 2.2 M | 351 MiB | MKL
Model 2 |  | |  | SplitODE

As a note for Model 1, 1 iteration over the linear part takes `0.000079 seconds (4 allocations: 160 bytes)`
 while the friction force takes `0.000094 seconds (46 allocations: 4.625 KiB)`

---

Additional Time required for Callbacks for actual simulations (simulation time = 20000 [s], 10x20 blocks)

Model |  CPU Time (avg)| No Alloc | Memory | Setup 
----- | --------- | -------- | ------ | -----
 Model 1 | 15.816 h | 415.73 M | 588.7 GiB | Std
 Model 1 | 19.271 h | 310.33 M | 580.9 GiB | + Pre-simulation of 10 [sec] 
 Model 2 | 17.638 h | 162.12 M | 28.4 GiB | -no Callbacks
### Further Checks to do
- [ ] Check [SplitODE](http://docs.juliadiffeq.org/latest/types/split_ode_types.html)'s to take advantage of (non-) linear separation and new developments i.e. [Exponential Krylov Integrator](http://juliadiffeq.org/2018/03/31/AdaptiveLowSDE.html)
- [ ] Check [Operators](http://docs.juliadiffeq.org/latest/features/diffeq_operator.html) for <a href="https://www.codecogs.com/eqnedit.php?latex=A\times&space;\textbf{x}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?A\times&space;\textbf{x}" title="A\times \textbf{x}" /></a> handling
- [ ] Check [GPU](http://juliadiffeq.org/2019/05/09/GPU.html) functionality
- [ ] Parallelism




