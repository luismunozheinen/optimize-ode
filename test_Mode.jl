using DifferentialEquations,ParameterizedFunctions,JLD2,LinearAlgebra,StaticArrays,BenchmarkTools,DiffEqOperators
#Load Function
#include("create_Amatrix.jl")
using create_Amatrix
#Initial Condition
@load "ic.jld2" u0
#Time domain
tspan = (0.0,100.0)
#import array containing Parameters
@load "parameter.jld2" p
#Load A matrix for linear interactions
const A = create_Amatrix.create_matrix(p)
#Various Implementations for the Friction Law
# 1.Standard Out-of-place-------------------------------------------------------
function friction(u,p)
    f=zeros(3*l*d)
    f[1:3:end] .= -p[4]
    f[2:3:end] .= p[5] .*p[8] .* asinh.(u[2:3:end] ./(2 .*p[4]) .* exp.(u[3:3:end] ./p[8]))
    f[3:3:end] .= p[9] .*p[4]/p[7] .*(exp.((p[6] .-u[3:3:end]) ./p[9]) .- u[2:3:end] ./p[4])
    return f
end

# 2.Standard In-place Version---------------------------------------------------
function friction!(f,u,p,t)
     f[1:3:end] .= -p[4]
     f[2:3:end] .= -1/p[1]*(p[5] .*p[8] .* asinh.(u[2:3:end] ./(2 .*p[4]) .* exp.(u[3:3:end] ./p[8])))
     f[3:3:end] .= p[9].*p[4]/p[7] .*(exp.((p[6] .-u[3:3:end]) ./p[9]) .- u[2:3:end] ./p[4])
    return f
end

# 3.Array Slicing Version-------------------------------------------------------
function friction_ArraySlicing!(f,u,p,t)
   pos = @view u[1:3:end]
   vel = @view u[2:3:end]
   state = @view u[3:3:end]
     f[1:3:end] .= -p[4]
     f[2:3:end] .= p[5] .*p[8] .* asinh.(vel ./(2 .*p[4]) .* exp.(state ./p[8]))
     f[3:3:end] .= p[9].*p[4]/p[7] .*(exp.((p[6] .-state) ./p[9]) .- vel ./p[4])
    return f
end

# 4. @. Version-----------------------------------------------------------------
function friction_noalloc!(f,u,p,t)
     f[1:3:end] .= @. -p[4]
     f[2:3:end] .= @. p[5]*p[8] .* asinh.(u[2:3:end] ./(2 .*p[4]) .* exp.(u[3:3:end] ./p[8]))
     f[3:3:end] .= @. p[9].*p[4]/p[7] .*(exp.((p[6] .-u[3:3:end]) ./p[9]) .- u[2:3:end] ./p[4])
     return f
end

# 5. @.+slicing Version-----------------------------------------------------------------
function friction_noalloc_slicing!(f,u,p,t)
    pos = @view u[1:3:end]
    vel = @view u[2:3:end]
    state = @view u[3:3:end]
    f[1:3:end] .= @. -p[4]
    f[2:3:end] .= @. p[5] .*p[8] .* asinh.(vel ./(2 .*p[4]) .* exp.(state ./p[8]))
    f[3:3:end] .= @. p[9].*p[4]/p[7] .*(exp.((p[6] .-state) ./p[9]) .- vel ./p[4])
    return f
end

# 6.Predefining u,v and θ Vectors-----------------------------------------------
const array_u = [(i+2) % 3 == 0 ? 1.0 : 0  for i=1:d*l*3]
const array_v = [(i+1) % 3 == 0 ? 1.0 : 0  for i=1:d*l*3]
const array_θ = [i % 3 == 0 ? 1 : 0.0  for i=1:d*l*3]

function friction_predef!(du,u,p,t)
    du = @. -p[4] .*array_u + p[5]*p[8] .* asinh.(u ./(2*p[4]) .* exp.(u ./p[8])) .*array_v +  (p[9]*p[4]/p[7] .*(exp.((p[6] .-u) ./p[9]) .- u ./p[4])) .*array_θ
end
#SOLVING SYSTEMS----------------------------------------------------------------
# 1. Static Array Version-------------------------------------------------------
u0S = @SArray [u0[i] for i=1:length(u0)]
A = SMatrix{300,300}(A1)

# 2. Split ODE------------------------------------------------------------------
f0 = create_Amatrix.friction(u0,p)  # allocate non-linear content into f0 (Initialization)
f1(du,u,p,t) = mul!(du,A,u) #linear part of the ode function
f2(du,u,p,t) = DiffEqArrayOperator(A) #linear part using DiffEqOperators
f2(du,u,p,t) = friction!(f0,u,p,t) #non-linear part of the ode function
tspan = (0.0,100.0)  #timespan over which to solve ode
prob = SplitODEProblem(f1,f2,u0,tspan,p) #set up ode problem to solve using initial conditions and timespan
sol = solve(prob,solver=Rodas5(autodiff=false),reltol=1e-6,abstol=1e-6,force_dtmin=true,maxiters=10^14,progress=true,save_everystep=false) #solve problem
@time f1(u0,u0,p,tspan); #check cpu requirements for single function call
#@time f2(u0,u0,p,tspan); #check cpu requirements for single function call

# 3. STANDARD-------------------------------------------------------------------
mf(du,u,p,t) = mul!(du,A,u) .+ friction!(du,u,p,t) #setup for ode function
tspan = (0.0,100.0) #timespan over which to solve ode
prob = ODEProblem(mf,u0,tspan,p)  #set up ode problem to solve using initial conditions and timespan
#@benchmark
sol = solve(prob,Rodas5(autodiff=false),reltol=1e-6,abstol=1e-6,save_everystep=true) #solve problem and benchmark cpu requirements
#@benchmark sol = solve(prob,Rodas5(autodiff=false,linsolve=LinSolveGMRES()),reltol=1e-6,abstol=1e-6,force_dtmin=true,maxiters=10^14,save_everystep=false)
@time mf(u0,u0,p,tspan);

#-------------------------------------------------------------------------------
#Plot Results to compare output-------------------------------------------------
using Plots
gr()
plot!(sol,label="")
