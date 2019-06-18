using DifferentialEquations,ParameterizedFunctions,JLD2,BenchmarkTools
# Load functions to setup the ODE system using ParameterizedFunctions.jl--------
#include("create_setup.jl")
#using Main.create_setup #
using create_setup

# Load ODE Problem Parameters---------------------------------------------------
@load "ic.jld2" u0  #initial Condition

tspan = (0.0,0.1) #time domain
@load "parameter.jld2" p #Load System Parameters
ν = 0.0001; m = 20;kp0 = 0.5; kc0 = 0.1; τ0 = 0.01; a0 = 0.012; b0 = 0.018; L = 0.005; σ0 = 1.0
p = (m,kp0,kc0,ν,σ0,τ0,L,a0,b0)
#Create and Solve ODE-----------------------------------------------------------
prob = ODEProblem(f,u0,tspan,p) #set up ode problem to solve using initial conditions and timespan
@time sol = solve(prob,Rodas5(autodiff=false),reltol=1e-6,abstol=1e-6,force_dtmin=true,maxiters=10^14,progress=true) #Solver
