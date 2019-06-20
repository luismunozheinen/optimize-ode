
Download the git into a local directory which is a recognized working directory by your Julia setup, i.e. push!(LOAD_PATH,"your_path/optimize-ode-master")

Files `setup_*` are modules, no need to run them! They will be imported into `optimize_ode`

`optimize_ode` is the main file to run. Make sure you first install all required packages listed in the first line
Subsequent code sets up the differential equation to solve and therefore:

 1. Loads the Initial Conditions
 1. Size of the Model/System
 1. Model Parameter
 1. 2 different Implementations for the Friction Forces (V1 and V2)
 1. Checks (`@time` and `@benchmark`) to verify 0 allocations for the friction functions (in-place) (OPTIONAL)  
 1. 2 different Implementations for the Setup of the full diff. eq. System (ODEProblem and SplitODEProblem)
 1. `@benchmark` sol_i To solve and benchmark the cpu cost for a given implementation + solver setup
 1. Options for Post-processing Results, i.e. extract motions and plot results
 
