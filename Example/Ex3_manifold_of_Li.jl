"""
Ex3_manifold_of_Li.jl
"""

using DifferentialEquations
using LinearAlgebra
using GLMakie
using Polynomials

# Function loading
functions_dir = normpath(joinpath(current_dir, ".", "Functions"))

if isdir(functions_dir)
    for file in readdir(functions_dir)
        if endswith(file, ".jl")
            full_path = joinpath(functions_dir, file)
            if isfile(full_path)
                include(full_path)
            else
                println("Warning: File not found - ", full_path)
            end
        end
    end
else
    println("Error: Functions directory not found at ", functions_dir)
end

# Compute stable and unstable manifolds
function compute_manifolds(L, mu, tf)
    sigma = (1 - mu) / abs(L[1] + mu)^3 + mu / abs(L[1] - 1 + mu)^3
    A = [0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1;
         2 * sigma + 1 0 0 0 2 0;
         0 1 - sigma 0 -2 0 0;
         0 0 -sigma 0 0 0]

    # Obtain eigenvalues and eigenvectors (Note the order: In MATLAB, `eig(A)` returns `[V, D] = eig(A)`. But in Julia, `D, V = eigen(A)`, where `D` is output as an array.)
    D, V = eigen(A)

    # 
    Vs = zeros(6) 
    Vu = zeros(6) 

    for i in 1:6
        # eigenvalues
        λ = D[i]
        # Extract only real eigenvalues
        if imag(λ) == 0 && real(λ) != 0
            if real(λ) < 0
                # Add eigenvectors corresponding to negative eigenvalues (Convert ComplexF64 → Float64)
                Vs = real(V[:, i])  
            elseif real(λ) > 0
                # Add eigenvectors corresponding to positive eigenvalues (Convert ComplexF64 → Float64)
                Vu = real(V[:, i])  
            end
        end
    end

    #
    x0_sp = vcat(L, zeros(3)) + 1e-10 * Vs / norm(Vs)
    x0_sm = vcat(L, zeros(3)) - 1e-10 * Vs / norm(Vs)
    x0_up = vcat(L, zeros(3)) + 1e-10 * Vu / norm(Vu)
    x0_um = vcat(L, zeros(3)) - 1e-10 * Vu / norm(Vu)

    # Specification integration time 
    tspan_s = (tf, 0)
    tspan_u = (0, tf)

    # Definition of ODEProblem
    prob_sp = ODEProblem(fun_cr3bp!, x0_sp, tspan_s, mu)
    prob_sm = ODEProblem(fun_cr3bp!, x0_sm, tspan_s, mu)
    prob_up = ODEProblem(fun_cr3bp!, x0_up, tspan_u, mu)
    prob_um = ODEProblem(fun_cr3bp!, x0_um, tspan_u, mu)

    # Execution of the ODE Solver
    sol_sp = solve(prob_sp, Vern7(), abstol=1e-14, reltol=1e-14)
    sol_sm = solve(prob_sm, Vern7(), abstol=1e-14, reltol=1e-14)
    sol_up = solve(prob_up, Vern7(), abstol=1e-14, reltol=1e-14)
    sol_um = solve(prob_um, Vern7(), abstol=1e-14, reltol=1e-14)

    return sol_sp, sol_sm, sol_up, sol_um
end

"""Zero_velocity_curve"""
# Retrieving Parameters for the Earth-Moon Circular Restricted Three-Body Problem
mu, a_1, a_s, w_1 = fun_cr3bp_parameter(2)
# Calculation of Lagrange points
L1, L2, L3, L4, L5 = fun_libration_points(mu)

# Mesh generation
x_range = -1.5:0.001:1.5
y_range = -1.5:0.001:1.5
x = repeat(collect(x_range), 1, length(y_range)) 
y = repeat(collect(y_range)', length(x_range), 1)  
z = zeros(size(x))  

# Potential function calculation
r1 = sqrt.((x .+ mu) .^ 2 .+ y .^ 2 .+ z .^ 2)
r2 = sqrt.((x .- (1 - mu)) .^ 2 .+ y .^ 2 .+ z .^ 2)
U = 0.5 .* (x .^ 2 .+ y .^ 2) .+ (1 - mu) ./ r1 .+ mu ./ r2
C = 2 .* U

"""manifolds"""
tf1, tf2, tf3 = 15.0, 25.0, 145.0
sol_sp_L1, sol_sm_L1, sol_up_L1, sol_um_L1 = compute_manifolds(L1, mu, tf1)
sol_sp_L2, sol_sm_L2, sol_up_L2, sol_um_L2 = compute_manifolds(L2, mu, tf2)
sol_sp_L3, sol_sm_L3, sol_up_L3, sol_um_L3 = compute_manifolds(L3, mu, tf3)

# Lagrange Points and Their Corresponding Jacobi Constants
C1 = fun_Jacobi_const(vcat(L1, [0.0, 0.0, 0.0]), mu)
C2 = fun_Jacobi_const(vcat(L2, [0.0, 0.0, 0.0]), mu)
C3 = fun_Jacobi_const(vcat(L3, [0.0, 0.0, 0.0]), mu)

# Save directory
save_dir = "./Example/Figure/"

# Generate each figure
for (Li, sol_sp, sol_sm, sol_up, sol_um, Ci, tf, label) in 
    zip([L1, L2, L3], [sol_sp_L1, sol_sp_L2, sol_sp_L3], [sol_sm_L1, sol_sm_L2, sol_sm_L3], 
        [sol_up_L1, sol_up_L2, sol_up_L3], [sol_um_L1, sol_um_L2, sol_um_L3], 
        [C1, C2, C3], [tf1, tf2, tf3], ["L1", "L2", "L3"])

    local fig = Figure(size=(800, 600))
    local ax = Axis(fig[1, 1], aspect = DataAspect(), xlabel = L"x\mathrm{[-]}", ylabel = L"y\mathrm{[-]}",
        xlabelsize = 16, ylabelsize = 16,
        xticklabelsize = 16, yticklabelsize = 16)

    p_primary = scatter!(ax, [-mu, 1 - mu], [0, 0], color = :black, markersize = 15)

    p_stable1 = lines!(ax, real(sol_sp[1, :]), real(sol_sp[2, :]), color = :blue, linewidth = 1.5)
    p_stable2 = lines!(ax, real(sol_sm[1, :]), real(sol_sm[2, :]), color = :blue, linewidth = 1.5)

    p_unstable1 = lines!(ax, real(sol_up[1, :]), real(sol_up[2, :]), color = :red, linewidth = 1.5)
    p_unstable2 = lines!(ax, real(sol_um[1, :]), real(sol_um[2, :]), color = :red, linewidth = 1.5)

    p_C = contour!(ax, x_range, y_range, C; levels = [Ci], linewidth = 1.5, color = :black)

    p_lagrange = scatter!(ax, [L1[1], L2[1], L3[1], L4[1], L5[1]],
                              [L1[2], L2[2], L3[2], L4[2], L5[2]],
                              color = :black, marker = :star6, markersize = 15)
    p_Li = scatter!(ax, [Li[1]], [Li[2]], color = :orange, marker = :star6, markersize = 15)

    fig[1, 2] = Legend(fig, [
            p_primary,
            p_stable1,
            p_unstable1,
            p_C,
            p_lagrange,
            p_Li
        ], [
            "Primary bodies",
            "Stable manifolds",
            "Unstable manifolds",
            "Zero-velocity curve",
            "Lagrange points",
            "Target L-point"
        ])

    # Save the figure
    local save_file_name = "Ex3_manifold_of_$(label)_mu=$(mu)_C=$(Ci)_tf=$(tf)"
    local save_file_name = replace(save_file_name, "." => ",")
    local save_path = joinpath(save_dir, save_file_name * ".png")
    save(save_path, fig)

    # Save trajectory data
    fun_save_data(sol_sp.t, sol_sp.u, save_file_name * "_stable_plus")
    fun_save_data(sol_sm.t, sol_sm.u, save_file_name * "_stable_minus")

    display(fig)
end

