using LinearAlgebra
using DifferentialEquations
using GLMakie

# Start of script
# closeall() # Close all figures

t0 = time()  # Start timer

"""Initial conditions and constants"""
# Mass of Sun [kg]
Ms = 1.9885e30
# Mass of Earth [kg]  
M1 = 5.9724e24
# Mass of Moon [kg]
M2 = 7.346e22
# Gravitational constant [m^3/(kg*s^2)]
G = 6.673e-11

"""Earth-Moon CR3BP(円制限3体問題)"""
# Characteristic length of Earth-Moon CR3BP [m] (地球-月系の特性長さは地球と月間距離を用いる)
chara_length_CR3BP = 3.8440e8  
# Characteristic mass of Earth-Moon CR3BP [kg] (地球-月系の特性質量は地球と月の合算質量を用いる)
chara_mass_CR3BP = M1 + M2
# Characteristic time of Earth-Moon CR3BP [s] (地球-月系の特性時間は地球と月の公転周期を用いる)
chara_time_CR3BP = sqrt(chara_length_CR3BP^3 / (G * chara_mass_CR3BP))
# Mass ratio of Earth-Moon CR3BP [-]
mu_EM = M2 / (M1 + M2)
# Mean-motion of Earth-Moon CR3BP [1/s]
N_CR3BP = sqrt(G * chara_mass_CR3BP / chara_length_CR3BP^3)

"""Earth-Moon ER3BP(楕円制限3体問題)"""
# Semimajor axis of Earth-Moon ER3BP [m]
a = 3.8440e8
# Mean-motion of Earth-Moon ER3BP [1/s]
N_ER3BP = sqrt(G * (M1 + M2) / (a^3))
# Eccentricity of Earth-Moon ER3BP [-]
e = 0.0549 

"""Earth-Moon-Sun BCR4BP(二重円制限4体問題)"""
# Characteristic length of the Sun-B1 [m] 
chara_length_SB1 = 1.4960e11
# Characteristic mass of the Sun-B1 [kg] 
chara_mass_SB1 = M1 + M2 + Ms 
# Characteristic time of the Sun-B1 [s] 
chara_time_SB1 = sqrt(chara_length_SB1^3 / (G * chara_mass_SB1))
# Nondimensional mass of Sun [s]
mS = Ms / (M1 + M2)
# Nondimensional Sun orbit radius [m]
aS = chara_length_SB1 / chara_length_CR3BP
# Nondimensional Sun angular velocity [1/s]
wS = sqrt((1 + mS) / (aS^3)) - 1
# Initial phase angle of the Sun-B1 [rad]（正しいかのちに確認）
thetaS0 = π 
# Mass ratio of the Sun-B1 [-]
mu_SB1 = (M1 + M2) / (M1 + M2 + Ms)
# Nondimensional distance of Earth-Moon-Sun BCR4BP [-]
aEM = chara_length_CR3BP / chara_length_SB1
# Characteristic time of the Sun-B1 [-]
wM = chara_time_SB1 / chara_time_CR3BP - 1
# Initial phase angle of Earth-Moon [rad]（正しいかのちに確認）
thetaM0 = 0

"""L2 Lyapunov orbit in Earth-Moon CR3BP"""
# Initial Condition （datファイルから読み取るようにする）
t_n_CR3BP = 3.467622949281189
x_n_CR3BP = [1.102866098413080, 0.0, 0.0, 0.0, 0.259155907029058, 0.0]

# Specification Parameter 
parameter = (mu_EM)
# Specification Integration Time 
tspan1 = (0.0, t_n_CR3BP)
# Setup ODE Problem 
prob = ODEProblem(fun_cr3bp!, x_n_CR3BP, tspan1, mu_EM)
# ODE Solver Execution
sol = solve(prob, Vern7(), abstol=1e-14, reltol=1e-14)

# Setup GLMakie Plot 
fig = Figure(size = (800, 600))
# Create an Axis for the Plot
ax = Axis(fig[1, 1], title = "L2 Lyapunov Orbit in Earth-Moon CR3BP", 
    xlabel = "x [-]", ylabel = "y [-]", aspect = DataAspect())

# Plot Elements
scatter!(ax, [1 - mu_EM], [0.0], color = "#EDB120", markersize = 20, label = "Moon")
scatter!(ax, [sol[1, 1]], [sol[2, 1]], color = "#0072BD", markersize = 10, label = "Initial Position")
lines!(ax, sol[1, :], sol[2, :], color = "#0072BD", linewidth = 2, label = "Orbit")

# Add a Legend
legend = Legend(fig, ax, "Legend", orientation = :vertical)
fig[1, 2] = legend

# Display the Figure
fig