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
parameter1 = (mu_EM)
# Specification Integration Time 
tspan1 = (0.0, t_n_CR3BP)
# Setup ODE Problem 
prob = ODEProblem(fun_cr3bp!, x_n_CR3BP, tspan1, parameter1)
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

"""trajectory in Earth-Moon ER3BP"""
# Orbital Period of ER3BP
t_ER3BP = t_n_CR3BP * chara_time_CR3BP * N_ER3BP
# transfomation from CR3BP reference frame to the inertial frame
X0_CR3BP = [x_n_CR3BP[1:3] .* chara_length_CR3BP; x_n_CR3BP[4:6] .* chara_length_CR3BP ./ chara_time_CR3BP]
C_CR3BP = [cos(0.0) -sin(0.0) 0.0;
           sin(0.0)  cos(0.0) 0.0;
                0.0       0.0 1.0]
dtheta_dt_CR3BP = N_CR3BP
# constructs a transformation matrix of CR3BP
rotating_matrix_CR3BP = fun_rotating_to_inertial_matrix(C_CR3BP, dtheta_dt_CR3BP)
# Convert the initial state vector X0_CR3BP expressed in the rotating frame of CR3BP to the inertial frame
X0_inertial = rotating_matrix_CR3BP * X0_CR3BP

# Transformation from inertial frame to ER3BP
C_ER3BP = [cos(0.0) -sin(0.0) 0.0;
           sin(0.0)  cos(0.0) 0.0;
                0.0       0.0 1.0]
# Calculation of angular velocity of ER3BP
dtheta_dt_ER3BP = sqrt(G * (M1 + M2) * (1 + e * cos(0.0))^4 / (a * (1 - e^2))^3)
# constructs a transformation matrix of ER3BP
rotating_matrix_ER3BP = fun_rotating_to_inertial_matrix(C_ER3BP, dtheta_dt_ER3BP)
# Convert the initial state vector X0_inertial expressed in the inertial frame to the rotating frame of ER3BP
X0_ER3BP = rotating_matrix_ER3BP \ X0_inertial # (inv(rotating_matrix_ER3BP)*X0_inertialと同じ意味)
x0_ER3BP = [X0_ER3BP[1:3] ./ (a * (1 - e)); X0_ER3BP[4:6] ./ (a * (1 - e)) ./ N_ER3BP] # 同じものを二回定義しているけど、後者を採用？

# Specification Parameter 
parameter2 = (mu_EM, e)
# Solve ODE in ER3BP
tspan2 = (0.0, t_ER3BP)
# Setup ODE Problem 
prob2 = ODEProblem(fun_er3bp!, x0_ER3BP, tspan2, parameter2)
# ODE Solver Execution
sol2 = solve(prob2, Vern7(), abstol=1e-14, reltol=1e-14)

# Visualization using GLMakie
fig2 = Figure(size = (800, 600))
ax2 = Axis(fig2[1, 1], title = "Trajectory in Earth-Moon ER3BP",
    xlabel = "x [-]", ylabel = "y [-]", aspect = DataAspect())

scatter!(ax2, [1 - mu_EM], [0.0], color = "#EDB120", markersize = 20, label = "Moon")
scatter!(ax2, [sol2[1, 1]], [sol2[2, 1]], color = "#0072BD", markersize = 10, label = "Initial Position")
lines!(ax2, sol2[1, :], sol2[2, :], color = "#0072BD", linewidth = 2, label = "Trajectory")

legend2 = Legend(fig2, ax2, "Legend", orientation = :vertical)
fig2[1, 2] = legend2

# Save figure
f2_name = "Ex1_ER3BP_x0=$(x0_ER3BP[1])_vy0=$(x0_ER3BP[5])_t=$(t_ER3BP)"
f2_name = replace(f2_name, "." => ",")
save("$(f2_name).png", fig2)

fig2