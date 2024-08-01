include("function_bundle.jl")

# * Parameter for simulated_annealing
N = 30
  # System with N x N spins
T_initial = 10.0
  # Initial temperature, kT units
T_final = 0.1
  # Final temperature, kT units
cooling_rate = 0.99
  # Next temperature = cooling_rate * current temperature
steps_per_temp = 5000
  # Number of Metropolis steps per temperature
perturbation_rate = 0.1
  # Intensity of random perturbation to spin state

# Simulated Annealing is done here.
spins = initialize_spins(N);
final_spins, final_energy = 
  simulated_annealing(spins, T_initial, T_final, cooling_rate, perturbation_rate, steps_per_temp)
# One more annealing at the final temperature
final_spins, final_energy = 
  simulated_annealing(final_spins, T_final, 0.01 * T_final, cooling_rate, perturbation_rate, steps_per_temp)


# * Energy and Spin Configuration
println("Final Spin Configuration: \n", final_spins)
println("Final Energy: ", final_energy)

# * With a Fourier Transform, we can extract the power spectrum of the spin configuration.
N = size(final_spins, 1)

Snx = [final_spins[i, j][1] for i in 1:N, j in 1:N];  Skx = fft(Snx);
Sny = [final_spins[i, j][2] for i in 1:N, j in 1:N];  Sky = fft(Sny);
Snz = [final_spins[i, j][3] for i in 1:N, j in 1:N];  Skz = fft(Snz);
Sk_Sk = real.(dot.(Skx, Skx) + dot.(Sky, Sky) + dot.(Skz, Skz))

recipBasis = [[√3/2, 1/2], [ 0.0, 1.0]];
x = [recipBasis[1][1]*i + recipBasis[2][1]*j for i in (-N+1):N, j in (-N+1):N];  x = x[:]
y = [recipBasis[1][2]*i + recipBasis[2][2]*j for i in (-N+1):N, j in (-N+1):N];  y = y[:]
z = Sk_Sk;  z = [z z; z z];                                                      z = z[:]

fig  = Figure();
axs = Axis(fig[1, 1], aspect = 1)
tr = tricontourf!(axs, x, y, z, colormap = :turbo)
xlims!(axs, -N*1/2, N*3/4);  ylims!(axs, -N*1/2, N*3/4);
fig

# * Visualize Spin configuration in a 3D plot with GLMakie

include("plot_spins_custom.jl")

fig = Figure(size = (1200, 400));  xs, ys, zs = define_lattice_point(final_spins);

ax1 = Axis3(fig[1, 1], azimuth = 3π/2, elevation = π/2, aspect = (1,√3/2,1));
hidexdecorations!(ax1);  hideydecorations!(ax1);  hidezdecorations!(ax1);
ff1 = scatter!(ax1, xs, ys, zs, markersize = 5, color = :grey)
gg1 = plot_spins_custom(ax1,final_spins,xs,ys);
xlims!(ax1,-5,10);  ylims!(ax1,0,6√3);

ax2 = Axis3(fig[1, 2]; aspect = (1,1,1), azimuth = π/4, elevation = π/12);
ff2 = draw_unit_sphere(ax2);  gg2 = reduce_on_unit_sphere(ax2, final_spins; frame=0);
hidexdecorations!(ax2);  hideydecorations!(ax2);  hidezdecorations!(ax2);  hidespines!(ax2);

fig