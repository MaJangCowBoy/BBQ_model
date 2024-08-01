include("spin_optimization.jl")

#* 파라미터 설정
N = 10  # 스핀 격자의 크기
T_initial = 10.0
T_final = 0.1
cooling_rate = 0.99
steps_per_temp = 1000
perturbation_rate = 0.1

# Simulated Annealing 수행
final_spins, final_energy = simulated_annealing(N, T_initial, T_final, cooling_rate, perturbation_rate, steps_per_temp)

# 최종 스핀 배열과 에너지 출력
println("Final Spin Configuration:")
println(final_spins)
println("Final Energy: ", final_energy)

# Fourier transform 을 통해, M1, M2, M3 성분이 얼마나 많이 존재하는지 확인
N = size(final_spins, 1)
Snx = [final_spins[i, j][1] for i in 1:N, j in 1:N]
Sny = [final_spins[i, j][2] for i in 1:N, j in 1:N]
Snz = [final_spins[i, j][3] for i in 1:N, j in 1:N]
Skx = fft(Snx)
Sky = fft(Sny)
Skz = fft(Snz)

argmax(abs.(Skx)), argmax(abs.(Sky)), argmax(abs.(Skz))

# GLMakie 패키지를 사용하여 스핀 배열을 시각화

include("plot_spins_custom.jl")

fig = Figure(size = (1200, 400));  xs, ys, zs = define_lattice_point(final_spins);

ax1 = Axis3(fig[1, 1], azimuth = 3π/2, elevation = π/2, aspect = (1,√3/2,1));
hidexdecorations!(ax1);  hideydecorations!(ax1);  hidezdecorations!(ax1);
ff1 = scatter!(ax1, xs, ys, zs, markersize = 5, color = :grey)
gg1 = plot_spins_custom(ax1,final_spins,xs,ys);
xlims!(ax1,-5,10);  ylims!(ax1,0,6√3);

ax2 = Axis3(fig[1, 2]; aspect = (1,1,1), azimuth = π/4, elevation = π/12);
ff2 = draw_unit_sphere(ax2);  gg2 = reduce_on_unit_sphere(ax2, final_spins; frame=0);
# hidexdecorations!(ax2);  hideydecorations!(ax2);  hidezdecorations!(ax2);  hidespines!(ax2);

fig