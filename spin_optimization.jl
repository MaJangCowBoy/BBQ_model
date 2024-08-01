using FFTW, LinearAlgebra, Random, ProgressBars

# 스핀 배열을 초기화하는 함수
function initialize_spins(N)
    spins = [normalize(randn(3)) for _ in 1:N, _ in 1:N]
    return reshape(spins, N, N)
end

# Fourier Transform을 사용하여 에너지를 계산하는 함수
function energy_fft(spins)
    N = size(spins, 1)
    
    # 각 스핀 성분을 분리
    Snx = [spins[i, j][1] for i in 1:N, j in 1:N]
    Sny = [spins[i, j][2] for i in 1:N, j in 1:N]
    Snz = [spins[i, j][3] for i in 1:N, j in 1:N]
    
    # 2D Fourier Transform 수행
    Skx = fft(Snx)
    Sky = fft(Sny)
    Skz = fft(Snz)
    
    # k1 = 1/2 a*, k2 = 1/2 b*, k3 = 1/2 (a* + b*) 성분 크기 추출
    idx = Int64(N/2+1)
    if N%2 != 0  error("N must be even")
    else
        Sk1_x = Skx[idx,1];   Sk1_y = Sky[idx,1];   Sk1_z = Skz[idx,1]
        Sk2_x = Skx[1,idx];   Sk2_y = Sky[1,idx];   Sk2_z = Skz[1,idx]
        Sk3_x = Skx[idx,idx]; Sk3_y = Sky[idx,idx]; Sk3_z = Skz[idx,idx]
    end

    Sk1 = [Sk1_x, Sk1_y, Sk1_z]
    Sk2 = [Sk2_x, Sk2_y, Sk2_z]
    Sk3 = [Sk3_x, Sk3_y, Sk3_z]
    Sks = [Sk1, Sk2, Sk3]

    # Fourier transform 된 성분으로부터 에너지 계산
    E = 0;
    J = -1.0;
    B = +0.10;
    Ntot = prod(size(spins));
    Norm = Ntot * Ntot;
    for Ski in Sks  
      E += J * dot(Ski, Ski)  
      E += B/Norm * (dot(Ski, Ski) * dot(Ski, Ski))
    end

    return E
end

# Metropolis 알고리즘을 사용하는 함수
function metropolis(spins, T, λ)
    N = size(spins, 1)
    i, j = rand(1:N), rand(1:N)
    # 기존 스핀에 작은 섭동을 더하여 새로운 스핀 제안
    old_spin = spins[i, j]
    new_spin = normalize(spins[i, j] + λ * randn(3))
    
    # 에너지 변화 계산
    old_energy = energy_fft(spins)
    spins[i, j] = new_spin
    new_energy = energy_fft(spins)
    ΔE = new_energy - old_energy
    
    # Metropolis 조건에 따라 새로운 스핀 상태 수락
    if real(ΔE) < 0 || rand() < real(exp(-ΔE / T))
        return spins
    else
        spins[i, j] = old_spin # 만약 거부하면 이전 상태로 되돌림
        return spins
    end
end

# Simulated Annealing을 사용하는 함수
function simulated_annealing(N, T_initial, T_final, cooling_rate, perturbation_rate, steps_per_temp)
    spins = initialize_spins(N)
    T = T_initial
    λ = perturbation_rate
    
    while T > T_final
        for _ in ProgressBar(1:steps_per_temp)
            spins = metropolis(spins, T, λ)
        end
        T *= cooling_rate
        # λ *= cooling_rate
    end
    
    return spins, energy_fft(spins)
end