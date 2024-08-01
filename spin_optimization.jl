using FFTW, LinearAlgebra, Random, ProgressBars

# Initialize spins, randomly distributed on the unit sphere
function initialize_spins(N)
    spins = [normalize(randn(3)) for _ in 1:N, _ in 1:N]
    return reshape(spins, N, N)
end

# BBQ model to calculate the energy of the spin configuration
function energy_fft(spins)
    N = size(spins, 1)
    
    # Fourier transform the spin structure.
    Snx = [spins[i, j][1] for i in 1:N, j in 1:N];    Skx = fft(Snx)
    Sny = [spins[i, j][2] for i in 1:N, j in 1:N];    Sky = fft(Sny)
    Snz = [spins[i, j][3] for i in 1:N, j in 1:N];    Skz = fft(Snz)
        
    # Extract the components of the Fourier transformed spins
    # symmetrically equivalent components should be all considered
    # pQ = (    N/3);  idx1 = Int64(pQ+1);
    # mQ = (N - N/3);  idx2 = Int64(mQ+1);
    pQ = (N * 10/30);  idx1 = Int64(pQ+1);
    mQ = (N * 20/30);  idx2 = Int64(mQ+1);
    

    if pQ%1 != 0  
        error("Ordering vector implementation is not correct")
    else
        Sk1p_x = Skx[idx1,   1];  Sk1p_y = Sky[idx1,   1];  Sk1p_z = Skz[idx1,   1];
        Sk1m_x = Skx[idx2,   1];  Sk1m_y = Sky[idx2,   1];  Sk1m_z = Skz[idx2,   1];
        Sk2p_x = Skx[   1,idx1];  Sk2p_y = Sky[   1,idx1];  Sk2p_z = Skz[   1,idx1];
        Sk2m_x = Skx[   1,idx2];  Sk2m_y = Sky[   1,idx2];  Sk2m_z = Skz[   1,idx2];
        Sk3p_x = Skx[idx1,idx2];  Sk3p_y = Sky[idx1,idx2];  Sk3p_z = Skz[idx1,idx2];
        Sk3m_x = Skx[idx2,idx1];  Sk3m_y = Sky[idx2,idx1];  Sk3m_z = Skz[idx2,idx1];
    end

    Sk1p = [Sk1p_x, Sk1p_y, Sk1p_z];    Sk1m = [Sk1m_x, Sk1m_y, Sk1m_z];
    Sk2p = [Sk2p_x, Sk2p_y, Sk2p_z];    Sk2m = [Sk2m_x, Sk2m_y, Sk2m_z];
    Sk3p = [Sk3p_x, Sk3p_y, Sk3p_z];    Sk3m = [Sk3m_x, Sk3m_y, Sk3m_z];

    Skp = [Sk1p, Sk2p, Sk3p];    Skm = [Sk1m, Sk2m, Sk3m];

    # Calculate the energy from the Fourier transformed spins
    E = 0;
    J = -1.00;
    B = +0.05;
    Ntot = prod(size(spins));
    Norm = Ntot * Ntot;
    for (S_plus_q,S_minus_q) in zip(Skp, Skm)  
      E += J * sum(S_plus_q .* S_minus_q)
      E += B/Norm * sum(S_plus_q .* S_minus_q) * sum(S_plus_q .* S_minus_q)
    end

    return E
end

# Metropolis algorithm to update the spin configuration
function metropolis(spins, T, λ)
    N = size(spins, 1)

    # Randomly select a spin to perturb
    i, j = rand(1:N), rand(1:N)
    old_spin = spins[i, j]
    new_spin = normalize(spins[i, j] + λ * randn(3))
    
    # Calculate the energy difference
    old_energy = energy_fft(spins)
    spins[i, j] = new_spin
    new_energy = energy_fft(spins)
    ΔE = new_energy - old_energy
    
    # Following Metropolis condition, 
    #   accept or reject the new spin configuration
    if real(ΔE) < 0 || rand() < real(exp(-ΔE / T))
        return spins
    else
        spins[i, j] = old_spin
        return spins
    end
end

# Simulated Annealing algorithm
function simulated_annealing(spins, T_initial, T_final, cooling_rate, perturbation_rate, steps_per_temp)
    T = T_initial
    λ = perturbation_rate
    
    while T > T_final
        for _ in ProgressBar(1:steps_per_temp)
            spins = metropolis(spins, T, λ)
        end
        T *= cooling_rate
    end
    
    return spins, energy_fft(spins)
end