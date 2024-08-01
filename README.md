# BBQ_model
custom BBQ model, bilinear-biquadratic model, for simulated annealing in momentum space

# How to use this code.

```
spins = initialize_spins(N);
```
N corresponds to system dimension.

```
final_spins, final_energy = 
  simulated_annealing(spins, T_initial, T_final, cooling_rate, perturbation_rate, steps_per_temp)
```
simulated annealing can make structure optimization.
