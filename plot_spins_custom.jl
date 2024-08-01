using GLMakie

function define_lattice_point(spins)

    lat = size(spins)
  
    bR = [ 1.0 -0.5;
           0.0 √3/2 ];
  
    xs = [];  ys = [];  zs = [];
  
    for i = 1:lat[1], j = 1:lat[2]
      push!(xs, bR[1,1]*i + bR[1,2]*j)
      push!(ys, bR[2,1]*i + bR[2,2]*j)
      push!(zs, 0.0)
    end
    
    xs = float.(xs);  ys = float.(ys);  zs = float.(zs);
  
    return xs, ys, zs;
  end
  
  
  function plot_spins_custom(ax,spins,xs,ys)
    lat = size(spins);  us = [];  vs = [];  ws = [];
    for i = 1:lat[1], j = 1:lat[2]
      push!(us, spins[i,j][1])
      push!(vs, spins[i,j][2])
      push!(ws, spins[i,j][3])
    end
    us = 0.15 * float.(us);  vs = 0.15 * float.(vs);  ws = 0.15 * float.(ws);
    fig = arrows!(ax, xs - 0.5 * us, ys - 0.5 * vs, zs - 0.5 * ws, us, vs, ws, color=ws)
  
    return fig;
  end
    
  function draw_unit_sphere(ax)
  
    θList = range(0,  π, 40);  φList = range(0, 2π, 40);
  
    eqx = [sin(θ)*cos(φ) for φ in φList, θ in [π/2]];
    eqy = [sin(θ)*sin(φ) for φ in φList, θ in [π/2]];
    eqz = [cos(θ)        for φ in φList, θ in [π/2]];
    gg = lines!(ax, eqx, eqy, eqz, color=:black, linewidth=1.5);
  
    for φ in φList
      xs = [sin(θ)*cos(φ) for θ in θList];
      ys = [sin(θ)*sin(φ) for θ in θList];
      zs = [cos(θ)        for θ in θList];
      hh = lines!(ax, xs, ys, zs, color=:grey, alpha = 0.75, transparency=true);
    end
  
    for θ in θList
      xs = [sin(θ)*cos(φ) for φ in φList];
      ys = [sin(θ)*sin(φ) for φ in φList];
      zs = [cos(θ)        for φ in φList];
      ii = lines!(ax, xs, ys, zs, color=:grey, alpha = 0.75, transparency=true);
    end
  
    return gg;
  end
    
  function reduce_on_unit_sphere(ax, spins; frame = 0)
  
    lat = size(spins);  us = [];  vs = [];  ws = [];
    for i = 1:lat[1], j = 1:lat[2]
      vec = [spins[i,j][1], spins[i,j][2], spins[i,j][3]];
      vec = vec/norm(vec);  u = vec[1];  v = vec[2];  w = vec[3];
      fig = scatter!(ax, u, v, w, markersize = 10, color = "green")
    end
  
    return fig;
  end