
tmod = now()
using Memento
using Plots
Plots.gr()
using ProgressMeter

dmodules = now() - tmod

const ρ = 0.4 # spins down ratio
const N = 120 # Lattice size
const N2 = N^2
const directions = CartesianIndex.([(0,1), (0,-1), (1,0), (-1,0)])

const SEED = 0x5e6c1f2e

const eq_steps = 10^3
const mc_steps = 9 * 10^3

function main()
  logger = Memento.config("debug"; fmt="[{level} | {name}]: {msg}")

  info(logger, "Seed: $SEED")
  prng = MersenneTwister(SEED)

  const tT = linspace(0.1, 6, 2)
  const outdir = "simulations/kawasaki_$SEED"
  mkpath(outdir)

  
  ω = init_lattice(prng; is_local=true)
  p0 = make_heatmap(ω)
  plot!(p0, title="Initial lattice (t=0)")
  savefig(p0, "$outdir/image0.svg")

  tE = similar(tT)
  tE2 = similar(tT)
  tc = similar(tT)

  @showprogress for (i_T,T) in enumerate(tT)
    subdir = "$outdir/$T"
    mkpath(subdir)

    totE = 0
    totE2 = 0
    t = 0
    for k in 1:eq_steps
      sweep(ω, T, prng)
      t = t + 1
    end

    for k in 1:mc_steps
      p = make_heatmap(ω)
      plot!(p, title="T=$T, t=$t")
      savefig(p, "$subdir/image$(dec(t, 9)).svg")
      sweep(ω, T, prng)
      t = t + 1

      E = H(ω)
      totE = totE + E
      totE2 = totE2 + E*E
    end

    mE = 1 / mc_steps * totE
    mE2 = 1 / mc_steps * totE2
    tE[iT] = mE
    tE2[iT] = mE2
    tc[iT] = k * β^2 * (mE2 - mE^2) / N2
  end

  l = @layout [  grid(1,2)]
  pE = plot(tT, tE)
  plot!(pE; xlabel="T", ylabel="E")
  plot!(pE; legend=false)

  pc = plot(tT, tc; xlabel="T", ylabel = "c")
  plot!(pc; legend=true)

  p = plot(pE,pclayout=l)
  savefig(p, "$outdir/mesures$N.svg") 
end


function H(ω)
  J = 1
  E = 0
  for I in CartesianRange(CartesianIndex(2,2), CartesianIndex(N-1,N-1))
    s = ω[I]
    for d in directions
      E += - J * s * ω[I + d]
    end
  end
  E
end






struct NonlocalModel
  ω::Array{Int,2}
  spins_down_coords::Array{ Tuple{Int,Int}, 1 }
  spins_up_coords::Array{ Tuple{Int,Int}, 1 }
end

function init_lattice(prng; is_local=true)
  mutable_sites_coordinates = [CartesianIndex(i,j) for i in 3:(N-2) for j in 3:(N-2)]
  #mutable_sites_coordinates = CartesianRange(CartesianIndex(3,3), CartesianIndex(N-2),(N-2))
  distrib = shuffle(prng, mutable_sites_coordinates)
  ω = ones(Int, (N, N))
  nb_spins_down = convert(Int, round(length(mutable_sites_coordinates) * ρ))
  spins_down_coords = collect(Iterators.take(distrib, nb_spins_down))
  spins_up_coords = collect(Iterators.drop(distrib, nb_spins_down))
  for X in spins_down_coords
    ω[X] = -1
  end
  if is_local
    ω
  else
    NonlocalModel(ω, spins_down_coords, spins_up_coords)
  end
end

function random_indexes(prng, n)
  # index in only one dimension
  # n: number of random indexes
  rand(prng, 3:(N - 2), n)
end

function random_coords(prng, n)
  CartesianIndex.(zip(random_indexes(prng, n), random_indexes(prng, n)))
end

function delta_neighbors(ω, X)
  s = 0
  for dX in directions
    s += ω[X + dX]
  end
  s
end

function sweep(ω, T, prng)
  for i in 1:N2
    kawasaki_local_move(ω, T, prng)
  end
end

function kawasaki_local_move(ω, T, prng)
  const k = 1
  const β = 1 / (k * T)
  const boltzmann_quotient = Dict([(i, exp(-i * β)) for i in -20:20])
  const RANDOM_COORDS = random_coords(model.lattice, N2)
  for X1 in RANDOM_COORDS
    spin1 = ω[X1]
    dX = rand(prng, directions)
    X2 = X1 + dX
    spin2 = ω[X2]
    if spin1 != spin2
      dneighbors = delta_neighbors(ω, X2) -
               delta_neighbors(ω, X1) -
               2 * spin1
      dE = - 2 * spin1 * dneighbors
      if dE <= 0 || get(boltzmann_quotient, dE, nothing) >= rand(prng)
        model.lattice.sites[X1] = spin2
        model.lattice.sites[X2] = spin1
      end
    end
  end
end

function make_heatmap(ω)
  Plots.heatmap(ω,
      color = :viridis,
      legend=false,
      grid=false,
      border=false,
      ticks=false,
      aspect_ratio=1)
end

main()
