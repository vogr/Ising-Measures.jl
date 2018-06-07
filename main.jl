
tmod = now()
using Memento
using Plots
Plots.gr()
using ProgressMeter
using DataStructures

dmodules = now() - tmod

const J = 1 # Neighbor interraction
const k = 1 # Boltzmann's constant
const dynamic = "k_local"  # in ["k_local", "k_nonlocal", "metropolis"]
const ρ = 0.4 # spins down ratio
const N = 50 # Lattice size
const N2 = N^2
const directions = CartesianIndex.([(0,1), (0,-1), (1,0), (-1,0)])


const SEED = 0x5e6c1f2e

const n_temperatures = 100
const eq_steps = 5 * 10^3
const mc_steps = 10^4
const plot_every = 10^3

const heat_bath_time = 10^5
const heat_bath_temp1 = 1.2
const heat_bath_temp2 = 0.1

function main()
  logger = Memento.config("debug"; fmt="[{level} | {name}]: {msg}")

  info(logger, "Loaded modules in $dmodules")

  info(logger, "Seed: $SEED")
  prng = MersenneTwister(SEED)

  const tT = linspace(0.1, 6, n_temperatures)
  const outdir = "simulations/$(dynamic)_$(SEED)_$(N)"
  info(logger, "Making output directory: $outdir")
  mkpath(outdir)

  
  model = init_model(prng; is_local=true)

  t0 = now() 
  p0 = make_heatmap(model.ω)
  plot!(p0, title="Initial lattice (t=0)")
  savefig(p0, "$outdir/image0.svg")
  dfirstplot = now() - t0
  info(logger, "Time for first plot: $dfirstplot")


  #First heat bath
  info(logger, "Beginning heat bath")
  hbath = now()
  β_hb1 = 1 / (k * heat_bath_temp1)
  for i in 1:heat_bath_time
    mcmove(model, β_hb1, prng)
  end
  p = make_heatmap(model.ω)
  plot!(p, title="After first heatbath")
  savefig(p, "$outdir/heatbath1.svg")
  β_hb2 = 1 / (k * heat_bath_temp2)
  for i in 1:heat_bath_time
    mcmove(model, β_hb2, prng)
  end
  p = make_heatmap(model.ω)
  plot!(p, title="After second heatbath")
  savefig(p, "$outdir/heatbath2.svg")
  dbath = now() - hbath
  info(logger, "Heat bath done in $dbath (2 * $heat_bath_time sweeps at T=$heat_bath_temp1 and $heat_bath_temp2")


  tE = similar(tT)
  tE2 = similar(tT)
  tc = similar(tT)

  h0 = now()
  for (i_T,T) in enumerate(tT)
    info(logger, "Starting T=$T")
    subdir = "$outdir/$T"
    mkpath(subdir)

    β = 1 / (k * T)
    totE = 0
    totE2 = 0
    t = 0

    qsize = 10
    durations = Queue(Int)

    p = make_heatmap(model.ω)
    plot!(p, title="T=$T, t=$t (eq.)")
    savefig(p, "$subdir/image$(dec(t, 9)).svg")

    # Reach a first equilibrium with eq_steps step then compute energy
    # on the next mc_steps configurations.
    h1 = now()
    for i in 1:(eq_steps + mc_steps)
      mcmove(model, β, prng)
      t = t + 1

      if i > eq_steps
        E = H(model.ω)
        totE = totE + E
        totE2 = totE2 + E*E
      end

      if i % plot_every == 0
        h2 = now()
        d = h2 - h1
        enqueue!(durations, Dates.value(d))
        if length(durations) > qsize
          dequeue!(durations)
        end

        hplot = now()
        p = make_heatmap(model.ω)
        if i > eq_steps
          plot!(p, title="T=$T, t=$t, E=$E")
        else
          plot!(p, title="T=$T, t=$t (eq.)")
        end
        savefig(p, "$subdir/image$(dec(t, 9)).svg")
        dplot = now() - hplot

        total_time = Dates.canonicalize(Dates.CompoundPeriod(now() - h0))
        sw_per_s = durations_to_swps(durations)
        message = """
        ($total_time) T=$(T) ; t=$(t)
        \t$plot_every sweeps in $d. On average (last $qsize): $sw_per_s sweeps/s.
        \t+$dplot to plot.
        """
        info(logger, message)
        h1 = now()
      end
    end

    mE = 1 / mc_steps * totE
    mE2 = 1 / mc_steps * totE2
    tE[i_T] = mE
    tE2[i_T] = mE2
    tc[i_T] = k * β^2 * (mE2 - mE^2) / N2
  end

  l = @layout [a b]
  pE = plot(tT, tE)
  plot!(pE; xlabel="T", ylabel="E")
  plot!(pE; legend=false)

  pc = plot(tT, tc)
  plot!(pc; xlabel="T", ylabel = "c")
  plot!(pc; legend=true)

  p = plot(pE,pc,layout=l)
  savefig(p, "$outdir/mesures.svg") 
end

function durations_to_swps(durations)
  sw_per_ms = (plot_every * length(durations)) / sum(durations)
  sw_per_s = round(sw_per_ms * 1000, 4)
	sw_per_s
end

function H(ω)
  E = 0
  for I in CartesianRange(CartesianIndex(2,2), CartesianIndex(N-1,N-1))
    s = ω[I]
    for d in directions
      E += - J * s * ω[I + d]
    end
  end
  E
end


struct LocalModel
  ω::Array{Int,2}
end


struct NonlocalModel
  ω::Array{Int,2}
  spins_down_coords::Array{CartesianIndex{2}, 1}
  spins_up_coords::Array{CartesianIndex{2}, 1}
end

struct MetropolisModel
  ω::Array{Int,2}
end

function init_model(prng; is_local=true)
  mutable_sites_coordinates = [CartesianIndex(i,j) for i in 3:(N-2) for j in 3:(N-2)]
  distrib = shuffle(prng, mutable_sites_coordinates)
  ω = ones(Int, (N, N))
  nb_spins_down = convert(Int, round(length(mutable_sites_coordinates) * ρ))
  spins_down_coords = collect(Iterators.take(distrib, nb_spins_down))
  spins_up_coords = collect(Iterators.drop(distrib, nb_spins_down))
  for X in spins_down_coords
    ω[X] = -1
  end
  if dynamic == "k_local"
    LocalModel(ω)
  elseif dynamic == "k_nonlocal"
    NonlocalModel(ω, spins_down_coords, spins_up_coords)
  elseif dynamic == "metropolis"
    MetropolisModel(ω)
  end
end

function random_indexes(prng, n)
  # index in only one dimension
  # n: number of random indexes
  rand(prng, 3:(N - 2), n)
end

function random_coords(prng, n)
  CartesianIndex.(collect(zip(random_indexes(prng, n), random_indexes(prng, n))))
end

function random_up_and_down_spins_indices(prng, model::NonlocalModel, n)
  zip(rand(prng, 1:length(model.spins_up_coords), n), rand(prng, 1:length(model.spins_down_coords), n))
end

function delta_neighbors(ω, X)
  s = 0
  for dX in directions
    s += ω[X + dX]
  end
  s
end


function mcmove(m::LocalModel, β, prng) # Kawasaki local move
  const boltzmann_quotient = Dict([(i, exp(-i * β)) for i in -20:20])
  const RANDOM_COORDS = random_coords(prng, N2)
  for X1 in RANDOM_COORDS
    spin1 = m.ω[X1]
    dX = rand(prng, directions)
    X2 = X1 + dX
    spin2 = m.ω[X2]
    if spin1 != spin2
      dneighbors = delta_neighbors(m.ω, X2) -
               delta_neighbors(m.ω, X1) -
               2 * spin1
      dE = - 2 * spin1 * dneighbors
      if dE <= 0 || get(boltzmann_quotient, dE, nothing) >= rand(prng)
        m.ω[X1] = spin2
        m.ω[X2] = spin1
      end
    end
  end
end

function mcmove(model::NonlocalModel, β, prng)
  const boltzmann_quotient = Dict([(i, exp(-i * β)) for i in -20:20])
  const RANDOM_UP_AND_DOWN_SPINS_INDICES = random_up_and_down_spins_indices(prng, model, N2)
  for (i_up, i_down) in RANDOM_UP_AND_DOWN_SPINS_INDICES
    # Coordinates on the graph
    Xd, Xu = model.spins_down_coords[i_down], model.spins_up_coords[i_up]
    dneighbors = delta_neighbors(model.ω, Xu)
          - delta_neighbors(model.ω, Xd)
    dE = 2 * dneighbors
    if dE <= 0 || get(boltzmann_quotient, dE, nothing) >= rand(prng)
      model.spins_down_coords[i_down] = Xu
      model.spins_up_coords[i_up] = Xd
      model.ω[Xd] = +1
      model.ω[Xu] = -1
    end
  end
end

function mcmove(model::MetropolisModel, β, prng)
  const boltzmann_quotient = Dict([(i, exp(-i * β)) for i in -20:20])
  const RANDOM_COORDS = random_coords(prng, N2)
  for X in RANDOM_COORDS
    spin = model.ω[X]
    dE = 2 * spin * delta_neighbors(model.ω, X)
    if dE <= 0 || get(boltzmann_quotient, dE, nothing) > rand(prng)
      model.ω[X] = -spin
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
