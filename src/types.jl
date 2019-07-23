mutable struct Genome
    n_loci::Int64
    n_haplo::Int64
    haplo::Array{Float64}
    Genome(n_loci::Int64, n_haplo::Int64) = new(n_loci, n_haplo, zeros(Float64, n_loci, n_haplo))
end

mutable struct Individual
    genome::Genome
    current_pop::Int64
    fitness::Float64
    Individual(current_pop::Int64, n_loci::Int64) = new(Genome(n_loci, 2), current_pop, 0)
end

mutable struct Population
    x::Float64
    y::Float64
    k::Float64
    mean_w::Float64
    efs::Array{Float64}
    Population(x::Float64, y::Float64, k::Float64, efs::Array{Float64}) = new(x,y,k,0.0,efs)
end

mutable struct GenomeDict
    n_loci::Int64
    mutation_rate::Float64
    selection_weights::Array{Float64}
    chromo_map::Array{Int64}
    chromo_lengths::Array{Float64}
    chromo_ends::Array{Int64}
    map_dist::Array{Float64}
    init_poly_ct::Array{Int64}
    ef_map::Array{Int64}
    allele_dict::Array{Array{Float64}}
end

mutable struct MPFMInstance
    params::Dict{String, Float64}
    indivs::Array{Individual};
    next_gen::Array{Individual};
    pops::Array{Population};
    genome_dict::GenomeDict
    dispersal_kernel::Array{Float64}
end
