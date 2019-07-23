function printsum(a)
    # summary generates a summary of an object
    println(summary(a), ": ", repr(a))
end

function init_default_populations(params::Dict{String, Float64})
    n_pops = params["num_populations"]
    pops::Array{Population} = []
    k::Float64 = params["num_individuals"] / n_pops
    efs::Array{Float64, 1} = [1.0]
    rng = MersenneTwister(1234);
    for p = 1:n_pops
        loc = rand(Uniform(), 2)
        pop = Population(loc[1], loc[2], k, efs)
        push!(pops, pop)
    end
    return pops
end

function init_default_genome_dict(params::Dict{String, Float64})
    n_loci = convert(Int64, params["n_loci"])
    n_chromo = convert(Int64, params["n_chromo"])
    n_loci_per_ef = convert(Int64,params["n_loci_per_ef"])
    n_ef = convert(Int64, params["n_ef"])
    genome_length::Float64 = params["genome_length"]
    n_sel_loci = n_ef * n_loci_per_ef
    n_loci_per_chromo::Int64 = n_loci / n_chromo

    mutation_rate = params["mutation_rate"]
    mean_poly_init = params["mean_poly_init"]

    fitness_loci = randperm(n_loci)[1:n_sel_loci]

    init_poly_ct::Array{Int64} = rand(Poisson(mean_poly_init), n_loci)
    chromo_map::Array{Int64} = []
    ef_map::Array{Int64} = zeros(Int64, n_loci)
    map_dist::Array{Float64} = []

    sel_weights::Array{Float64} = rand(Uniform(), n_sel_loci)
    fitness_locus_ct::Int64 = 1
    chromo_ct::Int64 = 1
    map_dist_ct::Float64 = 0.0

    map_dist_step::Float64 = genome_length / n_loci


    # length of each chromo array to easily pull num COs
    # uniform according to pdf of chromo length?
    chromo_length::Array{Float64} = zeros(n_chromo)

    chromo_ends::Array{Int64} = zeros(n_chromo)

    for l = 1:n_loci
        i = indexin(l,fitness_loci)
        if (l % n_loci_per_chromo == 0 && l > 1)
            chromo_ends[chromo_ct] = l
            if chromo_ct < n_chromo
                chromo_ct = chromo_ct + 1
            end
            map_dist_ct = 0.0
        end

        this_sel_weight::Float64 = 0.0
        ef_map[l] = 0
        if (i != nothing)
            ef_map[l] = rand(DiscreteUniform(1, n_ef))
            this_sel_weight = sel_weights[fitness_locus_ct]
            fitness_locus_ct = fitness_locus_ct + 1
        end
        this_chromo = chromo_ct
        chromo_length[chromo_ct] = chromo_length[chromo_ct] + map_dist_step
        this_map_dist = map_dist_ct + map_dist_step
        map_dist_ct = map_dist_ct + map_dist_step
        push!(sel_weights, this_sel_weight)
        push!(chromo_map, this_chromo)
        push!(map_dist, this_map_dist)
    end

    allele_dict::Array{Array{Float64}} = init_allele_dictionary(n_loci, init_poly_ct)
    genome_dict = GenomeDict(n_loci, mutation_rate, sel_weights, chromo_map, chromo_length, chromo_ends, map_dist, init_poly_ct, ef_map, allele_dict);
    return genome_dict
end

function init_indivs(pops::Array{Population}, genome_dict::GenomeDict, params::Dict{String, Float64})
    indivs::Array{Individual} = []
    n_pops::Int64 = convert(Int64, params["num_populations"])
    n_loci::Int64 = genome_dict.n_loci
    #for pop in pop, for i up to pop->k, etc.
    for p = 1:n_pops
        this_pop_k::Int64 = convert(Int64, pops[p].k)
        for i = 1:this_pop_k
            indiv::Individual = Individual(p, n_loci)
            push!(indivs, indiv)
        end
    end
    return indivs
end

function init_allele_dictionary(n_loci::Int64, init_poly_ct::Array{Int64})
    allele_vals::Array{Array{Float64}} = []

    for l = 1:n_loci
        push!(allele_vals, [])
        n_poly::Int64 = init_poly_ct[l]
        if (n_poly == 0)
            n_poly = 1
        end
        alleles = rand(Uniform(), n_poly)
        allele_vals[l] = alleles
    end
    return allele_vals
end

function init_genomes(mpfm::MPFMInstance)
    indivs::Array{Individual} = mpfm.indivs
    n_indivs = length(indivs)
    n_loci = mpfm.genome_dict.n_loci

    allele_vals = mpfm.genome_dict.allele_dict

    for ind = 1:n_indivs
        indiv::Individual = indivs[ind]
        for l = 1:n_loci
            n_poly = length(allele_vals[l])
            n_haplo = 2

            if (n_poly == 0)
                n_poly = 1
            end

            allele_indecies = rand(DiscreteUniform(1, n_poly), n_haplo)

            for h = 1:n_haplo
                ind = allele_indecies[h]
                set_haplo!(indiv, l, h, allele_vals[l][ind])
            end
        end
    end
end

function set_haplo!(indiv::Individual, locus::Int64, haplo::Int64, allele::Float64)
    indiv.genome.haplo[locus, haplo] = allele
end

function init_ibd_dispersal_kernel(pops::Array{Population}, params::Dict)
    n_pops = length(pops)
    dispersal_kernel::Array{Float64} = zeros(Float64, n_pops, n_pops)

    dispersal_decay::Float64 = params["dispersal_decay"]

    for i = 1:n_pops
        x1 = pops[i].x
        y1 = pops[i].y
        row_sum::Float64 = 0
        for j = 1:n_pops
            if (i != j)
                x2 = pops[j].x
                y2 = pops[j].y

                dist = sqrt((x2-x1)^2 + (y2-y1)^2)
                dispersal_kernel[i,j] = exp(-1*dispersal_decay*dist)
                row_sum += dispersal_kernel[i,j]
            end
        end

        for j = 1:n_pops
            if (i != j)
                dispersal_kernel[i,j] = dispersal_kernel[i,j] / row_sum
            end
        end
    end
    return dispersal_kernel
end

function init()
    desc = "
    Metapopulation Fragmentation Model (mpfm)
    v1.0

    A individual-based model of metapopopulation dynamics in spatiotemporally stochastic environments.

    University of Colorado at Boulder
    Dept. Ecology and Evolutionary Biology | Flaxman Lab
    Michael Catchen | michael.catchen@colorado.edu\n\n"
    println(desc)

    param_dict = Dict{String, Float64}("num_individuals" => 1000, "num_populations" => 20, "burn_in_length" => 500, "fragmentation_length" => 300, "mig_rate_burn_in" => 0.05, "mig_rate_fragmentation" => 0.005, "n_loci_per_ef" => 20, "n_ef" => 1, "n_loci" => 300, "n_chromo" => 5, "mean_poly_init" => 5, "dispersal_decay" => 5.0, "fecundity" => 3.0, "genome_length" => 100.0, "mutation_rate" => 0.00001, "sample" => -1)

    genome_dict::GenomeDict = init_default_genome_dict(param_dict)
    pops::Array{Population} = init_default_populations(param_dict)
    indivs::Array{Individual} = init_indivs(pops, genome_dict, param_dict)
    dispersal_kernel::Array{Float64} = init_ibd_dispersal_kernel(pops, param_dict)

    mpfm::MPFMInstance = MPFMInstance(param_dict, indivs, [], pops, genome_dict, dispersal_kernel)

    init_genomes(mpfm)

    return mpfm
end
