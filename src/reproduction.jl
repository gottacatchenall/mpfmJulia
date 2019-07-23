function reproduction(mpfm::MPFMInstance)
    pops::Array{Population} = mpfm.pops
    indivs::Array{Individual} = mpfm.indivs
    n_pops::Int64 = length(mpfm.pops)
    n_indivs::Int64 = length(mpfm.indivs)
    n_loci::Int64 = mpfm.genome_dict.n_loci
    b::Float64 = mpfm.params["fecundity"]
    genome_dict::GenomeDict = mpfm.genome_dict

    indivs_by_pop::Array{Array{Individual}} = split_indivs_by_pop(indivs, n_pops)

    for p = 1:n_pops
        this_pop_indivs::Array{Individual} = indivs_by_pop[p]

        n_indivs_this_pop::Int64 = length(this_pop_indivs)
        exp_num_offspring::Float64 = pops[p].mean_w * b * n_indivs_this_pop;
        n_off::Int64 = trunc(Int64, exp_num_offspring)
        if n_indivs_this_pop > 1
            for o = 1:n_off
                i1 = rand(DiscreteUniform(1,n_indivs_this_pop))
                i2 = i1
                while (i1 == i2)
                    i2 = rand(DiscreteUniform(1,n_indivs_this_pop))
                end

                parent1::Individual = this_pop_indivs[i1]
                parent2::Individual = this_pop_indivs[i2]

                offspring::Individual = Individual(p, n_loci)

                generate_haplo!(parent1, offspring, 1, genome_dict)
                generate_haplo!(parent2, offspring, 2, genome_dict)


                push!(mpfm.next_gen, offspring)
            end
        end
    end

#    println(mpfm.indivs)
#    println(mpfm.next_gen)

    mpfm.indivs = mpfm.next_gen
    mpfm.next_gen = []
#    println(mpfm.indivs)
#    println(mpfm.next_gen)

end

function generate_haplo!(parent::Individual, offspring::Individual, haplo_num::Int64, genome_dict::GenomeDict)
    n_loci::Int64 = genome_dict.n_loci

    mutation_sites = get_mutation_sites(genome_dict)
    crossing_over_sites::Array{Float64} = get_crossing_over_sites(genome_dict)
    chromo_map = genome_dict.chromo_map

    n_mutations::Int64 = (mutation_sites == nothing) ? 0 : length(mutation_sites)

    # could do this as blocks of chromosomes, but eh
    current_parent_haplotype::Int64 = rand(DiscreteUniform(1,2))

    co_pt::Int64 = 1
    this_chr::Int64 = 1
    mutation_ct::Int64 = 1

    for l = 1:n_loci
        # if there is a mutation this site, ignore haplotype changes etc...
        if n_mutations > 0
            if mutation_sites[mutation_ct] == l
                new_allele = rand(Uniform())
                set_haplo!(offspring, l, haplo_num, new_allele)
            end
        end

        if l > 1
            if chromo_map[l-1] != chromo_map[l]
                # new chr, pick parent haplo randomly
                current_parent_haplotype = rand(DiscreteUniform(1,2))
                this_chr = this_chr + 1
            end
            if co_pt < length(crossing_over_sites)
                # if co_put on [l-1, 1]
                if crossing_over_sites[co_pt] > (l-1) && crossing_over_sites[co_pt] < l
                    # flip haplo
                    (current_parent_haplotype == 1) ? (current_parent_haplotype = 2) : (current_parent_haplotype = 1)
                    co_pt = co_pt + 1
                end
            end
        end

        allele::Float64 = parent.genome.haplo[l, current_parent_haplotype]
        set_haplo!(offspring, l, haplo_num, allele)
    end

end

function get_mutation_sites(genome_dict::GenomeDict)
    mu::Float64 = genome_dict.mutation_rate
    n_loci::Int64 = genome_dict.n_loci
    mutation_ct::Int64 = rand(Binomial(n_loci, mu))
    if mutation_ct > 0
        return mutation_ct
    end

    return nothing
    # sorted array on indecies (preferablly unique, but hey)

end

function get_crossing_over_sites(genome_dict::GenomeDict)
    n_chromo = length(genome_dict.chromo_lengths)
    ret::Array{Float64} = []
    lo::Int64 = 0
    for c = 1:n_chromo
        exp_num_co::Float64 = genome_dict.chromo_lengths[c] / 100.0
        n_co::Int64 = rand(Poisson(exp_num_co))
        if n_co > 0
            if c > 1
                lo = genome_dict.chromo_ends[c-1]
                hi = genome_dict.chromo_ends[c]
            else
                lo = 0
                hi = genome_dict.chromo_ends[1]
            end

            co_locations::Array{Float64} = rand(Uniform(lo, hi), n_co)
            sort!(co_locations)
            ret = vcat(ret, co_locations)
            # sort this array and push it onto the return
        end

    end
    return ret
end
