function selection(mpfm::MPFMInstance)
    pops::Array{Population} = mpfm.pops
    indivs::Array{Individual} = mpfm.indivs

    n_pops::Int64 = length(mpfm.pops)
    n_indivs::Int64 = length(mpfm.indivs)
    n_loci::Int64 = length(mpfm.genome_dict.n_loci)
    indivs_by_pop::Array{Array{Individual}} = split_indivs_by_pop(indivs, n_pops)

    b = mpfm.params["fecundity"]

    rip_ct::Int64 = 0

    for p = 1:n_pops
        pop_sum::Float64 = 0.0
        pop_ct::Int64 = 0
        this_indivs = indivs_by_pop[p];
        n_indivs_this_pop::Int64 = length(this_indivs)
        efs::Array{Float64} = pops[p].efs
        for i = 1:n_indivs_this_pop
            indiv::Individual = this_indivs[i]
            fitness::Float64 = calc_fitness(indiv, efs, mpfm)
            indiv.fitness = fitness

            k_prime::Float64 = pops[p].k * fitness
            prop = beverton_holt_prob(k_prime, n_indivs_this_pop, b)

            surv::Bool = false;
            if (rand(Uniform()) < prop)
                surv = true;
            end


            if (surv)
                pop_sum += fitness
                pop_ct += 1
            else
                # rip :'(
                rip_ct += 1
                indiv.current_pop = -1
            end
        end
        this_pop_avg_fitness::Float64 = 0.0
        if (pop_ct > 0)
            this_pop_avg_fitness = pop_sum/pop_ct;
        end
        pops[p].mean_w = this_pop_avg_fitness
    end

#    println("n_rip = ", rip_ct)
#    println("n_total = ", n_indivs)
#    println("\n")

end

function beverton_holt_prob(k_prime::Float64, n_indivs_this_pop::Int64, b::Float64)
    prop_full::Float64 = n_indivs_this_pop/k_prime
    prob::Float64 = 1.0 / ((1.0 + (b - 1))*prop_full)
    return prob
end

function calc_fitness(indiv::Individual, efs::Array{Float64}, mpfm::MPFMInstance)
    n_loci::Int64 = indiv.genome.n_loci
    n_haplo = indiv.genome.n_haplo

    chromo_map::Array{Int64} = mpfm.genome_dict.chromo_map
    allele_dict::Array{Array{Float64}} = mpfm.genome_dict.allele_dict
    selection_weights::Array{Float64} = mpfm.genome_dict.selection_weights
    ef_map::Array{Int64} = mpfm.genome_dict.ef_map
    #println(allele_dict)

    w = 1.0
    s_max::Float64 = -0.004
    for l = 1:n_loci
        if (ef_map[l] != 0)
            if selection_weights[l] > 0
                for h = 1:n_haplo
                    x::Float64 = indiv.genome.haplo[l, h]
                    ef_val::Float64 = efs[ef_map[l]]
                    diff::Float64 = abs(x - ef_val)
                    @assert diff <= 1.0
                    gauss::Float64 = exp(-1*(diff)^2)
                    s_i::Float64 = s_max * selection_weights[l] * gauss
                    w_i::Float64 = 1.0 + s_i
                    w = w * w_i
                end
            end
        end
    end
    #println("W: ", w)
    return w
end
