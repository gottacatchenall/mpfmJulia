function dispersal(mpfm::MPFMInstance, mig_rate::Float64)
    pops::Array{Population} = mpfm.pops
    indivs::Array{Individual} = mpfm.indivs

    n_pops = length(mpfm.pops)
    n_indivs = length(mpfm.indivs)
    mean_migrants_per_gen::Float64 = 1.0

    indivs_by_pop::Array{Array{Individual}} = split_indivs_by_pop(indivs, n_pops)

    n_migrants_from_each = rand(Poisson(mean_migrants_per_gen), n_pops)
    for p = 1:n_pops
        this_pop_indivs::Array{Individual} = indivs_by_pop[p]
        n_indivs_this_pop::Int64 = length(this_pop_indivs)
        n_migrants_from_this_pop = n_migrants_from_each[p]
        if (n_migrants_from_this_pop > 0 && n_indivs_this_pop > 0)
            for m = 1:n_migrants_from_this_pop
                randindex = rand(DiscreteUniform(1, n_indivs_this_pop))
                this_pop_indivs[randindex].current_pop = p
            end
        end
    end

end

function split_indivs_by_pop(indivs::Array{Individual}, n_pops)
    n_indivs = length(indivs)

    indivs_by_pop::Array{Array{Individual}} = []
    for p = 1:n_pops
        push!(indivs_by_pop, [])
    end


    for i = 1:n_indivs
        indiv::Individual = indivs[i]
        this_pop = indiv.current_pop
        if (this_pop >= 0)
            push!(indivs_by_pop[this_pop], indiv)
        end
    end

    #println(ndims(indivs_by_pop))

    return indivs_by_pop
end
