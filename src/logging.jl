function logging(mpfm::MPFMInstance)
    # big poop heres the logging fcn
    n_pops::Int64 = length(mpfm.pops)
    # subsample indivs and genotype those suckers
    sample_size::Int64 = -1
    if (mpfm.params["sample"] != -1)
        sample_size = mpfm.params["sample"]
    end

    # get indivs by pop
    sampled_indivs_by_pop::Array{Array{Individual}} = []
    indivs_by_pop::Array{Array{Individual}} = split_indivs_by_pop(mpfm.indivs, n_pops)
    if (sample_size == -1)
        sampled_indivs_by_pop = indivs_by_pop
    else
        sampled_indivs_by_pop = sample_indivs_by_pop(indivs_by_pop, n_pops, sample_size)
    end


end

function sample_indivs_by_pop(indivs_by_pop::Array{Array{Individual}}, n_pops::Int64, sample_size::Int64)
    for p = 1:n_pops
        this_pops_indivs::Array{Individual} = indivs_by_pop[p]
        pick_n_random_indivs(this_pops_indivs, sample_size)
    end


end

function pick_n_random_indivs(indivs::Array{Population}, sample_size)
    
end


function log_pairwise_data(mpfm::MPFMInstance)

end


function log_population_data(mpfm::MPFMInstance)

end
