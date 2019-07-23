function run_mpfm(mpfm::MPFMInstance)
    params::Dict{String, Float64} = mpfm.params

    # Burn in
    n_gen_burn_in::Float64 = params["burn_in_length"]
    mig_rate_burn_in::Float64 = params["mig_rate_burn_in"]
    run_generations(mpfm, n_gen_burn_in, mig_rate_burn_in, false)

    # Fragmentation
    n_gen_fragmentation::Float64 = params["fragmentation_length"]
    mig_rate_fragmentation::Float64 = params["mig_rate_fragmentation"]
    run_generations(mpfm, n_gen_fragmentation, mig_rate_fragmentation, 20)
end

function run_generations(mpfm::MPFMInstance, n_gen::Number, mig_rate::Float64, log_freq)
    logging::Bool = false
    for gen = 1:n_gen
        if log_freq != false
            if (gen % 20 == 0)
                logging = gen % log_freq == 0 ? true : false
            end
        end
        run_gen(mpfm, mig_rate, logging)
    end
end

function run_gen(mpfm::MPFMInstance, mig_rate::Float64, log::Bool)
    dispersal(mpfm::MPFMInstance, mig_rate)
    selection(mpfm::MPFMInstance)

    if log
        logging(mpfm::MPFMInstance)
    end

    reproduction(mpfm::MPFMInstance)
end
