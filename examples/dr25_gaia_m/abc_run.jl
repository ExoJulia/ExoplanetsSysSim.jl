## ExoplanetsSysSim/examples/dr25_gaia_m/abc_run.jl
## (c) 2019 Danley C. Hsu & Eric B. Ford
# Script for producing DR25 M planet candidate occurrence rate estimates

include("abc_setup.jl")

using SysSimABC
using ExoplanetsSysSim
using JLD
using StatsBase

out2txt = false # Write occurrence rates & densities to text files
expandpart = false # Expand final generation for robust posteriors
prior_choice = "dirichlet"
bin_size_factor = 2.0
pop_num = 1

println("Setting up simulation...")
@time abc_plan = setup_abc(prior_choice = prior_choice, bin_size_factor = bin_size_factor)
println("")
println("Running simulation...")
@time output = run_abc(abc_plan)
# println("")
# println("Running simulation (part 2)...")
# @time abc_plan = setup_abc_p2(abc_plan)
# @time output = run_abc(abc_plan, output)
#@time abc_plan = change_distance()
#@time output = run_abc(abc_plan, output)
println("")

save(string("test-pop-out",pop_num,".jld"), "output", output, "ss_true", EvalSysSimModel.get_ss_obs())

if expandpart
    println("Expanding to large generation...")
    @time theta_largegen, weights_largegen = run_abc_largegen(abc_plan, output, EvalSysSimModel.get_ss_obs(), output.accept_log.epsilon[end-1], npart=1000)
    println("")

    save(string("test-pop-out",pop_num,".jld"), "output", output, "ss_true", EvalSysSimModel.get_ss_obs(), "theta_largegen", theta_largegen, "weights_largegen", weights_largegen)
end

if out2txt
    file_rate = open("rate_output.txt", "w")
    file_dens = open("dens_output.txt", "w")
end

limitP = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
limitR = get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1})
const r_dim = length(limitR)-1

if expandpart
    weight_vec = aweights(weights_largegen)
    #weight_vec = aweights(fill(1.0, length(weights_largegen)))
else
    weight_vec = aweights(output.weights)
    #weight_vec = aweights(fill(1.0, length(output.weights)))
end

for p_ind = 1:(length(limitP)-1)
    col_ind = (p_ind-1)*(r_dim)+1
    for r_ind = 1:r_dim
        bin_ind = (p_ind-1)*(r_dim)+r_ind
        dens_denom = 1.0/log(limitP[p_ind+1]/limitP[p_ind])/log(limitR[r_ind+1]/limitR[r_ind])

        if prior_choice == "dirichlet" && r_dim > 1
            col_ind = (p_ind-1)*(r_dim+1)+1
            bin_ind = (p_ind-1)*(r_dim+1)+r_ind+1
            if expandpart
                quant_arr = quantile(theta_largegen[bin_ind,:].*theta_largegen[col_ind,:], weight_vec, [0.1587, 0.5, 0.8413])
            else
                quant_arr = quantile(output.theta[bin_ind,:].*output.theta[col_ind,:], weight_vec, [0.1587, 0.5, 0.8413])
            end
        elseif prior_choice == "beta"
            col_lambda = bin_size_factor * 3 * log(limitP[p_ind+1]/limitP[p_ind])/log(2)
            if expandpart
                quant_arr = quantile(theta_largegen[bin_ind,:]*col_lambda, weight_vec, [0.1587, 0.5, 0.8413])
            else
                quant_arr = quantile(output.theta[bin_ind,:]*col_lambda, weight_vec, [0.1587, 0.5, 0.8413])
            end
        else
            if expandpart
                quant_arr = quantile(theta_largegen[bin_ind,:], weight_vec, [0.1587, 0.5, 0.8413])
            else
                quant_arr = quantile(output.theta[bin_ind,:], weight_vec, [0.1587, 0.5, 0.8413])
            end
        end

        println("-----------------------------")
        println("Orbital Period (day) = ", string(limitP[p_ind:p_ind+1]), " / Planet Radius (R_earth) = ", string(limitR[r_ind:r_ind+1]/ExoplanetsSysSim.earth_radius))
        println("")
        println("Rate = ", string(quant_arr[2], " + ", quant_arr[3]-quant_arr[2], " - ", quant_arr[2]-quant_arr[1]))
        println("Density = ", string(quant_arr[2]*dens_denom, " + ", (quant_arr[3]-quant_arr[2])*dens_denom, " - ", (quant_arr[2]-quant_arr[1])*dens_denom))

        if out2txt
            write(file_rate, string(quant_arr[2], " + ", quant_arr[3]-quant_arr[2], " - ", quant_arr[2]-quant_arr[1], "\n"))
            write(file_dens, string(quant_arr[2]*dens_denom, " + ", (quant_arr[3]-quant_arr[2])*dens_denom, " - ", (quant_arr[2]-quant_arr[1])*dens_denom, "\n"))
        end
    end
end

if out2txt
    close(file_rate)
    close(file_dens)
end

println("-----------------------------")
println("")
println(EvalSysSimModel.get_ss_obs())
