using ExoplanetsSysSim
using ApproximateBayesianComputing
const ABC = ApproximateBayesianComputing
using SpecialFunctions
using Statistics
import ApproximateBayesianComputing.CompositeDistributions.CompositeDist
import ApproximateBayesianComputing.TransformedBetaDistributions.LinearTransformedBeta
#import EvalSysSimModel

# https://en.wikipedia.org/wiki/Trigamma_function
function trigamma_x_gr_4(x::T) where T<: Real
   1/x + 0.5/x^2 + 1/(6*x^3) - 1/(30*x^5) + 1/(42*x^7) - 1/(30*x^9) + 5/(66*x^11) - 691/(2730*x^13) + 7/(6*x^15)
end

function trigamma_x_lt_4(x::T) where T<: Real
  n = floor(Int64,5-x)
  z = x+n 
  val = trigamma_x_gr_4(z)
  for i in 1:n
    z -= 1
    val += 1/z^2
  end
  val 
end

function trigamma(x::T) where T<: Real
   x >= 4 ? trigamma_x_gr_4(x) : trigamma_x_lt_4(x)
end

function var_weighted(x::AbstractArray{Float64,1}, w::AbstractArray{Float64,1} )
  #println("# size(x) = ",size(x), " size(w) = ", size(w)); flush(stdout)
  @assert(length(x)==length(w) )
  sumw = sum(w)
  @assert( sumw > 0. )
  if(sumw!= 1.0)
     w /= sum(w)
     sumw = 1.0
  end
  sumw2 = sum(w.*w)
  xbar =  sum(x.*w)
  covar = sum((x.-xbar).*(x.-xbar) .* w) * sumw/(sumw*sumw-sumw2)
end

function mom_alpha(x_bar::T, v_bar::T) where T<: Real 
    x_bar * (((x_bar * (1 - x_bar)) / v_bar) - 1)
end
function mom_beta(x_bar::T, v_bar::T) where T<: Real 
    (1 - x_bar) * (((x_bar * (1 - x_bar)) / v_bar) - 1)
end

# For algorithm, see https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=2613&context=etd
function fit_beta_mle(x::AbstractArray{T,1}; tol::T = 1e-6, max_it::Int64 = 10, init_guess::AbstractArray{T,1} = Array{T}(undef,0), w::AbstractArray{T,1} = Array{T}(undef,0), verbose::Bool = false ) where T<: Real
    lnxbar =   length(w)>1 ? Statistics.mean(log.(x),AnalyticWeights(w)) : Compat.Statistics.mean(log.(x))
    ln1mxbar = length(w)>1 ? Statistics.mean(log.(1.0.-x),AnalyticWeights(w)) : Compat.Statistics.mean(log.(1.0.-x))

    function itterate( mle_guess::Vector{T} ) where T<:Real
        (alpha, beta) = (mle_guess[1], mle_guess[2])
        dgab = digamma(alpha+beta)
        g1 = dgab - digamma(alpha) + lnxbar
        g2 = dgab - digamma(beta) + ln1mxbar
        tgab = trigamma(alpha+beta)
        G = [dgab-trigamma(alpha) tgab; tgab tgab-trigamma(beta)]
        mle_guess -= G \ [g1, g2]
    end 
  
    local mle_new 
    if length(init_guess) != 2
        xbar = length(w)>1 ? Compat.Statistics.mean(x,AnalyticWeights(w)) : Compat.Statistics.mean(x)
        vbar = length(w)>1 ? Compat.Statistics.varm(x,xbar,AnalyticWeights(w)) : Compat.Statistics.varm(x,xbar)
        mle_new = (vbar < xbar*(1.0-xbar)) ? [mom_alpha(xbar, vbar), mom_beta(xbar,vbar)] : ones(T,2)
    else
        mle_new = init_guess
    end
    if verbose
        println("it = 0: ", mle_new)
    end
    if any(mle_new.<=zero(T))
        println("# Warning: mean= ", xbar, " var= ",vbar," (alpha,beta)_init= ",mle_new," invalid, reinitializing to (1,1)")
        verbose = true
        mle_new = ones(T,2)
    end
    for i in 1:max_it
        mle_old = mle_new
        mle_new = itterate( mle_old )
        epsilon = max(abs.(mle_old.-mle_new))
        if verbose
            println("# it = ", i, ": ", mle_new, " max(Delta alpha, Delta beta)= ", epsilon)
        end
        if epsilon < tol
            break
        end
    end
    return mle_new
end

function make_beta(x::AbstractArray{T,1}, w::AbstractArray{T,1}; 
                   xmean::T = Compat.Statistics.mean(x,AnalyticWeights(w)), 
                   xvar::T = Compat.Statistics.varm(x,xmean,AnalyticWeights(w)), tau_factor::T=one(T) ) where T<:Real
    alpha_beta = (xvar < xmean*(1.0-xmean)) ? [mom_alpha(xmean, xvar), mom_beta(xmean,xvar)] : ones(T,2)
    if any(alpha_beta.<=zero(T))
        alpha_beta = fit_beta_mle(x, w=w, init_guess=alpha_beta, verbose=true)
    end
    if any(alpha_beta.<=zero(T))
        alpha_beta = ones(T,2)
    else 
        if minimum(alpha_beta)>1.5*tau_factor && sum(alpha_beta)>=20.0*tau_factor
            alpha_beta ./= tau_factor
        end
    end
    #println("Radius relative: a= ",alpha_beta[1], "  b= ",alpha_beta[2])
    Beta(alpha_beta[1], alpha_beta[2])
end

function make_beta_transformed(x::AbstractArray{T,1}, w::AbstractArray{T,1}; xmin::T=zero(T), xmax::T=one(T), xmean::T = Compat.Statistics.mean(x,AnalyticWeights(w)), xvar::T = Compat.Statistics.varm(x,xmean,AnalyticWeights(w)), tau_factor::T=one(T) ) where T<:Real
    alpha_beta = (xvar < xmean*(1.0-xmean)) ? [mom_alpha(xmean, xvar), mom_beta(xmean,xvar)] : ones(T,2)
    if any(alpha_beta.<=zero(T))
        alpha_beta = fit_beta_mle(x, w=w, init_guess=alpha_beta, verbose=true)
    end
    if any(alpha_beta.<=zero(T))
        alpha_beta = ones(T,2)
    else 
        if minimum(alpha_beta)>1.5*tau_factor && sum(alpha_beta)>=20.0*tau_factor
            alpha_beta ./= tau_factor
        end
    end
    #println("Total: a= ",alpha_beta[1], "  b= ",alpha_beta[2])
    LinearTransformedBeta(alpha_beta[1], alpha_beta[2], xmin=xmin, xmax=xmax)
end

function make_proposal_dist_multidim_beta(theta::AbstractArray{Float64,2}, weights::AbstractArray{Float64,1},  tau_factor::Float64; verbose::Bool = false)
    global sim_param_closure
    local limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    local p_dim = length(limitP)-1
    local r_dim = length(get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    
    theta_mean =  sum(theta.*weights',dims=2) # weighted mean for parameters
    theta_var = ABC.var_weighted(theta'.-theta_mean',weights)  # scaled, weighted covar for parameters
    tau_factor_indiv = fill(tau_factor,length(theta_var))

    dist_arr = ContinuousDistribution[]
    for j in 1:p_dim
        max_col_rate = 3*log(limitP[j+1]/limitP[j])/log(2)
        col_startidx = (j-1)*r_dim+1
        #tau_factor_indiv[col_startidx] = 2.0

        #=
        println("mean= ",theta_mean)
        println("var= ",theta_var)
        println("tau= ",tau_factor_indiv)
        for i in 1:length(theta_mean)
            println("a= ",alpha(theta_mean[i],tau_factor*theta_var[i]), "  b= ",beta(theta_mean[i],tau_factor*theta_var[i]))
        end
        =#
        
        dist_arr = vcat(dist_arr, ContinuousDistribution[make_beta_transformed(theta[i,:], weights, xmin=0.0, xmax=max_col_rate, xmean=theta_mean[i]/max_col_rate, xvar=theta_var[i]/max_col_rate^2, tau_factor=tau_factor_indiv[i]) for i in (col_startidx):(col_startidx+r_dim-1)])
    end

dist = CompositeDist(dist_arr)
end

function make_proposal_dist_multidim_beta(pop::abc_population_type, tau_factor::Float64; verbose::Bool = false)
    make_proposal_dist_multidim_beta(pop.theta, pop.weights, tau_factor, verbose=verbose)
end

function make_proposal_dist_multidim_beta_dirichlet(theta::AbstractArray{Float64,2}, weights::AbstractArray{Float64,1},  tau_factor::Float64; verbose::Bool = false)
    global sim_param_closure
    local limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    local p_dim = length(limitP)-1
    local r_dim = length(get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    
    theta_mean =  sum(theta.*weights',dims=2) # weighted mean for parameters
    theta_var = ABC.var_weighted(theta'.-theta_mean',weights)  # scaled, weighted covar for parameters
    tau_factor_indiv = fill(tau_factor,length(theta_var))

    dist_arr = ContinuousDistribution[]
    for j in 1:p_dim
        max_col_rate = 3*log(limitP[j+1]/limitP[j])/log(2)
        col_startidx = (j-1)*(r_dim+1)+1
        #tau_factor_indiv[col_startidx] = 2.0

        # if verbose
        #     println("total: ",theta_mean[col_startidx]," ",theta_var[col_startidx])
        # end
        # for i in (col_startidx+1):(col_startidx+r_dim)
        #     mean_ratio = sum(theta[col_startidx,:].*theta[i,:].*weights) /(theta_mean[col_startidx]*theta_mean[i]) # weighted mean for parameters
        #     var_ratio = var_weighted(vec(theta[col_startidx,:].*theta[i,:]).-(theta_mean[col_startidx]*theta_mean[i]),weights)/(2 * theta_mean[col_startidx] * theta_var[i]) # scaled, weighted covar for parameters
        #     if verbose
        #         println("i=",i,": ",theta_mean[i]," ",theta_var[i]," ratios: ",mean_ratio, " ",var_ratio)
        #     end
        #     var_ratio  = var_ratio  >= one(var_ratio)  ? var_ratio  : one(var_ratio)
        #     tau_factor_indiv[i] = tau_factor*var_ratio
        # end
        # if verbose
        #     flush(stdout)
        # end

        #=
        println("mean= ",theta_mean)
        println("var= ",theta_var)
        println("tau= ",tau_factor_indiv)
        for i in 1:length(theta_mean)
            println("a= ",alpha(theta_mean[i],tau_factor*theta_var[i]), "  b= ",beta(theta_mean[i],tau_factor*theta_var[i]))
        end
        =#
        
        dist_arr = vcat(dist_arr, make_beta_transformed(theta[col_startidx,:], weights, xmin=0.0, xmax=max_col_rate, xmean=theta_mean[col_startidx]/max_col_rate, xvar=theta_var[col_startidx]/max_col_rate^2, tau_factor=tau_factor_indiv[col_startidx]), ContinuousDistribution[ make_beta(theta[i,:], weights, xmean=theta_mean[i], xvar=theta_var[i], tau_factor=tau_factor_indiv[i]) for i in (col_startidx+1):(col_startidx+r_dim)])
    end

dist = CompositeDist(dist_arr)
end

function make_proposal_dist_multidim_beta_dirichlet(pop::abc_population_type, tau_factor::Float64; verbose::Bool = false)
    make_proposal_dist_multidim_beta_dirichlet(pop.theta, pop.weights, tau_factor, verbose=verbose)
end

function make_proposal_dist_multidim_beta_ratio(theta::AbstractArray{Float64,2}, weights::AbstractArray{Float64,1},  tau_factor::Float64; verbose::Bool = false)
    global sim_param_closure
    
    theta_mean =  sum(theta.*weights',dims=2) # weighted mean for parameters
    theta_var = ABC.var_weighted(theta'.-theta_mean',weights)  # scaled, weighted covar for parameters
    prior_max = 15.0

    dist_arr = ContinuousDistribution[make_beta_transformed(theta[1,:], weights, xmin=0.0, xmax=prior_max, xmean=theta_mean[1]/prior_max, xvar=theta_var[1]/prior_max^2, tau_factor=tau_factor)]

dist = CompositeDist(dist_arr)
end

function make_proposal_dist_multidim_beta_ratio(pop::abc_population_type, tau_factor::Float64; verbose::Bool = false)
    make_proposal_dist_multidim_beta_ratio(pop.theta, pop.weights, tau_factor, verbose=verbose)
end
