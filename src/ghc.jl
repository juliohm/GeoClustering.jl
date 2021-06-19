# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GHC(kernel,merge_crit;k=35,vars=nothing)

A method for partitioning spatial data into clusters 
using Geostatistical Hierarchical Clustering (GHC).

## Parameters

* `kernel`      - type of kernel used for estimator
* `merge_crit`  - strategy to specify difference between clusters
* `k`           - minimum number of observations contained in support 
                    of kernel function centered at each point. 
                    Defaults to 35.
* `vars`        - Variables (or features) to consider (default to all)

## References
* F. Fouedjio. 2016. [A hierarchical clustering method for multivariate geostatistical data]
  (https://www.sciencedirect.com/science/article/abs/pii/S2211675316300367)
"""

uniform(z;λ) = (z ≤ λ)
triangular(z;λ) = (λ - z)*(z ≤ λ)
Epanechnikov(z;λ) = (λ^2 - z^2)*(z ≤ λ )
Gaussian(z;λ) = exp((-z^2)/(2λ^2))

struct GHC <: PartitionMethod
  kernel::Float64 #not sure how to add type for functional input
  merge_crit::AbstractString
  k::Int
  λ::Float64
  vars::Union{Vector{Symbol},Nothing}

  kernels=[uniform,triangular,Epanechnikov,Gaussian]
  @assert kernel ⊆ kernels "Invalid kernel given"
  merge_crits=[:single,:average,:complete]
  @assert merge_crit ⊆ merge_crits "Invalid merge criterion given"
end

GHC(kernel::Float64, merge_crit::AbstractString; k=35,λ=10, vars=nothing) =
  GHC(kernel, merge_crit, k,λ, vars)

function partition(data::Data, method::GHC)
  # variables used for clustering
  dvars = Tables.schema(values(data)).names
  vars = isnothing(method.vars) ? dvars : method.vars

  @assert vars ⊆ dvars "GHC features not found in geospatial data"

  kernel = method.kernel
  merge_crit = method.merge_crit
  k = method.k
  λ = method.λ

  # view subset of variables
  ctor = constructor(typeof(data))
  dom  = domain(data)
  tab  = TableOperations.select(values(data), vars...)
  Ω    = ctor(dom, Dict(paramdim(dom) => tab))

  #not sure if we want to use Ω or data
  D = D_λ(data,kernel,λ)
  

  dendrogram = hclust(d,merge_crit)
  return cut(dendrogram,k)
end

function D_λ(data,kernel,λ)
    dom = domain(data)
    tab = values(data)
    n = nelements(data)
    p = length(tab[1])
    
    d_mat = zeros(p,p)
    #equation 1
    function γij(n,dom,tab,x,y,λ,kernel,i,j)
        
        if (x==y)
            return 0.0
        else
            num = 0.0
            for k in 1:n
                for l in 1:n
                    kernel_density=K_λ(λ,norm(x-centroid(dom,k)),kernel)*K_λ(λ,norm(y-centroid(dom,l)),kernel)
                    diff_prod=(tab[k][i]-tab[l][i])*(tab[k][j]-tab[l][j])
                    num += kernel_density*diff_prod
                end
            end
            denom = 0.0
            for k in 1:n
                for l in 1:n
                    kernel_density=K_λ(λ,norm(x-centroid(dom,k)),kernel)*K_λ(λ,norm(y-centroid(dom,l)),kernel)
                    denom += 2*kernel_density
                end
            end
            return num/denom
        end
    end

    #equation 2
    function d_λ(n,p,dom,tab,sₖ,sₗ,λ,kernel)
        
        d=0
        for i in 1:p
            for j in 1:p
                d += abs(γij(n,dom,tab,sₖ,sₗ,λ,kernel,i,j))
            end
        end
        return d
    end

    #equation 3, first pass
    for k in 1:p
        for l in 1:p
            sₖ,sₗ=centroid(dom,k),centroid(dom,l)
            d_mat[k,l]=d_λ(n,p,dom,tab,sₖ,sₗ,λ,kernel)
        end
    end
    return d_mat ./ max(d_mat) #normalized

end

