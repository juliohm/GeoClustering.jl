"""
    GHC(k, variogram)

Geostatistical Hierarchical Clustering.

Long description....

## References
...
"""
struct GHC <: PartitionMethod
  variogram::V
  k::Int
end

function partition(data::Data, method::GHC,kernel_type,merge_criterion)
  # find lambda
  𝛌 = choose_𝛌(data,35) #literature rule of thumb
  # initialize history to store cluster/S(q) data
  n = length(data) #number of points
  h = Array{partition_info, 1}(undef, n)
  C = 1:n #most granular partiton
  d = zeros(n,n)#trivial placeholder for d
  for q in 1:n
    # find cluster partition C
    C = merge_clusters(data,C,d,"avg")
    # compute S(q)
    S = compute_silhouette(h,q,C)
    # store info (C,S) in history
    h[q]=(C,S)
  end
  # select clustering partition with best S(q)
  best_C = h[argmax (x -> x.C,h)].C
  return best_C
end

function ghc_kernel(𝛌,x,s,type)
    #uniform
    #triangular
    #Epanechnikov
    #Gaussian
end

function kernel_estimator(x,y,𝛌,type) 
    #equation 1
end

function dissimilarity(s_k,s_l,𝛌)
    #equation 2
end

function find_d()
    #equation 3, first pass
end
function update_d(data,C,C_new,d,merge_criterion) #compute dissimilarity matrix
    #equation 3 modified by linkage rule
    #only update d for changed i,⋅ or ⋅,i where i is one of the two merged clusters
    #single linkage
    #complete linkage
    #average linkage
    return d
end

function update_C(data,C,d)
    # go through d matrix, find i,j with lowest d, and combine into one cluster
    return C
end

function merge_clusters(data,C,d,merge_criterion) #join most similar clusters
    #check if first pass
    if (length(data)!=length(unique(C)))
        d=find_d()
    end
    C_new = update_C(data,C,d)
    #recompute d matrix
    d = update_d(data,C,C_new,d,merge_criterion)
    return C_new,d
end

function choose_𝛌(data,k)
    #find maximum distance of kth neighbor
end

function compute_silhouette()
    #equation 4
end

function compute_importance()
    #equation 5
end

#for storing history
struct partition_info
    C::Vector{Int} #assigns each datum to cluster
    S::Float64
end


