
# using Pkg
using DataFrames
using Plots
using CSV
using Dates
using LinearAlgebra

function CreatePsi(vec)
    Nh = 24
    Ni = (365*41)
    Nc = 1
    matriz = zeros((Ni), Nc * Nh)
    m = 1
    a = 1
    for i in 1:length(vec)
        matriz[a, m] = vec[i]
        m = m + 1
        if m > Nh
            m = 1
            a = a + 1
        end
    end
    return matriz
end

############### wind ########################################################################################

#cd("/home/rafaela/Documents/PUC/IC/Wind Speed/dados/Floresta/")
path = pwd()

floresta = CSV.read("FLORESTA_RN_UCT_Julia.csv", DataFrame)
floresta = floresta[1:359160,:]

#florestaHourwindall.local_time = Date.(florestaHourwindall.local_time, "yyy-mm-dd H:M")

#wind = @df florestaHourwindall plot(:local_time, :wind_speed, title = "Wind Speed in Alegria 1 e 2 Wind Farms - RN", 
#xlab = "Hours", ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = :outerbottomright, label = "a")

#savefig(wind, "wind_all_anos_floresta.png")

vetor = Vector(floresta[!,4])

Ψ = CreatePsi(vetor)

#=n, m = size(Ψ)
function calcDist(Ψ, n)
    d = zeros(Float64, n, n)
    x = 0
    for p in 1:n
        for q in 1:n
            x = sum(((Ψ[p,h] - Ψ[q,h])^2) for h = 1:24)
            d[p,q] = sqrt(x)
        end
    end
    return d
end
D = calcDist(Ψ, n)

D[240,30]
=#
println(maximum(vetor)) #1100

#################################### k-medoids plots ################################################################

df_por_mes = CSV.read("porMes_floresta_agnes_4K.csv", DataFrame)

cluster = CSV.read("cluster_floresta_agnes_4K.csv", DataFrame)
cluster = cluster.x

a1 = []
a2 = []
a3 = []
a4 = []
for i in 1:length(cluster)
    if cluster[i] == 1
        a1 = push!(a1, i)
    elseif cluster[i] == 2
        a2 = push!(a2, i)
    elseif cluster[i] == 3
        a3 = push!(a3, i)
    elseif cluster[i] == 4
        a4 = push!(a4, i)
    end
end
aa1 = Array{Float64}(undef, length(a1), 24)
k = 1
aa2 = Array{Float64}(undef, length(a2), 24)
m = 1
aa3 = Array{Float64}(undef, length(a3), 24)
j = 1
aa4 = Array{Float64}(undef, length(a4), 24)
n = 1
for i in 1:length(cluster)
    if cluster[i] == 1
        aa1[k,:] = Ψ[i,:]
        k = k + 1
    elseif cluster[i] == 2
        aa2[m,:] = Ψ[i,:]
        m = m + 1
    elseif cluster[i] == 3
        aa3[j,:] = Ψ[i,:]
        j = j + 1
    elseif cluster[i] == 4
        aa4[n,:] = Ψ[i,:]
        n = n + 1
    end
end

plot1 = plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,1100), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot2 = plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,1100), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot3 = plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,1100), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot4 = plot(aa4', title = "Cluster 4", xlims = (1, 24), ylims = (0,1100), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)


savefig(plot1, "K4-K1-floresta-wind-agnes.png")
savefig(plot2, "K4-K2-floresta-wind-agnes.png")
savefig(plot3, "K4-K3-floresta-wind-agnes.png")
savefig(plot4, "K4-K4-floresta-wind-agnes.png")


plot(Ψ[630,:], title = "Chosen days", label = "day 630", ylims = (0,16), legend = :outertopright)
plot!(Ψ[1087,:], label = "day 1087")
plot!(Ψ[1250,:], label = "day 1250")
plot!(Ψ[1559,:], label = "day 1559")
plot!(Ψ[1679,:], label = "day 1679")
plot!(Ψ[1902,:], label = "day 1902")
plot!(Ψ[2607,:], label = "day 2607")
plot!(Ψ[2836,:], label = "day 2836")
plot!(Ψ[3311,:], label = "day 3311")
p = plot!(Ψ[3433,:], label = "day 3433")

savefig(p, "chosen-days-k10-floresta-winds.png")
