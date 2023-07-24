
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

#cd("/home/rafaela/Documents/PUC/IC/Wind Speed/dados/Janauba/")
path = pwd()

janauba = CSV.read("JANAUBA_MG_UCT.csv", DataFrame)
janauba = janauba[1:359160,:]

#janaubaHourwindall.local_time = Date.(janaubaHourwindall.local_time, "yyy-mm-dd H:M")

#wind = @df janaubaHourwindall plot(:local_time, :wind_speed, title = "Wind Speed in Alegria 1 e 2 Wind Farms - RN", 
#xlab = "Hours", ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = :outerbottomright, label = "a")

#savefig(wind, "wind_all_anos_janauba.png")

vetor = Vector(janauba[!,4])

Ψ = CreatePsi(vetor)

println(maximum(vetor)) #1200

#################################### k-medoids plots ################################################################

df_por_mes = CSV.read("porMes_janauba_agnes_3K.csv", DataFrame)

cluster = CSV.read("cluster_janauba_agnes_3K.csv", DataFrame)
cluster = cluster.x

a1 = []
a2 = []
a3 = []
for i in 1:length(cluster)
    if cluster[i] == 1
        a1 = push!(a1, i)
    elseif cluster[i] == 2
        a2 = push!(a2, i)
    elseif cluster[i] == 3
        a3 = push!(a3, i)
    end
end
aa1 = Array{Float64}(undef, length(a1), 24)
k = 1
aa2 = Array{Float64}(undef, length(a2), 24)
m = 1
aa3 = Array{Float64}(undef, length(a3), 24)
j = 1
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
    end
end

plot1 = plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot2 = plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot3 = plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

savefig(plot1, "K3-K1-janauba-wind-agnes.png")
savefig(plot2, "K3-K2-janauba-wind-agnes.png")
savefig(plot3, "K3-K3-janauba-wind-agnes.png")

plot(Ψ[630,:], title = "Chosen days", label = "day 630", ylims = (0,16), legend = :outertopright)
plot!(Ψ[1087,:], label = "day 1087")
p = plot!(Ψ[3433,:], label = "day 3433")

savefig(p, "chosen-days-k10-janauba-winds.png")
