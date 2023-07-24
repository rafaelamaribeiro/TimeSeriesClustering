
# using Pkg
using DataFrames
using Plots
using CSV
using Dates
using LinearAlgebra
using StatsPlots

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

path = pwd()

sertao = CSV.read("SERTAO_BA_UCT_Julia.csv", DataFrame)
sertao = sertao[1:359160,:]

vetor = Vector(sertao[!,4])

Ψ = CreatePsi(vetor)

println(maximum(vetor)) #1200

#################################### k-medoids plots ################################################################

df_por_mes = CSV.read("porMes_sertao_agnes_3K.csv", DataFrame)

cluster = CSV.read("cluster_sertao_agnes_3K.csv", DataFrame)
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

dfaa1 = DataFrame(aa1, :auto)
dfaa2 = DataFrame(aa2, :auto)
dfaa3 = DataFrame(aa3, :auto)
x = collect(1:24)

plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false, palette = :Blues_3)

plot1 = @df dfaa1 boxplot!(
    x', [:x1 :x2 :x3 :x4 :x5 :x6 :x7 :x8 :x9 :x10 :x11 :x12 :x13 :x14 :x15 :x16 :x17 :x18 :x19 :x20 :x21 :x22 :x23 :x24],
    xlims = (0, 25), ylims = (0,1200), xticks = 1:1:24,
    legend = false, color = "lightsteelblue1"
)

plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false, palette = :Blues_3)

plot2 = @df dfaa2 boxplot!(
    x', [:x1 :x2 :x3 :x4 :x5 :x6 :x7 :x8 :x9 :x10 :x11 :x12 :x13 :x14 :x15 :x16 :x17 :x18 :x19 :x20 :x21 :x22 :x23 :x24],
    xlims = (0, 25), ylims = (0,1200), xticks = 1:1:24, legend = false,
    color = "lightsteelblue1", title = "Cluster 2"
)

plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false, palette = :Blues_3)

plot3 = @df dfaa3 boxplot!(
    x', [:x1 :x2 :x3 :x4 :x5 :x6 :x7 :x8 :x9 :x10 :x11 :x12 :x13 :x14 :x15 :x16 :x17 :x18 :x19 :x20 :x21 :x22 :x23 :x24],
    xlims = (0, 25), ylims = (0,1200), xticks = 1:1:24, legend = false,
    color = "lightsteelblue1"
)

plot4 = plot(aa4', title = "Cluster 4", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

savefig(plot1, "K3-K1-sertao-sol-pam-boxplot-lines.png")
savefig(plot2, "K3-K2-sertao-sol-pam-boxplot-lines.png")
savefig(plot3, "K3-K3-sertao-sol-pam-boxplot-lines.png")
savefig(plot4, "K4-K4-sertao-wind-agnes.png")

plot(Ψ[630,:], title = "Chosen days", label = "day 630", ylims = (0,16), legend = :outertopright)
plot!(Ψ[1087,:], label = "day 1087")
p = plot!(Ψ[3433,:], label = "day 3433")

savefig(p, "chosen-days-k10-sertao-winds.png")
