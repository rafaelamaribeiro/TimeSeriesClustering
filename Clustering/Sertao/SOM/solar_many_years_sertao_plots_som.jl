
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

############### sol ########################################################################################

#cd("/home/rafaela/Documents/PUC/IC/Wind Speed/dados/Janauba/SOM")
#cd("./SOM")
path = pwd()

sertao = CSV.read("SERTAO_BA_UCT_Julia.csv", DataFrame)
sertao = sertao[1:359160,:]

vetor = Vector(sertao[!,4])

Ψ = CreatePsi(vetor)

println(maximum(vetor)) #1200

#################################### k-medoids plots ################################################################

cluster1 = CSV.read("sertao_clust1_som_12K_16n.csv", DataFrame)
cluster1 = cluster1.x
cluster2 = CSV.read("sertao_clust2_som_12K_16n.csv", DataFrame)
cluster2 = cluster2.x
cluster3 = CSV.read("sertao_clust3_som_12K_16n.csv", DataFrame)
cluster3 = cluster3.x
cluster4 = CSV.read("sertao_clust4_som_12K_16n.csv", DataFrame)
cluster4 = cluster4.x
cluster5 = CSV.read("sertao_clust5_som_12K_16n.csv", DataFrame)
cluster5 = cluster5.x
cluster6 = CSV.read("sertao_clust6_som_12K_16n.csv", DataFrame)
cluster6 = cluster6.x
cluster7 = CSV.read("sertao_clust7_som_12K_16n.csv", DataFrame)
cluster7 = cluster7.x
cluster8 = CSV.read("sertao_clust8_som_12K_16n.csv", DataFrame)
cluster8 = cluster8.x
cluster9 = CSV.read("sertao_clust9_som_12K_16n.csv", DataFrame)
cluster9 = cluster9.x
cluster10 = CSV.read("sertao_clust10_som_12K_16n.csv", DataFrame)
cluster10 = cluster10.x
cluster11 = CSV.read("sertao_clust11_som_12K_16n.csv", DataFrame)
cluster11 = cluster11.x
cluster12 = CSV.read("sertao_clust12_som_12K_16n.csv", DataFrame)
cluster12 = cluster12.x

aa1 = Array{Float64}(undef, length(cluster1), 24)
k = 1
aa2 = Array{Float64}(undef, length(cluster2), 24)
m = 1
aa3 = Array{Float64}(undef, length(cluster3), 24)
n = 1
aa4 = Array{Float64}(undef, length(cluster4), 24)
o = 1
aa5 = Array{Float64}(undef, length(cluster5), 24)
p = 1
aa6 = Array{Float64}(undef, length(cluster6), 24)
q = 1
aa7 = Array{Float64}(undef, length(cluster7), 24)
r = 1
aa8 = Array{Float64}(undef, length(cluster8), 24)
s = 1
aa9 = Array{Float64}(undef, length(cluster9), 24)
t = 1
aa10 = Array{Float64}(undef, length(cluster10), 24)
u = 1
aa11 = Array{Float64}(undef, length(cluster11), 24)
v = 1
aa12 = Array{Float64}(undef, length(cluster12), 24)
w = 1
#size(Ψ)[1]
for i in 1:size(Ψ)[1]
    if i in cluster1
        aa1[k,:] = Ψ[i,:]
        k = k + 1
    elseif i in cluster2
        aa2[m,:] = Ψ[i,:]
        m = m + 1
    elseif i in cluster3
        aa3[n,:] = Ψ[i,:]
        n = n + 1
    elseif i in cluster4
        aa4[o,:] = Ψ[i,:]
        o = o + 1
    elseif i in cluster5
        aa5[p,:] = Ψ[i,:]
        p = p + 1
    elseif i in cluster6
        aa6[q,:] = Ψ[i,:]
        q = q + 1
    elseif i in cluster7
        aa7[r,:] = Ψ[i,:]
        r = r + 1
    elseif i in cluster8
        aa8[s,:] = Ψ[i,:]
        s = s + 1
    elseif i in cluster9
        aa9[t,:] = Ψ[i,:]
        t = t + 1
    elseif i in cluster10
        aa10[u,:] = Ψ[i,:]
        u = u + 1
    elseif i in cluster11
        aa11[v,:] = Ψ[i,:]
        v = v + 1
    elseif i in cluster12
        aa12[w,:] = Ψ[i,:]
        w = w + 1
    end
end

plot1 = plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot2 = plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot3 = plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot4 = plot(aa4', title = "Cluster 4", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot5 = plot(aa5', title = "Cluster 5", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot6 = plot(aa6', title = "Cluster 6", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot7 = plot(aa7', title = "Cluster 7", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot8 = plot(aa8', title = "Cluster 8", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot9 = plot(aa9', title = "Cluster 9", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot10 = plot(aa10', title = "Cluster 10", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot11 = plot(aa11', title = "Cluster 11", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot12 = plot(aa12', title = "Cluster 12", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

savefig(plot1, "K12-K1-sertao-sol-som-16n.png")
savefig(plot2, "K12-K2-sertao-sol-som-16n.png")
savefig(plot3, "K12-K3-sertao-sol-som-16n.png")
savefig(plot4, "K12-K4-sertao-sol-som-16n.png")
savefig(plot5, "K12-K5-sertao-sol-som-16n.png")
savefig(plot6, "K12-K6-sertao-sol-som-16n.png")
savefig(plot7, "K12-K7-sertao-sol-som-16n.png")
savefig(plot8, "K12-K8-sertao-sol-som-16n.png")
savefig(plot9, "K12-K9-sertao-sol-som-16n.png")
savefig(plot10, "K12-K10-sertao-sol-som-16n.png")
savefig(plot11, "K12-K11-sertao-sol-som-16n.png")
savefig(plot12, "K12-K12-sertao-sol-som-16n.png")

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

savefig(p, "chosen-days-k10-sertao-sols.png")
