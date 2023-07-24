
# using Pkg
using DataFrames
using Plots
using CSV
using Dates
using LinearAlgebra

#=
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
=#
############### Sol ########################################################################################

#cd("/home/rafaela/Documents/PUC/IC/Global Horizontal Irradayrion/dados/Floresta/")
path = pwd()

#sertao = CSV.read("JANAUBA_MG_UCT.csv", DataFrame)
#sertao = sertao[1:359160,:]

#vetor = Vector(sertao[!,4])

#Ψ = CreatePsi(vetor)

#println(maximum(vetor)) #1200

#################################### k-medoids plots ################################################################

aa1 = CSV.read("profile1_sertao_k_means_3K.csv", DataFrame)
aa2 = CSV.read("profile2_sertao_k_means_3K.csv", DataFrame)
aa3 = CSV.read("profile3_sertao_k_means_3K.csv", DataFrame)

#df_por_mes = CSV.read("porMes_sertao_k_means_3K.csv", DataFrame)

#cluster = CSV.read("cluster_sertao_k_means_3K.csv", DataFrame)
#cluster = cluster.x

aa1 = Matrix(aa1)
aa2 = Matrix(aa2)
aa3 = Matrix(aa3)

plot1 = plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot2 = plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot3 = plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,1200), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

savefig(plot1, "K3-K1-sertao-sol-k_means.png")
savefig(plot2, "K3-K2-sertao-sol-k_means.png")
savefig(plot3, "K3-K3-sertao-sol-k_means.png")

plot(Ψ[630,:], title = "Chosen days", label = "day 630", ylims = (0,16), legend = :outertopright)
plot!(Ψ[1087,:], label = "day 1087")
p = plot!(Ψ[3433,:], label = "day 3433")

savefig(p, "chosen-days-k10-sertao-sols.png")
