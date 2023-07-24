
# using Pkg
using DataFrames
using Plots
using CSV
using Dates
using LinearAlgebra

#=
function CreatePsi(vec)
    Nh = 24
    Ni = (365*10)+1
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

#cd("/home/rafaela/Documents/PUC/IC/Global Horizontal Irradayrion/dados/florestaS/")
path = pwd()

#floresta = CSV.read("floresta_CE_UCT_Julia.csv", DataFrame)
#floresta = floresta[1:359160,:]

#vetor = Vector(floresta[!,4])

#Ψ = CreatePsi(vetor)

#println(maximum(vetor)) #1100

#################################### k-medoids plots ################################################################

aa1 = CSV.read("profile1_floresta_k_means_3K.csv", DataFrame)
aa2 = CSV.read("profile2_floresta_k_means_3K.csv", DataFrame)
aa3 = CSV.read("profile3_floresta_k_means_3K.csv", DataFrame)

df_por_mes = CSV.read("porMes_floresta_k_means_3K.csv", DataFrame)

cluster = CSV.read("cluster_floresta_k_means_3K.csv", DataFrame)
cluster = cluster.x

aa1 = Matrix(aa1)
aa2 = Matrix(aa2)
aa3 = Matrix(aa3)

plot1 = plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,1100), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot2 = plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,1100), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

plot3 = plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,1100), xticks = 0:1:24, xlab = "Hours", 
ylab = "Global Horizontal Irradayrion (Wh/m²)", legend = false)

savefig(plot1, "K3-K1-floresta-sol-k_means.png")
savefig(plot2, "K3-K2-floresta-sol-k_means.png")
savefig(plot3, "K3-K3-floresta-sol-k_means.png")

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

savefig(p, "chosen-days-k10-floresta-sols.png")
