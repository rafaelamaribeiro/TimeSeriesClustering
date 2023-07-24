
# using Pkg
using DataFrames
using Plots
using CSV
using Dates
using LinearAlgebra

#=function CreatePsi(vec)
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
end=#

############### wind ########################################################################################

#cd("/home/rafaela/Documents/PUC/IC/Wind Speed/dados/Cidreira/")
path = pwd()

#cidreira = CSV.read("hweol_CIDREIRA_RS_Julia.csv", DataFrame)
#cidreira = cidreira[1:359160,:]

#vetor = Vector(cidreira[!,2])

#Ψ = CreatePsi(vetor)

println(maximum(vetor)) #25

#################################### k-medoids plots ################################################################

aa1 = CSV.read("profile1_cidreira_pam_12K.csv", DataFrame)
aa2 = CSV.read("profile2_cidreira_pam_12K.csv", DataFrame)
aa3 = CSV.read("profile3_cidreira_pam_12K.csv", DataFrame)
aa4 = CSV.read("profile4_cidreira_pam_12K.csv", DataFrame)
aa5 = CSV.read("profile5_cidreira_pam_12K.csv", DataFrame)
aa6 = CSV.read("profile6_cidreira_pam_12K.csv", DataFrame)
aa7 = CSV.read("profile7_cidreira_pam_12K.csv", DataFrame)
aa8 = CSV.read("profile8_cidreira_pam_12K.csv", DataFrame)
aa9 = CSV.read("profile9_cidreira_pam_12K.csv", DataFrame)
aa10 = CSV.read("profile10_cidreira_pam_12K.csv", DataFrame)
aa11 = CSV.read("profile11_cidreira_pam_12K.csv", DataFrame)
aa12 = CSV.read("profile12_cidreira_pam_12K.csv", DataFrame)

df_por_mes = CSV.read("porMes_cidreira_pam_12K.csv", DataFrame)

cluster = CSV.read("cluster_cidreira_pam_12K.csv", DataFrame)
cluster = cluster.x

#=a11 = []
a12 = []
for i in 1:length(cluster)
    if cluster[i] == 11
        a11 = push!(a11, i)
    elseif cluster[i] == 12
        a12 = push!(a12, i)
    end
end
aa11 = Array{Float64}(undef, length(a11), 24)
k = 1
aa12 = Array{Float64}(undef, length(a12), 24)
m = 1
for i in 1:length(cluster)
    if cluster[i] == 11
        aa11[k,:] = Ψ[i,:]
        k = k + 1
    elseif cluster[i] == 12
        aa12[m,:] = Ψ[i,:]
        m = m + 1
    end
end=#

aa1 = Matrix(aa1)
aa2 = Matrix(aa2)
aa3 = Matrix(aa3)
aa4 = Matrix(aa4)
aa5 = Matrix(aa5)
aa6 = Matrix(aa6)
aa7 = Matrix(aa7)
aa8 = Matrix(aa8)
aa9 = Matrix(aa9)
aa10 = Matrix(aa10)
aa11 = Matrix(aa11)
aa12 = Matrix(aa12)


plot1 = plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot2 = plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot3 = plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot4 = plot(aa4', title = "Cluster 4", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot5 = plot(aa5', title = "Cluster 5", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot6 = plot(aa6', title = "Cluster 6", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot7 = plot(aa7', title = "Cluster 7", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot8 = plot(aa8', title = "Cluster 8", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot9 = plot(aa9', title = "Cluster 9", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot10 = plot(aa10', title = "Cluster 10", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot11 = plot(aa11', title = "Cluster 11", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot12 = plot(aa12', title = "Cluster 12", xlims = (1, 24), ylims = (0,25), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)


savefig(plot1, "K12-K1-cidreira-wind-pam.png")
savefig(plot2, "K12-K2-cidreira-wind-pam.png")
savefig(plot3, "K12-K3-cidreira-wind-pam.png")
savefig(plot4, "K12-K4-cidreira-wind-pam.png")
savefig(plot5, "K12-K5-cidreira-wind-pam.png")
savefig(plot6, "K12-K6-cidreira-wind-pam.png")
savefig(plot7, "K12-K7-cidreira-wind-pam.png")
savefig(plot8, "K12-K8-cidreira-wind-pam.png")
savefig(plot9, "K12-K9-cidreira-wind-pam.png")
savefig(plot10, "K12-K10-cidreira-wind-pam.png")
savefig(plot11, "K12-K11-cidreira-wind-pam.png")
savefig(plot12, "K12-K12-cidreira-wind-pam.png")


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

savefig(p, "chosen-days-k10-cidreira-winds.png")
