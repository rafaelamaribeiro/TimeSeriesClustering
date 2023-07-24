
# using Pkg
using DataFrames
using Plots
using CSV
using Dates
using LinearAlgebra
using StatsPlots

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

#cd("/home/rafaela/Documents/PUC/IC/Wind Speed/dados/Fernando/")
path = pwd()

#fernando = CSV.read("hweol_FERNANDO_RN_Julia.csv", DataFrame)
#fernando = fernando[1:359160,:]

#vetor = Vector(fernando[!,2])

#Ψ = CreatePsi(vetor)

#println(maximum(vetor)) #14

#################################### k-medoids plots ################################################################

aa1 = CSV.read("profile1_fernando_pam_2K.csv", DataFrame)
aa2 = CSV.read("profile2_fernando_pam_2K.csv", DataFrame)
aa3 = CSV.read("profile3_fernando_pam_12K.csv", DataFrame)
aa4 = CSV.read("profile4_fernando_pam_12K.csv", DataFrame)
aa5 = CSV.read("profile5_fernando_pam_12K.csv", DataFrame)
aa6 = CSV.read("profile6_fernando_pam_12K.csv", DataFrame)
aa7 = CSV.read("profile7_fernando_pam_12K.csv", DataFrame)
aa8 = CSV.read("profile8_fernando_pam_12K.csv", DataFrame)
aa9 = CSV.read("profile9_fernando_pam_12K.csv", DataFrame)
aa10 = CSV.read("profile10_fernando_pam_12K.csv", DataFrame)
aa11 = CSV.read("profile11_fernando_pam_12K.csv", DataFrame)
aa12 = CSV.read("profile12_fernando_pam_12K.csv", DataFrame)

df_por_mes = CSV.read("porMes_fernando_pam_12K.csv", DataFrame)

cluster = CSV.read("cluster_fernando_pam_12K.csv", DataFrame)
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

dfaa1 = DataFrame(aa1, :auto)
dfaa2 = DataFrame(aa2, :auto)
x = collect(1:24)

plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,14), xticks = 1:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false, palette = :Blues_3)

plot1 = @df dfaa1 boxplot!(
    x', [:x1 :x2 :x3 :x4 :x5 :x6 :x7 :x8 :x9 :x10 :x11 :x12 :x13 :x14 :x15 :x16 :x17 :x18 :x19 :x20 :x21 :x22 :x23 :x24],
    xlims = (0, 25), ylims = (0,14), xticks = 1:1:24,
    legend = false, color = "lightsteelblue1"
)

plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,14), xticks = 1:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false, palette = :Blues_3)

plot2 = @df dfaa2 boxplot!(
    x', [:x1 :x2 :x3 :x4 :x5 :x6 :x7 :x8 :x9 :x10 :x11 :x12 :x13 :x14 :x15 :x16 :x17 :x18 :x19 :x20 :x21 :x22 :x23 :x24],
    xlims = (0, 25), ylims = (0,14), xticks = 1:1:24, xlab = "Hours", ylab = "Wind Speed (m/s)", legend = false,
    color = "lightsteelblue1", title = "Cluster 2"
)

plot3 = plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot4 = plot(aa4', title = "Cluster 4", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot5 = plot(aa5', title = "Cluster 5", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot6 = plot(aa6', title = "Cluster 6", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot7 = plot(aa7', title = "Cluster 7", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot8 = plot(aa8', title = "Cluster 8", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot9 = plot(aa9', title = "Cluster 9", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot10 = plot(aa10', title = "Cluster 10", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot11 = plot(aa11', title = "Cluster 11", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot12 = plot(aa12', title = "Cluster 12", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)


savefig(plot1, "K2-K1-fernando-wind-pam-boxplot-lines.png")
savefig(plot2, "K2-K2-fernando-wind-pam-boxplot-lines.png")
savefig(plot3, "K12-K3-fernando-wind-pam.png")
savefig(plot4, "K12-K4-fernando-wind-pam.png")
savefig(plot5, "K12-K5-fernando-wind-pam.png")
savefig(plot6, "K12-K6-fernando-wind-pam.png")
savefig(plot7, "K12-K7-fernando-wind-pam.png")
savefig(plot8, "K12-K8-fernando-wind-pam.png")
savefig(plot9, "K12-K9-fernando-wind-pam.png")
savefig(plot10, "K12-K10-fernando-wind-pam.png")
savefig(plot11, "K12-K11-fernando-wind-pam.png")
savefig(plot12, "K12-K12-fernando-wind-pam.png")


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

savefig(p, "chosen-days-k10-fernando-winds.png")
