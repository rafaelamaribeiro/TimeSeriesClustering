
# using Pkg
using DataFrames
using Plots
using CSV
using Dates
using LinearAlgebra
#using StatsPlots

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

#cd("/home/rafaela/Documents/PUC/IC/Wind Speed/dados/Fernando/SOM")
cd("./SOM")
#path = pwd()

fernando = CSV.read("hweol_FERNANDO_RN_Julia.csv", DataFrame)
fernando = fernando[1:359160,:]

#fernandoHourwindall.local_time = Date.(fernandoHourwindall.local_time, "yyy-mm-dd H:M")

#wind = @df fernandoHourwindall plot(:local_time, :wind_speed, title = "Wind Speed in Alegria 1 e 2 Wind Farms - RN", 
#xlab = "Hours", ylab = "Wind Speed (m/s)", legend = :outerbottomright, label = "a")

#savefig(wind, "wind_all_anos_fernando.png")

vetor = Vector(fernando[!,2])

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
println(maximum(vetor)) #14

#################################### k-medoids plots ################################################################
#=
aa1 = CSV.read("profile1_fernando_som.csv", DataFrame)
aa2 = CSV.read("profile2_fernando_som.csv", DataFrame)
aa3 = CSV.read("profile3_fernando_som.csv", DataFrame)
aa4 = CSV.read("profile4_fernando_som.csv", DataFrame)
aa5 = CSV.read("profile5_fernando_som.csv", DataFrame)
aa6 = CSV.read("profile6_fernando_som.csv", DataFrame)
aa7 = CSV.read("profile7_fernando_som.csv", DataFrame)
aa8 = CSV.read("profile8_fernando_som.csv", DataFrame)
#aa9 = CSV.read("profile9_fernando_som.csv", DataFrame)
#aa10 = CSV.read("profile10_fernando_som.csv", DataFrame)
#aa11 = CSV.read("profile11_fernando_som.csv", DataFrame)
#aa12 = CSV.read("profile12_fernando_som.csv", DataFrame)=#

#df_por_mes = CSV.read("porMes_fernando_som.csv", DataFrame)

cluster1 = CSV.read("fernando_clust1_som.csv", DataFrame)
cluster1 = cluster1.x
cluster2 = CSV.read("fernando_clust2_som.csv", DataFrame)
cluster2 = cluster2.x

aa1 = Array{Float64}(undef, length(cluster1), 24)
k = 1
aa2 = Array{Float64}(undef, length(cluster2), 24)
m = 1
#size(Ψ)[1]
for i in 1:size(Ψ)[1]
    if i in cluster1
        aa1[k,:] = Ψ[i,:]
        k = k + 1
    elseif i in cluster2
        aa2[m,:] = Ψ[i,:]
        m = m + 1
    end
end

#fernando_3anos = vcat(fernandoHourwind_2017)#, fernandoHourwind_2018, fernandoHourwind_2019)

#plotMedoids = scatter(fernando_3anos.wind_speed, marker_z = cor.color, legend = false, title = "3 clusters fernando Nova", ylims = (0,16))

#savefig(plotMedoids, "fernando-10k-medoids-scatter-wind-all-anos")

#=aa1 = Matrix(aa1)
aa2 = Matrix(aa2)
aa3 = Matrix(aa3)
aa4 = Matrix(aa4)
aa5 = Matrix(aa5)
aa6 = Matrix(aa6)
aa7 = Matrix(aa7)
aa8 = Matrix(aa8)
aa9 = Matrix(aa9)
aa10 = Matrix(aa10)=#


plot1 = plot(aa1', title = "Cluster 1", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

plot2 = plot(aa2', title = "Cluster 2", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)

#=plot3 = plot(aa3', title = "Cluster 3", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
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
ylab = "Wind Speed (m/s)", legend = false)=#


savefig(plot1, "K2-K1-fernando-wind-som-20.png")
savefig(plot2, "K2-K2-fernando-wind-som-20.png")
#=savefig(plot3, "K12-K3-fernando-wind-som.png")
savefig(plot4, "K12-K4-fernando-wind-som.png")
savefig(plot5, "K12-K5-fernando-wind-som.png")
savefig(plot6, "K12-K6-fernando-wind-som.png")
savefig(plot7, "K12-K7-fernando-wind-som.png")
savefig(plot8, "K12-K8-fernando-wind-som.png")
savefig(plot9, "K12-K9-fernando-wind-som.png")
savefig(plot10, "K12-K10-fernando-wind-som.png")
savefig(plot11, "K12-K11-fernando-wind-som.png")
savefig(plot12, "K12-K12-fernando-wind-som.png")=#


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

cor = zeros(length(vetor))
m = 24
i = 1
while i < length(vetor)
    for k in 1:((365*10)+1)
        if cluster[k] == 1
            cor[i:i + m - 1] .= 1.0
        elseif cluster[k] == 2
            cor[i:i + m - 1] .= 2.0
        elseif cluster[k] == 3
            cor[i:i + m - 1] .= 3.0
        elseif cluster[k] == 4
            cor[i:i + m - 1] .= 4.0
        elseif cluster[k] == 5
            cor[i:i + m - 1] .= 5.0
        elseif cluster[k] == 6
            cor[i:i + m - 1] .= 6.0
        elseif cluster[k] == 7
            cor[i:i + m - 1] .= 7.0
        elseif cluster[k] == 8
            cor[i:i + m - 1] .= 8.0
        elseif cluster[k] == 9
            cor[i:i + m - 1] .= 9.0
        elseif cluster[k] == 10
            cor[i:i + m - 1] .= 10.0
        end
        i = i + 24
    end
end


v3 = zeros(length(vetor))
k = 1
i = 1
while i <= length(vetor)
    if cor[i] == 1
        if k < 24
            v3[i] = Ψ[1250,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[1250,k]
            k = 1
        end
    elseif cor[i] == 2
        if k < 24
            v3[i] = Ψ[2607,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[2607,k]
            k = 1
        end
    elseif cor[i] == 3
        if k < 24
            v3[i] = Ψ[3433,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[3433,k]
            k = 1
        end
    elseif cor[i] == 4
        if k < 24
            v3[i] = Ψ[3311,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[3311,k]
            k = 1
        end
    elseif cor[i] == 5
        if k < 24
            v3[i] = Ψ[1679,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[1679,k]
            k = 1
        end
    elseif cor[i] == 6
        if k < 24
            v3[i] = Ψ[1559,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[1559,k]
            k = 1
        end
    elseif cor[i] == 7
        if k < 24
            v3[i] = Ψ[1087,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[1087,k]
            k = 1
        end
    elseif cor[i] == 8
        if k < 24
            v3[i] = Ψ[1902,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[1902,k]
            k = 1
        end
    elseif cor[i] == 9
        if k < 24
            v3[i] = Ψ[630,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[630,k]
            k = 1
        end
    elseif cor[i] == 10
        if k < 24
            v3[i] = Ψ[2836,k]
            k = k + 1
        elseif k == 24
            v3[i] = Ψ[2836,k]
            k = 1
        end
    end
    i = i + 1
end



s = sum((v3[i]-vetor[i])^2 for i in 1:length(v3))
N = length(v3)
rmse = sqrt(s/N)



s = sum(abs(((v3[i]-vetor[i])^2)/vetor[i]) for i in 1:length(v3))
N = length(v3)
mape = (s/N)*100





@df fernandoHourwindall plot(:local_time, :wind_speed, title = "Wind Speed in fernando windar Farm - CE", 
xlab = "Hours", ylab = "Wind Speed (m/s)", legend = :outerbottomright, label = "Original Data")
fernandoHourwindall.cluster = v3
wind = @df fernandoHourwindall plot!(:local_time, :cluster, label = "Selected days")

savefig(wind, "comparacao_original_cluster_alegria.png")





s1 = scatter(vetor, v3, xlab = "Original data", ylab = "Representative days", legend = false)
s2 = scatter(v3, vetor, xlab = "Representative days", ylab = "Original data", legend = false)

residuals = vetor-v3
resid = plot(residuals, legend = false, title = "Residuals")
hist = histogram(residuals, legend = false, title = "Histogram of the residuals")


savefig(s1, "scatter1_alegria.png")
savefig(s2, "scatter2_alegria.png")
savefig(resid, "residuals_alegria.png")
savefig(hist, "hist_alegria.png")



dates = Vector(fernandoHourwindall.local_time)

dates_uniques = unique(dates)

vetor_dates = vcat(dates_uniques[630], dates_uniques[1087], dates_uniques[1250], dates_uniques[1559], 
dates_uniques[1679], dates_uniques[1902], dates_uniques[2607], dates_uniques[2836], 
dates_uniques[3311], dates_uniques[3433])



cluster[12]
v3[288:312]
for i in 1:5000
    if v3[i] == 0
        println(i)
    end
end
v3[5000]
