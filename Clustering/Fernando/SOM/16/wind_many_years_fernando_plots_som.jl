
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
#cd("./SOM")
path = pwd()

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

cluster1 = CSV.read("fernando_clust1_som_11K_16n.csv", DataFrame)
cluster1 = cluster1.x
cluster2 = CSV.read("fernando_clust2_som_11K_16n.csv", DataFrame)
cluster2 = cluster2.x
cluster3 = CSV.read("fernando_clust3_som_11K_16n.csv", DataFrame)
cluster3 = cluster3.x
cluster4 = CSV.read("fernando_clust4_som_11K_16n.csv", DataFrame)
cluster4 = cluster4.x
cluster5 = CSV.read("fernando_clust5_som_11K_16n.csv", DataFrame)
cluster5 = cluster5.x
cluster6 = CSV.read("fernando_clust6_som_11K_16n.csv", DataFrame)
cluster6 = cluster6.x
cluster7 = CSV.read("fernando_clust7_som_11K_16n.csv", DataFrame)
cluster7 = cluster7.x
cluster8 = CSV.read("fernando_clust8_som_11K_16n.csv", DataFrame)
cluster8 = cluster8.x
cluster9 = CSV.read("fernando_clust9_som_11K_16n.csv", DataFrame)
cluster9 = cluster9.x
cluster10 = CSV.read("fernando_clust10_som_11K_16n.csv", DataFrame)
cluster10 = cluster10.x
cluster11 = CSV.read("fernando_clust11_som_11K_16n.csv", DataFrame)
cluster11 = cluster11.x
#cluster12 = CSV.read("fernando_clust12_som_11K_16n.csv", DataFrame)
#cluster12 = cluster12.x

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

#=plot12 = plot(aa12', title = "Cluster 12", xlims = (1, 24), ylims = (0,14), xticks = 0:1:24, xlab = "Hours", 
ylab = "Wind Speed (m/s)", legend = false)=#


savefig(plot1, "K11-K1-fernando-wind-som-16n.png")
savefig(plot2, "K11-K2-fernando-wind-som-16n.png")
savefig(plot3, "K11-K3-fernando-wind-som.png")
savefig(plot4, "K11-K4-fernando-wind-som.png")
savefig(plot5, "K11-K5-fernando-wind-som.png")
savefig(plot6, "K11-K6-fernando-wind-som.png")
savefig(plot7, "K11-K7-fernando-wind-som.png")
savefig(plot8, "K11-K8-fernando-wind-som.png")
savefig(plot9, "K11-K9-fernando-wind-som.png")
savefig(plot10, "K11-K10-fernando-wind-som.png")
savefig(plot11, "K11-K11-fernando-wind-som.png")
#savefig(plot12, "K12-K12-fernando-wind-som.png")


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





# Each year

tail(fernandoHourwindall)
head(fernandoHourwindall)
fernandoHourwindall.wind_speed

fernandoHourwindall.YR = Dates.year.(fernandoHourwindall.local_time)
fernandoHourwindall.MON = Dates.month.(fernandoHourwindall.local_time)

Vetor = zeros(length(fernandoHourwindall.local_time))

for i in 1:length(fernandoHourwindall.wind_speed)
  Vetor[i] = fernandoHourwindall.wind_speed[length(fernandoHourwindall.wind_speed)+1-i]
end

#vec = Vetor

num_years = length(unique(fernandoHourwindall.YR))
num_month = length(unique(fernandoHourwindall.MON))

month1_2010 = 0
month1_2010_n = 0
month2_2010 = 0
month2_2010_n = 0
month3_2010 = 0
month3_2010_n = 0
month4_2010 = 0
month4_2010_n = 0
month5_2010 = 0
month5_2010_n = 0
month6_2010 = 0
month6_2010_n = 0
month7_2010 = 0
month7_2010_n = 0
month8_2010 = 0
month8_2010_n = 0
month9_2010 = 0
month9_2010_n = 0
month10_2010 = 0
month10_2010_n = 0
month11_2010 = 0
month11_2010_n = 0
month12_2010 = 0
month12_2010_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2010
        if fernandoHourwindall.MON[i] == 1
            month1_2010 = month1_2010 + fernandoHourwindall.wind_speed[i]
            month1_2010_n = month1_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2010 = month2_2010 + fernandoHourwindall.wind_speed[i]
            month2_2010_n = month2_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2010 = month3_2010 + fernandoHourwindall.wind_speed[i]
            month3_2010_n = month3_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2010 = month4_2010 + fernandoHourwindall.wind_speed[i]
            month4_2010_n = month4_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2010 = month5_2010 + fernandoHourwindall.wind_speed[i]
            month5_2010_n = month5_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2010 = month6_2010 + fernandoHourwindall.wind_speed[i]
            month6_2010_n = month6_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2010 = month7_2010 + fernandoHourwindall.wind_speed[i]
            month7_2010_n = month7_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2010 = month8_2010 + fernandoHourwindall.wind_speed[i]
            month8_2010_n = month8_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2010 = month9_2010 + fernandoHourwindall.wind_speed[i]
            month9_2010_n = month9_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2010 = month10_2010 + fernandoHourwindall.wind_speed[i]
            month10_2010_n = month10_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2010 = month11_2010 + fernandoHourwindall.wind_speed[i]
            month11_2010_n = month11_2010_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2010 = month12_2010 + fernandoHourwindall.wind_speed[i]
            month12_2010_n = month12_2010_n + 1
        end
    end
end

month1_2011 = 0
month1_2011_n = 0
month2_2011 = 0
month2_2011_n = 0
month3_2011 = 0
month3_2011_n = 0
month4_2011 = 0
month4_2011_n = 0
month5_2011 = 0
month5_2011_n = 0
month6_2011 = 0
month6_2011_n = 0
month7_2011 = 0
month7_2011_n = 0
month8_2011 = 0
month8_2011_n = 0
month9_2011 = 0
month9_2011_n = 0
month10_2011 = 0
month10_2011_n = 0
month11_2011 = 0
month11_2011_n = 0
month12_2011 = 0
month12_2011_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2011
        if fernandoHourwindall.MON[i] == 1
            month1_2011 = month1_2011 + fernandoHourwindall.wind_speed[i]
            month1_2011_n = month1_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2011 = month2_2011 + fernandoHourwindall.wind_speed[i]
            month2_2011_n = month2_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2011 = month3_2011 + fernandoHourwindall.wind_speed[i]
            month3_2011_n = month3_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2011 = month4_2011 + fernandoHourwindall.wind_speed[i]
            month4_2011_n = month4_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2011 = month5_2011 + fernandoHourwindall.wind_speed[i]
            month5_2011_n = month5_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2011 = month6_2011 + fernandoHourwindall.wind_speed[i]
            month6_2011_n = month6_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2011 = month7_2011 + fernandoHourwindall.wind_speed[i]
            month7_2011_n = month7_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2011 = month8_2011 + fernandoHourwindall.wind_speed[i]
            month8_2011_n = month8_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2011 = month9_2011 + fernandoHourwindall.wind_speed[i]
            month9_2011_n = month9_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2011 = month10_2011 + fernandoHourwindall.wind_speed[i]
            month10_2011_n = month10_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2011 = month11_2011 + fernandoHourwindall.wind_speed[i]
            month11_2011_n = month11_2011_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2011 = month12_2011 + fernandoHourwindall.wind_speed[i]
            month12_2011_n = month12_2011_n + 1
        end
    end
end

month1_2012 = 0
month1_2012_n = 0
month2_2012 = 0
month2_2012_n = 0
month3_2012 = 0
month3_2012_n = 0
month4_2012 = 0
month4_2012_n = 0
month5_2012 = 0
month5_2012_n = 0
month6_2012 = 0
month6_2012_n = 0
month7_2012 = 0
month7_2012_n = 0
month8_2012 = 0
month8_2012_n = 0
month9_2012 = 0
month9_2012_n = 0
month10_2012 = 0
month10_2012_n = 0
month11_2012 = 0
month11_2012_n = 0
month12_2012 = 0
month12_2012_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2012
        if fernandoHourwindall.MON[i] == 1
            month1_2012 = month1_2012 + fernandoHourwindall.wind_speed[i]
            month1_2012_n = month1_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2012 = month2_2012 + fernandoHourwindall.wind_speed[i]
            month2_2012_n = month2_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2012 = month3_2012 + fernandoHourwindall.wind_speed[i]
            month3_2012_n = month3_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2012 = month4_2012 + fernandoHourwindall.wind_speed[i]
            month4_2012_n = month4_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2012 = month5_2012 + fernandoHourwindall.wind_speed[i]
            month5_2012_n = month5_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2012 = month6_2012 + fernandoHourwindall.wind_speed[i]
            month6_2012_n = month6_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2012 = month7_2012 + fernandoHourwindall.wind_speed[i]
            month7_2012_n = month7_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2012 = month8_2012 + fernandoHourwindall.wind_speed[i]
            month8_2012_n = month8_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2012 = month9_2012 + fernandoHourwindall.wind_speed[i]
            month9_2012_n = month9_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2012 = month10_2012 + fernandoHourwindall.wind_speed[i]
            month10_2012_n = month10_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2012 = month11_2012 + fernandoHourwindall.wind_speed[i]
            month11_2012_n = month11_2012_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2012 = month12_2012 + fernandoHourwindall.wind_speed[i]
            month12_2012_n = month12_2012_n + 1
        end
    end
end

month1_2013 = 0
month1_2013_n = 0
month2_2013 = 0
month2_2013_n = 0
month3_2013 = 0
month3_2013_n = 0
month4_2013 = 0
month4_2013_n = 0
month5_2013 = 0
month5_2013_n = 0
month6_2013 = 0
month6_2013_n = 0
month7_2013 = 0
month7_2013_n = 0
month8_2013 = 0
month8_2013_n = 0
month9_2013 = 0
month9_2013_n = 0
month10_2013 = 0
month10_2013_n = 0
month11_2013 = 0
month11_2013_n = 0
month12_2013 = 0
month12_2013_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2013
        if fernandoHourwindall.MON[i] == 1
            month1_2013 = month1_2013 + fernandoHourwindall.wind_speed[i]
            month1_2013_n = month1_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2013 = month2_2013 + fernandoHourwindall.wind_speed[i]
            month2_2013_n = month2_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2013 = month3_2013 + fernandoHourwindall.wind_speed[i]
            month3_2013_n = month3_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2013 = month4_2013 + fernandoHourwindall.wind_speed[i]
            month4_2013_n = month4_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2013 = month5_2013 + fernandoHourwindall.wind_speed[i]
            month5_2013_n = month5_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2013 = month6_2013 + fernandoHourwindall.wind_speed[i]
            month6_2013_n = month6_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2013 = month7_2013 + fernandoHourwindall.wind_speed[i]
            month7_2013_n = month7_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2013 = month8_2013 + fernandoHourwindall.wind_speed[i]
            month8_2013_n = month8_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2013 = month9_2013 + fernandoHourwindall.wind_speed[i]
            month9_2013_n = month9_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2013 = month10_2013 + fernandoHourwindall.wind_speed[i]
            month10_2013_n = month10_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2013 = month11_2013 + fernandoHourwindall.wind_speed[i]
            month11_2013_n = month11_2013_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2013 = month12_2013 + fernandoHourwindall.wind_speed[i]
            month12_2013_n = month12_2013_n + 1
        end
    end
end

month1_2014 = 0
month1_2014_n = 0
month2_2014 = 0
month2_2014_n = 0
month3_2014 = 0
month3_2014_n = 0
month4_2014 = 0
month4_2014_n = 0
month5_2014 = 0
month5_2014_n = 0
month6_2014 = 0
month6_2014_n = 0
month7_2014 = 0
month7_2014_n = 0
month8_2014 = 0
month8_2014_n = 0
month9_2014 = 0
month9_2014_n = 0
month10_2014 = 0
month10_2014_n = 0
month11_2014 = 0
month11_2014_n = 0
month12_2014 = 0
month12_2014_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2014
        if fernandoHourwindall.MON[i] == 1
            month1_2014 = month1_2014 + fernandoHourwindall.wind_speed[i]
            month1_2014_n = month1_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2014 = month2_2014 + fernandoHourwindall.wind_speed[i]
            month2_2014_n = month2_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2014 = month3_2014 + fernandoHourwindall.wind_speed[i]
            month3_2014_n = month3_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2014 = month4_2014 + fernandoHourwindall.wind_speed[i]
            month4_2014_n = month4_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2014 = month5_2014 + fernandoHourwindall.wind_speed[i]
            month5_2014_n = month5_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2014 = month6_2014 + fernandoHourwindall.wind_speed[i]
            month6_2014_n = month6_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2014 = month7_2014 + fernandoHourwindall.wind_speed[i]
            month7_2014_n = month7_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2014 = month8_2014 + fernandoHourwindall.wind_speed[i]
            month8_2014_n = month8_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2014 = month9_2014 + fernandoHourwindall.wind_speed[i]
            month9_2014_n = month9_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2014 = month10_2014 + fernandoHourwindall.wind_speed[i]
            month10_2014_n = month10_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2014 = month11_2014 + fernandoHourwindall.wind_speed[i]
            month11_2014_n = month11_2014_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2014 = month12_2014 + fernandoHourwindall.wind_speed[i]
            month12_2014_n = month12_2014_n + 1
        end
    end
end

month1_2015 = 0
month1_2015_n = 0
month2_2015 = 0
month2_2015_n = 0
month3_2015 = 0
month3_2015_n = 0
month4_2015 = 0
month4_2015_n = 0
month5_2015 = 0
month5_2015_n = 0
month6_2015 = 0
month6_2015_n = 0
month7_2015 = 0
month7_2015_n = 0
month8_2015 = 0
month8_2015_n = 0
month9_2015 = 0
month9_2015_n = 0
month10_2015 = 0
month10_2015_n = 0
month11_2015 = 0
month11_2015_n = 0
month12_2015 = 0
month12_2015_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2015
        if fernandoHourwindall.MON[i] == 1
            month1_2015 = month1_2015 + fernandoHourwindall.wind_speed[i]
            month1_2015_n = month1_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2015 = month2_2015 + fernandoHourwindall.wind_speed[i]
            month2_2015_n = month2_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2015 = month3_2015 + fernandoHourwindall.wind_speed[i]
            month3_2015_n = month3_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2015 = month4_2015 + fernandoHourwindall.wind_speed[i]
            month4_2015_n = month4_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2015 = month5_2015 + fernandoHourwindall.wind_speed[i]
            month5_2015_n = month5_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2015 = month6_2015 + fernandoHourwindall.wind_speed[i]
            month6_2015_n = month6_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2015 = month7_2015 + fernandoHourwindall.wind_speed[i]
            month7_2015_n = month7_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2015 = month8_2015 + fernandoHourwindall.wind_speed[i]
            month8_2015_n = month8_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2015 = month9_2015 + fernandoHourwindall.wind_speed[i]
            month9_2015_n = month9_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2015 = month10_2015 + fernandoHourwindall.wind_speed[i]
            month10_2015_n = month10_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2015 = month11_2015 + fernandoHourwindall.wind_speed[i]
            month11_2015_n = month11_2015_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2015 = month12_2015 + fernandoHourwindall.wind_speed[i]
            month12_2015_n = month12_2015_n + 1
        end
    end
end


month1_2016 = 0
month1_2016_n = 0
month2_2016 = 0
month2_2016_n = 0
month3_2016 = 0
month3_2016_n = 0
month4_2016 = 0
month4_2016_n = 0
month5_2016 = 0
month5_2016_n = 0
month6_2016 = 0
month6_2016_n = 0
month7_2016 = 0
month7_2016_n = 0
month8_2016 = 0
month8_2016_n = 0
month9_2016 = 0
month9_2016_n = 0
month10_2016 = 0
month10_2016_n = 0
month11_2016 = 0
month11_2016_n = 0
month12_2016 = 0
month12_2016_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2016
        if fernandoHourwindall.MON[i] == 1
            month1_2016 = month1_2016 + fernandoHourwindall.wind_speed[i]
            month1_2016_n = month1_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2016 = month2_2016 + fernandoHourwindall.wind_speed[i]
            month2_2016_n = month2_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2016 = month3_2016 + fernandoHourwindall.wind_speed[i]
            month3_2016_n = month3_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2016 = month4_2016 + fernandoHourwindall.wind_speed[i]
            month4_2016_n = month4_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2016 = month5_2016 + fernandoHourwindall.wind_speed[i]
            month5_2016_n = month5_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2016 = month6_2016 + fernandoHourwindall.wind_speed[i]
            month6_2016_n = month6_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2016 = month7_2016 + fernandoHourwindall.wind_speed[i]
            month7_2016_n = month7_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2016 = month8_2016 + fernandoHourwindall.wind_speed[i]
            month8_2016_n = month8_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2016 = month9_2016 + fernandoHourwindall.wind_speed[i]
            month9_2016_n = month9_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2016 = month10_2016 + fernandoHourwindall.wind_speed[i]
            month10_2016_n = month10_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2016 = month11_2016 + fernandoHourwindall.wind_speed[i]
            month11_2016_n = month11_2016_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2016 = month12_2016 + fernandoHourwindall.wind_speed[i]
            month12_2016_n = month12_2016_n + 1
        end
    end
end


month1_2017 = 0
month1_2017_n = 0
month2_2017 = 0
month2_2017_n = 0
month3_2017 = 0
month3_2017_n = 0
month4_2017 = 0
month4_2017_n = 0
month5_2017 = 0
month5_2017_n = 0
month6_2017 = 0
month6_2017_n = 0
month7_2017 = 0
month7_2017_n = 0
month8_2017 = 0
month8_2017_n = 0
month9_2017 = 0
month9_2017_n = 0
month10_2017 = 0
month10_2017_n = 0
month11_2017 = 0
month11_2017_n = 0
month12_2017 = 0
month12_2017_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2017
        if fernandoHourwindall.MON[i] == 1
            month1_2017 = month1_2017 + fernandoHourwindall.wind_speed[i]
            month1_2017_n = month1_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2017 = month2_2017 + fernandoHourwindall.wind_speed[i]
            month2_2017_n = month2_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2017 = month3_2017 + fernandoHourwindall.wind_speed[i]
            month3_2017_n = month3_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2017 = month4_2017 + fernandoHourwindall.wind_speed[i]
            month4_2017_n = month4_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2017 = month5_2017 + fernandoHourwindall.wind_speed[i]
            month5_2017_n = month5_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2017 = month6_2017 + fernandoHourwindall.wind_speed[i]
            month6_2017_n = month6_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2017 = month7_2017 + fernandoHourwindall.wind_speed[i]
            month7_2017_n = month7_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2017 = month8_2017 + fernandoHourwindall.wind_speed[i]
            month8_2017_n = month8_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2017 = month9_2017 + fernandoHourwindall.wind_speed[i]
            month9_2017_n = month9_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2017 = month10_2017 + fernandoHourwindall.wind_speed[i]
            month10_2017_n = month10_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2017 = month11_2017 + fernandoHourwindall.wind_speed[i]
            month11_2017_n = month11_2017_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2017 = month12_2017 + fernandoHourwindall.wind_speed[i]
            month12_2017_n = month12_2017_n + 1
        end
    end
end

month1_2018 = 0
month1_2018_n = 0
month2_2018 = 0
month2_2018_n = 0
month3_2018 = 0
month3_2018_n = 0
month4_2018 = 0
month4_2018_n = 0
month5_2018 = 0
month5_2018_n = 0
month6_2018 = 0
month6_2018_n = 0
month7_2018 = 0
month7_2018_n = 0
month8_2018 = 0
month8_2018_n = 0
month9_2018 = 0
month9_2018_n = 0
month10_2018 = 0
month10_2018_n = 0
month11_2018 = 0
month11_2018_n = 0
month12_2018 = 0
month12_2018_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2018
        if fernandoHourwindall.MON[i] == 1
            month1_2018 = month1_2018 + fernandoHourwindall.wind_speed[i]
            month1_2018_n = month1_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2018 = month2_2018 + fernandoHourwindall.wind_speed[i]
            month2_2018_n = month2_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2018 = month3_2018 + fernandoHourwindall.wind_speed[i]
            month3_2018_n = month3_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2018 = month4_2018 + fernandoHourwindall.wind_speed[i]
            month4_2018_n = month4_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2018 = month5_2018 + fernandoHourwindall.wind_speed[i]
            month5_2018_n = month5_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2018 = month6_2018 + fernandoHourwindall.wind_speed[i]
            month6_2018_n = month6_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2018 = month7_2018 + fernandoHourwindall.wind_speed[i]
            month7_2018_n = month7_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2018 = month8_2018 + fernandoHourwindall.wind_speed[i]
            month8_2018_n = month8_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2018 = month9_2018 + fernandoHourwindall.wind_speed[i]
            month9_2018_n = month9_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2018 = month10_2018 + fernandoHourwindall.wind_speed[i]
            month10_2018_n = month10_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2018 = month11_2018 + fernandoHourwindall.wind_speed[i]
            month11_2018_n = month11_2018_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2018 = month12_2018 + fernandoHourwindall.wind_speed[i]
            month12_2018_n = month12_2018_n + 1
        end
    end
end

month1_2019 = 0
month1_2019_n = 0
month2_2019 = 0
month2_2019_n = 0
month3_2019 = 0
month3_2019_n = 0
month4_2019 = 0
month4_2019_n = 0
month5_2019 = 0
month5_2019_n = 0
month6_2019 = 0
month6_2019_n = 0
month7_2019 = 0
month7_2019_n = 0
month8_2019 = 0
month8_2019_n = 0
month9_2019 = 0
month9_2019_n = 0
month10_2019 = 0
month10_2019_n = 0
month11_2019 = 0
month11_2019_n = 0
month12_2019 = 0
month12_2019_n = 0
for i in 1:length(fernandoHourwindall.YR)
    if fernandoHourwindall.YR[i] == 2019
        if fernandoHourwindall.MON[i] == 1
            month1_2019 = month1_2019 + fernandoHourwindall.wind_speed[i]
            month1_2019_n = month1_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 2
            month2_2019 = month2_2019 + fernandoHourwindall.wind_speed[i]
            month2_2019_n = month2_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 3
            month3_2019 = month3_2019 + fernandoHourwindall.wind_speed[i]
            month3_2019_n = month3_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 4
            month4_2019 = month4_2019 + fernandoHourwindall.wind_speed[i]
            month4_2019_n = month4_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 5
            month5_2019 = month5_2019 + fernandoHourwindall.wind_speed[i]
            month5_2019_n = month5_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 6
            month6_2019 = month6_2019 + fernandoHourwindall.wind_speed[i]
            month6_2019_n = month6_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 7
            month7_2019 = month7_2019 + fernandoHourwindall.wind_speed[i]
            month7_2019_n = month7_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 8
            month8_2019 = month8_2019 + fernandoHourwindall.wind_speed[i]
            month8_2019_n = month8_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 9
            month9_2019 = month9_2019 + fernandoHourwindall.wind_speed[i]
            month9_2019_n = month9_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 10
            month10_2019 = month10_2019 + fernandoHourwindall.wind_speed[i]
            month10_2019_n = month10_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 11
            month11_2019 = month11_2019 + fernandoHourwindall.wind_speed[i]
            month11_2019_n = month11_2019_n + 1
        elseif fernandoHourwindall.MON[i] == 12
            month12_2019 = month12_2019 + fernandoHourwindall.wind_speed[i]
            month12_2019_n = month12_2019_n + 1
        end
    end
end

dfmes = DataFrame(mes = 
[month1_2010/month1_2010_n, month2_2010/month2_2010_n, month3_2010/month3_2010_n, month4_2010/month4_2010_n, 
month5_2010/month5_2010_n, month6_2010/month6_2010_n, month7_2010/month7_2010_n, month8_2010/month8_2010_n, 
month9_2010/month9_2010_n, month10_2010/month10_2010_n, month11_2010/month11_2010_n, month12_2010/month12_2010_n,
    
month1_2011/month1_2011_n, month2_2011/month2_2011_n, month3_2011/month3_2011_n, month4_2011/month4_2011_n, 
month5_2011/month5_2011_n, month6_2011/month6_2011_n, month7_2011/month7_2011_n, month8_2011/month8_2011_n, 
month9_2011/month9_2011_n, month10_2011/month10_2011_n, month11_2011/month11_2011_n, month12_2011/month12_2011_n,
    
month1_2012/month1_2012_n, month2_2012/month2_2012_n, month3_2012/month3_2012_n, month4_2012/month4_2012_n, 
month5_2012/month5_2012_n, month6_2012/month6_2012_n, month7_2012/month7_2012_n, month8_2012/month8_2012_n, 
month9_2012/month9_2012_n, month10_2012/month10_2012_n, month11_2012/month11_2012_n, month12_2012/month12_2012_n,
    
month1_2013/month1_2013_n, month2_2013/month2_2013_n, month3_2013/month3_2013_n, month4_2013/month4_2013_n, 
month5_2013/month5_2013_n, month6_2013/month6_2013_n, month7_2013/month7_2013_n, month8_2013/month8_2013_n, 
month9_2013/month9_2013_n, month10_2013/month10_2013_n, month11_2013/month11_2013_n, month12_2013/month12_2013_n,
    
month1_2014/month1_2014_n, month2_2014/month2_2014_n, month3_2014/month3_2014_n, month4_2014/month4_2014_n, 
month5_2014/month5_2014_n, month6_2014/month6_2014_n, month7_2014/month7_2014_n, month8_2014/month8_2014_n, 
month9_2014/month9_2014_n, month10_2014/month10_2014_n, month11_2014/month11_2014_n, month12_2014/month12_2014_n, 

month1_2015/month1_2015_n, month2_2015/month2_2015_n, month3_2015/month3_2015_n, month4_2015/month4_2015_n, 
month5_2015/month5_2015_n, month6_2015/month6_2015_n, month7_2015/month7_2015_n, month8_2015/month8_2015_n, 
month9_2015/month9_2015_n, month10_2015/month10_2015_n, month11_2015/month11_2015_n, month12_2015/month12_2015_n, 

month1_2016/month1_2016_n, month2_2016/month2_2016_n, month3_2016/month3_2016_n, month4_2016/month4_2016_n, 
month5_2016/month5_2016_n, month6_2016/month6_2016_n, month7_2016/month7_2016_n, month8_2016/month8_2016_n, 
month9_2016/month9_2016_n, month10_2016/month10_2016_n, month11_2016/month11_2016_n, month12_2016/month12_2016_n, 

month1_2017/month1_2017_n, month2_2017/month2_2017_n, month3_2017/month3_2017_n, month4_2017/month4_2017_n, 
month5_2017/month5_2017_n, month6_2017/month6_2017_n, month7_2017/month7_2017_n, month8_2017/month8_2017_n, 
month9_2017/month9_2017_n, month10_2017/month10_2017_n, month11_2017/month11_2017_n, month12_2017/month12_2017_n, 

month1_2018/month1_2018_n, month2_2018/month2_2018_n, month3_2018/month3_2018_n, month4_2018/month4_2018_n, 
month5_2018/month5_2018_n, month6_2018/month6_2018_n, month7_2018/month7_2018_n, month8_2018/month8_2018_n, 
month9_2018/month9_2018_n, month10_2018/month10_2018_n, month11_2018/month11_2018_n, month12_2018/month12_2018_n,

month1_2019/month1_2019_n, month2_2019/month2_2019_n, month3_2019/month3_2019_n, month4_2019/month4_2019_n, 
month5_2019/month5_2019_n, month6_2019/month6_2019_n, month7_2019/month7_2019_n, month8_2019/month8_2019_n, 
month9_2019/month9_2019_n, month10_2019/month10_2019_n, month11_2019/month11_2019_n, month12_2019/month12_2019_n])
vec = dfmes.mes

function test()
    matriz = zeros(num_years, num_month)
    m = 1
    a = 1
    for i in 1:length(vec)
        matriz[a, m] = vec[i]
        m = m + 1
        if m >= 13
            m = 1
            a = a + 1
        end
    end
    return matriz
end

m = test()

plot2 = @df fernandoHourwindall plot(m', title = "Mean speed by years in Alegria 1 e 2 - RN", 
xlims = (1, 15), xticks = 0:1:12, xlab = "Months", ylab = "Wind Speed (m/s)", 
palette = ["deeppink", "green", "blue", "coral", "yellow", "red", "purple", "black", "aquamarine", "darkorange2"],
label = ["2010" "2011" "2012" "2013" "2014" "2015" "2016" "2017" "2018" "2019"])

savefig(plot2, "por ano alegria.png")



# Filtro HP

using HPFilter
using Statistics

dates = fernandoHourwindall.local_time
intervalo_dates = dates[1:15000:length(dates)]

graf = plot(dates, fernandoHourwindall.wind_speed, xlabel = "Hours", ylabel = "Wind Speed (m/s)", 
title = "Power Generation in Alegria 1 and 2 Wind Farms - RN", label = "Original", 
xticks = (intervalo_dates))

T = length(fernandoHourwindall.wind_speed)

y = Vector(fernandoHourwindall.wind_speed)

hp = HP(y, 10_000_000_000)

grafQua10000000000 = plot(dates, fernandoHourwindall.wind_speed, xlabel = "Hours", 
ylabel = "Wind Speed (m/s)", title = "Serie with H-P Filter e lambda = 10000000000", 
label = "Original", xticks = (intervalo_dates))
plot!(dates, hp, lw = 2, color = "coral", label = "Filtered")
savefig(grafQua10000000000, "MacaubasQua10000000000.png")

unique(hp)

hppsi = CreatePsi(hp)

med_day = mean(hppsi, dims=2)

function CreatePsi2(vec)
    Nh = 24
    Ni = (365*10)+1
    Nc = 1
    matriz = zeros((Ni), Nc * Nh)
    m = 1
    a = 1
    for i in 1:length(vec)
        matriz[a, 1] = vec[i]
        matriz[a, 2] = vec[i]
        matriz[a, 3] = vec[i]
        matriz[a, 4] = vec[i]
        matriz[a, 5] = vec[i]
        matriz[a, 6] = vec[i]
        matriz[a, 7] = vec[i]
        matriz[a, 8] = vec[i]
        matriz[a, 9] = vec[i]
        matriz[a, 10] = vec[i]
        matriz[a, 11] = vec[i]
        matriz[a, 12] = vec[i]
        matriz[a, 13] = vec[i]
        matriz[a, 14] = vec[i]
        matriz[a, 15] = vec[i]
        matriz[a, 16] = vec[i]
        matriz[a, 17] = vec[i]
        matriz[a, 18] = vec[i]
        matriz[a, 19] = vec[i]
        matriz[a, 20] = vec[i]
        matriz[a, 21] = vec[i]
        matriz[a, 22] = vec[i]
        matriz[a, 23] = vec[i]
        matriz[a, 24] = vec[i]
        a = a + 1
    end
    return matriz
end

med_day_hour = CreatePsi2(med_day)
med_day_hour_vec = vec(med_day_hour')

plot!(dates, med_day_hour_vec, lw = 1, color = "deeppink", label = "Day Mean")



















cluster[12]
v3[288:312]
for i in 1:5000
    if v3[i] == 0
        println(i)
    end
end
v3[5000]
