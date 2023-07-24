
##################################### Importing the packages #####################################

using JuMP
using HiGHS #Gurobi
using CSV
using DataFrames

##################################### Opening the data #####################################

path = pwd()

buses = CSV.read("./data/bus.csv", DataFrame)
generators = CSV.read("./data/gen.csv", DataFrame)
branchs = CSV.read("./data/branch.csv", DataFrame)
cluster1 = CSV.read("./data/janauba_clust1_som_12K_16n.csv", DataFrame)
cluster2 = CSV.read("./data/janauba_clust2_som_12K_16n.csv", DataFrame)
cluster3 = CSV.read("./data/janauba_clust3_som_12K_16n.csv", DataFrame)
cluster4 = CSV.read("./data/janauba_clust4_som_12K_16n.csv", DataFrame)
cluster5 = CSV.read("./data/janauba_clust5_som_12K_16n.csv", DataFrame)
cluster6 = CSV.read("./data/janauba_clust6_som_12K_16n.csv", DataFrame)
cluster7 = CSV.read("./data/janauba_clust7_som_12K_16n.csv", DataFrame)
cluster8 = CSV.read("./data/janauba_clust8_som_12K_16n.csv", DataFrame)
cluster9 = CSV.read("./data/janauba_clust9_som_12K_16n.csv", DataFrame)
cluster10 = CSV.read("./data/janauba_clust10_som_12K_16n.csv", DataFrame)
cluster11 = CSV.read("./data/janauba_clust11_som_12K_16n.csv", DataFrame)
cluster12 = CSV.read("./data/janauba_clust12_som_12K_16n.csv", DataFrame)
cluster = CSV.read("./data/janauba_clus_som_12K_16n.csv", DataFrame)
ts = CSV.read("./data/DAY_AHEAD_regional_Load.csv", DataFrame) # 365 days -> no 29/02/2020
centres = CSV.read("./data/som_medoids_12.csv", DataFrame)
janauba = CSV.read("./data/JANAUBA_MG_UCT.csv", DataFrame)
z1 = CSV.read("./data/genZona1.csv", DataFrame)
z2 = CSV.read("./data/genZona2.csv", DataFrame)
z3 = CSV.read("./data/genZona3.csv", DataFrame)

##################################### Creating the zones #####################################

# How many % of MW each bus represent each zone total MW load
zone1 = buses[1:24,16]
zone2 = buses[25:48,16]
zone3 = buses[49:73,16]
sum(zone3)
busList = collect(buses[!,1])

##################################### Number of generators #####################################

unique(generators[!,7])
r = 0
ws = 0
reserve = 0
for i in 1:length(generators[!,7])
    if generators[i,7] == "Wind" || generators[i,7] == "Solar"
        global ws = ws + 1
    elseif generators[i,7] == "Hydro"
        global r = r + 1
    elseif generators[i,7] == "Storage" || generators[i,7] == "Sync_Cond"
        global reserve = reserve + 1
    end
end

##################################### Sets #####################################

Ng = length(generators[!,7])-r-reserve-ws # Number of conventional generators
Nn = r # Number of renewable generators
Nh = length(ts[!,1]) #24 # Number of hours of a day (24 hours)
Nb = length(buses[!,1]) # Number of nodes
CC = 4 # Set of candidate circuits
C = length(branchs[!,3]) # Set of existing circuits

##################################### Creating the centres #####################################

centres = Matrix(centres[:,2:end])

maxJanauba = maximum(janauba[!,4])
minJanauba = minimum(janauba[!,4])

janaubac1 = centres[1,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac2 = centres[2,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac3 = centres[3,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac4 = centres[4,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac5 = centres[5,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac6 = centres[6,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac7 = centres[7,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac8 = centres[8,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac9 = centres[9,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac10 = centres[10,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac11 = centres[11,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3
janaubac12 = centres[12,:] .* (maxJanauba - minJanauba) .+ minJanauba # Zone 3

janaubac1 = janaubac1 ./ maxJanauba # Zone 3
janaubac2 = janaubac2 ./ maxJanauba # Zone 3
janaubac3 = janaubac3 ./ maxJanauba # Zone 3
janaubac4 = janaubac4 ./ maxJanauba # Zone 3
janaubac5 = janaubac5 ./ maxJanauba # Zone 3
janaubac6 = janaubac6 ./ maxJanauba # Zone 3
janaubac7 = janaubac7 ./ maxJanauba # Zone 3
janaubac8 = janaubac8 ./ maxJanauba # Zone 3
janaubac9 = janaubac9 ./ maxJanauba # Zone 3
janaubac10 = janaubac10 ./ maxJanauba # Zone 3
janaubac11 = janaubac11 ./ maxJanauba # Zone 3
janaubac12 = janaubac12 ./ maxJanauba # Zone 3

println("zonas")

#Solar
solarZ3 = z3[1:34,8]
UniqZ3 = unique(z3[1:34,2])
solarZ3bus = Array{Float64}(undef, length(UniqZ3), 1)
for (s, sol) in enumerate(z3[1:34,2])
    if z3[s,2] == UniqZ3[1]
        solarZ3bus[1] = solarZ3bus[1] + solarZ3[s]
    elseif z3[s,2] == UniqZ3[2]
        solarZ3bus[2] = solarZ3bus[2] + solarZ3[s]
    elseif z3[s,2] == UniqZ3[3]
        solarZ3bus[3] = solarZ3bus[3] + solarZ3[s]
    elseif z3[s,2] == UniqZ3[4]
        solarZ3bus[4] = solarZ3bus[4] + solarZ3[s]
    elseif z3[s,2] == UniqZ3[5]
        solarZ3bus[5] = solarZ3bus[5] + solarZ3[s]
    elseif z3[s,2] == UniqZ3[6]
        solarZ3bus[6] = solarZ3bus[6] + solarZ3[s]
    elseif z3[s,2] == UniqZ3[7]
        solarZ3bus[7] = solarZ3bus[7] + solarZ3[s]
    elseif z3[s,2] == UniqZ3[8]
        solarZ3bus[8] = solarZ3bus[8] + solarZ3[s]
    end
end

solarZone3c1 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c1 = hcat(solarZone3c1, janaubac1.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c1, string.(UniqZ3))

solarZone3c2 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c2 = hcat(solarZone3c2, janaubac2.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c2, string.(UniqZ3))

solarZone3c3 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c3 = hcat(solarZone3c3, janaubac3.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c3, string.(UniqZ3))

solarZone3c4 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c4 = hcat(solarZone3c4, janaubac4.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c4, string.(UniqZ3))

solarZone3c5 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c5 = hcat(solarZone3c5, janaubac5.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c5, string.(UniqZ3))

solarZone3c6 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c6 = hcat(solarZone3c6, janaubac6.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c6, string.(UniqZ3))

solarZone3c7 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c7 = hcat(solarZone3c7, janaubac7.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c7, string.(UniqZ3))

solarZone3c8 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c8 = hcat(solarZone3c8, janaubac8.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c8, string.(UniqZ3))

solarZone3c9 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c9 = hcat(solarZone3c9, janaubac9.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c9, string.(UniqZ3))

solarZone3c10 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c10 = hcat(solarZone3c10, janaubac10.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c10, string.(UniqZ3))

solarZone3c11 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c11 = hcat(solarZone3c11, janaubac11.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c11, string.(UniqZ3))

solarZone3c12 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3c12 = hcat(solarZone3c12, janaubac12.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3c12, string.(UniqZ3))

Z3C1 = solarZone3c1
Z3C2 = solarZone3c2
Z3C3 = solarZone3c3
Z3C4 = solarZone3c4
Z3C5 = solarZone3c5
Z3C6 = solarZone3c6
Z3C7 = solarZone3c7
Z3C8 = solarZone3c8
Z3C9 = solarZone3c9
Z3C10 = solarZone3c10
Z3C11 = solarZone3c11
Z3C12 = solarZone3c12

println("cluster_zonas")

begin2020 = 14601
end2020 = 14965

YearZ3 = DataFrame()
for (i,k) in enumerate(cluster[begin2020:end2020,2])
    println(k)
    if k == 1
        global YearZ3 = vcat(YearZ3, Z3C1)
    elseif k == 2
        global YearZ3 = vcat(YearZ3, Z3C2)
    elseif k == 3
        global YearZ3 = vcat(YearZ3, Z3C3)
    elseif k == 4
        global YearZ3 = vcat(YearZ3, Z3C4)
    elseif k == 5
        global YearZ3 = vcat(YearZ3, Z3C5)
    elseif k == 6
        global YearZ3 = vcat(YearZ3, Z3C6)
    elseif k == 7
        global YearZ3 = vcat(YearZ3, Z3C7)
    elseif k == 8
        global YearZ3 = vcat(YearZ3, Z3C8)
    elseif k == 9
        global YearZ3 = vcat(YearZ3, Z3C9)
    elseif k == 10
        global YearZ3 = vcat(YearZ3, Z3C10)
    elseif k == 11
        global YearZ3 = vcat(YearZ3, Z3C11)
    elseif k == 12
        global YearZ3 = vcat(YearZ3, Z3C12)
    end
end

println("aaaaaaaaaa")

busZ3Uniq = unique(z3[!,2])
ΩZ3 = Array{Union{Nothing,Float64}}(nothing, 73, Nh)
for (i, index) in enumerate(busZ3Uniq)
    println(index)
    j = 1
    for k in 1:length(names(YearZ3))
        nome = names(YearZ3)[k][1:3]
        if parse(Int64, nome) == index
            row = findall(x->x==index, busList)[1]
            global ΩZ3[row,:] = YearZ3[!,k]
            #j = j + 1
        end
    end
end
for nc in 1:Nh
    for nr in 1:73
        if ΩZ3[nr,nc] == nothing
            ΩZ3[nr,nc] = 0
        end
    end
end

ΩsolEol = ΩZ3

println("dados")

##################################### Parameters #####################################

ρ = [] #[h,day]
for i in begin2020:end2020
    global ρ = vcat(ρ, repeat([1/(365*24)], outer = 24))
end
sum(ρ)

CG = generators[!,30]./0.293071 # The cost of generator g
Cens = 500 # The cost of non-served energy ($/MW)
α = 1
Cost = [5000 5000 2500 2500] * α # The anualized cost of investing on a circuit c (cost between Zones is 50% greater than cost only in Zone 3)
#Cost = [0 0 2500 2500]

d = DataFrame()
for (i, bus) in enumerate(zone1)
    global d = hcat(d, ts[!,5] .* bus, makeunique=true)
end
for (i, bus) in enumerate(zone2)
    global d = hcat(d, ts[!,6] .* bus, makeunique=true)
end
for (i, bus) in enumerate(zone3)
    global d = hcat(d, ts[!,7] .* bus, makeunique=true)
end # [h,b] # The demand at node b at hour h

d = Matrix(d)
d = d'

diff = d - ΩsolEol
minimum(sum(diff, dims = 1))

println("demand")

Ωend = branchs[!,3] # Existing circuits with b as an end node
Ωstart = branchs[!,2] # Existing circuits with b as an start node
Ωstartcc = [320, 218, 304, 322] # Candidate circuits with b as an start node
Ωendcc = [119, 319, 308, 314] # Candidate circuits with b as an end node

busIDrenew = []
genNameRenew = []
busIDconv = []
genNameConv = []
CGconv = []
CGrenew = []
Rmaxrenew = []
Rmaxconv = []
pnmax = [] # Maximum generation of renewable generator n
pgmax = [] # Maximum generation of conventional generator g
for i in 1:length(generators[!,2])
    if generators[i,7] == "Hydro" # || generators[i,7] == "Wind" || generators[i,7] == "Solar"
        global busIDrenew = vcat(busIDrenew, generators[i,2])
        global genNameRenew = vcat(genNameRenew, generators[i,3])
        global CGrenew = vcat(CGrenew, generators[i,30])
        global Rmaxrenew = vcat(Rmaxrenew, generators[i,17]*60) # Maximum ramp
        global pnmax = vcat(pnmax, generators[i,11])
    elseif generators[i,7] == "Oil" || generators[i,7] == "Coal" || generators[i,7] == "NG" || generators[i,7] == "Nuclear"
        global busIDconv = vcat(busIDconv, generators[i,2])
        global genNameConv = vcat(genNameConv, generators[i,3])
        global CGconv = vcat(CGconv, generators[i,30])
        global Rmaxconv = vcat(Rmaxconv, generators[i,17]*60) # Maximum ramp
        global pgmax = vcat(pgmax, generators[i,11])
    end
end
CGconv = (CGconv./0.293071) # The cost of conventional generator g
CGrenew = CGrenew./0.293071 # The cost of renewable generator g
ρ.*minimum(CGconv)
ρ.*maximum(CGconv)

busIDconvUniq = unique(busIDconv)
busIDrenewUniq = unique(busIDrenew)
ncols = 15 # count(x->x==mode(busIDrenew),busIDrenew)
ΩGrenew = Array{Union{Nothing,Int64}}(nothing, 73, ncols)
for (i, index) in enumerate(busIDrenewUniq)
    println(i)
    j = 1
    for k in 1:length(busIDrenew)
        if busIDrenew[k] == index
            row = findall(x->x==index, busList)[1]
            global ΩGrenew[row,j] = k # genNameRenew[k]
            j = j + 1
        end
    end
end
for nc in 1:ncols
    for nr in 1:73
        if ΩGrenew[nr,nc] == nothing
            ΩGrenew[nr,nc] = -1
        end
    end
end

ncols = 8 # count(x->x==mode(busIDconv),busIDconv)
ΩGconv = Array{Union{Nothing,Int64}}(nothing, 73, ncols)
for (i, index) in enumerate(busIDconvUniq)
    println(i)
    j = 1
    for k in 1:length(busIDconv)
        if busIDconv[k] == index
            row = findall(x->x==index, busList)[1]
            global ΩGconv[row,j] = k #genNameConv[k]
            j = j + 1
        end
    end
end
for nc in 1:ncols
    for nr in 1:73
        if ΩGconv[nr,nc] == nothing
            ΩGconv[nr,nc] = -1
        end
    end
end

Ωend
Ωstart
ncolsStart = 4 # count(x->x==mode(Ωstart),Ωstart)
ncolsEnd = 5 # count(x->x==mode(Ωend),Ωend)

ΩendUniq = unique(Ωend)
ΩstartUniq = unique(Ωstart)
ΩstartC = Array{Union{Nothing,Int64}}(nothing, 73, ncolsStart)
ΩstartCC = Array{Union{Nothing,Int64}}(nothing, 73, ncolsStart)
for (i, index) in enumerate(ΩstartUniq)
    println(index)
    j = 1
    for k in 1:length(Ωstart)
        if Ωstart[k] == index
            row = findall(x->x==index, busList)[1]
            global ΩstartC[row,j] = k #1
            j = j + 1
        end
    end
end
for nc in 1:ncolsStart
    for nr in 1:73
        if ΩstartC[nr,nc] == nothing
            ΩstartC[nr,nc] = -1
        end
    end
end
for (i, index) in enumerate(ΩstartUniq)
    println(index)
    j = 1
    for k in 1:length(Ωstartcc)
        if Ωstartcc[k] == index
            row = findall(x->x==index, busList)[1]
            global ΩstartCC[row,j] = k #1
            j = j + 1
        end
    end
end
for nc in 1:ncolsStart
    for nr in 1:73
        if ΩstartCC[nr,nc] == nothing
            ΩstartCC[nr,nc] = -1
        end
    end
end

ΩendC = Array{Union{Nothing,Int64}}(nothing, 73, ncolsEnd)
ΩendCC = Array{Union{Nothing,Int64}}(nothing, 73, ncolsEnd)
for (i, index) in enumerate(ΩendUniq)
    println(index)
    j = 1
    for k in 1:length(Ωend)
        if Ωend[k] == index
            row = findall(x->x==index, busList)[1]
            global ΩendC[row,j] = k #1
            j = j + 1
        end
    end
end
for nc in 1:ncolsEnd
    for nr in 1:73
        if ΩendC[nr,nc] == nothing
            ΩendC[nr,nc] = -1
        end
    end
end
for (i, index) in enumerate(ΩendUniq)
    println(index)
    j = 1
    for k in 1:length(Ωendcc)
        if Ωendcc[k] == index
            row = findall(x->x==index, busList)[1]
            global ΩendCC[row,j] = k #1
            j = j + 1
        end
    end
end
for nc in 1:ncolsEnd
    for nr in 1:73
        if ΩendCC[nr,nc] == nothing
            ΩendCC[nr,nc] = -1
        end
    end
end

S = ones(Nb) # Base power
Y = vcat(1 ./ branchs[!,5])
Yc = [1/branchs[119,5], 1/branchs[120,5], 1/branchs[99,5], 1/branchs[109,5]] #1/X # Admittance of circuit c
M = repeat([1000], CC) # Big M

fmax = vcat(branchs[!,9]) # Maximum capacity of an existing circuit c
maximum(fmax)
fcmax = [branchs[119,9], branchs[120,9], branchs[99,9], branchs[109,9]] # Maximum capacity of a candidate circuit c
#av = ones(length(generators[!,5]), Nh) # [g,h] # Availability state of a generator g during hour h (1 if yes, 0 otherwise)

##################################### Defining the model #####################################

TEP = Model(optimizer_with_attributes(HiGHS.Optimizer)) # Gurobi

println("modelo")

##################################### Variables #####################################

@variable(TEP, pg[1:Ng,1:Nh] >= 0) # Generation at conventional generator g at hour h must be greater than zero
@variable(TEP, pn[1:Nn,1:Nh] >= 0) # Generation at renewable generator n at hour h must be greater than zero
@variable(TEP, pns[1:Nb,1:Nh] >= 0) # Power non served (pns) at node b at hour h must be greater than zero
@variable(TEP, pns2[1:Nb,1:Nh] >= 0) # Power spilled (pns2) at node b at hour h must be greater than zero
@variable(TEP, x[1:CC], Bin) # Indicates if it will be invested in the candidate circuit (1 if yes, 0 otherwise)
@variable(TEP, f[1:C,1:Nh] >= 0) # Power flow through existing circuit c during hour h
@variable(TEP, fc[1:CC,1:Nh] >= 0) # Power flow through candidate circuit c during hour h
@variable(TEP, θ[1:Nb,1:Nh]) # Voltage angle at node b during hour h

# Objective Function

@objective(TEP, Min, sum(ρ[h] * sum(pg[g,h] * CG[g] for g in 1:Ng) for h in 1:Nh) + 
sum(ρ[h] * sum((pns[b,h] + pns2[b,h]) * Cens for b in 1:Nb) for h in 1:Nh) + 
sum(x[c] * Cost[c] for c in 1:CC))

println("OV")

# Constraints

for h in 1:Nh
    @constraint(TEP, [b = 1:Nb], 
    sum(pg[ΩGconv[b,g],h] for g in 1:8 if ΩGconv[b,g].>-1) + 
    sum(pn[ΩGrenew[b,n],h] for n in 1:15 if ΩGrenew[b,n].>-1) - 
    d[b,h] + pns[b,h] - pns2[b,h] - 
    sum(f[ΩstartC[b,c],h] for c in 1:ncolsStart if ΩstartC[b,c].>-1) + 
    sum(f[ΩendC[b,c],h] for c in 1:ncolsEnd if ΩendC[b,c].>-1) - 
    sum(fc[ΩstartCC[b,c],h] for c in 1:ncolsStart if ΩstartCC[b,c].>-1) + 
    sum(fc[ΩendCC[b,c],h] for c in 1:ncolsEnd if ΩendCC[b,c].>-1) + 
    ΩsolEol[b,h] == 0
    )
end

println("r1")

for h in 1:Nh
    for (b, buss) in enumerate(Nb)
        for c in 1:length(Ωstart)
            if buss == Ωstart[c]
                busstart = Ωstart[c]
                busend = Ωend[c]
                i = findall(x->x==busstart, busList)
                j = findall(x->x==busend, busList)
                @constraint(TEP, f[c,h] == S[b] * Y[c] * (θ[i,h] - θ[j,h]))
            end
        end
    end
end

println("r2")

for h in 1:Nh
    for (b, buss) in enumerate(Nb)
        for c in 1:length(Ωstartcc)
            if buss == Ωstartcc[c]
                busstart = Ωstartcc[c]
                busend = Ωendcc[c]
                i = findall(x->x==busstart, busList)
                j = findall(x->x==busend, busList)
                @constraint(TEP, fc[c,h] - (S[b] * Y[c] * (θ[i,h] - θ[j,h])) <= M[c] * (1 - x[c]))
            end
        end
    end
end

println("r3")

for h in 1:Nh
    for (b, buss) in enumerate(Nb)
        for c in 1:length(Ωstartcc)#MUDEI
            if buss == Ωstartcc[c]
                busstart = Ωstartcc[c]
                busend = Ωendcc[c]
                i = findall(x->x==busstart, busList)
                j = findall(x->x==busend, busList)
                @constraint(TEP, - fc[c,h] + (S[b] * Y[c] * (θ[i,h] - θ[j,h])) <= M[c] * (1 - x[c]))
            end
        end
    end
end

println("r4")

for h in 1:Nh
    @constraint(TEP, [c = 1:C], f[c,h] <= fmax[c])
end

for h in 1:Nh
    @constraint(TEP, [c = 1:C], -f[c,h] <= fmax[c])
end

for h in 1:Nh
    @constraint(TEP, [c = 1:CC], fc[c,h] <= x[c] * fcmax[c])
end

for h in 1:Nh
    @constraint(TEP, [c = 1:CC], -fc[c,h] <= x[c] * fcmax[c])
end

for h in 2:Nh
    @constraint(TEP, [g = 1:Ng], pg[g,h] - pg[g,h-1] <= Rmaxconv[g])
end

for h in 2:Nh
    @constraint(TEP, [g = 1:Ng], pg[g,h] - pg[g,h-1] >= -Rmaxconv[g])
end

for h in 2:Nh
    @constraint(TEP, [n = 1:Nn], pn[n,h] - pn[n,h-1] <= Rmaxrenew[n])
end

for h in 2:Nh
    @constraint(TEP, [n = 1:Nn], pn[n,h] - pn[n,h-1] >= -Rmaxrenew[n])
end

for h in 1:Nh
    @constraint(TEP, [b = 1:Nb], pns[b,h] <= d[b,h])
end

for h in 1:Nh
    @constraint(TEP, [g = 1:Ng], pg[g,h] <= pgmax[g])
end

for h in 1:Nh
    @constraint(TEP, [n = 1:Nn], pn[n,h] <= pnmax[n])
end

println("vai rodar")

optimize!(TEP)

println("rodou")

# Results

status = termination_status(TEP)
println(status)

obj_value = objective_value(TEP)
println(obj_value)

Time = solve_time(TEP)
println(Time)

line_built = value.(x)
println(line_built)

conv_gen = value.(pg)
#conv_gen1 = DataFrames(conv_gen, :auto)
conv_gen1 = convert(DataFrame, conv_gen)
CSV.write("conv_gen1.csv", conv_gen1, append=true)

renew_gen = value.(pn)
renew_gen1 = convert(DataFrame, renew_gen)
#renew_gen1 = DataFrames(renew_gen, :auto)
CSV.write("renew_gen1.csv", renew_gen1, append=true)

not_served = value.(pns)
not_served1 = convert(DataFrame, not_served)
CSV.write("not_served1.csv", not_served1, append=true)

#line_built1 = convert(DataFrame, line_built)
#line_built1 = DataFrames(line_built, :auto)
#CSV.write("line_built.csv", line_built1, append=true)

flows = value.(f)
flows1 = convert(DataFrame, flows)
CSV.write("flows1.csv", flows1, append=true)

fclows = value.(fc)
fclows1 = convert(DataFrame, fclows)
CSV.write("flowsCand1.csv", fclows1, append=true)

thetas = value.(θ)
thetas1 = convert(DataFrame, thetas)
CSV.write("thetas1.csv", thetas1, append=true)

#obj_value1 = convert(DataFrame, obj_value)
#CSV.write("obj_value.csv", obj_value1, append=true)

#status1 = convert(DataFrame, status)
#CSV.write("status.csv", status1, append=true)

println(size(renew_gen1))
println(size(conv_gen1))
