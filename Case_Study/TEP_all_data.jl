
##################################### Importing the packages #####################################

using JuMP
using Gurobi
using CSV
using DataFrames
using StatsBase

##################################### Opening the data #####################################

path = pwd()

buses = CSV.read("./Case_Study/data/bus.csv", DataFrame)
generators = CSV.read("./Case_Study/data/gen.csv", DataFrame)
branchs = CSV.read("./Case_Study/data/branch.csv", DataFrame)
ts = CSV.read("./Case_Study/data/DAY_AHEAD_regional_Load.csv", DataFrame) # 365 days -> no 29/02/2020
janauba = CSV.read("./Case_Study/data/JANAUBA_MG_UCT.csv", DataFrame)
z3 = CSV.read("./Case_Study/data/genZona3.csv", DataFrame)

##################################### Creating the zones #####################################

# How many % of MW each bus represent each zone total MW load
zone1 = buses[1:24,16]
zone2 = buses[25:48,16]
zone3 = buses[49:73,16]
sum(zone3)
busList = collect(buses[!,1])
begin2020 = 14601
end2020 = 14965
begin2020h = 350401
end2020h = 359160

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
Nh = length(ts[!,1]) # Number of hours of a day (24 hours)
Nb = length(buses[!,1]) # Number of nodes
CC = 4 # Set of candidate circuits
C = length(branchs[!,3]) # Set of existing circuits

##################################### Creating the centers #####################################

maxJanauba = maximum(janauba[!,4])
janauba = janauba[begin2020h:end2020h,4] # Zone 3
janauba = janauba ./ maxJanauba # Zone 3

println("zonas")

# Zone 3
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
solarZone3 = DataFrame()
for i in 1:length(solarZ3bus)
    global solarZone3 = hcat(solarZone3, janauba.*solarZ3bus[i], makeunique = true)
end
rename!(solarZone3, string.(UniqZ3))

Z3 = solarZone3

println("cluster_zonas")

busZ3Uniq = unique(z3[!,2])
ΩZ3 = Array{Union{Nothing,Float64}}(nothing, 73, Nh)
for (i, index) in enumerate(busZ3Uniq)
    j = 1
    for k in 1:length(names(Z3))
        nome = names(Z3)[k][1:3]
        if parse(Int64, nome) == index
            row = findall(x->x==index, busList)[1]
            global ΩZ3[row,:] = Z3[!,k]
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

Cens = 70 # The cost of non-served energy ($/MW)
ρ .* Cens
α = 1
Cost = [3000 3000 1500 1500] * α # The anualized cost of investing on a circuit c (cost between Zones is 50% greater than cost only in Zone 3)
#Cost = [0 0 0 0]

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
sumd = sum(d, dims=1)
col_max_d = findmax(sumd)[2][2]
demanda_dia = zeros(365)
for i in 1:365
    demanda_dia[i] = sum(sumd[(24*(i-1)+1):(24*i)])
end
col_max_dia = findmax(demanda_dia)[2]
demanda_dia[col_max_dia]
d2 = repeat(d[:,(24*(col_max_dia-1)+1):(24*col_max_dia)], 1, 365)
d = d2

diff = d - ΩsolEol
minimum(sum(diff, dims = 1))

println("demanda")

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
Cfix = CGconv.*10
ρ.*minimum(CGconv)
ρ.*maximum(CGconv)

busIDconvUniq = unique(busIDconv)
busIDrenewUniq = unique(busIDrenew)
ncolsrenew = count(x->x==StatsBase.mode(busIDrenew),busIDrenew)
ΩGrenew = Array{Union{Nothing,Int64}}(nothing, 73, ncolsrenew)
for (i, index) in enumerate(busIDrenewUniq)
    j = 1
    for k in 1:length(busIDrenew)
        if busIDrenew[k] == index
            row = findall(x->x==index, busList)[1]
            global ΩGrenew[row,j] = k # genNameRenew[k]
            j = j + 1
        end
    end
end
for nc in 1:ncolsrenew
    for nr in 1:73
        if ΩGrenew[nr,nc] == nothing
            ΩGrenew[nr,nc] = -1
        end
    end
end

ncolsconv = count(x->x==StatsBase.mode(busIDconv),busIDconv)
ΩGconv = Array{Union{Nothing,Int64}}(nothing, 73, ncolsconv)
for (i, index) in enumerate(busIDconvUniq)
    j = 1
    for k in 1:length(busIDconv)
        if busIDconv[k] == index
            row = findall(x->x==index, busList)[1]
            global ΩGconv[row,j] = k #genNameConv[k]
            j = j + 1
        end
    end
end
for nc in 1:ncolsconv
    for nr in 1:73
        if ΩGconv[nr,nc] == nothing
            ΩGconv[nr,nc] = -1
        end
    end
end

Ωend
Ωstart
ncolsStart = count(x->x==StatsBase.mode(Ωstart),Ωstart)
ncolsEnd = count(x->x==StatsBase.mode(Ωend),Ωend)

ΩendUniq = unique(Ωend)
ΩstartUniq = unique(Ωstart)
ΩendUniqcc = unique(Ωendcc)
ΩstartUniqcc = unique(Ωstartcc)

ΩstartC = Array{Union{Nothing,Int64}}(nothing, 73, ncolsStart)
ΩstartCC = Array{Union{Nothing,Int64}}(nothing, 73, ncolsStart)
for (i, index) in enumerate(ΩstartUniq)
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
for (i, index) in enumerate(ΩstartUniqcc)
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
for (i, index) in enumerate(ΩendUniqcc)
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
Y = 1 .+ 0 .*vcat(1 ./ branchs[!,5])
#Y = vcat(1 ./ branchs[!,5])
Yc = [1/branchs[119,5], 1/branchs[120,5], 1/branchs[99,5], 1/branchs[109,5]] #1/X # Admittance of circuit c
Yc = 1 .+ 0 .* Yc
M = repeat([1000], CC) # Big M

fmax = vcat(branchs[!,9]) # Maximum capacity of an existing circuit c
maximum(fmax)
fcmax = [branchs[119,9], branchs[120,9], branchs[99,9], branchs[109,9]] # Maximum capacity of a candidate circuit c
#av = ones(length(generators[!,5]), Nh) # [g,h] # Availability state of a generator g during hour h (1 if yes, 0 otherwise)

##################################### Defining the model #####################################

TEP = Model(optimizer_with_attributes(Gurobi.Optimizer))

time_lim_sec = 60

println("modelo")

##################################### Variables #####################################

@variable(TEP, pg[1:Ng,1:Nh] >= 0) # Generation at conventional generator g at hour h must be greater than zero
@variable(TEP, pn[1:Nn,1:Nh] >= 0) # Generation at renewable generator n at hour h must be greater than zero
@variable(TEP, pns[1:Nb,1:Nh] >= 0) # Power non served (pns) at node b at hour h must be greater than zero
@variable(TEP, pns2[1:Nb,1:Nh] >= 0) # Power spilled (pns2) at node b at hour h must be greater than zero
@variable(TEP, x[1:CC], Bin) # Indicates if it will be invested in the candidate circuit (1 if yes, 0 otherwise)
@variable(TEP, f[1:C,1:Nh]) # Power flow through existing circuit c during hour h
@variable(TEP, fc[1:CC,1:Nh]) # Power flow through candidate circuit c during hour h
@variable(TEP, θ[1:Nb,1:Nh]) # Voltage angle at node b during hour h
@variable(TEP, u[1:Ng,1:365], Bin) # Binary with number of days
@variable(TEP, psoleol[1:Nb,1:Nh] >= 0)

# Objective Function

@objective(TEP, Min, sum(ρ[h] * sum(pg[g,h] * CGconv[g] for g in 1:Ng) for h in 1:Nh) + 
sum(ρ[h] * sum((pns[b,h] + pns2[b,h]) * Cens for b in 1:Nb) for h in 1:Nh) + 
sum(x[c] * Cost[c] for c in 1:CC) +
sum(sum(Cfix[g] * u[g,k] for g in 1:Ng) for k in 1:365));

println("OV")

# Constraints

for h in 1:Nh
    @constraint(TEP, [b = 1:Nb], 
    sum(pg[ΩGconv[b,g],h] for g in 1:ncolsconv if ΩGconv[b,g].>-1) + 
    sum(pn[ΩGrenew[b,n],h] for n in 1:ncolsrenew if ΩGrenew[b,n].>-1) - 
    d[b,h] + pns[b,h] - pns2[b,h] - 
    sum(f[ΩstartC[b,c],h] for c in 1:ncolsStart if ΩstartC[b,c].>-1) + 
    sum(f[ΩendC[b,c],h] for c in 1:ncolsEnd if ΩendC[b,c].>-1) - 
    sum(fc[ΩstartCC[b,c],h] for c in 1:ncolsStart if ΩstartCC[b,c].>-1) + 
    sum(fc[ΩendCC[b,c],h] for c in 1:ncolsEnd if ΩendCC[b,c].>-1) + 
    psoleol[b,h] == 0
    )
end

for h in 1:Nh
    @constraint(TEP, [b = 1:Nb], psoleol[b,h] <= ΩsolEol[b,h])
end

println("r1")

for h in 1:Nh
    for (b, buss) in enumerate(buses[:,1])
        for c in 1:length(Ωstart)
            if buss == Ωstart[c]
                busstart = Ωstart[c]
                busend = Ωend[c]
                i = findall(x->x==busstart, busList)[1]
                j = findall(x->x==busend, busList)[1]
                @constraint(TEP, f[c,h] == S[b] * Y[c] * (θ[i,h] - θ[j,h]))
            end
        end
    end
end

println("r2")

for h in 1:Nh
    for (b, buss) in enumerate(buses[:,1])
        for c in 1:length(Ωstartcc)
            if buss == Ωstartcc[c]
                busstart = Ωstartcc[c]
                busend = Ωendcc[c]
                i = findall(x->x==busstart, busList)[1]
                j = findall(x->x==busend, busList)[1]
                @constraint(TEP, fc[c,h] - (S[b] * Yc[c] * (θ[i,h] - θ[j,h])) <= M[c] * (1 - x[c]))
            end
        end
    end
end

println("r3")

for h in 1:Nh
    for (b, buss) in enumerate(buses[:,1])
        for c in 1:length(Ωstartcc)
            if buss == Ωstartcc[c]
                busstart = Ωstartcc[c]
                busend = Ωendcc[c]
                i = findall(x->x==busstart, busList)[1]
                j = findall(x->x==busend, busList)[1]
                @constraint(TEP, - fc[c,h] + (S[b] * Yc[c] * (θ[i,h] - θ[j,h])) <= M[c] * (1 - x[c]))
            end
        end
    end
end

println("r4")

for h in 1:Nh
    @constraint(TEP, sum(f[c,h] for c in C) + sum(fc[c,h] for c in CC) == 0)
end

for h in 1:Nh
    @constraint(TEP, [c = 1:C], f[c,h] <= fmax[c])
end

for h in 1:Nh
    @constraint(TEP, [c = 1:C], f[c,h] >= -fmax[c])
end

for h in 1:Nh
    @constraint(TEP, [c = 1:CC], fc[c,h] <= x[c] * fcmax[c])
end

for h in 1:Nh
    @constraint(TEP, [c = 1:CC], fc[c,h] >= -x[c] * fcmax[c])
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

for k in 1:365
    @constraint(TEP, [g = 1:Ng, h = (24*(k-1)+1):(24*k)], pg[g,h] <= pgmax[g] * u[g,k])
end

println("vai rodar")

set_time_limit_sec(TEP, time_lim_sec)

optimize!(TEP)

println("rodou")

# Results

status = termination_status(TEP)
println("Status = $status")

obj_value = objective_value(TEP)
println("Objective Value = $obj_value")

cost_gen = sum(ρ[h] * sum(value(pg[g,h]) * CGconv[g] for g in 1:Ng) for h in 1:Nh)
println("Cost gen = $cost_gen")
geracao =  sum(sum(value(pg[g,h]) for g in 1:Ng) for h in 1:Nh)
println("Geração = $geracao")
cost_pns = sum(ρ[h] * sum((value(pns[b,h]) + value(pns2[b,h])) * Cens for b in 1:Nb) for h in 1:Nh)
total_imbalance_pos = sum(sum((value(pns[b,h])) for b in 1:Nb) for h in 1:Nh)
total_imbalance_neg = sum(sum((value(pns2[b,h])) for b in 1:Nb) for h in 1:Nh)
sum(sum(value(pg[g,h]) for g in 1:Ng) for h in 1:Nh)
println("Cost pns = $cost_pns")
println("Total Imbalance Positive = $total_imbalance_pos")
println("Total Imbalance Negative = $total_imbalance_neg")
cost_lines = sum(value(x[c]) * Cost[c] for c in 1:CC)
println("Cost lines = $cost_lines")
cost_days = sum(sum(Cfix[g] * value(u[g,k]) for g in 1:Ng) for k in 1:365)
println("Cost days = $cost_days")

Time = solve_time(TEP)
println("Time = $Time")

line_built = value.(x)
println("Lines = $line_built")

conv_gen = value.(pg)
conv_gen1 = DataFrame(conv_gen, :auto)
CSV.write("conv_gen_all_year.csv", conv_gen1, append=true)

renew_gen = value.(pn)
renew_gen1 = DataFrame(renew_gen, :auto)
CSV.write("renew_gen_all_year.csv", renew_gen1, append=true)

not_served = value.(pns)
not_served1 = DataFrame(not_served, :auto)
CSV.write("not_served_all_year.csv", not_served1, append=true)

spilled = value.(pns2)
spilled1 = DataFrame(spilled, :auto)
CSV.write("spilled_all_year.csv", spilled1, append=true)

flows = value.(f)
flows1 = DataFrame(flows, :auto)
CSV.write("flows_all_year.csv", flows1, append=true)

fclows = value.(fc)
fclows1 = DataFrame(fclows, :auto)
CSV.write("flowsCand_all_year.csv", fclows1, append=true)

thetas = value.(θ)
thetas1 = DataFrame(thetas, :auto)
CSV.write("thetas_all_year.csv", thetas1, append=true)

#obj_value1 = convert(DataFrame, obj_value)
#CSV.write("obj_value.csv", obj_value1, append=true)

#status1 = convert(DataFrame, status)
#CSV.write("status.csv", status1, append=true)

#line_built1 = convert(DataFrame, line_built)
#line_built1 = DataFrames(line_built, :auto)
#CSV.write("line_built.csv", line_built1, append=true)
