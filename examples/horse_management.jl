#  Copyright 2017, Eyob Zewdie, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=

    Inspired by

    R. McCardle,  Farm management optimization.,  Master’s thesis,  University
    of Louisville,  Louisville,  Kentucky, United States of America (2009)

    262.4213
=#

using SDDP, JuMP, Gurobi#Clp
using Base.Test

# cutting, stage
const S = [
    0 1 2;
    0 0 1;
    0 0 0
]

# days in period
const t = [60,60,245]
# demand
const D = [210,210,858]

# selling price per bale
const q = [
    [4.5 4.5 4.5; 4.5 4.5 4.5; 4.5 4.5 4.5],
    [5.5 5.5 5.5; 5.5 5.5 5.5; 5.5 5.5 5.5],
    [6.5 6.5 6.5; 6.5 6.5 6.5; 6.5 6.5 6.5]
]

# predicted yield (bales/acres) from cutting i in weather j
const b = [
    30 75 37.5;
    15 37.5 18.25;
    7.5 18.75 9.325
]
# max storage
const w = 3000
# cost to grow hay
const C = [50 50 50; 50 50 50; 50 50 50]
# Cost per bale of hay from cutting i during weather condition j;
const r = [
    [5 5 5; 5 5 5; 5 5 5],
    [6 6 6; 6 6 6; 6 6 6],
    [7 7 7; 7 7 7; 7 7 7]
]

# max acreage for planting
const M = 40.0#60.0
# initial inventory
const H = 0.0
# inventory cost
const V = [0.1, 0.1, 0.1]
# max demand for hay
const L = 3000.0

const transition = Array{Float64, 2}[
    [1.0]',
    [0.14 0.69 0.17],
    [1.0 1.0 1.0]',
    [0.14 0.69 0.17],
    [1.0 1.0 1.0]',
    [0.14 0.69 0.17]
]

m = SDDPModel(
             stages = 6,
    objective_bound = 0.0,
              sense = :Min,
  markov_transition = transition,
            #  solver = ClpSolver()
              solver = GurobiSolver(OutputFlag=0)
                        ) do sp, stage, weather
    @states(sp, begin
        # acres planted for each cutting
        0 <= acres <= M, acres0==M
        # bales from cutting i in storage
        bales[cutting=1:3] >= 0, bales0==[H,0,0][cutting]
        # because of the weird staging, we introduce additional state variables
        # quantity of bales to buy from cutting
        buy[cutting=1:3] >= 0, buy0==0
        # quantity of bales to sell from cutting
        sell[cutting=1:3] >= 0, sell0==0
        # quantity of bales to eat from cutting
        eat[cutting=1:3] >= 0, eat0==0
        # total quantity sold
        sold_bales >= 0, sold_bales0==0
    end)

    @variables(sp, begin
        pen_p[cutting=1:3] >= 0
        pen_n[cutting=1:3] >= 0
        pen_sell >= 0
    end)
    @expression(sp, total_penalties, sum(pen_p) + sum(pen_n) + pen_sell)

    if mod(stage, 2) == 1
        plan_cutting = round(Int, (stage+1)/2)
        # planning stage
        @constraints(sp, begin
            # plan planting for next stage
            acres <= acres0

            # plan to meet demand for next stage
            sum(eat) >= D[plan_cutting]
            # max sales over season
            sold_bales >= sold_bales0 + sum(sell)
            sold_bales <= L

            bales .== bales0

            # ensure bales fit in worst case (i.e. biggest yield)
            # sum(bales[c] + buy[c] - eat[c] - sell[c] for c in 1:3) + maximum(b[plan_cutting, :]) * acres <= w + max_bales_pen

            # ensure nonnegative bales in worst case (i.e. lowest yield)
            # bales[plan_cutting] + minimum(b[plan_cutting, :]) * acres + buy[plan_cutting] - eat[plan_cutting] - sell[plan_cutting] >= 0 - pen_n[plan_cutting]
            # [c=1:3;c!=plan_cutting], bales[c] + buy[c] - eat[c] - sell[c] >= 0 - pen_n[c]
        end)

        # planning objective
        stageobjective!(sp, 0.0)
    else
        op_cutting = round(Int, stage/2)
        @expression(sp, cut_ex[c=1:3], bales0[c] + buy0[c]  - eat0[c] - sell0[c] + pen_p[c] - pen_n[c])
        # action stage
        @constraints(sp, begin
            acres       == acres0
            sold_bales  == sold_bales0
            eat        .== eat0
            sell       .== sell0
            buy        .== buy0

            bales[op_cutting] == cut_ex[op_cutting] + acres0 * b[op_cutting, weather]
            # can buy and sell other cuttings
            [c=1:3;c!=op_cutting], bales[c] == cut_ex[c]
        end)

        stageobjective!(sp,
            1e2 * total_penalties +
            # cost of growing
            C[op_cutting, weather] * acres0 +
            sum(
                # inventory cost
                V[op_cutting] * bales0[cutting] * t[op_cutting] +
                # purchase cost
                r[cutting][op_cutting,weather] * buy0[cutting] +
                # feed cost
                S[cutting,op_cutting] * eat0[cutting] -
                # sell reward
                q[cutting][op_cutting,weather] * sell0[cutting]
            for cutting in 1:3)
        )
    end
end

srand(111)

solution = solve(m, max_iterations=50)

results = simulate(m,  # Simulate the policy
    100,               # number of monte carlo realisations
    [:acres,:acres0,:buy0,:sell0,:eat0,:pen_p,:pen_n,:sold_bales,:bales0, :pen_sell]       # variables to return
    )

@visualise(results, replication, stage, begin
	results[replication][:stageobjective][stage],              (title="Accumulated Profit", ylabel="Accumulated Profit (\$)", cumulative=true)
	results[replication][:stageobjective][stage],              (title="Weekly Income",      ylabel="Week Profit (\$)")
	results[replication][:acres][stage],    (title="acres")
    results[replication][:acres0][stage],    (title="acres0")
    sum(results[replication][:bales0][stage]),    (title="bales0")
    sum(results[replication][:buy0][stage]),    (title="buy0")
    sum(results[replication][:sell0][stage]),    (title="sell0")
    sum(results[replication][:eat0][stage]),    (title="eat0")
    sum(results[replication][:pen_p][stage]),    (title="pen_p")
    sum(results[replication][:pen_n][stage]),    (title="pen_n")
    results[replication][:markov][stage],    (title="weather")
    results[replication][:sold_bales][stage],    (title="sold bales")
    results[replication][:pen_sell][stage],    (title="pen_sell")
end)


# srand(111)
# @test isapprox(getbound(m), 262.4213, atol=1e-4)
