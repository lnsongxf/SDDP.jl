#=
  Copyright 2017, Oscar Dowson

  This file contins methods specific to the default value function
=#

function DefaultValueFunction(cutoracle=DefaultCutOracle())
    DefaultValueFunction(
        cutoracle,
        QuadExpr(0.0),
        JuMP.Variable(JuMP.Model(), 0)
    )
end

summarise{C}(::Type{DefaultValueFunction{C}}) = "Default Value Function"

## ==============================================================================
#   1. stageobjective!

function stageobjective!(vf::DefaultValueFunction, sp::JuMP.Model, obj)
    append!(vf.stageobjective, QuadExpr(obj))
    if ext(sp).finalstage
        JuMP.setobjective(sp, getsense(sp), obj)
    else
        JuMP.setobjective(sp, getsense(sp), obj + vf.theta)
    end
    JuMP.setobjectivesense(sp, getsense(sp))
end

# ==============================================================================
#   2. getstageobjective!

getstageobjective(vf::DefaultValueFunction, sp::JuMP.Model) = getvalue(vf.stageobjective)

# ==============================================================================
#   3. getstageobjective!

function init!{C}(vf::DefaultValueFunction{C}, m::JuMP.Model, sense, bound)
    vf.theta = futureobjective!(sense, m, bound)
    vf
end

# ==============================================================================
#   4. addcut!

function addcut!(vf::DefaultValueFunction, sp::JuMP.Model, cut::Cut)
    affexpr = cuttoaffexpr(sp, cut)
    _addcut!(ex.sense, sp, vf.theta, affexpr)
end
_addcut!(::Type{Min}, sp, theta, affexpr) = @constraint(sp, theta >= affexpr)
_addcut!(::Type{Max}, sp, theta, affexpr) = @constraint(sp, theta <= affexpr)

# ==============================================================================
#   5. rebuildsubproblem!

function rebuildsubproblem!{C<:AbstractCutOracle}(m::SDDPModel{DefaultValueFunction{C}}, sp::JuMP.Model)
    n = n_args(m.build!)
    ex = ext(sp)
    for i in 1:nstates(sp)
        pop!(ex.states)
    end
    for i in 1:length(ex.noises)
        pop!(ex.noises)
    end
    sp = Model(solver = m.lpsolver)

    vf.stageobjective = QuadExpr(0.0)
    vf.theta = futureobjective!(ex.sense, sp, ex.problembound)

    sp.ext[:SDDP] = ex
    if n == 2
        m.build!(sp, ex.stage)
    elseif n == 3
        m.build!(sp, ex.stage, ex.markovstate)
    end
    for cut in validcuts(vf.cutmanager)
        addcut!(vf, sp, cut)
    end
    m.stages[ex.stage].subproblems[ex.markovstate] = sp
end

rebuildsubproblem!(m::SDDPModel{DefaultValueFunction{DefaultCutOracle}}, sp::JuMP.Model) = nothing
