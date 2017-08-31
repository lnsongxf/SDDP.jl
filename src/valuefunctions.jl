#  Copyright 2017, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

stageobjective!(sp::JuMP.Model, obj...) = stageobjective!(valueoracle(sp), sp, obj...)
getstageobjective(sp::JuMP.Model) = getstageobjective(valueoracle(sp), sp)
modifyvaluefunction!(m::SDDPModel, settings::Settings, sp::JuMP.Model) = modifyvaluefunction!(valueoracle(sp), m, settings, sp)

type DefaultValueFunction{C<:AbstractCutOracle} <: AbstractValueFunction
    cutmanager::C
    stageobjective::QuadExpr
    theta::JuMP.Variable
end



function solvesubproblem!(direction, valuefunction, m::SDDPModel, sp::JuMP.Model)
    JuMPsolve(direction, m, sp)
end
solvesubproblem!(direction, m::SDDPModel, sp::JuMP.Model) = solvesubproblem!(direction, valueoracle(sp), m, sp)


# ==============================================================================
#   3. solvesubproblem!(ForwardPass, ...)

function writecut!(filename::String, cut::Cut, args...)
    open(filename, "a") do file
        for arg in args
            write("$(arg), ")
        end
        write(file, "$(cut.intercept)")
        for pi in cut.coefficients
            write(file, ",$(pi)")
        end
        write(file, "\n")
    end
end

function loadcuts!{C}(m::SDDPModel{DefaultValueFunction{C}}, filename::String)
    open(filename, "r") do file
        while true
            line      = readline(file)
            line == nothing || line == "" && break
            items     = split(line, ",")
            stage     = parse(Int, items[1])
            ms        = parse(Int, items[2])
            intercept = parse(Float64, items[3])
            coefficients = [parse(Float64, i) for i in items[4:end]]
            cut = Cut(intercept, coefficients)
            sp = getsubproblem(m, stage, ms)
            vf = valueoracle(sp)
            # Add cut to the cut manager
            storecut!(vf.cutmanager, m, sp, cut)
            # Add cut as a constraint to the subproblem
            addcut!(vf, sp, cut)
        end
    end
end

function modifyvaluefunction!(vf::DefaultValueFunction, m::SDDPModel, settings::Settings, sp::JuMP.Model)
    ex = ext(sp)
    I = 1:length(m.storage.objective)
    for i in I
        m.storage.probability[i] *= getstage(m, ex.stage+1).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
    modifyprobability!(ex.riskmeasure,
        view(m.storage.modifiedprobability.data, I),
        m.storage.probability.data[I],
        sp,
        getstage(m, ex.stage).state,
        m.storage.duals.data[I],
        m.storage.objective.data[I]
    )
    cut = constructcut(m, sp)

    if settings.cut_output_file != ""
        writecut!(settings.cut_output_file, cut, ex.stage, ex.markovstate)
    end

    storecut!(vf.cutmanager, m, sp, cut)
    addcut!(vf, sp, cut)
    storecut!(m, sp, cut)
    for i in I
        m.storage.probability[i] /= getstage(m, ex.stage+1).transitionprobabilities[ex.markovstate, m.storage.markov[i]]
    end
end
storecut!(m::SDDPModel, sp::JuMP.Model, cut::Cut, args...) = nothing

function solvesubproblem!(::Type{BackwardPass}, m::SDDPModel, sp::JuMP.Model, probability::Float64=1.0)
    ex = ext(sp)
    if hasnoises(sp)
        for (i, (rhs_noise, rhs_probability)) in enumerate(zip(ex.noises, ex.noiseprobability))
            setnoise!(sp, rhs_noise)
            JuMPsolve(BackwardPass, m, sp)
            push!(m.storage.objective, getobjectivevalue(sp))
            push!(m.storage.noise, i)
            push!(m.storage.probability, probability*rhs_probability)
            push!(m.storage.modifiedprobability, probability*rhs_probability)
            push!(m.storage.markov, ex.markovstate)
            push!(m.storage.duals, zeros(nstates(sp)))
            saveduals!(m.storage.duals[end], sp)
        end
    else
        JuMPsolve(BackwardPass, m, sp)
        push!(m.storage.objective, getobjectivevalue(sp))
        push!(m.storage.noise, 0)
        push!(m.storage.probability, probability)
        push!(m.storage.modifiedprobability, probability)
        push!(m.storage.markov, ex.markovstate)
        push!(m.storage.duals, zeros(nstates(sp)))
        saveduals!(m.storage.duals[end], sp)
    end
end
