module StructuralTransformations

const UNVISITED = 0
const UNASSIGNED = 0

using ModelingToolkit
using ModelingToolkit: ODESystem, var_from_nested_derivative, Differential,
                       states, equations, vars, Symbolic, diff2term, value,
                       operation, arguments, Sym, Term, simplify

include("utils.jl")

# Equation-variable bipartite matching
"""
    find_augmenting_path(edges, eq, assign, active, vcolor=falses(length(active)), ecolor=falses(length(edges))) -> path_found::Bool

Try to find augmenting paths.
"""
function find_augmenting_path(edges, eq, assign, active, vcolor=falses(length(active)), ecolor=falses(length(edges)))
    # Note: `edge[active]` are active edges
    ecolor[eq] = true

    # if a `var` is unassigned and the edge `eq <=> var` exists
    for var in edges[eq]
        if active[var] && assign[var] == UNASSIGNED
            assign[var] = eq
            return true
        end
    end

    # for every `var` such that edge `eq <=> var` exists and `var` is uncolored
    for var in edges[eq]
        (active[var] && !vcolor[var]) || continue
        vcolor[var] = true
        if find_augmenting_path(edges, assign[var], assign, active, vcolor, ecolor)
            assign[var] = eq
            return true
        end
    end
    return false
end

"""
    find_augmenting_path(edges, nvars, active=trues(nvars)) -> assign

Find equation-variable bipartite matching. `edges` is a bipartite graph. `nvars`
is the number of variables. Matched variables and equations are stored in like
`assign[var] => eq`.
"""
function matching(edges, nvars, active=trues(nvars))
    assign = fill(UNASSIGNED, nvars)
    for eq in 1:length(edges)
        find_augmenting_path(edges, eq, assign, active)
    end
    return assign
end

###
### Reassemble: structural information -> system
###

# Naive subtree matching, we can make the dict have levels
# Going forward we should look into storing the depth in Terms
function walk_and_substitute(expr, substitution_dict)
    if haskey(substitution_dict, expr)
        return substitution_dict[expr]
    elseif expr isa Term
        return Term(walk_and_substitute(expr.op, substitution_dict),
                         map(x->walk_and_substitute(x, substitution_dict), expr.args))
    else
        @assert !(expr isa Equation) "RHS cannot contain equations"
        return expr
    end
end

function pantelides_reassemble(sys, vars, vars_asso, eqs_asso, assign)
    # Step 1: write derivative equations
    in_eqs = equations(sys)
    out_eqs = Vector{Any}(undef, length(eqs_asso))
    fill!(out_eqs, nothing)
    out_eqs[1:length(in_eqs)] .= in_eqs

    out_vars = Vector{Any}(undef, length(vars_asso))
    fill!(out_vars, nothing)
    out_vars[1:length(vars)] .= vars

    D = ModelingToolkit.Differential(sys.iv)

    for (i, v) in enumerate(vars_asso)
        # vars[v] = D(vars[i])
        v == 0 && continue
        vi = out_vars[i]
        @assert vi !== nothing "Something went wrong on reconstructing states from variable association list"
        # `vars[i]` needs to be not a `D(...)`, because we want the DAE to be
        # first-order.
        if vi isa Term && operation(vi) isa Differential
            vi = out_vars[i] = diff2term(vi)
        end
        out_vars[v] = D(vi)
    end

    d_dict = Dict(zip(vars, 1:length(vars)))
    lhss = Set{Any}([x.lhs for x in in_eqs if x.lhs isa Term && x.lhs.op isa Differential])
    for (i, e) in enumerate(eqs_asso)
        if e === 0
            continue
        end
        # LHS variable is looked up from vars_asso
        # the vars_asso[i]-th variable is the differentiated version of var at i
        eq = out_eqs[i]
        lhs = if !(eq.lhs isa Symbolic)
            0
        elseif eq.lhs isa Term && eq.lhs.op isa Differential
            # look up the variable that represents D(lhs)
            @assert !(eq.lhs.args[1] isa Differential) "The equation $eq is not first order"
            i = get(d_dict, eq.lhs.args[1], nothing)
            if i !== nothing
                lhs = D(out_vars[vars_asso[i]])
                if lhs in lhss
                    # check only trivial equations are removed
                    @assert isequal(diff2term(D(eq.rhs)), diff2term(lhs)) "The duplicate equation is not trivial: $eq"
                    lhs = Num(nothing)
                end
                lhs
            else
                D(eq.lhs)
            end
        else
            D(eq.lhs)
        end
        rhs = ModelingToolkit.expand_derivatives(D(eq.rhs))
        substitution_dict = Dict(x.lhs => x.rhs for x in out_eqs if x !== nothing && x.lhs isa Symbolic)
        sub_rhs = walk_and_substitute(rhs, substitution_dict)
        out_eqs[e] = lhs ~ sub_rhs
    end

    final_vars = unique(filter(x->!(x.op isa Differential), vars))
    final_eqs = map(identity, filter(x->value(x.lhs) !== nothing, out_eqs[sort(filter(!iszero, assign))]))

    # remove clashing equations (from order lowering vs index reduction)
    return ODESystem(final_eqs, sys.iv, final_vars, sys.ps)
end

"""
    pantelides(sys::ODESystem; kwargs...)

Perform graph based Pantelides algorithm.
"""
function pantelides(sys::ODESystem; kwargs...)
    edges, fullvars, vars_asso = sys2bigraph(sys)
    return pantelides!(edges, vars_asso, sys.iv; kwargs...)
end

"""
    pantelides!(edges, vars, vars_asso, iv; maxiters = 8000)

Perform Pantelides algorithm.
"""
function pantelides!(edges, vars_asso, iv; maxiters = 8000)
    neqs = length(edges)
    nvars = length(vars_asso)
    vcolor = falses(nvars)
    ecolor = falses(neqs)
    assign = fill(UNASSIGNED, nvars)
    eqs_asso = fill(0, neqs)
    neqs′ = neqs
    D = Differential(iv)
    for k in 1:neqs′
        eq′ = k
        pathfound = false
        # In practice, `maxiters=8000` should never be reached, otherwise, the
        # index would be on the order of thousands.
        for iii in 1:maxiters
            # run matching on (dx, y) variables
            #
            # the derivatives and algebraic variables are zeros in the variable
            # association list
            active = vars_asso .== 0
            resize!(vcolor, nvars)
            fill!(vcolor, false)
            resize!(ecolor, neqs)
            fill!(ecolor, false)
            pathfound = find_augmenting_path(edges, eq′, assign, active, vcolor, ecolor)
            pathfound && break # terminating condition
            for var in eachindex(vcolor); vcolor[var] || continue
                # introduce a new variable
                nvars += 1
                # the new variable is the derivative of `var`
                vars_asso[var] = nvars
                push!(vars_asso, 0)
                push!(assign, UNASSIGNED)
            end

            for eq in eachindex(ecolor); ecolor[eq] || continue
                # introduce a new equation
                neqs += 1
                push!(edges, copy(edges[eq]))
                # the new equation is created by differentiating `eq`
                eqs_asso[eq] = neqs
                for var in edges[eq]
                    if !(vars_asso[var] in edges[neqs])
                        push!(edges[neqs], vars_asso[var])
                    end
                end
                push!(eqs_asso, 0)
            end

            for var in eachindex(vcolor); vcolor[var] || continue
                # the newly introduced `var`s and `eq`s have the inherits
                # assignment
                assign[vars_asso[var]] = eqs_asso[assign[var]]
            end
            eq′ = eqs_asso[eq′]
        end # for _ in 1:maxiters
        pathfound || error("maxiters=$maxiters reached! File a bug report if your system has a reasonable index (<100), and you are using the default `maxiters`. Try to increase the maxiters by `pantelides(sys::ODESystem; maxiters=1_000_000)` if your system has an incredibly high index and it is truly extremely large.")
    end # for k in 1:neqs′
    return edges, assign, vars_asso, eqs_asso
end

"""
    dae_index_lowering(sys::ODESystem) -> ODESystem

Perform the Pantelides algorithm to transform a higher index DAE to an index 1
DAE.
"""
function dae_index_lowering(sys::ODESystem; kwargs...)
    edges, fullvars, vars_asso = sys2bigraph(sys)
    edges, assign, vars_asso, eqs_asso = pantelides!(edges, vars_asso, sys.iv; kwargs...)
    return pantelides_reassemble(sys, fullvars, vars_asso, eqs_asso, assign)
end

###
### BLT ordering
###

"""
    find_scc(edges, assign=nothing)

Find strongly connected components of the graph defined by `edges`. When `assign
=== nothing`, we assume that the ``i``-th variable is assigned to the ``i``-th
equation.
"""
function find_scc(edges, assign=nothing)
    id = 0
    stack = Int[]
    components = Vector{Int}[]
    n = length(edges)
    onstack = falses(n)
    lowlink = zeros(Int, n)
    ids = fill(UNVISITED, n)

    for eq in 1:length(edges)
        if ids[eq] == UNVISITED
            id = strongly_connected!(stack, onstack, components, lowlink, ids, edges, assign, eq, id)
        end
    end
    return components
end

"""
    strongly_connected!(stack, onstack, components, lowlink, ids, edges, assign, eq, id)

Use Tarjan's algorithm to find strongly connected components.
"""
function strongly_connected!(stack, onstack, components, lowlink, ids, edges, assign, eq, id)
    id += 1
    lowlink[eq] = ids[eq] = id

    # add `eq` to the stack
    push!(stack, eq)
    onstack[eq] = true

    # for `adjeq` in the adjacency list of `eq`
    for var in edges[eq]
        if assign === nothing
            adjeq = var
        else
            # assign[var] => the equation that's assigned to var
            adjeq = assign[var]
            # skip equations that are not assigned
            adjeq == UNASSIGNED && continue
        end

        # if `adjeq` is not yet idsed
        if ids[adjeq] == UNVISITED # visit unvisited nodes
            id = strongly_connected!(stack, onstack, components, lowlink, ids, edges, assign, adjeq, id)
        end
        # at the callback of the DFS
        if onstack[adjeq]
            lowlink[eq] = min(lowlink[eq], lowlink[adjeq])
        end
    end

    # if we are at a start of a strongly connected component
    if lowlink[eq] == ids[eq]
        component = Int[]
        repeat = true
        # pop until we are at the start of the strongly connected component
        while repeat
            w = pop!(stack)
            onstack[w] = false
            lowlink[w] = ids[eq]
            # put `w` in current component
            push!(component, w)
            repeat = w != eq
        end
        push!(components, component)
    end
    return id
end

end # module
