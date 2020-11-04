###
### System -> structural information
###

# V-nodes `[x_1, x_2, x_3, ..., dx_1, dx_2, ..., y_1, y_2, ...]` where `x`s are
# differential variables and `y`s are algebraic variables.
function get_vnodes(sys)
    dxvars = []
    eqs = equations(sys)
    edges = map(_->Int[], 1:length(eqs))
    for (i, eq) in enumerate(eqs)
        if eq.lhs isa Symbolic
            # Make sure that the LHS is a first order derivative of a var.
            @assert eq.lhs.op isa Differential "The equation $eq is not in the form of `D(...) ~ ...`"
            @assert !(eq.lhs.args[1] isa Differential) "The equation $eq is not first order"

            push!(dxvars, eq.lhs)
            # For efficiency we note down the diff edges here
            push!(edges[i], length(dxvars))
        end
    end

    xvars = (first ∘ var_from_nested_derivative).(dxvars)
    algvars  = setdiff(states(sys), xvars)
    return xvars, dxvars, edges, algvars
end

function sys2bigraph(sys; find_solvables = false)
    xvars, dxvars, edges, algvars = get_vnodes(sys)
    xvar_offset = length(xvars)
    algvar_offset = 2xvar_offset
    for edge in edges
        isempty(edge) || (edge .+= xvar_offset)
    end

    eqs = equations(sys)
    for (i, eq) in enumerate(eqs)
        # T or D(x):
        # We assume no derivatives appear on the RHS at this point
        vs = vars(eq.rhs)
        for v in vs
            for (j, target_v) in enumerate(xvars)
                if isequal(v, target_v)
                    push!(edges[i], j)
                end
            end
            for (j, target_v) in enumerate(algvars)
                if isequal(v, target_v)
                    push!(edges[i], j+algvar_offset)
                end
            end
        end
    end

    fullvars = [xvars; dxvars; algvars] # full list of variables
    vars_asso = Int[(1:xvar_offset) .+ xvar_offset; zeros(Int, length(fullvars) - xvar_offset)] # variable association list

    if find_solvables
        solvable_edges = map(_->Int[], 1:length(eqs))
        J = ModelingToolkit.jacobian(map(eq->eq.rhs-eq.lhs, eqs), fullvars)
        for eq in 1:length(eqs), var in 1:length(fullvars)
            v = value(J[eq, var])
            if !(v isa Symbolic) && v isa Integer && v != 0
                push!(solvable_edges[eq], var)
            end
        end
        return edges, solvable_edges, fullvars, vars_asso
    else
        return edges, fullvars, vars_asso
    end
end

#=
TODO: optimized solvable variable detection
function add_coeff!(dict, term, coeff, sign)
    @show term, coeff
    if coeff === nothing
        dict[term] = nothing
        return dict
    end
    iszero(coeff) && return dict
    coeff = sign ? coeff : -coeff
    if haskey(dict, term)
        if dict[term] !== nothing
            dict[term] += coeff
        end
    else
        dict[term] = coeff
    end
    return dict
end

issym(x::Sym) = true
function issym(x::Term)
    op = operation(x)
    op isa Sym || op isa Differential
end
issym(::Number) = false

find_linear_coeffs(eq::Equation) = find_linear_coeffs(eq.rhs - eq.lhs)
function find_linear_coeffs(t)
    dict = Dict()
    find_linear_coeffs!(dict, t)
    filter(x->x[2]!==nothing, dict)
end
function find_linear_coeffs!(dict, t, sign=true)
    if issym(t)
        add_coeff!(dict, t, 1, sign)
    elseif t isa Number
        add_coeff!(dict, 1, t, sign)
    elseif t isa Term
        op = operation(t)
        args = arguments(t)
        if op === (*) && length(args) == 2
            if (issym(args[1]) ⊻ issym(args[2]))
                if !(args[1] isa Symbolic)
                    add_coeff!(dict, args[2], args[1], sign)
                elseif !(args[2] isa Symbolic)
                    add_coeff!(dict, args[1], args[2], sign)
                else
                    for x in vars(t)
                        add_coeff!(dict, x, nothing, sign)
                    end
                end
            else
                for x in vars(t)
                    add_coeff!(dict, x, nothing, sign)
                end
            end
        elseif op === (+)
            foreach(x->find_linear_coeffs!(dict, x, sign), args)
        elseif op === (-)
            foreach(x->find_linear_coeffs!(dict, x, !sign), args)
        else
            for x in vars(t)
                add_coeff!(dict, x, nothing, sign)
            end
        end
    end
    return dict
end
find_linear_coeffs!(dict, t::Num) = find_linear_coeffs!(dict, value(t))
=#
