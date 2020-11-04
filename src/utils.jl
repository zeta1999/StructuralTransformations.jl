function add_coeff!(dict, term, coeff)
    if coeff === nothing
        dict[term] = nothing
        return dict
    end
    iszero(coeff) && return dict
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
issym(x::Term) = operation(x) isa Sym
issym(::Number) = false

find_linear_coeffs(eq::Equation) = find_linear_coeffs(simplify(eq.rhs - eq.lhs; polynorm=true))
function find_linear_coeffs(t)
    dict = Dict()
    find_linear_coeffs!(dict, t)
    filter(x->x[2]!==nothing, dict)
end
function find_linear_coeffs!(dict, t)
    if issym(t)
        add_coeff!(dict, t, 1)
    elseif t isa Number
        add_coeff!(dict, 1, t)
    elseif t isa Term
        op = operation(t)
        args = arguments(t)
        if op === (*) && length(args) == 2
            if (issym(args[1]) ⊻ issym(args[2]))
                if !(args[1] isa Symbolic)
                    add_coeff!(dict, args[2], args[1])
                elseif !(args[2] isa Symbolic)
                    add_coeff!(dict, args[1], args[2])
                else
                    for x in vars(t)
                        add_coeff!(dict, x, nothing)
                    end
                end
            else
                for x in vars(t)
                    add_coeff!(dict, x, nothing)
                end
            end
        elseif op === (+)
            foreach(x->find_linear_coeffs!(dict, x), args)
        else
            for x in vars(t)
                add_coeff!(dict, x, nothing)
            end
        end
    end
    return dict
end
find_linear_coeffs!(dict, t::Num) = find_linear_coeffs!(dict, value(t))
