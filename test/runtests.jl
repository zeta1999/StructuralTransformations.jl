using SafeTestsets

# Debug use only
print_bigraph(sys, vars, edges) = print_bigraph(stdout, sys, vars, edges)
function print_bigraph(io::IO, sys, vars, edges)
    println(io, "Equations:")
    eqs = equations(sys)
    foreach(x->println(io, x), [i => eqs[i] for i in 1:length(eqs)])
    for (i, edge) in enumerate(edges)
        println(io, "\nEq $i has:")
        print(io, '[')
        for e in edge
            print(io, "$(vars[e]), ")
        end
        print(io, ']')
    end
    return nothing
end

@safetestset "Index Reduction & SCC" begin include("index_reduction.jl") end
