using Test
using StructuralTransformations
using ModelingToolkit

# Define some variables
@parameters t L g
@variables x(t) y(t) w(t) z(t) T(t)
@derivatives D'~t

# Simple pendulum in cartesian coordinates
eqs = [D(x) ~ w,
       D(y) ~ z,
       D(w) ~ T*x,
       D(z) ~ T*y - g,
       0 ~ x^2 + y^2 - L^2]
pendulum = ODESystem(eqs, t, [x, y, w, z, T], [L, g], name=:pendulum)
edges, solvable_edges, fullvars, vars_asso = StructuralTransformations.sys2bigraph(pendulum; find_solvables=true)
@test edges == [[5, 3], [6, 4], [7, 1, 9], [8, 2, 9], [2, 1]]
@test solvable_edges == [[3, 5], [4, 6], [7], [8], Int64[]]
@test isequal(fullvars, [x, y, w, z, D(x), D(y), D(w), D(z), T])
@test vars_asso == [5, 6, 7, 8, 0, 0, 0, 0, 0]

#=
using StructuralTransformations: find_linear_coeffs

@parameters t
@derivatives D'~t
@variables x y z
eq = sin(z) + z - x ~ x + y + z + 2x + z
dict = find_linear_coeffs(eq)
@test dict == Dict(x=>4, y=>1)

eq2 = sin(z) + z - x ~ x + y + z + 2x + z + x*y
dict = find_linear_coeffs(eq2)
@test isempty(dict)

eq3 = x ~ x*y + y
dict = find_linear_coeffs(eq3)
@test isempty(dict)

eq3 = x ~ z + x*y + y + 3z
dict = find_linear_coeffs(eq3)
@test dict == Dict(z => 4)

eq4 = -15 + y + 2x ~ 5 + sin(x)
dict = find_linear_coeffs(eq4)
@test dict == Dict(1 => 20, y => -1)

eq5 = D(x) - 15 + y + 2x ~ 5 + sin(x+y) + z
dict = find_linear_coeffs(eq5)
@test dict == Dict(1 => 20, z => 1, D(x) => -1)
=#
