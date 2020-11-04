using Test
using ModelingToolkit
using StructuralTransformations: find_linear_coeffs

@parameters t
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

eq5 = -15 + y + 2x ~ 5 + sin(x+y) + z
dict = find_linear_coeffs(eq5)
@test dict == Dict(1 => 20, z => 1)
