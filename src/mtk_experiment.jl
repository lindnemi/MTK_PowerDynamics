using ModelingToolkit, OrdinaryDiffEq, Plots

@parameters t P D σ
@variables θ(t) ω(t) c(t)
@derivatives Dt'~t

eqs1 = [Dt(θ) ~ ω,
       Dt(ω) ~ P - D * ω - σ * sin(θ - c) ]

kuramoto1 = ODESystem(eqs1, pins=[c], name=:kuramoto1)


eqs2 = [Dt(θ) ~ ω,
       Dt(ω) ~ P - D * ω - σ * sin(θ - c) ]

kuramoto2 = ODESystem(eqs2, pins=[c], name=:kuramoto2)


connections = [ kuramoto2.c ~ kuramoto1.θ; kuramoto1.c ~ kuramoto2.θ]

sys = ODESystem(Array{Equation, 1}(), t, [], [], observed = connections, systems=[kuramoto1, kuramoto2])

flattened_system = ModelingToolkit.flatten(sys) # Turns the composite system into one
afs = alias_elimination(flattened_system) # Connects the pins and the observations, and also removes other equalities.

u0 = [kuramoto1.θ => 0.,
      kuramoto1.ω => 1.,
      kuramoto2.θ => .5,
      kuramoto2.ω => -1.]

p = [kuramoto1.P => 1.,
     kuramoto2.P => -1.,
     kuramoto1.D => .75,
     kuramoto2.D => .75,
     kuramoto1.σ => 5,
     kuramoto2.σ => 5]

tspan = (0., 15.)

prob = ODEProblem(afs, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, vars=[2,4])
