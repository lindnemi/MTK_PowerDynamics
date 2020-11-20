## Setup

using ModelingToolkit, OrdinaryDiffEq

## Our package provides a standard time that everything is in reference to:

struct SystemBase
  D
  t
  function SystemBase()
    @parameters t
    @derivatives D'~t
    new(D, t)
  end       
end

standard_base = SystemBase()

## and components:

function PhaseNodeBase(;b::SystemBase=standard_base)
  @variables ϕ(b.t), ω(b.t), c(b.t)
  ODESystem([b.D(ϕ) ~ ω + c],
  b.t,[ϕ, ω, c],[],pins=[ω],name=:PNBase)
end

function ConstantVelocity(;velocity = 1., b::SystemBase=standard_base, name=:ConstVelo)
  @variables ω(b.t)
  ODESystem([ω ~ velocity],
  b.t,[ω],[],name=name)
end

## The user then can initialize the components:

pnb = PhaseNodeBase()
cv = ConstantVelocity()
cv2 = ConstantVelocity(name = :CV2)

## And connect them (is there a better way to call ODESystem here?)

connected = ODESystem(Array{Equation, 1}(),standard_base.t,[],[],
  observed=[pnb.ω ~ cv.ω, pnb.c ~ cv2.ω],
  systems=[pnb, cv, cv2])

flattened_system = ModelingToolkit.flatten(connected) # Turns the composite system into one
aliased_flattened_system = alias_elimination(flattened_system) # Connects the pins and the observations, and also removes other equalities.

##

println(aliased_flattened_system.states)
println(aliased_flattened_system.eqs)
println(aliased_flattened_system.observed)
println(aliased_flattened_system.pins)

##

#=
Questions and observations:
- Is there a difference for flatten and alias_elimination whether an equation
is in eqs or in observed?
- Pins can be states, or not.
  - It seems like connecting to a pin when its not a state doesn't work.
  - Should pins always be states?
- It seems like making a variable a pin doesn't change the behaviour of anything (?).
- What happens to unconnected pins when we flatten/alias remove?
  - Answer: removing the aliases fails with a BoundsError error.
  - Actually the pins are already eliminated for the flattened_system
- It's tough to programmatically handle arrays of terms
(e.g. you can't filter them as comparing terms provides
you a comparison term rather than a bool), is it documented
somewhere how to acutally work with terms?
- Can we make sure the connected system is sane? E.g. only observed variables connected to pins?
- The terms are namespaced when the systems are put into a joined system.
What happens if two systems with the same name are connected?
  - This fails with a misleading error: The differential variable ... is not unique. (misleading because in my test case the non-unique variable is not differential)
- If we have equalities, how is it decided which term ends up as state and which as alias?
=#

## Below is a sketch for how nake_vertex could look.

function make_vertex(vertex_system, vertex_variables, edge_variables, aggregator)
  """Make an ODEVertex based on vertex_system, with the following convention:
  f(dx, x, e, p, t) where for dx and x the vertex variables are in
  order at the beginning of the vector, and the edge variables are in the order
  given by edge_variables. The e is calculated using the aggregator: e = aggregator(e_s, e_d).
  """
  flattened_system = ModelingToolkit.flatten(vertex_system)
  afs = alias_elimination(flattened_system)
  # order states so vertex_variables come first
  bf = build_function(afs, states, edge_variables, afs.ps, afs.ivs)
  f(dx, x, e_s, e_d, p, t) = bf(dx, x, aggregator(e_s, e_d), p, t)
  # get mass matrix
  # get dimensions
  # ODEVertex(f, dim, mm)
end
