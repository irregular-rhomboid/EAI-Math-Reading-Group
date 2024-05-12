using DynamicalSystems
using OrdinaryDiffEq
using GLMakie



# the second range is a convenience for intermittency example of logistic
rrange = 1:0.001:4.0
# rrange = (rc = 1 + sqrt(8); [rc, rc - 1e-5, rc - 1e-3])

lo = Systems.logistic(0.4; r = rrange[1])
fig1 = interactive_cobweb(lo, rrange, 5)

fig2 = interactive_orbitdiagram(lo, 1, 0.0, 4.0; u0 = [0.5])



lorenz = Systems.lorenz()
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)

ds = CoupledODEs(lorenz, diffeq)


u0s = [
    [0.0,1.0,1.0],
    [0.0,1.0,1.01],
]

fig3, dsob = interactive_trajectory_timeseries(
    lorenz,
    [1,2,3],
    u0s;
    statespace_axis = true, Δt = 0.01,
    tail = 1000,
    figure = (size = (1100, 650),)
)

fig3

u0 = [10 .*rand(3) for i in 1:10]

trs = [trajectory(ds, 50000, u0)[1][:, SVector(1,2,3)] for u0 ∈ u0s]
j = 2 # the dimension of the plane

interactive_poincaresos_scan(trs, j;
    linekw = (transparency = true,),
    scatterkw = (markersize=3,)
    )

tent = Systems.tentmap(0.2,1.0)

function tent_rule(u, p, n)
    x = u[1]
    μ = p[1]
    y = x < 0.5 ? μ*x : μ*(1-x)
    return SVector(y)
end

tent = DeterministicIteratedMap(tent_rule, 0.25, [2.0])

fig4 = interactive_cobweb(tent, 0.0:0.01:4.0, 5)
