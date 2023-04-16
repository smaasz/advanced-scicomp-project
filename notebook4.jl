### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 25c86226-be8b-11ed-0611-dd2967a6d7d8
begin
	using VoronoiFVM
	using SimplexGridFactory
	using ExtendableGrids
	using Triangulate
	using GridVisualize
	using PlutoVista
	default_plotter!(PlutoVista)
	using PyPlot
	using PlutoUI
	using BenchmarkTools
	using Printf
	using VoronoiFVMDiffEq
	using DifferentialEquations
end

# ╔═╡ 778994ac-d9ee-4482-80f1-ea76eec6c6d5
using LinearAlgebra

# ╔═╡ 40dc1e16-ecf9-435b-8272-c94bc303bec5
md"""
## Grid
"""

# ╔═╡ ceebef24-a96c-4253-955e-6ed7b6f9f070
begin
	L = 10.0
	H = 1.0
end;

# ╔═╡ 1edfff44-8085-4d25-99cd-0a65d28e441d
;

# ╔═╡ 8206d5f8-b953-4bc2-b9ce-758d62ab22cf
function get_grid(; 
				isstructured 	= true,
				isrefined 		= false,
				xstep 			= 0.2,
				ystep 			= 0.1,
				maxvolume 		= 0.01)
	
	if isstructured
		X 		= 0.0:xstep:L
		Y 		= 0.0:ystep:H
		grid 	= simplexgrid(X, Y)
		
		bfacemask!(grid, [0, H], [L, H], 1)
		bfacemask!(grid, [0.2 * L, 0.0], [0.8 * L, 0.0], 3)
	else
		builder = let
			b = SimplexGridBuilder(Generator=Triangulate)
		
			# nodes
			p1 = point!(b, 0.0, 0.0)
			p2 = point!(b, L, 0.0)
			p3 = point!(b, L, H)
			p4 = point!(b, 0.0, H)
			p5 = point!(b, 0.2 * L, 0.0)
			p6 = point!(b, 0.8 * L, 0.0)
		
			# cell
			cellregion!(b, 1)
			regionpoint!(b, 0.5 * L, 0.5 * H)
		
			# faces
			facetregion!(b, 1)
			facet!(b, p1, p5)
			facet!(b, p6, p2)
			facet!(b, p3, p4)
		
			facetregion!(b, 2)
			facet!(b, p2, p3)
		
			facetregion!(b, 3)
			facet!(b, p5, p6)
		
			facetregion!(b, 4)
			facet!(b, p4, p1)

			if isrefined
				function unsuitable(x1, y1, x2, y2, x3, y3, area)
			        bary_y 	= (y1 + y2 + y3) / 3.0
			        dy 		= bary_y
			        area 	> min(maxvolume, 0.1 * max(0.01, dy^2))
				end
				options!(b, unsuitable = unsuitable)
			else
				options!(b, maxvolume = maxvolume)
			end
			
			b
		end
		grid = simplexgrid(builder; refine = false)
	end
	return grid
end;

# ╔═╡ f071e323-1e76-4774-933f-7178990e9331


# ╔═╡ 4f0b3a90-9c86-4602-84f7-49fd2407cd8e


# ╔═╡ 456af0e4-b93c-4741-a54a-b4454a6bffee
md"""
- ``\Gamma_{\text{cat}}`` has the region index 3
- ``\Gamma_{\text{out}}`` has the region index 2
- ``\Gamma_{\text{in}}`` has the region index 4
"""

# ╔═╡ b3d068b7-13dc-4bf6-b33f-87cd70efd50a
md"""
Apply local refinement at the bottom: $(@bind isrefined PlutoUI.CheckBox())
"""

# ╔═╡ 2950bf88-76d3-4bce-bfe7-054eee34e9d6
md"""
Use structured grid: $(@bind isstructured PlutoUI.CheckBox())
"""

# ╔═╡ 8ec02aea-816f-4500-942e-9f2a4b7c6712
grid = get_grid(; isstructured = isstructured, isrefined = isrefined)

# ╔═╡ d5748251-d585-445e-b64a-eaf2ea795243
gridplot(grid, resolution=(600, 400), legend=:rt, linewidth=0.5, zoom=2.5)

# ╔═╡ d21a4876-144a-4fb9-bbb1-6d3b210012a1
bgrid = subgrid(grid, [3], boundary=true);

# ╔═╡ 839219ec-fbf9-43c6-88fe-b2122fa0f9ce
grid2= uniform_refine(grid)

# ╔═╡ 4028515a-3a7f-4381-831c-61d2a3ae8290


# ╔═╡ 96bcf4da-5562-4b58-a48e-e63a1d0f1e23
gridplot(grid2, resolution=(600, 400), legend=:rt, linewidth=0.5, zoom=2.5)

# ╔═╡ 92d2c851-86d4-411d-8908-d6c026dd3755


# ╔═╡ 77f61c64-39bb-46d0-afdb-615d51ccf653


# ╔═╡ b88df152-2ef2-4547-be80-d320197b46ac


# ╔═╡ ac3dc830-20b6-405e-9584-55f6bee0eecb


# ╔═╡ ac48c632-bd37-4faf-a30f-54a5f0d3dbc5


# ╔═╡ d30b26b7-6c3d-4d15-8263-0bf0373214d2


# ╔═╡ 763004cf-39e1-49f0-abe8-e4b207969f86
md"""
## Physics
"""

# ╔═╡ 567d36a0-d082-4fed-8385-5267962c6e71
begin
	const dA = 1.0
	const dB = 0.001
	
	const kp_AC = 100.0
	const km_AC = 1.0
	const kp_BC = 0.1
	const km_BC = 1.0
end;

# ╔═╡ 5fbe7493-9b59-471d-8807-518e95880526
begin
	idA = 1
	idB = 2
	idC = 3
end;

# ╔═╡ 6200ebb3-74ed-418e-82d2-1d1539c85808
function get_physics(grid;
			velo 		= "none",
			vx 			= 0.5,
			fluxscheme 	= "upwind",
			S 			= 0.1)

	# velocity field
	if velo == "Hagen-Poisseuille"
		v = (x, y) -> [6 * vx * (1 - y / H) * y / H, 0.0]
	elseif velo == "uniform"
		v = (x, y) -> [vx, 0.0]
	else
		v = (x, y) -> [0.0, 0.0]
	end

	# flux discretization
	# think there is a error in vfvm_formfactors.jl line 378 -> multiply with -1
	evelo = -edgevelocities(grid, v) 
	
	function flux_expfit!(f, u, edge)
		vh = evelo[edge.index]
		f[idA] = dA * (fbernoulli(-vh / dA) * u[idA, 1] - fbernoulli(vh / dA) * u[idA, 2])
		f[idB] = dB * (fbernoulli(-vh / dB) * u[idB, 1] - fbernoulli(vh / dB) * u[idB, 2])
	end
	function flux_upwind!(f, u, edge)
		vh = evelo[edge.index]
		f[idA] = dA * (u[idA, 1] - u[idA, 2]) + vh * (vh > 0.0 ? u[idA, 1] : u[idA, 2])
		f[idB] = dB * (u[idB, 1] - u[idB, 2]) + vh * (vh > 0.0 ? u[idB, 1] : u[idB, 2])
	end
	function flux_central!(f, u, edge)
		vh = evelo[edge.index]
		f[idA] = dA * (u[idA, 1] - u[idA, 2]) + vh * 0.5 * (u[idA, 1] + u[idA, 2])
		f[idB] = dB * (u[idB, 1] - u[idB, 2]) + vh * 0.5 * (u[idB, 1] + u[idB, 2])
	end
	
	if fluxscheme == "expfit"
		flux! = flux_expfit!
	elseif fluxscheme == "central"
		flux! = flux_central!
	else
		flux! = flux_upwind!
	end

	# storage
	function storage!(f, u, node)
		f[idA] = u[idA]
		f[idB] = u[idB]
	end

	# bcondition
	R_AC(uA, uC) = kp_AC * uA * (1 - uC) - km_AC * uC
	R_BC(uB, uC) = kp_BC * uB * (1 - uC) - km_BC * uC
	function bcondition!(f, u, bnode)
		boundary_dirichlet!(f, u, bnode, species=idA, region=2, value=0.0)
		v = ramp(bnode.time, du = (0, 1.0), dt = (0, 1.0e-2))
		boundary_dirichlet!(f, u, bnode, species=idA, region=4, value=v)
		#boundary_dirichlet!(f, u, bnode, species=idA, region=4, value=1.0)
		#boundary_neumann!(f, u, bnode, species=idA, region=1, value=0.0)
	
		boundary_dirichlet!(f, u, bnode, species=idB, region=2, value=0.0)
		boundary_dirichlet!(f, u, bnode, species=idB, region=4, value=0.0)
		#boundary_neumann!(f, u, bnode, species=idB, region=1, value=0.0)
	
		if bnode.region == 3
			f[idA] = S * R_AC(u[idA], u[idC])
			f[idB] = S * R_BC(u[idB], u[idC])
			f[idC] = S *(- R_AC(u[idA], u[idC]) - R_BC(u[idB], u[idC]))
		end
	end

	# bstorage
	function bstorage!(f, u, bnode)
		if bnode.region == 3
			f[idC] = u[idC]
		end
	end

	return flux!, storage!, bcondition!, bstorage!
end;

# ╔═╡ 3a0d7dcb-f85e-451b-95db-cc06599b276a
md"""
__Select velocity field__: $(@bind velo PlutoUI.Select(["none", "uniform", "Hagen-Poisseuille"], default="uniform"))
"""

# ╔═╡ e072661b-db1e-4aed-85ac-38432be7ae7d
md"""
__Choose horizontal velocity__: $(@bind vx PlutoUI.Slider(0.0:0.1:1.0, show_value=true, default=0.5))
"""

# ╔═╡ 957b2e21-7a44-44a4-b42e-40573288e7ca
md"""
__Choose a discretization of the flux__: $(@bind fluxscheme PlutoUI.Select(["upwind", "expfit", "central"]))
"""

# ╔═╡ 10077a0a-6c1f-4e73-873a-3acb95d4c918
md"""
__S__: $(@bind S PlutoUI.Slider(0.0:0.1:1.0, show_value=true, default=0.3))
"""

# ╔═╡ 6e4a8cb7-c231-4f29-bf2f-9cfa28cd1d1d
flux!, storage!, bcondition!, bstorage! = get_physics(grid; velo = velo, vx = vx, fluxscheme = fluxscheme, S = S);

# ╔═╡ d81771d6-3d0d-4c07-ba3e-d8a390e999bd
md"""
## Solving the System
"""

# ╔═╡ 34e337cc-75d8-4a8b-8f4e-1662b1045a9c
function get_system(grid, flux!, storage!, bcondition!, bstorage!)
	sys = VoronoiFVM.System(grid; unknown_storage = :sparse,
                                  bstorage = bstorage!,
                                  flux = flux!,
                                  storage = storage!,
								  bcondition=bcondition!)

	enable_species!(sys, idA, [1])
	enable_species!(sys, idB, [1])
	enable_boundary_species!(sys, idC, [3])
	return sys
end;

# ╔═╡ 713d7370-f821-4740-a4f0-6e08e50eff7a
function heterogeneous_reaction(; 	isstructured 	= false,
									isrefined 		= false,
									maxvolume       = 0.01,
								  	velo 			= "none",
								  	vx 				= 0.5,
									fluxscheme 		= "upwind",
									S 				= 0.1,
									solvermethod 	= nothing,
									Δt 				= 0.00001,
									Δu_opt 			= 1.0,
									damp_initial 	= 0.5,
									tend 			= 100)

	grid 									= get_grid(; isstructured = isstructured, 													 isrefined = isrefined, maxvolume = maxvolume)
	flux!, storage!, bcondition!, bstorage! = get_physics(grid; 
																	velo = "uniform", 
																	vx = vx, 
																	fluxscheme = fluxscheme, 
																	S = S)
	sys 									= get_system(grid, flux!, storage!, bcondition!, bstorage!)
	inival 									= unknowns(sys; inival=0.0)
	with_terminal() do
		println("length of initial value= $(length(inival))")
	end
	
	if isnothing(solvermethod)
		# solver control
		control 				= VoronoiFVM.SolverControl()
		control.Δt_min 			= 0.01 * Δt
		control.Δt 				= Δt
		control.Δt_max 			= 0.1 * tend
		control.Δu_opt 			= Δu_opt
		control.damp_initial 	= damp_initial
		control.log 			= true
		
		tsol 			= solve(sys; inival=inival, times=[0, tend], control=control)
		#ssol 			= solve(sys; inival=inival)
	else
		problem 		= ODEProblem(sys, inival, (0, tend))
		odesolution 	= DifferentialEquations.solve(problem, solvermethod())
		tsol 			= reshape(odesolution, sys)
	end
	return sys, tsol#, ssol
end;

# ╔═╡ 9b01f154-7cc4-455f-9636-2bc337458134
function heterogeneous_reaction_grid(grid; 	isstructured 	= false,
									isrefined 		= false,
									maxvolume       = 0.01,
								  	velo 			= "none",
								  	vx 				= 0.5,
									fluxscheme 		= "upwind",
									S 				= 0.1,
									solvermethod 	= nothing,
									Δt 				= 0.00001,
									Δu_opt 			= 1.0,
									damp_initial 	= 0.5,
									tend 			= 100)
	flux!, storage!, bcondition!, bstorage! = get_physics(grid; 
																	velo = "uniform", 
																	vx = vx, 
																	fluxscheme = fluxscheme, 
																	S = S)
	sys 									= get_system(grid, flux!, storage!, bcondition!, bstorage!)
	inival 									= unknowns(sys; inival=0.0)
	with_terminal() do
		println("length of initial value= $(length(inival))")
	end
	
	if isnothing(solvermethod)
		# solver control
		control 				= VoronoiFVM.SolverControl()
		control.Δt_min 			= 0.01 * Δt
		control.Δt 				= Δt
		control.Δt_max 			= 0.1 * tend
		control.Δu_opt 			= Δu_opt
		control.damp_initial 	= damp_initial
		control.log 			= true
		
		tsol 			= solve(sys; inival=inival, times=[0, tend], control=control)
		#ssol 			= solve(sys; inival=inival)
	else
		problem 		= ODEProblem(sys, inival, (0, tend))
		odesolution 	= DifferentialEquations.solve(problem, solvermethod())
		tsol 			= reshape(odesolution, sys)
	end
	return sys, tsol#, ssol
end;

# ╔═╡ 17a230d4-b536-428e-a29b-7d53e0dae5f4
# ╠═╡ show_logs = false
sys, tsol = heterogeneous_reaction( isstructured 	= isstructured,
									isrefined 		= isrefined,
								  	velo 			= velo,
								  	vx 				= vx,
									fluxscheme 		= fluxscheme,
									S 				= S,
									#solvermethod 	= QNDF1,
									Δt 				= 0.00001,
									Δu_opt 			= 1.0,
									damp_initial 	= 0.5,
									tend 			= 100);

# ╔═╡ 0d6da6c6-63a9-4269-a401-08a41587c282
md"""
## Visualization
"""

# ╔═╡ c84d6e53-885c-48ba-9bc6-e85cc87bdefd
md"""
### Steady state solution
"""

# ╔═╡ fe737c36-b08b-465a-b6c7-e8656badc5a1
md"""
### Transient Solution
"""

# ╔═╡ 64f8b862-fa71-4479-bbfe-4d300388a6ba
@bind step PlutoUI.Slider(1:length(tsol.t), default=length(tsol.t), show_value=true)

# ╔═╡ 50637189-a467-4f1b-bf22-606d06ea0981
begin
	p = GridVisualizer(; Plotter = PlutoVista, layout = (3, 1), legend=:lt)
	scalarplot!(p[1,1], grid, tsol[idA, :, 1]; limits=(0.0, 1.0), zoom = 4.5, markerevery=2, label="Distribution of species A")
	scalarplot!(p[2,1], grid, tsol[idA, :, 1]; limits=(0.0, 5.0), zoom = 4.5, title="Distribution of the species B")
	scalarplot!(p[3, 1], bgrid, view(tsol[idA, :, 1], bgrid); limits=(0.0, 1.0), title="Distribution of the surface species C", ylabel="u_C")
	reveal(p);
end

# ╔═╡ 2b9f5ab8-2c99-43c0-9189-92f3de7f7da6
begin
	scalarplot!(p[1,1], grid, tsol[idA, :, step]; clear = true, limits=(0.0, 1.0), label="Distribution of species A")
	scalarplot!(p[2,1], grid, tsol[idB, :, step]; clear = true, limits=(0.0, 5.0), label="Distribution of the species B")
	scalarplot!(p[3, 1], bgrid, view(tsol[idC, :, step], bgrid); clear=true, limits=(0.0, 1.0), title="Distribution of species A")
	reveal(p);
end

# ╔═╡ 3cc2d23e-45b8-4725-99a2-8a980b239ccc
max(tsol[idB, :, step]...)

# ╔═╡ 13bef003-732f-48af-b7b0-aba678189ab9
sys.grid

# ╔═╡ 4405f2d9-69f5-4331-aafa-27991ec52a77
tsol[idB, :, end]

# ╔═╡ 0c7b040f-6217-44be-8a3f-629b715e8a8c
history(sys)

# ╔═╡ a48a0989-cd6e-4644-a6ab-e204dd7a403b
md"""
## Production rate of ``B``
"""

# ╔═╡ 8afe233d-1bce-435d-a689-aeddfb93578d
bfluxesB = VoronoiFVM.integrate(sys, bcondition!, tsol.u[end]; boundary=true)[idB, :]

# ╔═╡ 19f07991-ff8a-499e-9d63-b73cde9fbf74
pb = -bfluxesB[3]

# ╔═╡ f95ed529-fe7e-443f-9578-f1b8da779dad
Markdown.parse("__The production rate of the species B is ``p_B`` = $(@sprintf "%.4f" pb) in the stationary state for horizontal velocity ``v_x = ``$(vx) and fraction of available catalytic sites ``S = ``$(S).__")

# ╔═╡ 8cc09ef6-fd13-47b8-a6d3-f3a49c643913
begin 
	tf 		= TestFunctionFactory(sys)
	Tout 	= testfunction(tf, [2], [4])
	Iout 	= integrate(sys, Tout, tsol[:, :, end], tsol[:, :, end - 1], tsol.t[end] - tsol.t[end - 1])
	#Jout = integrate(sys, Tout, ssol)
	#Jin = integrate(sys, Tin, ssol)
	Iout, Iout[1]+Iout[2] #Jout[idB]
end

# ╔═╡ 25d16219-af4a-434a-a116-236363494a69
begin
	function mayb!(f, u, bnode)
		if bnode.region == 3
			f[idB] = Tout[bnode.index] * (-S) * (kp_BC * u[idB] * (1 - u[idC]) - km_BC * u[idC])
		end
	end
	integrate(sys, mayb!, tsol[:, :, end]; boundary = true)[idB, 3]
end

# ╔═╡ f3ab8364-8f10-4c02-a7d6-854f055d11c5
md"""
## Tests
"""

# ╔═╡ 14c04f9b-27b2-4638-a6ee-414cb6bd993e
function test_productionrate(vxs, Ss; 	isstructured 	= false,
										isrefined 		= false,
										velo 			= "uniform",
										fluxscheme 		= "upwind",
										Δt 				= 0.00001,
										Δu_opt 			= 1.0,
										damp_initial 	= 0.5,
										tend 			= 100)
	
	function bflux!(f, u, bnode)
		if bnode.region == 3
				f[idB] = -(kp_BC * u[idB] * (1 - u[idC]) - km_BC * u[idC])
		end
	end	

	fluxesB = Matrix(undef, length(vxs), 0)
	for S in Ss
		col = []
		for vx in vxs
			sys, tsol = heterogeneous_reaction( isstructured 	= isstructured,
												isrefined 		= isrefined,
											  	velo 			= velo,
											  	vx 				= vx,
												fluxscheme 		= fluxscheme,
												S 				= S,
												Δt 				= Δt,
												Δu_opt 			= Δu_opt,
												damp_initial 	= damp_initial,
												tend 			= tend)
			fluxB = integrate(sys, bflux!, tsol.u[end]; boundary = true)[idB, 3]
			push!(col, fluxB)
		end
		fluxesB = hcat(fluxesB, col)
	end
	return fluxesB
end;

# ╔═╡ e9b07b6f-2b47-478a-a65e-3473e6554ced
md"""
__Run tests__: $(@bind runtests PlutoUI.CheckBox())
"""

# ╔═╡ fa4a65d8-08d1-4704-ae33-77dfb4e80268
# ╠═╡ show_logs = false
if runtests
	vxs 		= collect(0.0:0.2:1.0)
	Ss 			= [0.1, 0.3, 0.6]
	fluxesB 	= test_productionrate(vxs, Ss)
	test_plot 	= PlutoVista.PlutoVistaPlot(; resolution = (300, 300), 
							title 	= "Production rate of B in dependence of horizontal velocity",
							xlabel 	= "horizontal velocity vx",
							ylabel 	= "Production rate of B",
							legend 	= :lt)

	# plot the results
	for (i, col) in enumerate(eachcol(fluxesB))
		PlutoVista.plot!(test_plot, vxs, col; 
			markertype = :square, 
			limits = (0.0, 6.0), 
			label="S = $(Ss[i])")
	end
	test_plot
end

# ╔═╡ 9866ffec-e48f-45de-93c3-bba2d4dcbaf8
md"""### Convergence Test"""

# ╔═╡ 3ec7c7ee-a353-44c8-8f8d-2064020542d1
function test_convergence(n_grid, initial_volume; 
										norm_idx        = 2,
										isstructured 	= false,
										isrefined 		= false,
										velo 			= "uniform",
										fluxscheme 		= "upwind",
										Δt 				= 0.00001,
										Δu_opt 			= 1.0,
										damp_initial 	= 0.5,
										tend 			= 100)

	grid_range = range(1,n_grid)
	col=[]
	volume = initial_volume

	function bflux!(f, u, bnode)
		if bnode.region == 3
				f[idB] = -(kp_BC * u[idB] * (1 - u[idC]) - km_BC * u[idC])
		end
	end	
	
	for n in grid_range
		# evaluation at each grid point
		sys, tsol = heterogeneous_reaction( isstructured 	= isstructured,
												isrefined 		= isrefined,
												maxvolume 		= volume,
											  	velo 			= velo,
											  	vx 				= vx,
												fluxscheme 		= fluxscheme,
												S 				= S,
												Δt 				= Δt,
												Δu_opt 			= Δu_opt,
												damp_initial 	= damp_initial,
												tend 			= tend)
		# norm vector
		grid = sys.grid
		bgrid = subgrid(grid, [3], boundary=true);
		# evalA = tsol[idA, :, end]
		# evalB = tsol[idB, :, end]
		# evalC = view(tsol[idA, :, end], bgrid)

		# normA = norm(evalA,norm_idx)
		# normB = norm(evalB,norm_idx)
		# normC = norm(evalC,norm_idx)
		# values = [normA, normB, normC]
		
		fluxB = integrate(sys, bflux!, tsol.u[end]; boundary = true)[idB, 3]
		
		# save value
		# push!(norm_vector, values)
		push!(col, fluxB)
		
		# volume berechnen 
		volume = volume*0.75

		# prod rate
	end 
	return col	
end;

# ╔═╡ 64160199-56a2-4a07-91c0-2a7f890f7f6b
function test_convergence2(n_grid, initial_volume;
										norm_idx        = 2,
										isstructured 	= false,
										isrefined 		= false,
										velo 			= "uniform",
										fluxscheme 		= "upwind",
										Δt 				= 0.00001,
										Δu_opt 			= 1.0,
										damp_initial 	= 0.5,
										tend 			= 100)
	#initial system
	init_sys, init_tsol = heterogeneous_reaction( isstructured 	= isstructured,
												isrefined 		= isrefined,
												maxvolume 		= initial_volume,
											  	velo 			= velo,
											  	vx 				= vx,
												fluxscheme 		= fluxscheme,
												S 				= S,
												Δt 				= Δt,
												Δu_opt 			= Δu_opt,
												damp_initial 	= damp_initial,
												tend 			= tend)
	init_grid = init_sys.grid #intial grid with 860 points
	grid_range = range(1, n_grid) #how many times do we refine the grid
	idx_tsol = num_nodes(init_grid) #first values of the solutions that are common 
	values = [] #save the values
	for k in grid_range
		# refine the grid
		k_grid = uniform_refine(init_grid)
		# compute the new system 
		k_sys, k_tsol = heterogeneous_reaction_grid(k_grid; isrefined 		= false,
									maxvolume       = 0.1,
								  	velo 			= "none",
								  	vx 				= 0.5,
									fluxscheme 		= "upwind",
									S 				= 0.1,
									solvermethod 	= nothing,
									Δt 				= 0.00001,
									Δu_opt 			= 1.0,
									damp_initial 	= 0.5,
									tend 			= 100)
		
		# find the interesting values
		value = k_tsol[:,1:idx_tsol,end]
		# save the values
		push!(values, value)
		# update init_grid
		init_grid = k_grid
		init_tsol = k_tsol
	end
	error = values[1:n_grid-1] - values[2:n_grid]
end

# ╔═╡ 76761ecc-9ca4-4b8f-849c-4192d9d6cf3a
# ╠═╡ disabled = true
#=╠═╡
function evaluate_convergence(n_grid, initial_volume)
	values, idx_tsol = test_convergence2(n_grid, initial_volume)
	diff = values[1:idx_tsol-1] - values[2:idx_tsol]
	return diff
end
  ╠═╡ =#

# ╔═╡ fe0ed7cb-5671-49ae-9003-cfefdf0936d0
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
error_vector = test_convergence2(6, 1)
  ╠═╡ =#

# ╔═╡ 7790a1f9-df1e-484e-a401-c8ae27eed721
function compare_grids()
	# compute convergence test for structured
	struc_error = test_convergence2(7,1; isstructured = true)
	unstruc_error = test_convergence2(7,1; isstructured = false)
	struc_normA = []
	struc_normB = []
	unstruc_normA = []
	unstruc_normB = []
	# compute convergence test for unstructured
	for idx in range(1,5)
		strucA = norm(struc_error[idx][idA], Inf)
		unstrucA = norm(unstruc_error[idx], Inf)
		push!(struc_normA, strucA)
		push!(unstruc_normA, unstrucA)

		strucB = norm(struc_error[idx][idB], Inf)
		unstrucB = norm(unstruc_error[idB], Inf)
		push!(struc_normB, strucB)
		push!(unstruc_normB, unstrucB)
		
	# plot to compare the errors
	end
	x_axis = collect(range(1,5))
	p = GridVisualizer(; Plotter = PlutoVista, layout = (2, 2), legend=:lt)
	scalarplot!(p[1,1], x_axis, struc_normA; label="Structured grid: Error for species A")
	scalarplot!(p[1,2], x_axis, unstruc_normA; label="Unstructured grid: Error for species A")
	scalarplot!(p[2,1], x_axis, struc_normB; label="Structured grid: Error for species B")
	scalarplot!(p[2,2], x_axis, unstruc_normB; label="Unstructured grid: Error for species B")
	reveal(p)
	return struc_normA, unstruc_normA, struc_normB, unstruc_normB, struc_error, unstruc_error
end

# ╔═╡ d82a2da3-144c-4780-8e59-dac86620c45c


# ╔═╡ acb36288-ecd5-4ab6-8c85-e07b9f20bec9


# ╔═╡ 400a739f-6cc2-46ba-b2aa-9873bf109f15
# ╠═╡ show_logs = false
result = compare_grids()

# ╔═╡ bfad73fe-eaa5-4b68-a8ed-e70a75753c03


# ╔═╡ 4e4a8844-5548-4b33-b00c-2d423979cd30
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
sys2, tsol2 = heterogeneous_reaction_grid(grid2; isstructured 	= isstructured,
									isrefined 		= isrefined,
								  	velo 			= velo,
								  	vx 				= vx,
									fluxscheme 		= fluxscheme,
									S 				= S,
									#solvermethod 	= QNDF1,
									Δt 				= 0.00001,
									Δu_opt 			= 1.0,
									damp_initial 	= 0.5,
									tend 			= 100)
  ╠═╡ =#

# ╔═╡ c1ed66a0-7d86-4c2d-a8a7-03d22008bc9e


# ╔═╡ b0692cf3-8eed-4b85-abab-ecf139c3b2c3
# ╠═╡ disabled = true
#=╠═╡
tsol[idB]
  ╠═╡ =#

# ╔═╡ c659d686-2917-4035-9b12-d1f0327fa752
# ╠═╡ disabled = true
#=╠═╡
grid[Coordinates]
  ╠═╡ =#

# ╔═╡ 00f22fd5-d247-4c88-a9ae-c2ef8ee770bd
# ╠═╡ disabled = true
#=╠═╡
grid2[Coordinates][1:2,1:860] == grid[Coordinates]
  ╠═╡ =#

# ╔═╡ 1c5aef0a-9f5c-4471-9956-3ef457eef49b
# ╠═╡ disabled = true
#=╠═╡
begin
	run_convergencetests = true
	n_grid = 15
	initial_volume = 1
end
  ╠═╡ =#

# ╔═╡ a3478c2b-28ba-4d3b-ba73-4772137b5359
# ╠═╡ disabled = true
#=╠═╡
abs([-3, 2])
  ╠═╡ =#

# ╔═╡ 89ad2f8d-72bf-4dc2-9c5b-d3d19124ba2b
# ╠═╡ disabled = true
#=╠═╡
begin 
a = []
push!(a, 2)
push!(a,-3)
end
  ╠═╡ =#

# ╔═╡ 12c02070-4c24-4038-8e39-8781466c4e0c
# ╠═╡ disabled = true
#=╠═╡
collect(range(1,5))
  ╠═╡ =#

# ╔═╡ 642f9235-714a-404b-acdb-1f702412eba8
# ╠═╡ disabled = true
#=╠═╡
error_vector
  ╠═╡ =#

# ╔═╡ 0a2bbfb4-ce90-4442-9583-bfd76fdb2010
# ╠═╡ disabled = true
#=╠═╡
begin 
	t1 = tsol2[:,1:50, end]
	t2 = tsol2[:,51:100, end]
	e = t1-t2
	normA = norm(e[idA])
	normB = norm(e[idB])
	normC = norm(e[idC])
end
  ╠═╡ =#

# ╔═╡ bc7bfea9-68fd-4831-b262-82e6d38c632a
# ╠═╡ disabled = true
#=╠═╡
normC
  ╠═╡ =#

# ╔═╡ fb622e06-d9b4-4f58-a099-721d8cc4bc9d


# ╔═╡ 734993cd-062b-4883-ab1b-6c85f8302fa8


# ╔═╡ e6701190-0d1f-4c8f-99a2-5708bb32794a


# ╔═╡ 87019342-1089-4fd9-b088-6f57b56bec53


# ╔═╡ b26344e1-4229-4cff-a4a0-bb5bf40e5500
# ╠═╡ disabled = true
#=╠═╡
e
  ╠═╡ =#

# ╔═╡ 7bec601b-ac99-4446-a582-1415b83a754a
# ╠═╡ disabled = true
#=╠═╡

  ╠═╡ =#

# ╔═╡ a5a43c84-3671-4567-a040-f5c4c1ed2508
# ╠═╡ disabled = true
#=╠═╡
t1
  ╠═╡ =#

# ╔═╡ 0bfc2e2d-8009-4b03-b7f8-ace74d99013e
# ╠═╡ disabled = true
#=╠═╡
t2
  ╠═╡ =#

# ╔═╡ c4c860a5-e01f-4d9e-bbc3-0700cda36d1f


# ╔═╡ 9f5339aa-c96d-4b16-8f4c-c94a4a2e4a53


# ╔═╡ dbe648ea-f2b2-418b-a5b9-f2a9558473d1


# ╔═╡ 3d145c69-54d5-4175-abaf-c9f673ec2a88


# ╔═╡ c8d29792-23b1-4afe-891d-bf3a27063656


# ╔═╡ d4f79819-d23b-4683-9eae-691abae28697


# ╔═╡ c2688388-c1af-433d-bc90-844cdb32cf08


# ╔═╡ 76c2df03-3e5a-4fb0-9ea7-dd2c04cc9131


# ╔═╡ 5adcd720-0b19-402e-b638-9857b5a0db27


# ╔═╡ c96bd5c1-5f75-4412-a90c-4820454fcf48


# ╔═╡ c0a5ea95-01db-4428-bd7a-89af1d75752a


# ╔═╡ e0126cb0-3ae9-49d9-8753-235baf94b6a8


# ╔═╡ 023d053b-60e8-4e45-837b-18b754764f40


# ╔═╡ 6b3f828e-541d-4ee6-a5f9-acd4fcd2d790


# ╔═╡ 155ce2a2-28ad-4cdd-9230-cbf12f7cdc7d


# ╔═╡ 7ced64c4-5339-4e6f-998d-e8ea915695a4


# ╔═╡ 40f0251a-9ac2-4ce9-a180-e763b4380716


# ╔═╡ 3378bb7b-936d-44e6-8529-314bc33ed67f


# ╔═╡ 9d3d52e0-d64e-40d0-b149-3a4163b44070


# ╔═╡ 7896788e-4f37-4c99-96ce-a8006d3f6de6


# ╔═╡ 96ad3dbe-950d-4ecd-b919-211f0f36e549


# ╔═╡ a5621f61-7853-4fc9-9f46-aae0f1dbe360


# ╔═╡ ae9cce49-0303-4fc6-a5e0-eea616ab82de


# ╔═╡ b60b8e7a-24f8-49f0-9d01-fa5f3c53cdaa


# ╔═╡ fb401ab1-878a-41c8-b44c-806c409697da


# ╔═╡ 2800de62-c6a7-47e1-a977-d263f52624d5
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
if run_convergencetests
	p_rate = test_convergence(n_grid, initial_volume)
	range_grid = collect(range(1,n_grid))
	# A = [i[1] for i in norm_vec]
	# B = [i[2] for i in norm_vec]
	# C = [i[3] for i in norm_vec]
	q = PlutoVista.plot(range_grid, p_rate; limits = (0,4), title="Production rate vs. Grid size", ylabel="Production rate")
	# q = PlutoVista.plot(range_grid, A; title="Norm A", ylabel="norm")
	# PlutoVista.plot!(q, range_grid, B; title="Norm A", ylabel="norm")
	# PlutoVista.plot!(q, range_grid, C; title="Norm A", ylabel="norm")
	
end
  ╠═╡ =#

# ╔═╡ 41e2de15-dcbf-4340-9469-3ddc5184af49
md"""
## Benchmarks
"""

# ╔═╡ a0c6ec1a-bb4d-41e5-8689-04a9691796aa
# ╠═╡ disabled = true
#=╠═╡
function benchmark_solve(; 	velo 			= "uniform",
							vx 				= 0.5,
							S 				= 0.1,
							Δt 				= 0.00001,
							Δu_opt 			= 1.0,
							damp_initial 	= 0.5,
							tend 			= 100)
	
	suite = BenchmarkGroup()

	# solver control
	control 				= VoronoiFVM.SolverControl()
	control.Δt_min 			= 0.01 * Δt
	control.Δt 				= Δt
	control.Δt_max 			= 0.1 * tend
	control.Δu_opt 			= Δu_opt
	control.damp_initial 	= damp_initial
	
	# different flux functions
	suite["flux function"] = BenchmarkGroup(["spatial discretization", "flux"])
	for fluxscheme in ["upwind", "expfit", "central"]
		grid 									= get_grid()
		flux!, storage!, bcondition!, bstorage! = get_physics(grid; velo = velo, vx = vx, fluxscheme = fluxscheme, S = S)
		sys 									= get_system(grid, flux!, storage!, bcondition!, bstorage!)
		inival 									= unknowns(sys; inival=0.0)

		# create benchmark
		suite["flux function"][fluxscheme] 	= @benchmarkable solve($(sys), inival=$(inival), times=[0, $(tend)], control=$(control)) samples=5
	end

	run(suite, verbose = true, seconds = 60)
end
  ╠═╡ =#

# ╔═╡ 33cb6da8-d2c9-46a9-8667-19c06df0d7f7
md"""
__Run benchmarks__: $(@bind runbenchmarks PlutoUI.CheckBox())
"""

# ╔═╡ 62968b3b-f4be-44c2-a2bb-a5c858bbe97b
# ╠═╡ disabled = true
#=╠═╡
if runbenchmarks
	suite = benchmark_solve()
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
SimplexGridFactory = "57bfcd06-606e-45d6-baf4-4ba06da0efd5"
Triangulate = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
VoronoiFVMDiffEq = "c4d4532c-19f6-4375-9ced-3808a145cd37"

[compat]
BenchmarkTools = "~1.3.2"
DifferentialEquations = "~7.7.0"
ExtendableGrids = "~0.9.17"
GridVisualize = "~1.0.3"
PlutoUI = "~0.7.50"
PlutoVista = "~0.8.24"
PyPlot = "~2.11.1"
SimplexGridFactory = "~0.5.19"
Triangulate = "~2.2.0"
VoronoiFVM = "~1.2.1"
VoronoiFVMDiffEq = "~0.1.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "7a886aaf0ff352c3ebb5c479ebb5dd82c1bb4dc0"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "a69dbe3b376ace7d9eebe2db43216e8b52ba6da9"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.29.2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38911c7737e123b28182d89027f4216cfc8a9da7"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.3"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4aff5fa660eb95c2e0deb6bcdabe4d9a96bc4667"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.18"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SnoopPrecompile", "SparseArrays"]
git-tree-sha1 = "6ef8fc1d77b60f41041d59ce61ef9eb41ed97a83"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.18"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase", "SparseArrays"]
git-tree-sha1 = "ed8e837bfb3d1e3157022c9636ec1c722b637318"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.11.0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "2c144ddb46b552f72d7eafe7cc2f50746e41ea21"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "70232f82ffaab9dc52585e0dd043b5e0c6b714f1"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.12"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "e32a90da027ca45d84678b826fffd3110bb3fc90"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.8.0"

[[deps.Configurations]]
deps = ["ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "62a7c76dbad02fdfdaa53608104edf760938c4ca"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.17.4"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack"]
git-tree-sha1 = "89f3fbfe78f9d116d1ed0721d65b0b2cf9b36169"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.42.0"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DataStructures", "Distributions", "DocStringExtensions", "EnumX", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "Static", "StaticArraysCore", "Statistics", "Tricks", "TruncatedStacktraces", "ZygoteRules"]
git-tree-sha1 = "117b2d02e737aeefd58cd4a4803abecadd37c8cc"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.122.2"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "63b6be7b396ad395825f3cc48c56b53bfaf7e69d"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.26.1"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "2c4ed3eedb87579bfe9f20ecc2440de06b9f3b89"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.16.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "ac145e3d718157c679fc4febf2fcef73ec77b067"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.7.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "13027f188d26206b9e7b863036f87d2f2e7d013a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.87"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "698124109da77b6914f64edd696be8dccf90229e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "8b84876e31fa39479050e2d3395c4b3b210db8b0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "e1c40d78de68e9a2be565f0202693a158ec9ad85"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.11"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "Printf", "SnoopPrecompile", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "fb7dbef7d2631e2d02c49e2750f7447648b0ec9b"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.24.0"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

[[deps.ExproniconLite]]
deps = ["Pkg", "TOML"]
git-tree-sha1 = "c2eb763acf6e13e75595e0737a07a0bec0ce2147"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.7.11"

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "StaticArrays", "Test", "WriteVTK"]
git-tree-sha1 = "2921bf0ffab4c8b7eda6a36c7b06a0dde6df0137"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.9.17"

[[deps.ExtendableSparse]]
deps = ["DocStringExtensions", "ILUZero", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "Sparspak", "SuiteSparse", "Test"]
git-tree-sha1 = "bbb16c582df45544612cd703fe1b8179339ccd2b"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "1.0.1"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "d1248fceea0b26493fd33e8e9e8c553270da03bd"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c1293a93193f0ae94be7cf338d33e162c39d8788"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.9"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7072f1e3e5a8be51d525d64f63d3ec1287ff2790"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.11"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "03fcb1c42ec905d15b305359603888ec3e65f886"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.19.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.FuzzyCompletions]]
deps = ["REPL"]
git-tree-sha1 = "e16dd964b4dfaebcded16b2af32f05e235b354be"
uuid = "fb4132e2-a121-4a70-b8a1-d5b831dcdcc2"
version = "0.5.1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "0eb6de0b312688f852f347171aba888658e29f20"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "303202358e38d2b01ba46844b92e48a3c238fd9e"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.6"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "GridVisualizeTools", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "e806d85549e112f306002e8ad030a0c77efe9f4f"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "1.0.3"

[[deps.GridVisualizeTools]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "StaticArraysCore"]
git-tree-sha1 = "7c892c426f8d03a180366411566d0f6ac1790f6c"
uuid = "5573ae12-3b76-41d9-b48c-81d0b6e61cc5"
version = "0.3.0"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "827f29c95676735719f8d6acbf0a3aaf73b3c9e5"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.3.2"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "734fd90dd2f920a2f1921d5388dcebe805b262dc"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.14"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "432b5b03176f8182bd6841fbfc42c718506a2d5f"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.15"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.ILUZero]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b007cfc7f9bee9a958992d2301e9c5b63f332a90"
uuid = "88f59080-6952-5380-9ea5-54057fb9a43f"
version = "0.2.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "42c17b18ced77ff0be65957a591d34f4ed57c631"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.31"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "740c685ba3d7f218663436b2152041563c19db6e"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.6.1"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "764164ed65c30738750965d55652db9c94c59bfe"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "dd90aacbfb622f898a97c2a4411ac49101ebab8a"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.0"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "88b8f66b604da079a627b6fb2860d3704a6729a1"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.14"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "410fe4739a4b092f2ffe36fcb0dcc3ab12648ce1"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.2.1"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "EnumX", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "4a4f8cc7a59fadbb02d1852d1e0cef5dca3a9460"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.42.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "SpecialFunctions", "Static", "StaticArrayInterface", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "defbfba8ddbccdc8ca3edb4a96a6d6fd3cd33ebd"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.157"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MeshIO]]
deps = ["ColorTypes", "FileIO", "GeometryBasics", "Printf"]
git-tree-sha1 = "8be09d84a2d597c7c0c34d7d604c039c9763e48c"
uuid = "7269a6da-0436-5bbc-96c2-40638cbb6118"
version = "0.4.10"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "fc8c15ca848b902015bd4a745d350f02cf791c2a"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.2.0"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "eaa98afe2033ffc0629f9d0d83961d66a021dfcc"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3295d296288ab1a0a2528feb424b854418acff57"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.2.3"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "EnumX", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "a6000c813371cd3cd9cbbdf8a356fc3a97138d92"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "1.6.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "1903afc76b7d01719d9c30d3c7d501b61db96721"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.4"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "IfElse", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLNLSolve", "SimpleNonlinearSolve", "SimpleUnPack", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "9fb1f72106bfa1370006b90771cfbcce6c7468b6"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.49.4"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.Pluto]]
deps = ["Base64", "Configurations", "Dates", "Distributed", "FileWatching", "FuzzyCompletions", "HTTP", "HypertextLiteral", "InteractiveUtils", "Logging", "LoggingExtras", "MIMEs", "Markdown", "MsgPack", "Pkg", "PrecompileSignatures", "REPL", "RegistryInstances", "RelocatableFolders", "SnoopPrecompile", "Sockets", "TOML", "Tables", "URIs", "UUIDs"]
git-tree-sha1 = "cb1e1b21261e6aa3cf5547caebff0baa84403451"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.19.24"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PlutoVista]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "GridVisualizeTools", "HypertextLiteral", "Pluto", "UUIDs"]
git-tree-sha1 = "30675d4a579f50e60e14a72e36cc453610af7b76"
uuid = "646e1f28-b900-46d7-9d87-d554eb38a413"
version = "0.8.24"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "a0f1159c33f846aa77c3f30ebbc69795e5327152"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.4"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "0fe4e7c4d8ff4c70bfa507f0dd96fa161b115777"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.3"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

[[deps.PrecompileSignatures]]
git-tree-sha1 = "18ef344185f25ee9d51d80e179f8dad33dc48eb1"
uuid = "91cefc8d-f054-46dc-8f8c-26e11d7c5411"
version = "3.0.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "62f417f6ad727987c755549e9cd88c46578da562"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.95.1"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "92e7ca803b579b8b817f004e74b205a706d9a974"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.11.1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "7a1a306b72cfa60634f03a911405f4e64d1b718b"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.0"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "140cddd2c457e4ebb0cdc7c2fd14a7fbfbdf206e"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.3"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "9088515ad915c99026beb5436d0a09cd8c18163e"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.18"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "f139e81a81e6c29c40f1971c9e5309b09c03f2c3"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.6"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "cda0aece8080e992f6370491b08ef3909d1c04e7"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.38"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SnoopPrecompile", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "392d3e28b05984496af37100ded94dc46fa6c8de"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.91.7"

[[deps.SciMLNLSolve]]
deps = ["DiffEqBase", "LineSearches", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "2e1606c282fae6bd9aed4f159695774a44b9c75f"
uuid = "e9a6253c-8580-4d32-9898-8661bb511710"
version = "0.1.4"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "e61e48ef909375203092a6e83508c8416df55a83"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.2.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterface", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Reexport", "Requires", "SciMLBase", "SnoopPrecompile", "StaticArraysCore"]
git-tree-sha1 = "54c78ac3cc0343a16785adabe5bbf4063c737967"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.14"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.SimplexGridFactory]]
deps = ["DocStringExtensions", "ElasticArrays", "ExtendableGrids", "FileIO", "GridVisualize", "LinearAlgebra", "MeshIO", "Printf", "Test"]
git-tree-sha1 = "a02fa2a95e241c4145285475d55d640497b922ad"
uuid = "57bfcd06-606e-45d6-baf4-4ba06da0efd5"
version = "0.5.19"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "e19ac47477c9a8fcca06dab5e5471417d5d9d723"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.31.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "342cf4b449c299d8d1ceaf00b7a49f4fbc7940e7"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.9"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "08be5ee09a7632c32695d954a602df96a877bf0d"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.6"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "fd5f417fd7e103c121b0a0b4a6902f03991111f4"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.3.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b8d897fe7fa688e93aef573711cb207c08c9e11e"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.19"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "04a7d0bb1c824857ba0bb0c17bc5950dccbfdd5d"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.14.0"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "073da86200349ddf4ef8bc3e3f3acd62e1d554f7"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.60.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "e2d60a1cd52d0583471f83bd5d2dcefa626d271f"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.10"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SciMLBase", "SnoopPrecompile", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "a4e8491c163d74fb92905c6443e59558f5e609a4"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.16.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "5cb1f963f82e7b81305102dd69472fcd3e0e1483"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.5"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "7ecd651e3829d2957478516e92f693f12d5b4781"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.2.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "c97f60dd4f2331e1a495527f80d242501d2f9865"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.1"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f2fd3f288dfc6f507b0c3a2eb3bac009251e548b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.22"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.Triangle_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "fe28e9a4684f6f54e868b9136afb8fd11f1734a7"
uuid = "5639c1d2-226c-5e70-8d55-b3095415a16a"
version = "1.6.2+0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Triangulate]]
deps = ["DocStringExtensions", "Libdl", "Printf", "Test", "Triangle_jll"]
git-tree-sha1 = "bbca6ec35426334d615f58859ad40c96d3a4a1f9"
uuid = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
version = "2.2.0"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "7bc1632a4eafbe9bd94cf1a784a9a4eb5e040a91"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

[[deps.VTKBase]]
git-tree-sha1 = "c2d0db3ef09f1942d08ea455a9e252594be5f3b6"
uuid = "4004b06d-e244-455f-a6ce-a5f9919cc534"
version = "1.0.1"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "b182207d4af54ac64cbc71797765068fdeff475d"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.64"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.VoronoiFVM]]
deps = ["BandedMatrices", "CommonSolve", "DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "InteractiveUtils", "JLD2", "LinearAlgebra", "LinearSolve", "Printf", "Random", "RecursiveArrayTools", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrays", "Statistics", "SuiteSparse", "Symbolics", "Test"]
git-tree-sha1 = "33974fbb5978a37188e13f793a3e9740b9003f0d"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "1.2.1"

[[deps.VoronoiFVMDiffEq]]
deps = ["DifferentialEquations", "LinearAlgebra", "RecursiveArrayTools", "Reexport", "SciMLBase", "VoronoiFVM"]
git-tree-sha1 = "ab60171c6238859d1502212890675ebc2b55fd97"
uuid = "c4d4532c-19f6-4375-9ced-3808a145cd37"
version = "0.1.2"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams", "VTKBase"]
git-tree-sha1 = "7b46936613e41cfe1c6a5897d243ddcab8feabec"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.18.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═25c86226-be8b-11ed-0611-dd2967a6d7d8
# ╟─40dc1e16-ecf9-435b-8272-c94bc303bec5
# ╠═ceebef24-a96c-4253-955e-6ed7b6f9f070
# ╠═1edfff44-8085-4d25-99cd-0a65d28e441d
# ╠═8206d5f8-b953-4bc2-b9ce-758d62ab22cf
# ╠═f071e323-1e76-4774-933f-7178990e9331
# ╠═4f0b3a90-9c86-4602-84f7-49fd2407cd8e
# ╟─456af0e4-b93c-4741-a54a-b4454a6bffee
# ╟─b3d068b7-13dc-4bf6-b33f-87cd70efd50a
# ╟─2950bf88-76d3-4bce-bfe7-054eee34e9d6
# ╠═8ec02aea-816f-4500-942e-9f2a4b7c6712
# ╠═d5748251-d585-445e-b64a-eaf2ea795243
# ╠═d21a4876-144a-4fb9-bbb1-6d3b210012a1
# ╠═839219ec-fbf9-43c6-88fe-b2122fa0f9ce
# ╠═4028515a-3a7f-4381-831c-61d2a3ae8290
# ╠═96bcf4da-5562-4b58-a48e-e63a1d0f1e23
# ╠═92d2c851-86d4-411d-8908-d6c026dd3755
# ╠═77f61c64-39bb-46d0-afdb-615d51ccf653
# ╠═b88df152-2ef2-4547-be80-d320197b46ac
# ╠═ac3dc830-20b6-405e-9584-55f6bee0eecb
# ╠═ac48c632-bd37-4faf-a30f-54a5f0d3dbc5
# ╠═d30b26b7-6c3d-4d15-8263-0bf0373214d2
# ╟─763004cf-39e1-49f0-abe8-e4b207969f86
# ╠═567d36a0-d082-4fed-8385-5267962c6e71
# ╠═5fbe7493-9b59-471d-8807-518e95880526
# ╠═6200ebb3-74ed-418e-82d2-1d1539c85808
# ╟─3a0d7dcb-f85e-451b-95db-cc06599b276a
# ╟─e072661b-db1e-4aed-85ac-38432be7ae7d
# ╟─957b2e21-7a44-44a4-b42e-40573288e7ca
# ╟─10077a0a-6c1f-4e73-873a-3acb95d4c918
# ╠═6e4a8cb7-c231-4f29-bf2f-9cfa28cd1d1d
# ╟─d81771d6-3d0d-4c07-ba3e-d8a390e999bd
# ╠═34e337cc-75d8-4a8b-8f4e-1662b1045a9c
# ╠═713d7370-f821-4740-a4f0-6e08e50eff7a
# ╠═9b01f154-7cc4-455f-9636-2bc337458134
# ╠═17a230d4-b536-428e-a29b-7d53e0dae5f4
# ╟─0d6da6c6-63a9-4269-a401-08a41587c282
# ╟─c84d6e53-885c-48ba-9bc6-e85cc87bdefd
# ╟─fe737c36-b08b-465a-b6c7-e8656badc5a1
# ╠═64f8b862-fa71-4479-bbfe-4d300388a6ba
# ╠═50637189-a467-4f1b-bf22-606d06ea0981
# ╟─2b9f5ab8-2c99-43c0-9189-92f3de7f7da6
# ╠═3cc2d23e-45b8-4725-99a2-8a980b239ccc
# ╠═13bef003-732f-48af-b7b0-aba678189ab9
# ╠═4405f2d9-69f5-4331-aafa-27991ec52a77
# ╠═0c7b040f-6217-44be-8a3f-629b715e8a8c
# ╟─a48a0989-cd6e-4644-a6ab-e204dd7a403b
# ╠═8afe233d-1bce-435d-a689-aeddfb93578d
# ╠═19f07991-ff8a-499e-9d63-b73cde9fbf74
# ╟─f95ed529-fe7e-443f-9578-f1b8da779dad
# ╠═8cc09ef6-fd13-47b8-a6d3-f3a49c643913
# ╠═25d16219-af4a-434a-a116-236363494a69
# ╟─f3ab8364-8f10-4c02-a7d6-854f055d11c5
# ╠═14c04f9b-27b2-4638-a6ee-414cb6bd993e
# ╟─e9b07b6f-2b47-478a-a65e-3473e6554ced
# ╠═fa4a65d8-08d1-4704-ae33-77dfb4e80268
# ╟─9866ffec-e48f-45de-93c3-bba2d4dcbaf8
# ╠═3ec7c7ee-a353-44c8-8f8d-2064020542d1
# ╠═64160199-56a2-4a07-91c0-2a7f890f7f6b
# ╠═76761ecc-9ca4-4b8f-849c-4192d9d6cf3a
# ╠═fe0ed7cb-5671-49ae-9003-cfefdf0936d0
# ╠═7790a1f9-df1e-484e-a401-c8ae27eed721
# ╠═d82a2da3-144c-4780-8e59-dac86620c45c
# ╠═acb36288-ecd5-4ab6-8c85-e07b9f20bec9
# ╠═400a739f-6cc2-46ba-b2aa-9873bf109f15
# ╠═bfad73fe-eaa5-4b68-a8ed-e70a75753c03
# ╠═4e4a8844-5548-4b33-b00c-2d423979cd30
# ╠═c1ed66a0-7d86-4c2d-a8a7-03d22008bc9e
# ╠═b0692cf3-8eed-4b85-abab-ecf139c3b2c3
# ╠═c659d686-2917-4035-9b12-d1f0327fa752
# ╠═00f22fd5-d247-4c88-a9ae-c2ef8ee770bd
# ╠═1c5aef0a-9f5c-4471-9956-3ef457eef49b
# ╠═a3478c2b-28ba-4d3b-ba73-4772137b5359
# ╠═89ad2f8d-72bf-4dc2-9c5b-d3d19124ba2b
# ╠═12c02070-4c24-4038-8e39-8781466c4e0c
# ╠═642f9235-714a-404b-acdb-1f702412eba8
# ╠═0a2bbfb4-ce90-4442-9583-bfd76fdb2010
# ╠═bc7bfea9-68fd-4831-b262-82e6d38c632a
# ╠═fb622e06-d9b4-4f58-a099-721d8cc4bc9d
# ╠═734993cd-062b-4883-ab1b-6c85f8302fa8
# ╠═e6701190-0d1f-4c8f-99a2-5708bb32794a
# ╠═87019342-1089-4fd9-b088-6f57b56bec53
# ╠═b26344e1-4229-4cff-a4a0-bb5bf40e5500
# ╠═7bec601b-ac99-4446-a582-1415b83a754a
# ╠═a5a43c84-3671-4567-a040-f5c4c1ed2508
# ╠═0bfc2e2d-8009-4b03-b7f8-ace74d99013e
# ╠═c4c860a5-e01f-4d9e-bbc3-0700cda36d1f
# ╠═9f5339aa-c96d-4b16-8f4c-c94a4a2e4a53
# ╠═dbe648ea-f2b2-418b-a5b9-f2a9558473d1
# ╠═3d145c69-54d5-4175-abaf-c9f673ec2a88
# ╠═c8d29792-23b1-4afe-891d-bf3a27063656
# ╠═d4f79819-d23b-4683-9eae-691abae28697
# ╠═c2688388-c1af-433d-bc90-844cdb32cf08
# ╠═76c2df03-3e5a-4fb0-9ea7-dd2c04cc9131
# ╠═5adcd720-0b19-402e-b638-9857b5a0db27
# ╠═c96bd5c1-5f75-4412-a90c-4820454fcf48
# ╠═c0a5ea95-01db-4428-bd7a-89af1d75752a
# ╠═e0126cb0-3ae9-49d9-8753-235baf94b6a8
# ╠═023d053b-60e8-4e45-837b-18b754764f40
# ╠═6b3f828e-541d-4ee6-a5f9-acd4fcd2d790
# ╠═155ce2a2-28ad-4cdd-9230-cbf12f7cdc7d
# ╠═7ced64c4-5339-4e6f-998d-e8ea915695a4
# ╠═40f0251a-9ac2-4ce9-a180-e763b4380716
# ╠═3378bb7b-936d-44e6-8529-314bc33ed67f
# ╠═9d3d52e0-d64e-40d0-b149-3a4163b44070
# ╠═7896788e-4f37-4c99-96ce-a8006d3f6de6
# ╠═96ad3dbe-950d-4ecd-b919-211f0f36e549
# ╠═a5621f61-7853-4fc9-9f46-aae0f1dbe360
# ╠═ae9cce49-0303-4fc6-a5e0-eea616ab82de
# ╠═b60b8e7a-24f8-49f0-9d01-fa5f3c53cdaa
# ╠═fb401ab1-878a-41c8-b44c-806c409697da
# ╠═2800de62-c6a7-47e1-a977-d263f52624d5
# ╠═778994ac-d9ee-4482-80f1-ea76eec6c6d5
# ╟─41e2de15-dcbf-4340-9469-3ddc5184af49
# ╠═a0c6ec1a-bb4d-41e5-8689-04a9691796aa
# ╟─33cb6da8-d2c9-46a9-8667-19c06df0d7f7
# ╠═62968b3b-f4be-44c2-a2bb-a5c858bbe97b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
