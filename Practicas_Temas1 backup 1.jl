### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 69c7e4b0-d636-11ee-2408-9baffb96902a
begin
	using Pkg
	Pkg.add("DifferentialEquations")
	Pkg.add("Plots")
	
end

# ╔═╡ a05ac043-3be0-4296-8649-ae3e0082c5a4
begin
	using DifferentialEquations
	using Plots
	
	# Definir el sistema de ecuaciones diferenciales
	function system!(du, u, p, t)
	    x, v = u
	    du[1] = v
	    du[2] = -0.6 * v - 3 * x - x^2
	end
	
	# Condiciones iniciales
	u0 = [1.0, 0.0]  # Cambia las condiciones iniciales si es necesario
	
	# Definir el rango de tiempo
	t_span = (0.0, 10.0)
	
	# Resolver el sistema de ecuaciones
	prob = ODEProblem(system!, u0, t_span)
	sol = solve(prob, Vern7(), dt=0.01)
	
	# Graficar el plano fase
	plot(sol, xlabel="x", ylabel="x'", legend=false)
	scatter!([0, -3], [0, 0], color=:red, label="Puntos Singulares", markersize=5)
	
end

# ╔═╡ 3703cc87-98da-4f83-96fd-efe7b7623f88
begin
	
	using LinearAlgebra
	
	# Definir la ecuación diferencial
	function system(x, t)
	    return [x[2], -0.6*x[2] - 3*x[1] - x[1]^2]
	end
	
	# Definir el espacio de fase
	function phase_diagram()
	    x_range = -5:0.1:5
	    v_range = -5:0.1:5
	
	    plot(size=(500, 500), xlabel="x", ylabel="x'", legend=false)
	
	    for x in x_range, v in v_range
	        u = [x, v]
	        result = system(u, 0)
	        quiver!([x], [v], quiver=([result[1]], [result[2]]), color=:blue)
	    end
	
	    scatter!([0, -3], [0, 0], color=:red, label="Puntos Singulares", markersize=5)
	end
	
	# Llamar a la función para generar el diagrama de fases
	phase_diagram()

		# Definir la matriz jacobiana
	A = [0 1; -3 -0.6]
	
	# Calcular los autovalores
	autovalores,_=eigen(A)
	
	md"""
	
	
	autovalores= $(autovalores)"""

end

# ╔═╡ Cell order:
# ╠═69c7e4b0-d636-11ee-2408-9baffb96902a
# ╠═a05ac043-3be0-4296-8649-ae3e0082c5a4
# ╠═3703cc87-98da-4f83-96fd-efe7b7623f88
# ╠═75948c33-474f-4e72-ab77-f729009cb657
