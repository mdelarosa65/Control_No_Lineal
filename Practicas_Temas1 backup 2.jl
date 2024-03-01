### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 69c7e4b0-d636-11ee-2408-9baffb96902a
begin
	using Pkg
	Pkg.add("DifferentialEquations")
	Pkg.add("Plots")
	Pkg.add("GR")
	
end

# ╔═╡ 59dd98b0-0009-40c6-baeb-3cf79d407071
function oscilador_lineal(z, p, t)
	dz=[0.0 ; 0.0]
    dz[1] = z[2]
    dz[2] = p[1]*z[1]+p[2]*z[2]
	return dz
end

# ╔═╡ 27cb7947-f899-43f5-804f-2cefaa7ee874
scatter(real(autovalores),imag(autovalores),
	grid=true,legend=false, framestyle = :origin,
	xaxis=[-1,1], yaxis=[-1,1] )

# ╔═╡ 47381f13-b2c1-49e6-a611-e6b50a380db7
function flechas!(figura,f,p,rangos;N=10)
	#modifica una figura añadiendo las flechas por eso tiene una exclamación !
	# N es un parámetro opcional, si no se lo paso el valor es 10
	# los otros son posicionales
	xmin,xmax,ymin,ymax=rangos #desestructurar un iterable (explicar)
	longitud=max(abs(xmax-xmin),abs(ymax-ymin))/N
	
    #hago una malla de x e ys y pongo las flechas en cada punto
	x1s=xmin:((xmax-xmin)/N):xmax
	x2s=ymin:((ymax-ymin)/N):ymax

    #barrido para las flechas
	X1=[];
	X2=[];
	U=[];
	V=[];
	for x1 in x1s #barrido para calcular las flechas
		for x2 in x2s
			push!(X1,x1); #Push añade un componente al final de un vector
			push!(X2,x2);
			F=f([x1,x2],p,0)
			X=0.1*F/(0.1+sqrt(F[1]^2+F[2]^2)) #Saturo el tamaño de la flecha
			push!(U,X[1]);
			push!(V,X[2]);
		end
	end
	#Pinta las flechas
	quiver!(figura,X1,X2,quiver=(U,V), arrow=arrow(0.1*longitud, 0.1*longitud))
	return figura  
end

# ╔═╡ 54f9d537-c158-42b7-8124-337668677940
#= Es un euler modificado que integra hacia delante y hacia atrás
    congtinúa integrando hasta que se sale o completa un número de pasos. 
	   -hacia delante se detectan bien los PE y ciclos límite estables
	   -Hacia atrás permite ver también los inestables
	Ajusta dt para que el paso tenga un tamaño que sea aproximadamente constante=#

function trayectoria(x0,f,p,rangos;N=10,pasosMax=100)
	xmin,xmax,ymin,ymax=rangos
	long=max(abs(xmax-xmin),abs(ymax-ymin))/N
	x,y=x0
	pasos=0
	
	X1=[x]
	Y1=[y]
	X2=[x] #parte de la trayectoria hacia atrás en el tiempo
	Y2=[y]
	
	while (xmin<x<xmax) & (ymin<y<ymax) & (pasos<pasosMax)
		derivada=f([x,y],p,0)
		dt=0.05*long/(0.1+sqrt(derivada[1]^2+derivada[2]^2))
		Xn=[x;y]+derivada*dt
		x,y=Xn
		pasos=pasos+1
		push!(X1,x)
		push!(Y1,y)
	end
	pasos=0
	x,y=x0
	while (xmin<x<xmax) & (ymin<y<ymax) & (pasos<pasosMax)
		derivada=f([x,y],p,0)
		dt=-0.05*long/(0.1+sqrt(derivada[1]^2+derivada[2]^2))
		Xn=[x;y]+derivada*dt
		x,y=Xn
		pasos=pasos+1
		push!(X2,x)
		push!(Y2,y)
	end
	X=append!(reverse(X2),X1) #junto ambas partes
	Y=append!(reverse(Y2),Y1)
	return (X,Y)
end

# ╔═╡ 249db41c-d93b-4207-b697-e3ef4dcf83a7
function trayectorias!(figura,f,p,rangos;N=10)
	#modifica una figura añadiendo las trayectorias por eso tiene una exclamación !
	#barrido para las soluciones
	xmin,xmax,ymin,ymax=rangos
	N=round(N/2);
	x1s=xmin:((xmax-xmin)/N):xmax
	x2s=ymin:((ymax-ymin)/N):ymax
	T=10;
	for x1 in x1s
		for x2 in x2s 
			X,Y=trayectoria([x1;x2],f,p,rangos,N=N,pasosMax=100)
			figura=plot!(figura,X,Y,legend=false)
		end
	end 
	return figura
end

# ╔═╡ 02a2de67-5aa5-4a38-b57e-0b7ada77fd58
function fases(f,p;rangos=[-1,1,-1,1],N=10)
	xmin,xmax,ymin,ymax=rangos
	x1s=xmin:((xmax-xmin)/N):xmax
	x2s=ymin:((ymax-ymin)/N):ymax
	
	p1=plot()
	
	p1=flechas!(p1,f,p,rangos,N=N) #argumento opcional
	p1=trayectorias!(p1,f,p,rangos,N=N)
	
	xlims!(p1,(xmin, xmax))
	ylims!(p1,(ymin, ymax))
	return p1 #si lo quisieseis usar fuera de pluto hay que usar display(figura)
end

# ╔═╡ fb48d1a6-4b0b-4c05-b7ef-da6084e25f75
begin
	using LinearAlgebra
	using GR
    
	omega = 1.0
gr()
	for alpha in -1.0:0.5:1.0
		
	    k1 = -omega^2
	    k2 = -2 * alpha * omega
	    A = [0 1; k1 k2]
	    autovalores, _ = eigen(A)  
	    println("Para alpha = $alpha, los autovalores son: $autovalores")
		println()  # Agrega otro salto de línea
		display(fases(oscilador_lineal, [k1, k2]))

	end


end

# ╔═╡ a4be0a88-8c0b-4f05-bc3a-3fec921e5816
using LinearAlgebra

# ╔═╡ Cell order:
# ╠═69c7e4b0-d636-11ee-2408-9baffb96902a
# ╠═a05ac043-3be0-4296-8649-ae3e0082c5a4
# ╠═3703cc87-98da-4f83-96fd-efe7b7623f88
# ╠═59dd98b0-0009-40c6-baeb-3cf79d407071
# ╠═a4be0a88-8c0b-4f05-bc3a-3fec921e5816
# ╠═27cb7947-f899-43f5-804f-2cefaa7ee874
# ╠═10e51a1d-d9f1-408f-8d46-49906b34adcb
# ╠═fb48d1a6-4b0b-4c05-b7ef-da6084e25f75
# ╟─47381f13-b2c1-49e6-a611-e6b50a380db7
# ╟─54f9d537-c158-42b7-8124-337668677940
# ╟─249db41c-d93b-4207-b697-e3ef4dcf83a7
# ╟─02a2de67-5aa5-4a38-b57e-0b7ada77fd58
