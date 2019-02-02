module graficar
#este módulo sirve para graficar mapeos.
export IterarMap,GraficarEst


function IterarMap(f, x_ini, n)   #Definimos una función para iterar el mapeo

    x = x_ini[1]                      #                                     #Damos condiciones iniciales
    y = x_ini[2]                          #

    iteradosMapX = [x]
    iteradosMapY = [y]               #Definimos dos listas que tendran los valores de cada par ordenado de theta y P, y agregamos las condiciones iniciales

    for i=0:n              #iniciamos un ciclo de iteraciones donde se calculan x_n, y_n y se agregan a lalista correspodiente
        F = f(x,y)

        push!(iteradosMapX,F[1])
        push!(iteradosMapY,F[2])

        x = F[1]
        y = F[2]

    end

    return iteradosMapX, iteradosMapY  #La funcion iterados regresa las listas que corresponden a la trayectoria del

end

function GraficarEstandarMap(f)
    n = 50
    s=2pi/14.
    i=0.
    for p_i=0:s:2pi
        for x_i =0:s:2pi
            i=i+1
            x_ini=[p_i,x_i]
            a,b = IterarMap(f,x_ini,n)
            p = scatter(a,b,marker=".",s=0.1,color="gray")
        end
    end

end

function Evaluar(Tiempo, paso, Pol, PuntoFijo, col)
    ValX = Float64[]
    ValY = Float64[]

    push!(ValX, PuntoFijo[1])
    push!(ValY, PuntoFijo[2])
    #@show(ValX)
    for t = Tiempo[1]:paso:Tiempo[2]

        x = Float64(Pol[1](t))
        y = Float64(Pol[2](t))

        push!(ValX,-x)
        push!(ValY,-y)

    end
    @show(ValX)

    p = plot(ValX,ValY,linestyle = "-",marker=",",color=col)


#=
    title("Espacio Fase")
    xlabel(L"$x$")
    ylabel(L"$y$")
    xlim(0.,2*pi)
    =#

end

function GraficarEst(Tiempo,paso,f,Pol,PuntoFijo,col)
    GraficarEstandarMap(f)
    Evaluar(Tiempo,paso,Pol,PuntoFijo,col)

end



end
