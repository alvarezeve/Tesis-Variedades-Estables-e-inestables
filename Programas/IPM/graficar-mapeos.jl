#este módulo sirve para graficar mapeos.


function IterarMap(f,x_ini,n)   #Definimos una función para iterar el mapeo

    x = x_ini[1]                      #
                                     #Damos condiciones iniciales
    y = x_ini[2]                          #

    iteradosMapX = [x]

    iteradosMapY = [y]               #Definimos dos listas que tendran los valores de cada par ordenado de theta y P, y agregamos las condiciones iniciales

    for i=0:n              #iniciamos un ciclo de iteraciones donde se calculan x_n, y_n y se agregan a lalista correspodiente

        F = f(x,y)

        push!(iteradosMX,F[1])

        push!(iteradosMY,F[2])


        x = F[1]

        p = F[2]
    end

    return iteradosMX, iteradosMY  #La funcion iterados regresa las listas que corresponden a la trayectoria del
end


function GraficarEstandarMap(k)
    n = 50
    s=2pi/14
    i=0
    for p_i=0:s:2pi
        for x_i =0:s:2pi
            i=i+1
            t,p = IterarMap(f,x_ini,n)
            p = scatter(t,p,marker=".",s=0.1,color="gray")
        end
    end

end

function Evaluar(Tiempo,paso,A,B,PuntoFijo,col)
    ValX = Float64[]
    ValY = Float64[]

    push!(ValX, PuntoFijo[1])
    push!(ValY, PuntoFijo[2])
    for t = Tiempo[1]:-paso:Tiempo[2]

        x = Float64(evaluate(A,t))
        y = Float64(evaluate(B,t))

        push!(ValX,x)
        push!(ValY,y)

    end

    p = plot(ValX,ValY,linestyle = "-",marker=",",color=col)


    title("Espacio Fase")
    xlabel("$\theta$")
    ylabel("$p$")
    xlim(0.,2*pi)
    ylim(0.,2*pi)
    legend(loc="upper right",fancybox="true")



end


function Graficar(Tiempo,paso,k,A,B,PuntoFijo,col)
    GraficarEstandarMap(k)
    Evaluar(Tiempo,paso,A,B,PuntoFijo,col)
end
