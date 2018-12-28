

#..........................................................................................................................................
"""

PolinomioTaylor1(g,CX,CY)

   Es una función cuyo objetivo es recibir dos listas con valores iniciales y crear dos polinomios de grado g.


Argumentos:





   - g       : Grado del polinomio
   - CX,CY  :  Arreglo que contiene los coeficientes iniciales, es del tipo Array{TaylorSeries.TaylorN{T}}




"""
function PolinomioTaylor1(g::Int64,CX::Array{TaylorN{T}},CY::Array{TaylorN{T}}) where T<:Real
    #=
    Creamos x,y como variables tipo TaylorN de orden 2
    =#
    x,y = set_variables(T, "x y", order = 2)
    #especificamos que LX,LY son de arreglos que guardarán coeficientes del polinomio,
    # sólo son auxiliares en esta función
    LX = Array{TaylorN{T},1}
    LY = Array{TaylorN{T},1}

    # usamos un condicional para separar el caso 1 del resto, en el caso uno sólo se creal el polinomio
    # mientras que para g>1 hay que actualizar las lista y luego crear el polinomio
    if g == 1
        Taylor = [Taylor1([x], g),Taylor1([y], g)]
    #en el caso en que g>1 entonces usamos las listas que van guardando los coeficientes
    else
    #como CX,CY están guardando los coeficientes necesitamos primero agregar el último término que será una variable
        LX = push!(CX, x)
        LY = push!(CY, y)
    # luego crear el polinomio
        Taylor = [Taylor1(LX, g), Taylor1(LY, g)]

    end
    return Taylor
end


"""
InvarianciaR(λ_v,g,CX,CY)
Es una función que calcula la parte derecha de la ecuación comohológica, es decir la parte que involucra
el valor propio.
Regresa un arreglo de tipo TaylorSeries.TaylorN{T}


Argumentos:



- g      :  grado del polinomio.
- λ_v    :  Arreglo de dos dimensiones que contiene el valor propio y sus potencias.
- CX,CY  :  Los arreglos con los polinomios que se calculan en PolinomioTaylor1.

"""

function InvarianciaR(λ_v, g, CX, CY,λ,Polinomio)
    # El arreglo λ_v contiene las potencias del valor propio en orden ascendente .
    # Es importante hacer la distinción de una lista para X y otra para Y puesto que dependiendo del
    # punto fijo donde se esté calculando el primer valor de λ en X será diferente del primer valor de λ en Y
    xλt = Taylor1(λ_v[1].*CX, g)
    yλt = Taylor1(λ_v[2].*CY, g)


    variable_t = Taylor1([0.,1.])

    #= por alguna razón esto no le gusta ...
    xλt1 = Polinomio[1](λ*variable_t)
    yλt1 = Polinomio[2](λ*variable_t)
    resta  = yλt1-yλt
    @show(resta)
    =#
    InvarR = [xλt,yλt]

    return InvarR
end


#=
_______________________________________________________________________________________________________________________________________________

    El siguiente grupo de funciones es usado para calcular los valores y
    vectores propios en el caso que las variables sean del tipo BigFloat
=#

"""
EigenValores(M)


Función que calcula los valores y vectores propios para una matriz de elementos de tipo BigFloat.
El método usado es Inverse Power Method


Argumentos:


-M: Matriz tipo BigFloat (2-element Array{BigFloat,1})

"""
function EigenValores(M::Array{T,2}) where T<:BigFloat
    # MM sirve como matriz auxiliar para el cálculo de eigenvalores con la función eig()
    MM = [Float64(M[1]) Float64(M[3]); Float64(M[2]) Float64(M[4])]
    vfloat = [0. 0. ;0. 0.]
    h = 0.
    λfloat_aux,vfloat_aux = eigen(MM)
    imag(λfloat_aux[1]) == 0. ? h=1. :   error("Error: el valor propio 1 es complejo, es decir es un punto elíptico")
    imag(λfloat_aux[2]) == 0. ?  h=1. :  error("Error: el valor propio 2 es complejo, es decir es un punto elíptico")
    if λfloat_aux[1]>λfloat_aux[2]
        λfloat = [λfloat_aux[2],λfloat_aux[1]]
        vfloat[:,1] = vfloat_aux[:,2]
        vfloat[:,2] = vfloat_aux[:,1]
    else
        λfloat,vfloat = eigen(MM)
    end

    #Por otra parte calculamos los valores propios de M usando :
    Vals = real(GenericLinearAlgebra._eigvals!(M))
    #Creamos matrices que nos ayudarán en el método
    Id = [big(1.) big(0.); big(0.) big(1.)] #matriz identidad
    EigVec = [big(0.) big(0.);big(0.) big(0.)] #Matriz que guardará los vectores propios
    ValsOrden = [big(0.0),big(0.)] #guardará los valores propios ordenados



    #k recorre los dos valores propios
    for k in [1,2]
        #la tolerancia sirve para saber cuando paramos el método
        tol = 0.0

        vec_i = big.(rand(2,1))         # el vector inicial será un vector random

        # comienza el método


        while abs(tol)<1e-180

            H = inv(M - Vals[k]*Id)*vec_i                          #calculamos H

            vec_n = H/norm(H)                                      #normalizamos
            vec_i = vec_n                                          #actualizamos
            tol = norm(MM*vec_i-Vals[k]*vec_i, Inf)  #calculamos la diferencia
        end

        #los ordenamos de acuerdo a el orden de la función eig
        if abs(λfloat[1]-Vals[k])<1e-10

            ValsOrden[1] = sign(λfloat[1])*abs(Vals[k])
            EigVec[1] = sign(vfloat[1])*abs(vec_i[1])
            EigVec[2] = sign(vfloat[2])*abs(vec_i[2])

        else

            ValsOrden[2] = sign(λfloat[2])*abs(Vals[k])
            EigVec[3] = sign(vfloat[3])*abs(vec_i[1])
            EigVec[4]= sign(vfloat[4])*abs(vec_i[2])

        end


    end

    return ValsOrden,EigVec

end
#_________________________________________________________________________________________________________________________________________________

function EigenValores(M::Array{T,2}) where T<:Float64
    EigVec=[0. 0.;0. 0.]
    ValsOrden = T[0.,0.]
    #@show(M)

    ValsOrden_aux,EigVec_aux = eigen(M)
    #@show(ValsOrden_aux)
    #@show(EigVec_aux)
    if ValsOrden_aux[1] > ValsOrden_aux[2]

        ValsOrden = [ValsOrden_aux[2],ValsOrden_aux[1]]

        EigVec[:,1]=EigVec_aux[:,2]
        EigVec[:,2]=EigVec_aux[:,1]
    #    @show(EigVec_aux[:,2])

    else


        ValsOrden,EigVec = eigen(M)


    end


	return ValsOrden,EigVec

end


#=
_________________________________________________________________________________________________________________________________________________




 							Comienza el método de parametrización
La función Orden uno hace la inicialización del método.
=#
function Orden1(CX::Array{TaylorN{T}}, CY::Array{TaylorN{T}}, Mapeo, PuntoFijo::Array{T,1}, tipo_v, λarrayX::Array{T}, λarrayY::Array{T}) where T<:Real

            t = PolinomioTaylor1(1, CX, CY)

            #Aplicamos el mapeo a los polinomios que resultan de la función anterior.
            Or1 = Mapeo(PuntoFijo[1]+t[1], PuntoFijo[2]+t[2])
            AuxOr1 = [Or1[1][0], Or1[2][0]]

            #Calculamos el jacobiano del Orden 1 para obtener sus valores y vectores propios.
            JPO = jacobian(AuxOr1)
            JPO = [JPO[1] JPO[3];JPO[2] JPO[4]]

            #Calculamos los valores y vectores propios
            eigval,eigvec = EigenValores(JPO)

            #escogemos el tipo de variedad que queremos calcular. Como se ordenan de menor a mayor la inestable es la
            # segunda
            λ = eigval[tipo_v]

            # preguntamos si hay parte imaginaria del valor propio, en tal caso muestra un error
            tt = imag(λ)
            tt == 0. ?  Coef = eigvec[:,tipo_v] : error("Error: el valor propio es complejo, es decir es un punto elíptico")

            #Ponemos los coeficientes en una variable nueva cada uno y los agregamos a las listas CX,CP,λ
            push!(CX, Coef[1])
            push!(CY, Coef[2])
            push!(λarrayX, λ)
            push!(λarrayY, λ)
            λ_v = [λarrayX,λarrayY]



    return CX, CY, λarrayX, λarrayY, λ_v,λ, eigval,eigvec
end

#=
............................................................................................................................................
     							El resto del método se hace con la función Variedades

Creamos una función que reciba el orden del polinomio , el punto fijo, el parámetro k y
el tipo de varidad que queremos(estable=1, inestable=2)
=#

"""
Variedades(Mapeo,orden, PuntoFijo,tipo_v)
Es una función que calcula las variedades del Mapeo. Usa las funciones de PolinomioTaylor1 e InvarianciaR para
calcular los polinomios de cada lado de la ecuación cohomológica.



Argumentos:



- Mapeo : Mapeo de dos dimensiones, debe recibir al menos dos parámetros que son los polinomios antes calculados.
- orden : se trata del orden del polinomio.
- PuntoFijo : ES el punto fijo donde queremos calcular la variedad.
- tipo_v : 1 si la variedad es inestable, 2 si es estable.



"""
function Variedades(Mapeo, orden::Int64, PuntoFijo::Array{T,1},tipo_v) where T<:Real

    # definimos unas listas donde se guardarán los coeficientes  de todo el polinomio, tales deben ser
    # de tipo "Array{TaylorSeries.TaylorN{T},1}" dado que los términos que se van agregando
    # en cada orden son de tipo TaylorN.


    a=T(PuntoFijo[1])
    b=T(PuntoFijo[2])

    CX = [a+TaylorN(HomogeneousPolynomial([0.,0.]))]
    CY = [b+TaylorN(HomogeneousPolynomial([0.,0.]))]


    #λarray es la lista que contiene a los coeficientes del polinomio de λ
    λarrayX = [a^0]
    λarrayY = [b^0]


    #Calculamos el orden uno
    CX, CY, λarayX, λarrayY, λ_v,λ,eigval,eigvec = Orden1(CX, CY, Mapeo, PuntoFijo,tipo_v, λarrayX, λarrayY)

    for g in 2:orden

            #Creamos los polinomios con las listas correspondientes
            t = PolinomioTaylor1(g, CX, CY)

            # aplicamos el mapeo estándar y al resultado le llamamos OrG por Orden g.
            OrG = Mapeo(t[1], t[2])

            push!(λarrayX, λarrayX[2]^g)
            push!(λarrayY, λarrayY[2]^g)
            λ_v = [λarrayX, λarrayY]

            #agregamos el término correspondiente a λ
            λ_vec = InvarianciaR(λ_v, g, CX, CY,λ,t)


            # ahora ya tenemos las dos partes de la ecuación hay que igualarlas para resolver.
            EcuaCohomo = OrG-λ_vec


            # de esta ecuación necesitamos solo los de orden g, así que los extraemos manualmente
            X_g = EcuaCohomo[1].coeffs[g+1]
            Y_g = EcuaCohomo[2].coeffs[g+1]
            vec_orden_g = [X_g, Y_g]



            #calculamos el término independiene en la ecuación
            X_ind = EcuaCohomo[1].coeffs[g+1].coeffs[1].coeffs[1]
            Y_ind = EcuaCohomo[2].coeffs[g+1].coeffs[1].coeffs[1]
            vec_ind = [-X_ind, -Y_ind]


            #calculamos el jacobiano
            JacOrdenG = jacobian(vec_orden_g)




            #Con esta información podemos evaluar lo siguiente:
            # Si el vector de términos independientes es cero y el determinante del jacobiano es distinto de cero
            # entonces la solución a la ecuación cohomológica es la trivial
            if norm(vec_ind) == 0.
                if det(JacOrdenG) != 0.

                    CX[g+1] = 0.
                    CY[g+1] = 0.
                end
            else
                # Si el vector de términos independientes es distinto de ceroentonces necesitamos
                #resolver la ecuación JacOrdenG[x_g,p_g]*[x,p]**=vec_ind[x_g,p_g]
                # entonces solo se trata de invertir el jacobiano y multiplicar con el vector del lado izquierdo
                TermG = JacOrdenG \ vec_ind


                CX[g+1] = TermG[1]
                CY[g+1] = TermG[2]

            end


    end
    return CX, CY, λarrayX, λarrayY, eigval,eigvec
end
#=
......................................................................................................................................
					El método sigue con el calculo del error
=#

"""

Esta función calcula la ecuación cohomológica con los polinomios que ya se calcularon. Regresa un arreglo de dos
elementos que son los valores de x,y del mapeo.


Argumentos:




- Mapeo : función o mapeo del cual calculamos las variedades.Debe ser una función que reciba cuatro parámetros
 que son dos de sus variables y dos constantes del mapeo. Como salida debe tener un arreglo de dos elementos.
- Pol_vec : Es un arreglo de dos elementos que son los polinomios calculados con anterioridad.

- λvec : Arreglo de las potencias del valor propio
-modulo : si el mapeo tiene algún modulo, ej. 2pi

"""
function Invariancia(Mapeo, Pol_vec::Array{Taylor1{T},1}, λvec::Array{Taylor1{T},1} ,modulo::T) where T<:Real
    #@show(Pol_vec)
    Map_vec = Mapeo(Pol_vec[1],Pol_vec[2])
    #@show(Map_vec)
    Ec_Invariancia = mod(Map_vec-λvec, modulo)

    return Ec_Invariancia
end

function Invariancia(Mapeo, Pol_vec::Array{Taylor1{T},1}, λvec::Array{Taylor1{T},1} ) where T<:Real
    #@show(Pol_vec)
    #@show(λvec)
    Map_vec = Mapeo(Pol_vec[1],Pol_vec[2])
    #@show(Map_vec)

    Ec_Invariancia = Map_vec-λvec
    return Ec_Invariancia
end
#...................................................................................................................................

"""
EvaluarPol(Ec_2var,Tiempo,paso)

Es una función que toma un arreglo de dos dimensiones que contiene polinomios y los evalúa en el tiempo dado en los pasos deseados




Argumentos:




- Ec_2var : Arreglo de dos dimensiones que contiene polinomios en cada una de ellas.
- Tiempo  : Arreglo que dice el intervalo de tiempos a calcular
- paso    : es el paso que se considera en cada evaluación del polinomio.

"""
function EvaluarPol(Ec_2var,Tiempo::Array{T}, paso::T) where T<:Real

    Val = T[]
    Tiem = T[]


    for t = Tiempo[1]:paso:Tiempo[2]
        x = evaluate(Ec_2var[1], t)
        y = evaluate(Ec_2var[2], t)


        norma = norm([x,y], Inf)


        push!(Val, norma)
        push!(Tiem, t)

    end


    return Tiem, Val
end




#=
.....................................................................................................................................
 					El método termina con la creación de los polinomios
=#

"""
CreaPol es una función que dadas dos listas y un grado crea  un arreglo de dos entradas , en cada una de ellas se encuentra
el polinomio de grado g con los coeficientes de las listas.


Argumentos:



- A,B : arreglos que contienen lo que serán los coeficientes del polinomio.
- orden : grado del polinomio
"""
function CreaPol(A::Array{T}, B::Array{T}, orden::Int64) where T<:Real

    Taylor = [Taylor1(A, orden), Taylor1(B, orden)]

    return Taylor
end

#=
....................................................................................................................................

_______________________________Convergencia_________________________________________________________
para evaluar la convergencia de los polinomios se hizo la siguiente función que
 evalúa el cociente entre el coeficiente n+1 y el n.
=#


function Convergencia(A::Taylor1{T}, B::Taylor1{T}) where T<:Real
    A_aux = []
    B_aux = []
    Con_x = []
    Con_y = []
    suma_A = 0.
    suma_B = 0.

    for i in 1:length(A)
        if abs(A.coeffs[i]) >= 1e-180
            push!(A_aux,A.coeffs[i])
        end
    end
    for i in 1:length(B)
        if abs(B.coeffs[i]) >= 1e-180
            push!(B_aux,B.coeffs[i])
        end
    end

    for i in 1:(length(A_aux)-1)
        push!(Con_x,A_aux[i+1]/A_aux[i])
    end
    for i in 1:(length(B_aux)-1)
        push!(Con_y,B_aux[i+1]/B_aux[i])
    end




    return Con_x,Con_y
end
#___________________________________________ Test de tres términos_____________________________________
#Para analizar la convergencia de las soluciones tenemos el test de tres términos

function Convergencia3(A,B)
    A_aux = []
    B_aux = []
    Con_x = []
    Con_y = []
    suma_A = 0.
    suma_B = 0.

    for i in 1:length(A)
        if abs(A.coeffs[i]) > 1e-180
            push!(A_aux,A.coeffs[i])
        end
    end
    for i in 1:length(B)
        if abs(B.coeffs[i]) > 1e-180
            push!(B_aux,B.coeffs[i])
        end
    end
    for i in 2:(length(A_aux)-1)
        L = i*(A_aux[i+1]/A_aux[i])-(i-1)*(A_aux[i]/A_aux[i-1])
        L1 = abs(L)
        push!(Con_x,L1)
    end
    for i in 2:(length(B_aux)-1)
        L = i*(B_aux[i+1]/B_aux[i])-(i-1)*(B_aux[i]/B_aux[i-1])
        L1 = abs(L)
        push!(Con_y,L1)
    end
    return Con_x,Con_y
end






function CalculoError(Mapeo,modulo::Real,Taylor::Array{Taylor1{T}}, λ_vec::Array{Taylor1{T}},Tiempo::Array{T},paso::Real) where T<:Real


    Ecua_Invariancia = Invariancia(Mapeo, Taylor, λ_vec, modulo)

    Valor_t , Error = EvaluarPol(Ecua_Invariancia,Tiempo,paso)
    ErrorV = [Valor_t, Error]
    return ErrorV
end

function CalculoError(Mapeo,Taylor::Array{Taylor1{T}}, λ_vec::Array{Taylor1{T}},Tiempo::Array{T},paso::Real) where T<:Real


    Ecua_Invariancia = Invariancia(Mapeo, Taylor, λ_vec)

    Valor_t , Error = EvaluarPol(Ecua_Invariancia,Tiempo,paso)
    ErrorV = [Valor_t, Error]
    return ErrorV
end

function Error2(Mapeo,Taylor::Array{Taylor1{T}}, λ_vec::Array{Taylor1{T}},Tiempo::Array{T},paso::Real) where T<:Real

    Ecua_Invariancia = Invariancia(Mapeo, Taylor, λ_vec)

    norma_vec = norm.([Ecua_Invariancia[1].coeffs,Ecua_Invariancia[2].coeffs],Inf)

    Val = norm(norma_vec,Inf)

    return [Tiempo, Val]

end

function Error3(PuntoFijo::Array{T,1},Taylor::Array{Taylor1{T}},eigval, eigvec,Tiempo::Array{T}, paso::Real,orden::Real) where T<:Real


    var_aux = Taylor1([0.,1.])
    U = eigvec
    cero = Taylor1([0.])

    prim = 1-eigval[1]
    cuar = 1+1.5-eigval[1]
    matriz_a = [prim*var_aux 1.5*var_aux; 1.0*var_aux cuar*var_aux ]


    U_inv =inv(U)

    Mat = U*matriz_a*U_inv

    h=transpose(PuntoFijo)*matriz_a


    Invariancia_O =[abs(Taylor[1](0.))^orden,abs(Taylor[2](0.))^orden]

    norma = norm([Invariancia_O,Inf])

    return norma
end

function Error4(Mapeo,Taylor::Array{Taylor1{T}}, λ_vec::Array{Taylor1{T}},Tiempo::Array{T},paso::Real) where T<:Real
    DW_X = derivative(Taylor[1])
    DW_Y = derivative(Taylor[2])
    Ecua_derecha = [DW_X*λ_vec[1],DW_Y*λ_vec[2]]

    Ecua_Invariancia = Mapeo(Taylor[1],Taylor[2])-Ecua_derecha

    norma_vec = norm.([Ecua_Invariancia[1].coeffs,Ecua_Invariancia[2].coeffs],Inf)

    Val = norm(norma_vec,Inf)

    return [Tiempo, Val]
end






function Evaluar(Tiempo, paso,A,B)
    ValX=Float64[]
    ValY=Float64[]

    for t = Tiempo[1]:paso:Tiempo[2]

        x = evaluate(A,t)
        y = evaluate(B,t)

        push!(ValX,x)
        push!(ValY,y)

    end
    return ValX,ValY

end






#				 Las funciones que  engloban todo lo anterior son la siguientes : Inestable, Estable





function Inestable(Mapeo , orden::Int64, PuntoFijo::Array{T}, Tiempo::Array{T}, paso::T) where T<:Real
    CoeficienteX, CoeficienteY, λarrayX, λarrayY,eigval,eigvec = Variedades(Mapeo, orden, PuntoFijo,1)

    X = T[]
    Y = T[]

    for i in 1:orden+1

        push!(X,CoeficienteX[i].coeffs[1].coeffs[1])
        push!(Y,CoeficienteY[i].coeffs[1].coeffs[1])

    end

    Taylor = CreaPol(X, Y, orden+100)

    Xn = [X[i] for i in 1:orden+1]
    Yn = [Y[i] for i in 1:orden+1]

    λ_vec = CreaPol(Xn.*λarrayX, Yn.*λarrayY, orden+100)



    return Taylor,λ_vec , eigval,eigvec

end





#..............................................................................................................................................
function Estable(Mapeo, orden::Int64, PuntoFijo::Array{T}, Tiempo::Array{T}, paso::T) where T<:Real
    CoeficienteX, CoeficienteY, λarrayX, λarrayY,eigval,eigvec = Variedades(Mapeo, orden, PuntoFijo,2)


    X = T[]
    Y = T[]
    for i in 1:orden+1

        push!(X,CoeficienteX[i].coeffs[1].coeffs[1])
        push!(Y,CoeficienteY[i].coeffs[1].coeffs[1])

    end

    Taylor = CreaPol(X, Y, orden+100)
    Xn = [X[i] for i in 1:orden+1]
    Yn = [Y[i] for i in 1:orden+1]

    λ_vec = CreaPol(Xn.*λarrayX, Yn.*λarrayY, orden+100)



    return Taylor, λ_vec,eigval,eigvec

end
