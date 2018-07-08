"""

PolinomioTaylor(Cθ,CP, TipoVariable)
   
   Es una función cuyo objetivo es recibir dos listas con valores iniciales y crear dos polinomios de grado 1.

#.......................Método de parametrización usando factorización SVD ........................

-------------------------------------------------------------------------------------------------------------------------------------


Argumentos:


   - Cθ,CP  :  Arreglo que contiene los coeficientes iniciales, es del tipo Array{TaylorSeries.TaylorN{Float64}}
   - TipoVariable :  es el tipo de número que queremos uasr: Real, Float64, BigFloat, etc.
   

    Esta función regresa dos arreglos que contienen elementos del tipo Taylor1.TaylorN, creados con los valores propios del sistema.


"""
function PolinomioTaylor(Cθ,CP,TipoVariable)
    #=  
    Separamos el primer orden puesto que es diferente de los otros.
    Creamos θ,p como variables tipo TaylorN de orden 2, recordemos que aquí es importante el orden puesto que si se usan
    más adelante para otras operaciones el polinomio se puede cortar si no se define adecuadamente
    =#
    
    
    
    θ,p = set_variables(TipoVariable, "θ p", order=2)
    
    #especificamos que Lθ,LP son de arreglos que guardarán coeficientes del polinomio, sólo son auxiliares en esta función
    
   
   # Lθ = Array{TaylorSeries.TaylorN{TipoVariable}}(1)
   # LP = Array{TaylorSeries.TaylorN{TipoVariable}}(1)
        
    T = [Taylor1([θ], 1),Taylor1([p], 1)]
    return T
end
    

---------------------------------------------------------------------------------------------------------------------------------------------


"""

PolinomioTaylor1(g,Cθ,CP, TipoVariable)
   
   Es una función cuyo objetivo es recibir dos listas con valores iniciales y crear dos polinomios de grado g.




Argumentos:





   - g       : grado del polinomio
   - Cθ,CP  :  Arreglo que contiene los coeficientes iniciales, es del tipo Array{TaylorSeries.TaylorN{Float64}}
   - TipoVariable :  es el tipo de número que queremos uasr: Real, Float64, BigFloat, etc.
   

    Esta función regresa dos arreglos que contienen elementos del tipo Taylor1.TaylorN, creados con las listas iniciales. 
    Funciona específicamente para orden g>1.


"""
function PolinomioTaylor1(g::Real,Cθ,CP,TipoVariable)
    #=
    g es el grado del polinomio, CX y CP son
    arreglos que contienen los coeficientes que forman la variedad.
    
    
    Creamos x,p como variables tipo TaylorN de orden 2
    =#
    
    
    
    θ,p = set_variables(TipoVariable, "θ p", order=2)
    
    #especificamos que LX,LP son de arreglos que guardarán coeficientes del polinomio, sólo son auxiliares en esta función
    
   
    Lθ = Array{TaylorSeries.TaylorN{TipoVariable}}(1)
    LP = Array{TaylorSeries.TaylorN{TipoVariable}}(1)
    
    
    
    
        #como CX,CP están guardando los coeficientes pero necesitamos agregar el último término que será una variable 
        
    Lθ = push!(Cθ, θ) 
    LP = push!(CP, p)

    T=[Taylor1(Lθ, g),Taylor1(LP, g)]
        
    
    return T
end
#Esta función regresa tθ,tp que son Taylor1.TaylorN
-------------------------------------------------------------------------------------------------------------------------------------------------


#=Esta función toma el arreglo que contiene las lambdas que se van calculando, los coeficientes de los polinomios
y el orden de los mismos, lo que hace es generar el lado derecho de la ecuación cohomológica, multiplicando a_n*λ^n
y generando un polinomio de gradno g con estos coeficientes
=#
"""
Vecλ(λ_v,g,Cθ,CP)
Es una función que calcula la parte derecha de la ecuación comohológica, es decir la parte que involucra el valor propio.
Regresa un arreglo de tipo TaylorSeries.TaylorN{Float64}


Sus argumentos son:
-g      :  grado del polinomio.
-λ_v    :  Arreglo de dos dimensiones que contiene el valor propio y sus potencias. 
-Cθ,CP  :  Los arreglos con los polinomios que se calculan en PolinomioTaylor1.

"""
function Vecλ(λ_v,g,Cθ,CP)
   # el arreglo de λ_v contiene los arreglos que corresponden a la parte derecha de la ecuación cohomológica
    # en θ,p. Es importante hacer la distinción puesto que dependiendo del punto fijo donde se esté calculando
    # el primer valor de λ en θ serpa diferente del primer valor de λ en P
    θλt=Taylor1(λ_v[1].*Cθ,g)
    pλt=Taylor1(λ_v[2].*CP,g)
    
    λvec=[θλt,pλt]
    
    return λvec
end
---------------------------------------------------------------------------------------------------------------------------------------------

function EigenValores(M)
    M_facto = svdfact(M)
    @show(M_facto)
    Vals = real(LinearAlgebra.EigenGeneral.eigvals!(M))
    @show(Vals)
    EigVec=[[big(0.),big(0.)],[big(0.), big(0.)]]
    M_aux = M_facto.Vt
    for i in 1:2
        if Vals[i]== M_facto.S[i]
            EigVec[i] = M_aux[i,:]
        else
            EigVec[i] = -1*M_aux[i,:]
        end
    end
    if Vals[1]>Vals[2]
        Vals = [Vals[1],Vals[2]]
        EigVectors = [EigVec[1][1] EigVec[1][2]; EigVec[2][1] EigVec[2][2]]
    else 
        Vals = [Vals[2],Vals[1]]
        EigVectors = [EigVec[1][2] EigVec[2][2]; EigVec[1][1] EigVec[2][1]]
    end
        
    return (Vals,EigVectors)
        
        
        
    
end
-------------------------------------------------------------------------------------------------------------------------------------------


function Orden1(Cθ,CP,TipoVariable,Mapeo,k,l,PuntoFijo,tipo_v,λarrayθ,λarrayP)
            #usamos la función PolinomioTaylor para crear el polinomio tipo Taylor1.TaylorN{T}
            t = PolinomioTaylor(Cθ,CP,TipoVariable)
            @show(typeof(t))
            #Aplicamos el mapeo a los polinomios que resultan de la función anterior.
            Or1 = Mapeo(t[1],t[2],k,l)
            @show(typeof(Or1))
            AuxOr1=[Or1[1][1],Or1[2][1]]
            @show(typeof(AuxOr1))
            #Calculamos el jacobiano del Orden 1 para obtener sus valores y vectores propios.
            JPO = jacobian(AuxOr1,[PuntoFijo[1],PuntoFijo[2]])
            @show(typeof(JPO))
            
            
            #Calculamos los valores y vectores propios
            if TipoVariable == BigFloat
                eigval,eigvec = EigenValores(JPO)
            else
                eigval,eigvec = eig(JPO)
            end
            #escogemos el tipo de variedad que queremos calcular. Como se ordenan de menor a mayor la inestable es la segunda
            λ = eigval[tipo_v]
            
            #Ponemos los coeficientes en una variable nueva cada uno y los agregamos a las listas CX,CP,λ
            Coefθ,CoefP = eigvec[:,tipo_v]
            @show(typeof(Coefθ))
            #@show(Coefθ,CoefP)
            push!(Cθ, Coefθ)
            push!(CP, CoefP)
            push!(λarrayθ, λ)
            push!(λarrayP, λ)
            λ_v=[λarrayθ,λarrayP]
            
            #@show(λarray)  
            @show(typeof(Cθ))
            @show(typeof(λarrayθ))
    return Cθ, CP,λarrayθ,λarrayP, λ_v
end
            
-------------------------------------------------------------------------------------------------------------------------------------------


#Creamos una función que reciba el orden del polinomio , el punto fijo, el parámetro k y 
#el tipo de varidad que queremos(estable=1, inestable=2)
"""
Variedades(Mapeo,orden, PuntoFijo,k,tipo_v,TipoVariable)
Es una función que calcula las variedades de cierto mapeo. Usa las funciones de PolinomioTaylor1 y Vecλ para calcular los
polinomios de cada lado de la ecuación cohomológica y les aplica el mapeo dado. 



Argumentos:


- Mapeo : Mapeo de dos dimensiones, debe recibir al menos dos parámetros que son los polinomios antes calculados.
- orden : se trata del orden del polinomio.
- PuntoFijo : ES el punto fijo donde queremos calcular la variedad.
- k     : Es la constante del mapeo.
- tipo_v : 1 si la variedad es estable, 2 si es inestable.
- TipoVariable :  Float64,BigFloat, Integer,etc.


"""
function Variedades(Mapeo,orden, PuntoFijo,k,l,tipo_v, TipoVariable)
   
    #definimos unas listas donde se guardarán los coeficientes  de todo el polinomio, tales deben ser
    # de tipo "Array{TaylorSeries.TaylorN{Int64},1}" dado que los términos que se van agregando 
    # en cada orden son de tipo TaylorN.
    
    a=TipoVariable(PuntoFijo[1])
    b=TipoVariable(PuntoFijo[2])
    Cθ = [a+TaylorN(0.)]
    CP = [b+TaylorN(0.)]
    
    
    #λarray es la lista que contiene a los coeficientes del polinomio de λ
    λarrayθ = [a^0]
    λarrayP = [b^0]
    
    #definimos un vector que contiene el punto en el que se evalúa el jacobiano que se calcula después
    #dado que sólo lo usamos para obtener los valores que resultaron en el mapeo evaluamos siempre en [1.,1.]
    
    
    
    
    Cθ,CP,λarayθ, λarrayP,λ_v = Orden1(Cθ,CP,TipoVariable,Mapeo,k,l,PuntoFijo,tipo_v,λarrayθ,λarrayP)
    


    for g in 2:orden
        
            #Creamos los polinomios con las listas correspondientes 
            t = PolinomioTaylor1(g,Cθ,CP,TipoVariable)
            
            # aplicamos el mapeo estándar y al resultado le llamamos OrG por Orden g.
            OrG = Mapeo(t[1],t[2],k,l)
            
            push!(λarrayθ,λarrayθ[2]^g)
            push!(λarrayP,λarrayP[2]^g)
            λ_v=[λarrayθ,λarrayP]
            
            #agregamos el término correspondiente a λ 
            λ_vec=Vecλ(λ_v,g,Cθ,CP)
            
            
            #@show(λvec)
            
            # ahora ya tengo las dos partes de la ecuación y debo igualarlas para resolver.
            EcuaCohomo=OrG-λ_vec
            
            
            # de esta ecuación necesitamos solo los de orden g, así que los extraemos manualmente 
            θ_g=EcuaCohomo[1].coeffs[g+1]
            p_g=EcuaCohomo[2].coeffs[g+1]
            vec_orden_g=[θ_g,p_g]
            
            
            #calculamos el término independiene en la ecuación
            θ_ind=EcuaCohomo[1].coeffs[g+1].coeffs[1].coeffs[1]
            p_ind=EcuaCohomo[2].coeffs[g+1].coeffs[1].coeffs[1]
            vec_ind=[-θ_ind,-p_ind]
            
            #calculamos el jacobiano
            JacOrdenG = jacobian(vec_orden_g)
            
            
            
            
            #Con esta información podemos evaluar lo siguiente:
            # Si el vector de términos independientes es cero y el determinante del jacobiano es distinto de cero
            # entonces la solución a la ecuación cohomológica es la trivial
            if norm(vec_ind)==0.
                if det(JacOrdenG)!=0.
                    
                    Cθ[g+1]=0.
                    CP[g+1]=0.
                end
            else
                # Si el vector de términos independientes es distinto de ceroentonces necesitamos 
                #resolver la ecuación JacOrdenG[x_g,p_g]*[x,p]**=vec_ind[x_g,p_g]
                # entonces solo se trata de invertir el jacobiano y multiplicar con el vector del lado izquierdo
                TermG=JacOrdenG \ vec_ind
                
                Cθ[g+1]=TermG[1]
                CP[g+1]=TermG[2]
            
            end
            

    end
    return Cθ,CP,λarrayθ, λarrayP
end
-----------------------------------------------------------------------------------------------------------------------------------------


"""
PolinomioCohomo(Mapeo,Pol_vec,λvec, k)
Esta función calcula la ecuación cohomológica con los polinomios que ya se calcularon. Regresa un arreglo de dos 
elementos que son los valores de x,θ del mapeo.


Argumentos:
-Mapeo : función o mapeo del cual calculamos las variedades.Debe ser una función que reciba tres parámetros
 que son dos de sus variables y la constante del mapeo. Como salida debe tener un arreglo de dos elementos. 
-Pol_vec : Es un arreglo de dos elementos que son los polinomios calculados con anterioridad. 
-k     : es el valor de la constante del mapeo 
-λvec : 

"""
function PolinomioCohomo(Mapeo,Pol_vec,λvec, k,l ,PuntoFijo,modulo)
    Map_vec=Mapeo(Pol_vec[1],Pol_vec[2],k,l)
    if modulo==2*pi
        Ec_Cohomo = mod(Map_vec-λvec,modulo)
    else
    Ec_Cohomo = Map_vec-λvec
    end
    return Ec_Cohomo
end
---------------------------------------------------------------------------------------------------------------------------------------


"""
EvaluarPol(Ec_2var,Tiempo,paso)

Es una función que toma un arreglo de dos dimensiones que contiene polinomios y los evalúa en el tiempo dado en los pasos deseados




Argumentos:

-Ec_2var : Arreglo de dos dimensiones que contiene polinomios en cada una de ellas. 
-Tiempo  : Valor hasta el cual se quiere evaluar cada polinomio
-paso    : es el paso que se considera en cada evaluación del polinomio. 

"""
function EvaluarPol(Ec_2var,Tiempo,paso,TipoVariable)
    
    
    
    
    
    Val=TipoVariable[]
    Tiem=TipoVariable[]
    
    
    for t = 0:paso:Tiempo
        x = evaluate(Ec_2var[1], t)
        y = evaluate(Ec_2var[2], t)

        
        norma = norm([x,y],Inf)
        push!(Val,norma)
        push!(Tiem,t)
    
    end
    return Tiem,Val
end
-----------------------------------------------------------------------------------------------------------------------------------

#....................................................Error......................................................................
"""
CreaPol es una función que dadas dos listas y un grado crea  un arreglo de dos entradas , en cada una de ellas se encuentra 
el polinomio de grado g con los coeficientes de las listas. 


Argumentos:


- A,B : arreglos que contienen lo que serán los coeficientes del polinomio.
- orden : grado del polinomio
"""
function CreaPol(A,B,orden)
    Taylor = [Taylor1(A,orden),Taylor1(B,orden)]
    return Taylor
end

------------------------------------------------------------------------------------------------------------------------------------------

function MetParametrización(Mapeo,modulo,orden,PuntoFijo,k,l,tipo_v,Tiempo,paso, TipoVariable)
    Coeficienteθ,CoeficienteP,λarrayθ,λarrayP = Variedades(Mapeo,orden,PuntoFijo,k,l,tipo_v,TipoVariable)
        
    
    θ = TipoVariable[]
    P = TipoVariable[]
    
    for i in 1:orden+1
            
        push!(θ,Coeficienteθ[i].coeffs[1].coeffs[1])
        push!(P,CoeficienteP[i].coeffs[1].coeffs[1])
        
    end
    
    Taylor=CreaPol(θ,P,orden)
    
    λ_vec=CreaPol(θ.*λarrayθ,P.*λarrayP,orden)
    
    
    
    Ecua_Cohomo = PolinomioCohomo(Mapeo,Taylor,λ_vec, k,l,PuntoFijo,modulo)
    Valor_t , Error = EvaluarPol(Ecua_Cohomo,Tiempo,paso, TipoVariable)
    ErrorV = [Valor_t,Error]
    
    
    return Taylor,ErrorV,λ_vec
   
end
------------------------------------------------------------------------------------------------------------------------------------------




