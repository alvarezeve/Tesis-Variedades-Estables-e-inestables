

#..........................................................................................................................................
"""

PolinomioTaylor1(g,CX,CY,TipoVariable)
   
   Es una función cuyo objetivo es recibir dos listas con valores iniciales y crear dos polinomios de grado g.




Argumentos:





   - g       : Grado del polinomio
   - CX,CY  :  Arreglo que contiene los coeficientes iniciales, es del tipo Array{TaylorSeries.TaylorN{}}
   - TipoVariable :  es el tipo del que serán los coeficientes ( Real, Float64, BigFloat)
   

Esta función regresa dos arreglos que contienen elementos del tipo Taylor1.TaylorN, creados con las listas iniciales. 


"""
function PolinomioTaylor1(g::Real,CX,CY,TipoVariable)
    #=
    g es el grado del polinomio, CX y CY son
    arreglos que contienen los coeficientes que forman la variedad.
    
    
    Creamos x,y como variables tipo TaylorN de orden 2
    =#
    x,y = set_variables(TipoVariable, "x y", order=2)
    
    #especificamos que LX,LY son de arreglos que guardarán coeficientes del polinomio, sólo son auxiliares en esta función
    LX = Array{TaylorSeries.TaylorN{TipoVariable}}(1)
    LY = Array{TaylorSeries.TaylorN{TipoVariable}}(1)
    
    #usamos un condicional para separar el caso 1 del resto
    if g == 1  
        T = [Taylor1([x], g),Taylor1([y], g)]
    #en el caso en que g>1 entonces usamos las listas que van guardando los coeficientes
    else
        #como CX,CY están guardando los coeficientes pero necesitamos agregar el último término que será una variable 
        LX = push!(CX, x) 
        LY = push!(CY, y)

        T=[Taylor1(LX, g),Taylor1(LY, g)]
        
    end
    return T
end
#Esta función regresa tx,ty que son Taylor1.TaylorN




#..................................................................................................................................

#Esta función toma el arreglo que contiene las lambdas que se van calculando, los coeficientes de los polinomios
#y el orden de los mismos, lo que hace es generar el lado derecho de la ecuación cohomológica, multiplicando a_n*λ^n
#y generando un polinomio de gradno g con estos coeficientes

"""
Vecλ(λ_v,g,CX,CY)
Es una función que calcula la parte derecha de la ecuación comohológica, es decir la parte que involucra el valor propio.
Regresa un arreglo de tipo TaylorSeries.TaylorN{T}


Sus argumentos son:



- g      :  grado del polinomio.
- λ_v    :  Arreglo de dos dimensiones que contiene el valor propio y sus potencias. 
- CX,CY  :  Los arreglos con los polinomios que se calculan en PolinomioTaylor1.

"""
function Vecλ(λ_v,g,CX,CY)
   # el arreglo de λ_v contiene los arreglos que corresponden a la parte derecha de la ecuación cohomológica
    # en θ,p. Es importante hacer la distinción puesto que dependiendo del punto fijo donde se esté calculando
    # el primer valor de λ en θ serpa diferente del primer valor de λ en P
    xλt=Taylor1(λ_v[1].*CX,g)
    yλt=Taylor1(λ_v[2].*CY,g)
    
    λvec=[xλt,yλt]
    
    return λvec
end
#_______________________________________________________________________________________________________________________________________________

#    El siguiente grupo de funciones es usado para calcular los valores y vectores propios en el caso que las variables sean del tipo BigFloat

"""
EigenValores(M)


Función que calcula los valores y vectores propios para una matriz de elementos de tipo BogFloat.
El método usado es el de Inverse Power Method 




Argumentos:
-M: Matriz tipo BigFloat (2-element Array{BigFloat,1})
"""
function EigenValores(M)
    # MM sirve como matriz auxiliar para el cálculo de eigenvalores con la función eig()
    
    MM = [Float64(M[1]) Float64(M[3]); Float64(M[2]) Float64(M[4])]
    λfloat,vfloat = eig(MM)
    #Por otra parte calculamos los valores propios de M usando :
    Vals = real(LinearAlgebra.EigenGeneral.eigvals!(M))
    #Creamos matrices que nos ayudarán en el método
    Id = [big(1.) big(0.); big(0.) big(1.)] #matriz identidad
    EigVec=[big(0.) big(0.);big(0.) big(0.)] #MAtriz que guardará los vectores propios
    ValsOrden=[big(0.0),big(0.)] #guardará los valores propios ordenados 
    
    #vec_i =([0.6,0.3],[1.,0.]) #debemos dar un vector inicial para el cálculo
    
    #k recorre los dos valores propios
    for k in [1,2]
        #la tolerancia sirve para saber cuando paramos el método
        tol = 0.0
        
        vec_i = big(rand(2,1)) # el vector inicial será un vector random
        # comienza el método
        #for i in[1:10]
            
        while abs(tol)<1e-180
            
            H = inv(M-Vals[k]*Id)*vec_i #calculamos H
            
            vec_n = H/norm(H) #normalizamos
            vec_i =vec_n #actualizamos
            #vec_i[k][1] = vec_n[1]
            #vec_i[k][2] = vec_n[2]
            tol = norm(MM*vec_i-Vals[k]*vec_i,Inf)  #calculamos el error o discrepancia
        end
        
    
        if abs(λfloat[1]-Vals[k])<1e-10

            ValsOrden[1]=sign(λfloat[1])*abs(Vals[k])
            EigVec[1]=sign(vfloat[1])*abs(vec_i[1])
            EigVec[2]=sign(vfloat[2])*abs(vec_i[2])
            
        else
            
            ValsOrden[2]=sign(λfloat[2])*abs(Vals[k])
            EigVec[3]=sign(vfloat[3])*abs(vec_i[1])
            EigVec[4]=sign(vfloat[4])*abs(vec_i[2])
    
        end
            
        
    end
    
    #V=(ValsOrden,EigVec)
    
    return ValsOrden,EigVec
        
        
        
    
end     
#_________________________________________________________________________________________________________________________________________________

# 							Comienza el método de parametrización
#La función Orden uno hace la inicialización del método.

function Orden1(CX,CY,TipoVariable,Mapeo,k,l,PuntoFijo,tipo_v,λarrayX,λarrayY)
            #usamos la función PolinomioTaylor para crear el polinomio tipo Taylor1.TaylorN{T}
            t = PolinomioTaylor1(1,CX,CY,TipoVariable)
            
            #Aplicamos el mapeo a los polinomios que resultan de la función anterior.
            Or1 = Mapeo(PuntoFijo[1]+t[1],PuntoFijo[2]+t[2],k,l)
            
            AuxOr1=[Or1[1][1],Or1[2][1]]
            
            #Calculamos el jacobiano del Orden 1 para obtener sus valores y vectores propios.
            JPO = jacobian(AuxOr1)
            
            
            
            #Calculamos los valores y vectores propios
            if TipoVariable == BigFloat
                eigval,eigvec = EigenValores(JPO)
            else
                eigval,eigvec = eig(JPO)
            end
            #escogemos el tipo de variedad que queremos calcular. Como se ordenan de menor a mayor la inestable es la segunda
            λ = eigval[tipo_v]
    
            tt = imag(λ)
            
            #Ponemos los coeficientes en una variable nueva cada uno y los agregamos a las listas CX,CP,λ
            tt == 0.?  Coef = eigvec[:,tipo_v] : error("Error: el valor propio es complejo, es decir es un punto elíptico")
            
            
            push!(CX, Coef[1])
            push!(CY, Coef[2])
            push!(λarrayX, λ)
            push!(λarrayY, λ)
            λ_v=[λarrayX,λarrayY]
            
            
    return CX, CY,λarrayX,λarrayY, λ_v
end
            
            
#............................................................................................................................................
# 							El resto del método se hace con la función Variedades    

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
    CX = [a+TaylorN(0.)]
    CY = [b+TaylorN(0.)]
    
    
    #λarray es la lista que contiene a los coeficientes del polinomio de λ
    λarrayX = [a^0]
    λarrayY = [b^0]
    
    #definimos un vector que contiene el punto en el que se evalúa el jacobiano que se calcula después
    #dado que sólo lo usamos para obtener los valores que resultaron en el mapeo evaluamos siempre en [1.,1.]
    
    
    
    
    CX,CY,λarayX, λarrayY,λ_v = Orden1(CX,CY,TipoVariable,Mapeo,k,l,PuntoFijo,tipo_v,λarrayX,λarrayY)
    


    for g in 2:orden
        
            #Creamos los polinomios con las listas correspondientes 
            t = PolinomioTaylor1(g,CX,CY,TipoVariable)
            
            # aplicamos el mapeo estándar y al resultado le llamamos OrG por Orden g.
            OrG = Mapeo(t[1],t[2],k,l)
            
            push!(λarrayX,λarrayX[2]^g)
            push!(λarrayY,λarrayY[2]^g)
            λ_v=[λarrayX,λarrayY]
            
            #agregamos el término correspondiente a λ 
            λ_vec=Vecλ(λ_v,g,CX,CY)
            
            
            
            # ahora ya tengo las dos partes de la ecuación y debo igualarlas para resolver.
            EcuaCohomo=OrG-λ_vec
            
            
            # de esta ecuación necesitamos solo los de orden g, así que los extraemos manualmente 
            X_g=EcuaCohomo[1].coeffs[g+1]
            Y_g=EcuaCohomo[2].coeffs[g+1]
            vec_orden_g=[X_g,Y_g]
            
            
            #calculamos el término independiene en la ecuación
            X_ind=EcuaCohomo[1].coeffs[g+1].coeffs[1].coeffs[1]
            Y_ind=EcuaCohomo[2].coeffs[g+1].coeffs[1].coeffs[1]
            vec_ind=[-X_ind,-Y_ind]
            
            #calculamos el jacobiano
            JacOrdenG = jacobian(vec_orden_g)
            
            
            
            
            #Con esta información podemos evaluar lo siguiente:
            # Si el vector de términos independientes es cero y el determinante del jacobiano es distinto de cero
            # entonces la solución a la ecuación cohomológica es la trivial
            if norm(vec_ind)==0.
                if det(JacOrdenG)!=0.
                    
                    CX[g+1]=0.
                    CY[g+1]=0.
                end
            else
                # Si el vector de términos independientes es distinto de ceroentonces necesitamos 
                #resolver la ecuación JacOrdenG[x_g,p_g]*[x,p]**=vec_ind[x_g,p_g]
                # entonces solo se trata de invertir el jacobiano y multiplicar con el vector del lado izquierdo
                TermG=JacOrdenG \ vec_ind
                
                CX[g+1]=TermG[1]
                CY[g+1]=TermG[2]
            
            end
            

    end
    return CX,CY,λarrayX, λarrayY
end

#......................................................................................................................................
#					El método sigue con el calculo del error

"""
PolinomioCohomo(Mapeo,Pol_vec,λvec, k)
Esta función calcula la ecuación cohomológica con los polinomios que ya se calcularon. Regresa un arreglo de dos 
elementos que son los valores de x,θ del mapeo.


Argumentos:




- Mapeo : función o mapeo del cual calculamos las variedades.Debe ser una función que reciba tres parámetros
 que son dos de sus variables y la constante del mapeo. Como salida debe tener un arreglo de dos elementos. 
- Pol_vec : Es un arreglo de dos elementos que son los polinomios calculados con anterioridad. 
- k     : es el valor de la constante del mapeo 
- λvec : 

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
#...................................................................................................................................

"""
EvaluarPol(Ec_2var,Tiempo,paso)

Es una función que toma un arreglo de dos dimensiones que contiene polinomios y los evalúa en el tiempo dado en los pasos deseados




Argumentos:




- Ec_2var : Arreglo de dos dimensiones que contiene polinomios en cada una de ellas. 
- Tiempo  : Valor hasta el cual se quiere evaluar cada polinomio
- paso    : es el paso que se considera en cada evaluación del polinomio. 

"""
function EvaluarPol(Ec_2var,Tiempo,paso,TipoVariable)
    
    
    
    
    
    Val=TipoVariable[]
    Tiem=TipoVariable[]
    
    
    for t = Tiempo[1]:paso:Tiempo[2]
        x = evaluate(Ec_2var[1], t)
        y = evaluate(Ec_2var[2], t)

        
        norma = norm([x,y],Inf)
        push!(Val,norma)
        push!(Tiem,t)
    
    end
    return Tiem,Val
end
#.....................................................................................................................................
# 					El método termina con la creación de los polinomios 


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
#....................................................................................................................................

#____________________________________________________________Convergencia_________________________________________________________
#para evaluar la convergencia de los polinomios se hizo la siguiente función que evalúa el cociente entre el coeficiente n+1 y el n


function Convergencia(A,B)
    A_aux = []
    B_aux = []
    Con_x = []
    Con_y = []
    suma_A = 0.
    suma_B = 0.
    
    for i in 1:length(A)-1
        if A[i] != 0.
            push!(A_aux,A[i])
        end
        if B[i] != 0.
            push!(B_aux,B[i])
        end
    end
    for i in 1:(length(A_aux)-2)
        push!(Con_x,A_aux[i+1]/A_aux[i])
        push!(Con_y,B_aux[i+1]/B_aux[i])
    end
        
   
    return Con_x,Con_y
end
        







#				 La función que engloba todo lo anterior es MetParametrización

function MetParametrización(Mapeo,modulo,orden,PuntoFijo,k,l,tipo_v,Tiempo,paso, TipoVariable)
    CoeficienteX,CoeficienteY,λarrayX,λarrayY = Variedades(Mapeo,orden,PuntoFijo,k,l,tipo_v,TipoVariable)
    Conver_X,Conver_Y=Convergencia(CoeficienteX,CoeficienteY)    
    
    X = TipoVariable[]
    Y = TipoVariable[]
    
    for i in 1:orden+1
            
        push!(X,CoeficienteX[i].coeffs[1].coeffs[1])
        push!(Y,CoeficienteY[i].coeffs[1].coeffs[1])
        
    end
    
    Taylor=CreaPol(X,Y,orden)
    
    λ_vec=CreaPol(X.*λarrayX,Y.*λarrayY,orden)
    
    
    
    Ecua_Cohomo = PolinomioCohomo(Mapeo,Taylor,λ_vec, k,l,PuntoFijo,modulo)
    Valor_t , Error = EvaluarPol(Ecua_Cohomo,Tiempo,paso, TipoVariable)
    ErrorV = [Valor_t,Error]
    
    
    return Taylor,ErrorV,λ_vec,Conver_X,Conver_Y
   
end
#..............................................................................................................................................






