"""
La función Estandar(x::Any, y::Any k::Real)

El mapeo estándar con módulo 2pi
"""
function Estandar(x,y,k)

    x_n = mod(x+k*sin(y),2.*pi)
    y_n = mod(x+y+k*sin(y),2.*pi)

    return [x_n,y_n]
end

"""
La función EstandarI(x::Any, y::Any k::Real)

El mapeo estándar inverso con módulo 2pi
"""

function EstandarI(x,y,k)

    x_n = mod(x-y+k*sin(x),2.*pi)
    y_n = mod(y-k*sin(x),2.*pi)

    return [x_n,y_n]
end
"""
La función Jung(x::Any, y::Any, a::Real )

El mapeo de Jung(exponencial)
"""

function Jung(x,y,a)

    x_n = x+y
    y_n = y+a*x_n*(x_n-1.)*e^(-x_n)

    return [x_n,y_n]
end


"""
La función Henon(x::Any,y::Any,a::Real,b::Real)

EL mapeo de Hénon con dos parámetros.
"""
function Henon(x,y,a,b)

    x_n = a-b*y-x^2
    y_n = x

    return [x_n,y_n]
end


"""
HenonI(x::Any, y::Any,a::Real,b::Real)
Función que implementa el mapeo de Hénon inverso
"""
function HenonI(x,y,a,b)

    x_n = y
    y_n = (x+y^2-a)/(-b)

    return [x_n,y_n]
end
