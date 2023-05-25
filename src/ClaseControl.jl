module ClaseControl

export bode, fdos

using Plots, LaTeXStrings, DSP, Roots, Interpolations
plotly()

struct Bode
    fig
    wco
    RAco
    w1
    phi1
end

"""
**bode(Gol; wmin=1e-1, wmax=1e1, points=100, co=false, ra1=false, RAlabel="RA")**

Representación del diagrama de Bode de la función de transferencia _Gol_.

Parámetros:

- `Gol`: Función de transferencia dependiente de _s_
- `wmin`: Frecuencia angular mínima a representar
- `wmax`: Frecuencia angular máxima a representar
- `points`: Número de puntos a representar
- `co`: Calcula y representa la frecuencia de cruce
- `ra1`: Calcula y representa la frecuencia que hace que _RA_ = 1
- `RAlabel`: Etiqueta del gráfico de RA en el diagrama de Bode

_Ejemplo:_

`salida = bode(s->20*(s+1)/s/(s+5)/(s^2+2s+10); wmax=10, co=true, ra1=true)`

Para los resultados:

- `salida.fig`: Diagrama de Bode
- `salida.wco` y `salida.RAco`: Frecuencia de cruce y razón de amplitudes
    para la frecuencia de cruce
- `salida.w1` y `salida.phi1`: Frecuencia para RA = 1 y desfase para RA=1
"""
function bode(Gol; wmin=1e-1, wmax=1e1, points=100, co=false, ra1=false, RAlabel="RA")
    wco = nothing
    RAco = nothing
    w1 = nothing
    phi1 = nothing

    RAol(w) = sqrt(real(Gol(im*w))^2+imag(Gol(im*w))^2)
    phiol(w) = atan(imag(Gol(im*w)), real(Gol(im*w)))
    
    # Hace una lista equiespaciada para los valores del log10 del desfase
    wlog = range(log10(wmin), log10(wmax); length=points)
    # Genera la versión lineal
    wlin = 10 .^wlog
    
    # Calcula los datos de desfase para los valores de frecuencia
    phidata = phiol.(wlin)
    # Tiene en cuanta las "vueltas" del desfase
    DSP.unwrap!(phidata)
    # Calcularemos los valores del desfase mediante interpolación
    # lineal
    itp = LinearInterpolation(wlin,phidata)
    
    # Calculo de la frecuencia de cruce (desfase = -180°) 
    # y la RAco, si es necesario
    if co
        f(x) = itp(x)+pi
        wco = find_zero(f, (wmin, wmax))
        RAco = RAol(wco)
        phico = itp(wco)
    end
    
    #Cálculo de la frecuencia para RA=1 y su desfase
    if ra1    
        w1 = find_zero(w->RAol(w)-1, (wmin, wmax))
        RA1 = RAol(w1)
        phi1 = itp(w1)
    end

    # Representación del diagrama de Bode
    #l = @layout [a; b]
    RAplot = plot(RAol, wmin, wmax, xscale=:log10, yscale=:log10,
        legend=false, lw=2, xlabel="", ylabel=RAlabel,
        minorticks=:auto)
    # Encuentra el valor mínimo de la escala del eje y 
    RAminscale = ylims(RAplot)[1]
    if co
        plot!([wmin, wco, wco],[RAco, RAco, RAminscale], color=:red, 
            annotations= (wco, float(RAco), text(L"w_{co}, RA_{co}", pointsize=10, :left)))
    end
    if ra1
        plot!([wmin, w1, w1], [RA1, RA1, RAminscale], color=:lime, 
            annotations= (w1, float(RA1), text(L"w_1, 1", pointsize=10, :left)))
    end
    phiplot = plot(wlin, phidata*180/pi, xscale=:log10,
        legend=false, lw=2, xlabel="ω", ylabel="φ",
        minorticks=:auto)
    phiminscale = ylims(phiplot)[1]
    if co
        plot!([wmin, wco, wco],[phico*180/pi, phico*180/pi, phiminscale],
            color=:red, annotations= (wco ,float(phico*180/pi), text(L"w_{co}, -180\degree",
                    pointsize=10, :left)))
    end
    if ra1
        plot!([wmin, w1, w1], [phi1*180/pi, phi1*180/pi, phiminscale], color=:lime,
        annotations = (w1, float(phi1*180/pi), text(L"w_1, \phi_1", pointsize=10, :left)))
    end
    
    fig = plot(RAplot, phiplot, layout=grid(2,1), show=true)
    output = Bode(fig, wco, RAco, w1, phi1)
    return output
end



end  # fin del módulo
