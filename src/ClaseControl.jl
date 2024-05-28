module ClaseControl

export bode

using DSP, Roots, Interpolations, PlutoPlotly

struct Bode
    fig
    wco
    RAco
    w1
    phi1
end

"""
bode(Gol; wmin=1e-1, wmax=1e1, points=100, co=false, ra1=false, RAlabel="RA")

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
function bode(Gol::Function; wmin=1e-1, wmax=1e1, points=100, co=false, ra1=false, RAlabel="RA")
    wco = nothing
    RAco = nothing
    w1 = nothing
    phi1 = nothing
    
    trace_RAco = scatter()
    trace_RAra1 = scatter()
    trace_phico = scatter()
    trace_phira1 = scatter()

    RAol(w) = abs(Gol(im*w))
    phiol(w) = angle(Gol(im*w))
    
    # Hace una lista equiespaciada para los valores del log10 del desfase
    wlog = range(log10(wmin), log10(wmax); length=points)
    # Genera la versión lineal
    wlin = 10 .^wlog
    
    # Cálculo de los datos de RA para los valores de frecuencia
    RAdata = RAol.(wlin)
    
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
    trace_RA = scatter(; x=wlin, y=RAdata)
    # Encuentra el valor mínimo de la escala del eje y 
    RAminscale = minimum(RAdata)

    if co
        trace_RAco = scatter(; x=[wmin, wco, wco], y=[RAco, RAco, RAminscale],
                            mode="markers+lines+text",
                            text=["RA<sub>co</sub>", "", "ω<sub>co</sub>"],
                            textposition="left",
                            hoverinfo="text",
                            hovertext=["RA<sub>co</sub> = $(round(RAco; sigdigits=3))", "ω<sub>co</sub> = $(round(wco; sigdigits=3)), RA<sub>co</sub> = $(round(RAco; sigdigits=3))", "ω<sub>co</sub> = $(round(wco; sigdigits=3))"])
        trace_RAco["maker.color"] = "red"
    end
    if ra1
        trace_RAra1 = scatter(; x=[wmin, w1, w1], y=[RA1, RA1, RAminscale],
                        mode="markers+lines+text",
                        text=["1", "", "ω₁"],
                        textposition="left",
                        hoverinfo="text",
                        hovertext=["RA = 1", "RA = 1, ω₁ = $(round(w1; sigdigits=3))", "ω₁ = $(round(w1; sigdigits=3))"])
        trace_RAra1["marker.color"] = "lime"
    end
    trace_phi = scatter(; x=wlin, y=phidata*180/pi)
    phiminscale = minimum(phidata)
    if co
        trace_phico = scatter(; x=[wmin, wco, wco],
                        y= [phico*180/pi, phico*180/pi, phiminscale],
                        mode="markers+lines+text",
                        text=["-180°", "", "ω<sub>co</sub>"],
                        textposition=["left", "", "right"],
                        hoverinfo="text",
                        hovertext=["φ<sub>co</sub> = -180°", "φ<sub>co</sub> = -180°, ω<sub>co</sub> = $(round(wco; sigdigits=3))", "ω<sub>co</sub> = $(round(wco; sigdigits=3))"])
        trace_phico["marker.color"] = "red"
    end
    if ra1
        trace_phira1 = scatter(; x=[wmin, w1, w1],
                        y=[phi1*180/pi, phi1*180/pi, phiminscale],
                        mode="markers+lines+text",
                        text=["φ₁", "", "ω₁"],
                        textposition="left",
                        hoverinfo="text",
                        hovertext=["φ₁ = $(round(phi1*180/pi; sigdigits=3))°", "φ₁ = $(round(phi1*180/pi; sigdigits=3))°, ω₁ = $(round(w1; sigdigits=3))", "ω₁ = $(round(w1; sigdigits=3))"])
        trace_phira1["marker.color"] = "lime"
    end
    
    RAplot = plot([trace_RA, trace_RAco, trace_RAra1], Layout(xaxis_type="log",
                yaxis_type="log", xaxis_title="ω", yaxis_title=RAlabel))
    phiplot = plot([trace_phi, trace_phico, trace_phira1], Layout(xaxis_type="log",
                xaxis_title="ω", yaxis_title="φ"))
    
    fig = [RAplot; phiplot]
    fig = relayout(fig, showlegend=false)
    output = Bode(fig, wco, RAco, w1, phi1)
    return output
end

#function bode(f::Sym; kwargs...)
#	if length(f.free_symbols) > 1
#		error("Demasiadas variables: $(f.free_symbols)")
#	else
#		bode(lambdify(f); kwargs...)
#	end
#end

end  # fin del módulo
