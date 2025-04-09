# TODO ADD DOCSTRING

function lookangles(svpos::ECEF, recpos::ECEF; unit=:rad) 
    los = lineofsight(svpos, recpos)
    lla = LLA(recpos, wgs84)
    lat = lla.lat * pi/180
    lon = lla.lon * pi/180
    h = lla.alt
    ecef2tangent = [ -sin(lon) -sin(lat)*cos(lon) cos(lat)*cos(lon);
                     cos(lon) -sin(lat)*sin(lon) cos(lat)*sin(lon);
                     0 cos(lat) sin(lat)]'
    los = ecef2tangent * los
    az = 1/2*pi - atan(los[2], los[1])
    el = asin(los[3])
    if unit == :deg
        az = az * 180/pi
        el = el * 180/pi
    elseif unit == :rad
        nothing
    else
        error("unit must be :rad or :deg")
    end
    return az, el
end


function lookangles(svpos::Vector{ECEF{T}}, recpos::ECEF; unit=:rad) where T <: Number
    azel = lookangles.(svpos, Ref(recpos); unit=unit)
    az = [a[1] for a in azel]
    el = [e[2] for e in azel]
    return az, el
end

lookangles(svpos::ECEF, recpos::LLA; datum=:wgs84, unit=:rad) = lookangles(svpos, ECEF(recpos, datum); unit=unit)

function lookangles(svpos::Vector{ECEF{T}}, recpos::LLA; unit=:rad)  where T <: Number 
    return lookangles(svpos, ECEF(recpos, wgs84); unit=unit)
end