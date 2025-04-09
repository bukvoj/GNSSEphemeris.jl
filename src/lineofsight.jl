"""
    lineofsight(svpos::ECEF, recpos::ECEF)
    lineofsight(svpos::ECEF, recpos::LLA; datum=wgs84)

Computes the unit line of sight vector from the receiver position to the satellite position.
"""
function lineofsight(svpos::ECEF, recpos::ECEF)
    los = svpos - recpos
    los = los / norm(los)
    return los
end

function lineofsight(svpos::ECEF, recpos::LLA; datum=wgs84)
    return lineofsight(svpos, ECEF(recpos, datum))
end