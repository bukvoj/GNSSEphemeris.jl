"""
    getsvpos(
            time, 
            toe::Number,
            gpsweek::Int,
            svclockbias::Real, svclockdrift::Real, svclockdriftrate::Real,
            sqrtA::Real,
            M₀::Real,
            Δn::Real,
            e::Real,
            ω::Real,
            cus::Real, cuc::Real, crs::Real, crc::Real, cis::Real, cic::Real,
            i0::Real, 
            idot::Real,
            Ω₀::Real, Ωdot::Real,
            rho::Real,
            timeoffset::Number = 0
            )

    
Returns satellite position and velocity in ECEF frame. The return type is ECEF object from Geodesy package.

Arguments:
    time: Time of interest in GPS time. Can be a DateTime, TimeDate or Real (seconds since epoch).
    toe: Time of ephemeris in seconds since epoch.
    gpsweek: GPS week number.
    svclockbias: Satellite clock bias in seconds.
    svclockdrift: Satellite clock drift in seconds per second.
    svclockdriftrate: Satellite clock drift rate in seconds per second squared.
    sqrtA: Square root of the semi-major axis in meters.
    M₀: Mean anomaly at reference time in radians.
    Δn: Mean motion difference from computed value in radians per second.
    e: Eccentricity of the satellite orbit.
    ω: Argument of perigee in radians.
    cus, cuc, crs, crc, cis, cic: Harmonic coefficients for the satellite orbit.
    i0: Inclination angle at reference time in radians.
    idot: Rate of inclination angle in radians per second.
    Ω₀: Longitude of ascending node at reference time in radians.
    Ωdot: Rate of right ascension in radians per second.
    rho: Optional parameter to account for the time delay due to the signal travel time. Default is 0.0.
    timeoffset: When different time than GPS is used, this is the offset in seconds.
"""
function getsvpos(
                time::Real, 
                toe::Number,
                gpsweek::Int,
                svclockbias::Real, svclockdrift::Real, svclockdriftrate::Real,
                sqrtA::Real,
                M₀::Real,
                Δn::Real,
                e::Real,
                ω::Real,
                cus::Real, cuc::Real, crs::Real, crc::Real, cis::Real, cic::Real,
                i0::Real, 
                idot::Real,
                Ω₀::Real, Ωdot::Real,
                rho::Real = 0.0,
                timeoffset::Number = 0
                )
    # Define constants
    μ = 3.986005e14             # m^3/s^2 ::: Earth's gravitational constant
    Ωₑdot = 7.2921151467e-5     # rad/s ::: Earth's rotation rate
    ω_e = 7.2921151467e-5
    c = 299792458
    
    # Compute the time of ephemeris
    epochtoe = weeksecondstoutc(gpsweek, toe)
    epochtoe = datetime2unix(epochtoe)

    # Compute the time of transmission
    time = time - timeoffset - rho/c
    dt = svclockbias + (time-epochtoe)*svclockdrift + 1/2*svclockdriftrate*(time-epochtoe).^2
    time = time + dt

    tₖ = time - epochtoe
    if tₖ < -302400
        tₖ += 604800
    elseif tₖ > 302400
        tₖ -= 604800
    end

    # A₀ = sqrtA^2      # sometimes different calculation is used, hence A0 and then Aₖ. I only use Aₖ
    Aₖ = sqrtA^2        #A₀

    n = sqrt(μ)/(sqrtA^3) + Δn
    Mₖ = M₀ + n * tₖ

    # Newton-Raphson iteration for eccentric anomaly
    Eₖ = Mₖ
    for i in 1:25 # Would it be better to use a while loop with a tolerance?
        Eₖ = Eₖ + (Mₖ - Eₖ + e * sin(Eₖ))/(1-e*cos(Eₖ))
    end
    
    vₖ =  2*atan(sqrt((1+e)/(1-e)) * tan(Eₖ/2))
    # vₖ = atan((sqrt(1 - e^2) * sin(Eₖ))/(1 - e* cos(Eₖ)), (cos(Eₖ) - e)/ (1 - e* cos(Eₖ))) # This is the same as above, but less readable

    Φₖ = vₖ + ω
    δuₖ = cus * sin(2*Φₖ) + cuc * cos(2*Φₖ)
    δrₖ = crs * sin(2*Φₖ) + crc * cos(2*Φₖ)
    δiₖ = cis * sin(2*Φₖ) + cic * cos(2*Φₖ)

    uₖ = Φₖ + δuₖ
    rₖ = Aₖ * (1 - e * cos(Eₖ)) + δrₖ 
    iₖ = i0 + idot * tₖ + δiₖ    

    Xk′ = rₖ * cos(uₖ)
    Yk′ = rₖ * sin(uₖ)

    Ωₖ = Ω₀ + (Ωdot - Ωₑdot) * tₖ - Ωₑdot * toe
    
    Xk = Xk′ * cos(Ωₖ) - Yk′ * cos(iₖ) * sin(Ωₖ)
    Yk = Xk′ * sin(Ωₖ) + Yk′ * cos(iₖ) * cos(Ωₖ)
    Zk = Yk′ * sin(iₖ)


    # velocities:
    Eₖdot = n / (1 - e * cos(Eₖ))

    vₖdot = Eₖdot * sqrt(1-e^2) / (1 - e * cos(Eₖ))
    iₖdot = idot + 2*vₖdot*(cis*cos(2*Φₖ) - cic*sin(2*Φₖ))
    uₖdot = vₖdot + 2 * vₖdot * (cus*cos(2*Φₖ) - cuc*sin(2*Φₖ))
    rₖdot = Aₖ*e*sin(Eₖ)*Eₖdot + 2*vₖdot*(crs*cos(2*Φₖ) - crc*sin(2*Φₖ))

    Ωₖdot = Ωdot - Ωₑdot

    Xkdot′ = rₖdot * cos(uₖ) - rₖ * uₖdot * sin(uₖ)
    Ykdot′ = rₖdot * sin(uₖ) + rₖ * uₖdot * cos(uₖ)

    Xkdot = (-Xk′*Ωₖdot*sin(Ωₖ) + Xkdot′*cos(Ωₖ) - Ykdot′*cos(iₖ)*sin(Ωₖ)
        -Yk′*(Ωₖdot*cos(Ωₖ)*cos(iₖ)-iₖdot*sin(iₖ)*sin(Ωₖ)))
    Ykdot = (Xk′*Ωₖdot*cos(Ωₖ) + Xkdot′*sin(Ωₖ) + Ykdot′*cos(iₖ)*cos(Ωₖ)
        -Yk′*(Ωₖdot*sin(Ωₖ)*cos(iₖ)+iₖdot*sin(iₖ)*cos(Ωₖ)))
    Zkdot = Ykdot′*sin(iₖ) + Yk′*iₖdot*cos(iₖ)
    
    # Rotate the position and velocity vectors to account for the signal travel time
    α = ω_e .* rho ./ c
    R = [cos(α) sin(α) 0;-sin(α) cos(α) 0;0 0 1]
    
    svpos = ECEF((R*[Xk, Yk, Zk])...)
    svvel = ECEF((R*[Xkdot,Ykdot,Zkdot])...)
    return svpos, svvel
end


# Alternative methods to get the satellite position and velocity
function getsvpos(
    time::DateTime, 
    toe::Number,
    gpsweek::Int,
    svclockbias::Real, svclockdrift::Real, svclockdriftrate::Real,
    sqrtA::Real,
    M₀::Real,
    Δn::Real,
    e::Real,
    ω::Real,
    cus::Real, cuc::Real, crs::Real, crc::Real, cis::Real, cic::Real,
    i0::Real, 
    idot::Real,
    Ω₀::Real, Ωdot::Real,
    rho::Real = 0.0,
    timeoffset::Number = 0
    )
    time = datetime2unix(time)
    return getsvpos(time, toe, gpsweek, svclockbias, svclockdrift, svclockdriftrate, sqrtA, M₀, Δn, e, ω, cus, cuc, crs, crc, cis, cic, i0, idot, Ω₀, Ωdot, rho, timeoffset)
end

function getsvpos(
    time::TimeDate, 
    toe::Number,
    gpsweek::Int,
    svclockbias::Real, svclockdrift::Real, svclockdriftrate::Real,
    sqrtA::Real,
    M₀::Real,
    Δn::Real,
    e::Real,
    ω::Real,
    cus::Real, cuc::Real, crs::Real, crc::Real, cis::Real, cic::Real,
    i0::Real, 
    idot::Real,
    Ω₀::Real, Ωdot::Real,
    rho::Real = 0.0,
    timeoffset::Number = 0
    )
    time = timedates2unix(time)
    return getsvpos(time, toe, gpsweek, svclockbias, svclockdrift, svclockdriftrate, sqrtA, M₀, Δn, e, ω, cus, cuc, crs, crc, cis, cic, i0, idot, Ω₀, Ωdot, rho, timeoffset)
end

# Methods that work with RinexContent navigation data (RinexRead.jl)
"""
    getsvpos(
        time, 
        svid::Int, 
        constellation::Char, 
        navdata::RinexContent, 
        rho::Real = 0.0
    )

Returns satellite position and velocity in ECEF frame. The return type is ECEF object from Geodesy package.
Uses the navigation data from the RinexContent object obtained using RinexRead.jl.

Arguments:
    time: Time of interest in GPS time. Can be a DateTime, TimeDate or Real (seconds since epoch).
    svid: Satellite ID.
    constellation: Constellation character ('G', 'R', 'E', 'C', 'J', 'I').
    navdata: RinexContent object containing the navigation data.
    rho: Optional parameter to account for the time delay due to the signal travel time. Default is 0.0.

    Constellation characters:
        'G' - GPS
        'R' - GLONASS   NOT YET IMPLEMENTED
        'E' - GALILEO
        'C' - BEIDOU
        'J' - QZSS  NOT YET IMPLEMENTED
        'I' - IRNSS NOT YET IMPLEMENTED
"""
function getsvpos(time::TimeDate, svid::Int, constellation::Char, navdata::RinexContent, rho::Real = 0.0)
    time = timedate2unix(time)
    svpos, svvel = getsvpos(time, svid, constellation, navdata, rho)
    return svpos, svvel
end

function getsvpos(time::DateTime, svid::Int, constellation::Char, navdata::RinexContent, rho::Real = 0.0)
    time = datetime2unix(time)
    svpos, svvel = getsvpos(time, svid, constellation, navdata, rho)
    return svpos, svvel
end

function getsvpos(time::Real, svid::Int, constellation::Char, navdata::RinexContent, rho::Real = 0.0)
    nav, timeoffset = extractnav(navdata, constellation, svid)
    ephid = argmin(abs.(timedate2unix.(nav.Time) .- time))

    toe = nav.Toe[ephid]
    gpsweek = nav.GPSWeek[ephid]
    svclockbias = nav.SVClockBias[ephid]
    svclockdrift = nav.SVClockDrift[ephid]
    svclockdriftrate = nav.SVClockDriftRate[ephid]
    sqrtA = nav.sqrtA[ephid]
    M₀ = nav.M0[ephid]
    Δn = nav.Delta_n[ephid]
    e = nav.Eccentricity[ephid]
    ω = nav.omega[ephid]
    cus = nav.Cus[ephid]
    cuc = nav.Cuc[ephid]
    crs = nav.Crs[ephid]
    crc = nav.Crc[ephid]
    cis = nav.Cis[ephid]
    cic = nav.Cic[ephid]
    i0 = nav.i0[ephid]
    idot = nav.IDOT[ephid]
    Ω₀ = nav.OMEGA0[ephid]
    Ωdot = nav.OMEGA_DOT[ephid]

    return getsvpos(time, toe, gpsweek, svclockbias, svclockdrift, svclockdriftrate, sqrtA, M₀, Δn, e, ω, cus, cuc, crs, crc, cis, cic, i0, idot, Ω₀, Ωdot, rho, timeoffset)
end


# Helper functions
function extractnav(navdata::RinexContent, constellation::Char, svid::Int)
    if constellation == 'G'
        return navdata.data.GPS[navdata.data.GPS.SatelliteID .== svid, :], 0
    elseif constellation == 'R'
        return navdata.data.GLONASS[navdata.data.GLONASS.SatelliteID .== svid, :], Nan
    elseif constellation == 'E'
        out = navdata.data.GALILEO[navdata.data.GALILEO.SatelliteID .== svid, :]
        out.GPSWeek = out.GALWeek
        return out, 0
    elseif constellation == 'C'
        out = navdata.data.BEIDOU[navdata.data.BEIDOU.SatelliteID .== svid, :]
        out.GPSWeek = out.BDTWeek .+ 1356   #GPS-BDT week offset
        return out, 14                      #BDT GPS second offset
    elseif constellation == 'J'
        return navdata.data.QZSS[navdata.data.QZSS.SatelliteID .== svid, :], Nan
    elseif constellation == 'I'
        return navdata.data.IRNSS[navdata.data.IRNSS.SatelliteID .== svid, :], Nan
    end
end

