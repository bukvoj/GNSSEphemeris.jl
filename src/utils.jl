function weeksecondstoutc(week::Number, tow::Number)
    epoch = DateTime(1980,1,6)
    weekseconds = week * 604800 + tow
    return epoch + Dates.Second(weekseconds)
end

function timedate2unix(timedate::TimeDate)
    return datetime2unix(DateTime(timedate)) + Microsecond(timedate).value/1e6 + Nanosecond(timedate).value/1e9
end