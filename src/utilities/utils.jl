function is_leap_feb29(date)
    y = year(date)
    return (isleapyear(y) && month(date) == 2 && day(date) == 29)
end

# Calculate the day of the year assuming a fixed 365-day year (non-leap year assumption)
function day_of_year_fixed(date)
    # Days in each month (fixed, non-leap year)
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    m = month(date)
    d = day(date)
    # Sum the days in all previous months and add the current day
    return sum(days_in_month[1:m-1]) + d
end

# Return the day index counting from 2003-01-01 (which is considered day 1)
function day_index(date, 
                   base_year::Integer
)
    # Calculate the number of full years since 2003 (each year assumed to be 365 days)
    full_years = year(date) - base_year
    return full_years * 365 + day_of_year_fixed(date)
end
