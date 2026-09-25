"""Static crop-calendar configuration."""
mutable struct CropCalendarAuxiliary{I}
    sowing_date::I            # Current daily sowing trigger (1–365; zero otherwise).
    prescribed_sowing_date::I # Externally supplied sowing date retained across dynamic resolution.
end

"""Discrete crop events that occur during the current day."""
mutable struct CropEvents{I}
    sowing::I # Current-day sowing/cultivation event flag (0/1).
    harvest::I # Current-day harvest event flag (0/1).
end

function init_crop_calendar_auxiliary(cell_size::Int, device)
    int_state() = device(zeros(Int32, cell_size))
    return CropCalendarAuxiliary(int_state(), int_state())
end

function init_crop_events(cell_size::Int, device)
    int_event() = device(zeros(Int32, cell_size))
    return CropEvents(ntuple(_ -> int_event(), 2)...)
end
