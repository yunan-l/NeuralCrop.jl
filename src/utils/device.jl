_adapt_to_device(::typeof(identity), value) = value
_adapt_to_device(device, value) = Adapt.adapt(device, value)
