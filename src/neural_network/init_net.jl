function get_mlp((input_size, output_size)::Pair,
                 hidden_layers,
                 hidden_width,
                 activation,
                 net_struct,
                 rng #rng = Xoshiro(0), random number generator
)
    nn = Chain(
        Dense(input_size => hidden_width, activation),                       # Input layer
        [Dense(hidden_width => hidden_width, activation) for _ = 1:hidden_layers]...,  # Remaining hidden layers
        Dense(hidden_width => output_size),                                  # Output layer
    )
    mlp = net_struct(nn)
    ps, st = Lux.setup(rng, mlp)
    return mlp, ps, st
end

function get_mlp((input_size, output_size)::Pair,
                 hidden_layers,
                 hidden_width,
                 activation,
                 net_struct,
                 rng, #rng = Xoshiro(0), random number generator
                 constrain_f
)
    nn = Chain(
        Dense(input_size => hidden_width, activation),                       # Input layer
        [Dense(hidden_width => hidden_width, activation) for _ = 1:hidden_layers]...,  # Remaining hidden layers
        Dense(hidden_width => output_size, constrain_f),                                  # Output layer
    )
    mlp = net_struct(nn)
    ps, st = Lux.setup(rng, mlp)
    return mlp, ps, st
end

function get_node((input_size, output_size)::Pair,
                  hidden_layers,
                  hidden_width,
                  activation,
                  net_struct,
                  solver,
                  rng #rng = Xoshiro(0), random number generator
)
    nn = Chain(
        Dense(input_size => hidden_width, activation),                       # Input layer
        [Dense(hidden_width => hidden_width, activation) for _ = 1:hidden_layers]...,  # Remaining hidden layers
        Dense(hidden_width => output_size),                                  # Output layer
    )
    node = net_struct(nn; solver = solver)
    ps, st = Lux.setup(rng, node)
    return node, ps, st
end

function get_node((input_size, output_size)::Pair,
                  hidden_layers,
                  hidden_width,
                  activation,
                  net_struct,
                  solver,
                  rng, #rng = Xoshiro(0), random number generator
                  constrain_f
)
    nn = Chain(
        Dense(input_size => hidden_width, activation),                       # Input layer
        [Dense(hidden_width => hidden_width, activation) for _ = 1:hidden_layers]...,  # Remaining hidden layers
        Dense(hidden_width => output_size, constrain_f),                                  # Output layer
    )
    node = net_struct(nn; solver = solver)
    ps, st = Lux.setup(rng, node)
    return node, ps, st
end


# function residual_block(neurons, activation, layers)
#     hidden_layers = layers - 2
#     f = Chain(
#         Dense(neurons => neurons, activation),                       # Input layer
#         [Dense(neurons => neurons, activation) for _ = 1:hidden_layers]...,  # Remaining hidden layers
#         Dense(neurons => neurons),                                  # Output layer
#     )
#     return SkipConnection(f, +)
# end

# function epd_mlp(input_dim, neurons, output_dim, activation, layers, blocks, net_struct, rng)
    
#     # Step 1: Encode
#     encoder = Dense(input_dim, neurons)

#     # Step 2: Process (residual-style update)
#     process = Chain([
#         residual_block(neurons, activation, layers) for _ in 1:blocks
#     ]...)

#     # Step 3: Decode
#     decoder = Dense(neurons, output_dim)

#     model = Chain(encoder, process, decoder)

#     mlp = net_struct(model)
#     ps, st = Lux.setup(rng, mlp)
#     return mlp, ps, st
# end


# function epd_mlp(input_dim, neurons, output_dim, activation, layers, blocks, net_struct, rng, constrain_f)
    
#     # Step 1: Encode
#     encoder = Dense(input_dim, neurons)

#     # Step 2: Process (residual-style update)
#     process = Chain([
#         residual_block(neurons, activation, layers) for _ in 1:blocks
#     ]...)

#     # Step 3: Decode
#     decoder = Dense(neurons, output_dim, constrain_f)

#     model = Chain(encoder, process, decoder)

#     mlp = net_struct(model)
#     ps, st = Lux.setup(rng, mlp)
#     return mlp, ps, st
# end

# function epd_node(input_dim, neurons, output_dim, activation, layers, blocks, net_struct, solver, rng)
    
#     # Step 1: Encode
#     encoder = Dense(input_dim, neurons)

#     # Step 2: Process (residual-style update)
#     process = Chain([
#         residual_block(neurons, activation, layers) for _ in 1:blocks
#     ]...)

#     # Step 3: Decode
#     decoder = Dense(neurons, output_dim)
    
#     model = Chain(encoder, process, decoder)
    
#     node = net_struct(model; solver = solver)
#     ps, st = Lux.setup(rng, node)
#     return node, ps, st
    
# end

# function epd_node(input_dim, neurons, output_dim, activation, layers, blocks, net_struct, solver, rng, constrain_f)
    
#         # Step 1: Encode
#     encoder = Dense(input_dim, neurons)

#     # Step 2: Process (residual-style update)
#     process = Chain([
#         residual_block(neurons, activation, layers) for _ in 1:blocks
#     ]...)

#     # Step 3: Decode
#     decoder = Dense(neurons, output_dim, constrain_f)
    
#     model = Chain(encoder, process, decoder)
    
#     node = net_struct(model; solver = solver)
#     ps, st = Lux.setup(rng, node)
#     return node, ps, st
    
# end