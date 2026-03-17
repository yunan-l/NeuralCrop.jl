function plot_loss_curve(epoch::AbstractArray{S},
                         train_loss::AbstractArray{T},
                         val_loss::AbstractArray{T},
                         fig_path::String
) where {T <: AbstractFloat, S <: Integer}

    plt = plot(xtickfont = font(12), ytickfont = font(12), guidefont = font(14), legendfont = font(12), titlefont = font(16), size = (1800, 600), bottom_margin=5Plots.mm)

    plot!(plt, epoch, train_loss; xlabel = "Epoch", ylabel = "MSE", label = "Training", title = "Loss curve", linewidth = 2)
    plot!(plt, epoch, val_loss; xlabel = "Epoch", ylabel = "MSE", label = "Validation", title = "Loss curve", linewidth = 2)

    savefig(plt, fig_path)

end