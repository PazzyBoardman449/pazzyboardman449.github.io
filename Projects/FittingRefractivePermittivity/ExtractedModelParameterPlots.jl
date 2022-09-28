using Statistics
using DataFrames
using CSV
#using Plots
using LaTeXStrings
using CairoMakie
using ColorSchemes

#Code to plot the model outputs

Data = CSV.read(
   "/Users/patricboardman/PhD Local Files/CODE/DATA/18 May/Collated Data/ThickSampleFittedParams.csv",
   DataFrame,
)

t = Data[:, 1]
L = Data[:, 2]

N_Scans = length(t)

p_w_WA = Data[:, 3]
p_w_MG = Data[:, 4]
p_w_BR = Data[:, 5]

ε_gel_WA = Data[:, 6] .+ im * Data[:, 9]
ε_gel_WA = Data[:, 7] .+ im * Data[:, 10]
ε_gel_BR = Data[:, 8] .+ im * Data[:, 11]

n_140GHz_WA = Data[:, 12]
n_140GHz_MG = Data[:, 13]
n_140GHz_BR = Data[:, 14]

κ_140GHz_WA = Data[:, 15]
κ_140GHz_MG = Data[:, 16]
κ_140GHz_BR = Data[:, 17]


#Plot of refractive index vs water concentration

fontsize_theme = Theme(fontsize = 20, font = "Garamond")
set_theme!(fontsize_theme)



fig = Figure()
supertitle =
   Label(fig[1, 1:2], "Phantom sample extrapolation to 140GHz", textsize = 30)
scatter(
   fig[2, 1],
   p_w_WA,
   n_140GHz_WA,
   color = 1:N_Scans,
   colormap = :haline,
   rev = true,
   marker = :circle,
   markersize = 13,
   label = "Weighted Average",
   axis = (
      aspect = 1,
      xlabel = L"$p_{w}$ %",
      title = "Refractive Index",
      ylabel = L"$n'(p_{w})$",
   ),
)
scatter!(
   fig[2, 1],
   p_w_MG,
   n_140GHz_MG,
   color = 1:N_Scans,
   colormap = :haline,
   rev = true,
   marker = :utriangle,
   markersize = 13,
   label = "Maxwell-Garnett",
)
scatter!(
   fig[2, 1],
   p_w_BR,
   n_140GHz_BR,
   color = 1:N_Scans,
   colormap = :haline,
   rev = true,
   marker = :rect,
   markersize = 13,
   label = "Bruggemann",
)
axislegend(position = (5.1, 0))
axislegend(position = (1.0, 0))
scatter(
   fig[2, 2],
   p_w_WA,
   κ_140GHz_WA,
   color = 1:N_Scans,
   colormap = :haline,
   rev = true,
   marker = :circle,
   markersize = 13,
   label = "Weighted Average",
   axis = (
      aspect = 1,
      xlabel = L"$p_{w}$ %",
      title = "Extinction Coefficient",
      ylabel = L"$\kappa'(p_{w})$",
   ),
)
scatter!(
   fig[2, 2],
   p_w_MG,
   κ_140GHz_MG,
   color = 1:N_Scans,
   colormap = :haline,
   rev = true,
   marker = :utriangle,
   markersize = 13,
   label = "Maxwell-Garnett",
)
scatter!(
   fig[2, 2],
   p_w_BR,
   κ_140GHz_BR,
   color = 1:N_Scans,
   colormap = :haline,
   rev = true,
   marker = :rect,
   markersize = 13,
   label = "Bruggemann",
)
Colorbar(
   fig[3, 1:2],
   limits = (0, 300),
   colormap = :haline,
   label = "Time/ mins",
   justification = :left,
   vertical = false,
)
fig



save(
   "/Users/patricboardman/PhD Local Files/CODE/Outputs/2022-05-30/18MayExtrapolatedParams.pdf", fig
)


#Plot of Phantom index at 140GHz over time:
