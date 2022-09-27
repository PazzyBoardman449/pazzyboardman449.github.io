using DelimitedFiles
using DataFrames
using CSV
using NLsolve
using LaTeXStrings
using LsqFit
using CairoMakie

#Debye Model Parameters (Free water)
ε_∞_FW = 3.5 #Assume this to be constant across all frequencies
Δε_1_FW = 73.5 #Asume this to be constant across all freqiencies
Δε_2_FW = 1.4 #Asume this to be constant across all freqiencies
τ_D_1_FW = 8.2E-12
τ_D_2_FW = 0.18E-12

#Debye Model Parameters (Bound water)
ε_∞_BW = 3.2 #Assume this to be constant across all frequencies
Δε_1_BW = 37.45 #Asume this to be constant across all freqiencies
Δε_2_BW = 1.67 #Asume this to be constant across all freqiencies
τ_D_1_BW = 284.1E-12
τ_D_2_BW = 0.248E-15

#Reading in the Data
Data = CSV.read(
   "/Users/patricboardman/PhD Local Files/CODE/DATA/29 April/Data/Collated Data/ThickSampleResults.csv",
   DataFrame,
)

L_Data = CSV.read(
   "/Users/patricboardman/PhD Local Files/CODE/DATA/29 April/Data/Collated Data/ThickThicknessData.csv",
   DataFrame,
)

#For the real part of the index:
function EffectiveMedium(f, p)
   #Three component effective medium theory, treating free water, bound water and collagen separately.

   ε_FW =
      ε_∞_FW .+ Δε_1_FW ./ (1 .+ 2 * pi * im * τ_D_1_FW .* f) .+
      Δε_2_FW ./ (1 .+ 2 * pi * im * f * τ_D_2_FW)

   ε_eff = p[1] * ε_FW .+ (p[2] + im * p[3]) .* (1 - p[1])

   #=
   ε_eff_new =
      p[1] * ε_BW .+ p[2] * (p[3] + im * p[4]) .+ (1 - p[1] - p[2])* ε_FW
   =#

   n_eff = sqrt.(conj.(ε_eff))
   #n_eff_new = sqrt.(conj.(ε_eff_new))

   return n_eff

end

#For the real part of the index:

function EffectiveMediumMaxwell(f, p)
   ε_FW =
      ε_∞_FW .+ Δε_1_FW ./ (1 .+ 2 * pi * im * τ_D_1_FW .* f) .+
      Δε_2_FW ./ (1 .+ 2 * pi * im * f * τ_D_2_FW)

      #=

   ε_BW =
      ε_∞_BW .+ Δε_1_BW ./ (1 .+ 2 * pi * im * τ_D_1_BW .* f) .+
      Δε_2_BW ./ (1 .+ 2 * pi * im * f * τ_D_2_BW)

      =#

   ε_eff =
      ε_FW .*
      (ε_FW .+ ((3 - 2 .* p[1]) ./ 3) .* ((p[2] + im * p[3]) .- ε_FW)) ./
      (ε_FW .+ (p[1] / 3) .* ((p[2] + im * p[3]) .- ε_FW))

   n_eff = sqrt.(conj.(ε_eff))

   return n_eff
end



function EffectiveMediumBruggeman(f, p)
   ε_FW =
      ε_∞_FW .+ Δε_1_FW ./ (1 .+ 2 * pi * im * τ_D_1_FW .* f) .+
      Δε_2_FW ./ (1 .+ 2 * pi * im * f * τ_D_2_FW)

   ε_eff =
      (1 / 12) * (
         6 .* ε_FW .- 3 .* (p[2] .+ im * p[3]) .+
         9 .* ((p[2] + im * p[3]) .- ε_FW) .* (1 - p[1]) .+
         sqrt.(
            72 .* ε_FW .* (p[2] + im * p[3]) +
            (
               6 .* ε_FW .- 3 .* (p[2] + im * p[3]) .+
               9 .* ((p[2] .+ im * p[3]) .- ε_FW) .* (1 - p[1])
            ) .^ 2,
         )
      )

   n_eff = sqrt.(conj.(ε_eff))

   return n_eff
end


N_Freq = size(Data)[1]
N_Scans = size(Data)[2] - 1

freq_arr = real.(parse.(ComplexF64, Data[:, 1]))
n_results = zeros(Float64, N_Freq, N_Scans)
κ_results = zeros(Float64, N_Freq, N_Scans)
Results = zeros(ComplexF64, N_Freq, N_Scans)


for i = 1:N_Freq
   for j = 2:N_Scans+1
      n_results[i, j-1] = real(parse.(ComplexF64, Data[i, j]))
      κ_results[i, j-1] = imag(parse.(ComplexF64, Data[i, j]))
      Results[i, j-1] = parse.(ComplexF64, Data[i, j])
   end
end

scan_num = 10

#p0 = [0.3, 3, 0.2]
p0 = [0.5, 3.5, 0.5]

TheFit = curve_fit(
   EffectiveMedium,
   freq_arr[31:150],
   Results[31:150, scan_num],
   p0,
   lower = float.([0, 1, -4]),
   upper = float.([1, 5, 4]),
)
TheFitMaxwell = curve_fit(
   EffectiveMediumMaxwell,
   freq_arr[31:150],
   Results[31:150, scan_num],
   p0,
   lower = float.([0, 1, -4]),
   upper = float.([1, 5, 4]),
)
TheFitBruggeman = curve_fit(
   EffectiveMediumBruggeman,
   freq_arr[31:150],
   Results[31:150, scan_num],
   p0,
   lower = float.([0, 1, -4]),
   upper = float.([1, 5, 4]),
)

n_fit = [
   real(EffectiveMedium.(freq_arr, Ref(TheFit.param))),
   real(EffectiveMediumMaxwell.(freq_arr, Ref(TheFitMaxwell.param))),
   real(EffectiveMediumBruggeman.(freq_arr, Ref(TheFitBruggeman.param))),
]
κ_fit = [
   imag(EffectiveMedium.(freq_arr, Ref(TheFit.param))),
   imag(EffectiveMediumMaxwell.(freq_arr, Ref(TheFitMaxwell.param))),
   imag(EffectiveMediumBruggeman.(freq_arr, Ref(TheFitBruggeman.param))),
]


fitted_params = zeros(Float64, N_Scans, 17)

for i = 1:N_Scans
   TheFit = curve_fit(
      EffectiveMedium,
      freq_arr[30:150],
      Results[30:150, i],
      p0,
      lower = float.([0, 1, -4]),
      upper = float.([1, 5, 4]),
   )
   TheFitMaxwell = curve_fit(
      EffectiveMediumMaxwell,
      freq_arr[30:150],
      Results[30:150, i],
      p0,
      lower = float.([0, 1, -4]),
      upper = float.([1, 5, 4]),
   )
   TheFitBruggeman = curve_fit(
      EffectiveMediumBruggeman,
      freq_arr[30:150],
      Results[30:150, i],
      p0,
      lower = float.([0, 1, -4]),
      upper = float.([1, 5, 4]),
   )

   fitted_params[i, 1] = L_Data[i, 1]
   fitted_params[i, 2] = L_Data[i, 2]
   fitted_params[i, 3] = round(TheFit.param[1] * 100, digits = 1)
   fitted_params[i, 4] = round(TheFitMaxwell.param[1] * 100, digits = 1)
   fitted_params[i, 5] = round(TheFitBruggeman.param[1] * 100, digits = 1)
   fitted_params[i, 6] = round(TheFit.param[2], digits = 3)
   fitted_params[i, 7] = round(TheFitMaxwell.param[2], digits = 3)
   fitted_params[i, 8] = round(TheFitBruggeman.param[2], digits = 3)
   fitted_params[i, 9] = round(TheFit.param[3], digits = 3)
   fitted_params[i, 10] = round(TheFitMaxwell.param[3], digits = 3)
   fitted_params[i, 11] = round(TheFitBruggeman.param[3], digits = 3)
   fitted_params[i, 12] =
      real(EffectiveMedium.(freq_arr, Ref(TheFit.param)))[14]
   fitted_params[i, 13] =
      real(EffectiveMedium.(freq_arr, Ref(TheFitMaxwell.param)))[14]
   fitted_params[i, 14] =
      real(EffectiveMedium.(freq_arr, Ref(TheFitBruggeman.param)))[14]
   fitted_params[i, 15] =
      imag(EffectiveMedium.(freq_arr, Ref(TheFit.param)))[14]
   fitted_params[i, 16] =
      imag(EffectiveMedium.(freq_arr, Ref(TheFitMaxwell.param)))[14]
   fitted_params[i, 17] =
      imag(EffectiveMedium.(freq_arr, Ref(TheFitBruggeman.param)))[14]
end


fontsize_theme = Theme(fontsize = 16, font = "Garamond")
set_theme!(fontsize_theme)

ThePlot = Figure()
supertitle = Label(
   ThePlot[1, 1:2],
   "Thick fresh phantom sample with 3 parameter fits",
   textsize = 24,
)
ho = scatter(
   ThePlot[2, 1],
   freq_arr / 1E12,
   n_results[:, scan_num],
   markersize = 3,
   label = "Data",
   axis = (
      limits = (0.3, 1.5, 1.5, 2.5),
      aspect = 1,
      xlabel = L"$f$ THz",
      title = "Refractive Index",
      ylabel = L"n(f)",
   ),
)
hi = lines!(ThePlot[2, 1], freq_arr / 1E12, n_fit[1], linewidth = 1)
pazzylol = lines!(ThePlot[2, 1], freq_arr / 1E12, n_fit[2],linewidth = 1)
pazzy = lines!(ThePlot[2, 1], freq_arr / 1E12, n_fit[3], linewidth = 1)
scatter(
   ThePlot[2, 2],
   freq_arr / 1E12,
   κ_results[:, scan_num],
   markersize = 3,
   axis = (
      limits = (0.3, 1.5, 0, 1),
      aspect = 1,
      xlabel = L"$f$ THz",
      title = "Extinction Coefficient",
      ylabel = L"$\kappa(f)$",
   ),
)
lines!(ThePlot[2, 2], freq_arr / 1E12, κ_fit[1], label = "Weighted Average", linewidth = 1)
lines!(ThePlot[2, 2], freq_arr / 1E12, κ_fit[2], label = "Maxwell-Garnett", linewidth = 1)
lines!(ThePlot[2, 2], freq_arr / 1E12, κ_fit[3], label = "Bruggemann", linewidth = 1)
Legend(
   ThePlot[3, 1:2],
   [hi, pazzylol, pazzy],
   [
      "Debye/Weighted Average: Gel Fraction: = $(100 - fitted_params[scan_num,3])% and Gel Permittivity = $(fitted_params[scan_num,6] - im*fitted_params[scan_num,9])",
      "Debye/Maxwell Garnett: Gel Fraction: = $(100 - fitted_params[scan_num,4])% and Gel Permittivity = $(fitted_params[scan_num,7] - im*fitted_params[scan_num,10])",
      "Debye/Bruggeman: Gel Fraction: = $(100 - fitted_params[scan_num,5])% and Gel Permittivity = $(fitted_params[scan_num,8] - im*fitted_params[scan_num,11])",
   ],
)
#trim!(ThePlot.layout)
ThePlot


save(
   "/Users/patricboardman/PhD Local Files/CODE/Outputs/2022-05-30/ThickFit.pdf", ThePlot
)


#CSV.write("/Users/patricboardman/PhD Local Files/CODE/DATA/18 May/Collated Data/ThickSampleFittedParams.csv", Tables.table(fitted_params))
