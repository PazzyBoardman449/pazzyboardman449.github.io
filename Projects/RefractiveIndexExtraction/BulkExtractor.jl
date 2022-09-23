using DelimitedFiles
using DataFrames
using CSV
#using Plots
using FFTW
using LaTeXStrings
using ColorSchemes
using Peaks
using Images
using StatsBase
using NaNStatistics
using CairoMakie

#Parameters and Constants:
const c = 299792458 #Speed of Light in a vacuum
Δt = 0.05E-12 #Time between measurements (ps)
n_sub = 1.5
n_q = 1.966

function RatioMinimizer(L, ω, n_ph, κ_ph, R_data)

      F =

      R_data -

      (
         (2 / (1 + (n_ph + im * κ_ph))) *
         ((2 * (n_ph + im * κ_ph)) / (n_q + (n_ph + im * κ_ph))) *
         exp(im * (n_ph + im * κ_ph - 1) * ω * L / c)
      ) / (
         1 +
         ((1 - (n_ph + im * κ_ph)) / (1 + (n_ph + im * κ_ph))) *
         (((n_ph + im * κ_ph) - n_q) / ((n_ph + im * κ_ph) + n_q)) *
         exp(2 * im * (n_ph + im * κ_ph) * ω * L / c)

      )


   return F
end

#Collated DataFrames
Ref_Data_Coll = CSV.read(
   "/Users/patricboardman/PhD Local Files/CODE/DATA/18 May/Collated Data/Reference.csv",
   DataFrame,
)

Samp_1_Data_Coll = CSV.read(
   "/Users/patricboardman/PhD Local Files/CODE/DATA/18 May/Collated Data/ThickSample.csv",
   DataFrame,
)

L_Data = CSV.read(
   "/Users/patricboardman/PhD Local Files/CODE/DATA/18 May/Collated Data/Thickness.csv",
   DataFrame,
)

t_array = L_Data[:,1]
L_1_array = 1E-3*L_Data[:,2]

function Extractor(freq, Ref_Data, Samp_Data, L_array)

   #Assuming that number of scans is the same:
   N_Scans = size(Ref_Data)[2]
   E_Ref_Data = zeros(Float64,size(Ref_Data)[1],N_Scans)
   E_Samp_Data = zeros(Float64,size(Samp_Data)[1],N_Scans)

   for i in 1:N_Scans
      E_Ref_Data[:,i] = Ref_Data[:,i]
      E_Samp_Data[:,i] = Samp_Data[:,i]
   end

   E_Ref_Data_New = cat(E_Ref_Data, fill(0, 2000 - size(E_Ref_Data)[1], N_Scans); dims = 1)
   E_Samp_Data_New = cat(E_Samp_Data, fill(0, 2000 - size(E_Samp_Data)[1], N_Scans); dims = 1)

   N_Steps = size(E_Ref_Data_New)[1]
   Δf = (1 / Δt) / N_Steps  #Frequency step (Inverse of the time step)
   freq_index = Int(round(freq/(Δf)))

   #Initialise arrays to store the FFT of the field data
   E_Ref_Freq_Data = zeros(ComplexF64,N_Steps,N_Scans)
   E_Samp_Freq_Data = zeros(ComplexF64,N_Steps,N_Scans)
   Sample_Ratio = zeros(ComplexF64,Int(ceil(N_Steps / 10)),N_Scans)

   for i in 1:N_Scans
      E_Ref_Freq_Data[:,i] = fft(E_Ref_Data_New[:,i])
      E_Samp_Freq_Data[:,i] = fft(E_Samp_Data_New[:,i])

      Sample_Ratio[:,i] = (E_Samp_Freq_Data[:,i]./E_Ref_Freq_Data[:,i])[1:Int(N_Steps / 10)]
      Sample_Ratio[:,i] = abs.(Sample_Ratio[:,i]).*(exp.(-im* angle.(Sample_Ratio[:,i])))
   end

   #Extraction of the Index:ß
   N_S = 500 #Number of Steps
   n = 1+2/N_S:2/N_S:3 #Guess range for n
   κ = 0+1/N_S:1/N_S:1 #Guess range for κ

   f_step = 2E12 / size(Sample_Ratio)[1]
   f = 0+f_step:f_step:2E12 #Frequency Range for the refractive index
   ω = 2π * f #Angular Frequency

   F_1_abs = zeros(Float64, length(n), length(κ), length(ω), length(L_array))

   n_1_stored_minima = Array{Vector}(undef, length(κ), length(ω), length(L_array))
   n_1_min_index = zeros(Int64, length(κ), length(ω), length(L_array))
   n_1_min_index_alt = zeros(Int64, length(ω), length(L_array))
   n_1_sol = zeros(Float64, length(ω), length(L_array))

   κ_1_stored_minima = Array{Vector}(undef, length(κ), length(ω), length(L_array))
   #κ_1_min_index = zeros(Int64, length(κ), length(ω), length(L_array))
   κ_1_min_index = zeros(Int64, length(ω), length(L_array))
   κ_1_sol = zeros(Float64, length(ω), length(L_array))

   for l = 1:length(L_array)
      for i = 1:length(κ)
         for j = 1:length(n)
            for k = 1:length(ω)
               F_1_abs[i, j, k, l] = abs(RatioMinimizer(L_array[l], ω[k], n[j], κ[i], Sample_Ratio[k,l]))
            end
         end
      end
   end


   for l = 1:length(L_array)
      for j = 1:length(n)
         n_1_stored_minima[j, 1, l] = [N_S/2,N_S/2]
         κ_1_stored_minima[j, 1, l] = [N_S/2,N_S/2]
         for k = 2:length(ω)
            n_1_stored_minima[j, k, l] = findminima(F_1_abs[j, :, k, l])[1]
            κ_1_stored_minima[j, k, l] = findminima(F_1_abs[j, :, k, l])[1]
            if length(n_1_stored_minima[j, k, l]) > 10 ||
               length(n_1_stored_minima[j, k, l]) == 0
               n_1_stored_minima[j, k, l] = n_1_stored_minima[j, k-1, l]
               κ_1_stored_minima[j, k, l] = κ_1_stored_minima[j, k-1, l]
            else
            end
         end
      end
   end

   for l = 1:length(L_array)
      for j = 1:length(n)
         n_1_min_index[j, 1, l] = n_1_stored_minima[j, 1, l][1]

         n_1_sol[1, l] = first(n_1_min_index[j, 1, l])

         for k = 2:length(ω)
            n_1_min_index[j, k, l] =
               argmin(x -> abs(n_1_min_index[j, k-1, l] - x), n_1_stored_minima[j, k, l])

            n_1_sol[k, l] = n[mode(n_1_min_index[j, k, l])]
         end

      end
   end


   for l = 1:length(L_array)
      for k = 1:length(ω)
         κ_1_min_index[k,l] = findmin(F_1_abs[:, :, k, l])[2][1]
         κ_1_sol[k,l] = κ[κ_1_min_index[k,l]]
      end
   end

   solution = zeros(ComplexF64, length(ω),length(L_array))

   for l = 1:length(L_array)
      for k = 1:length(ω)
         solution[k,l] = n_1_sol[k, l] + im * κ_1_sol[k, l]
      end
   end

   Sample_1_Solution_Data = zeros(ComplexF64,length(ω),length(L_array)+1)

   for k in 1:length(ω)
      Sample_1_Solution_Data[k,1] = f[k]
      for l in 2:size(L_Data)[1] + 1
         Sample_1_Solution_Data[k,l] = solution[k,l-1]
      end
   end

   #CSV.write("/Users/patricboardman/PhD Local Files/CODE/DATA/18 May/Collated Data/Results.csv", Tables.table(Sample_1_Solution_Data))

   return solution, freq_index, Sample_1_Solution_Data
end

freq = 600E9
Solution1 = Extractor(freq, Ref_Data_Coll, Samp_1_Data_Coll, L_1_array)

fontsize_theme = Theme(fontsize = 16, font = "Garamond")
set_theme!(fontsize_theme)


IndexTimePlot = Figure()
supertitle = Label(IndexTimePlot[1, 1:2], "18 May Experiment", textsize = 30)
scatter(IndexTimePlot[2,1:2], t_array, real(Solution1[1])[Solution1[2],:],markersize = 5, color =:blue, axis = (aspect = 3, xlabel = L"$t$ mins", title = "Refractive Index", ylabel = L"$n(t)$"))
lines!(IndexTimePlot[2,1:2], t_array, real(Solution1[1])[Solution1[2],:], label =  "Thick Sample, $(freq/1E9) GHz", color =:blue)
axislegend(position=(0,0))
scatter(IndexTimePlot[3,1:2], t_array, imag(Solution1[1])[Solution1[2],:],markersize = 5, color =:purple, axis = (aspect = 3, xlabel = L"$t$ mins", title = "Extinction Coefficient", ylabel = L"$\kappa(t)$"))
lines!(IndexTimePlot[3,1:2], t_array, imag(Solution1[1])[Solution1[2],:], label = "Thick Sample, $(freq/1E9) GHz",color =:purple)
axislegend(position=(0,0))
IndexTimePlot


#=


ExtinctionCoeffThicknessPlot = Figure()
supertitle = Label(ExtinctionCoeffThicknessPlot[1, 1:2], "29 April Thin Sample", textsize = 30)
scatter(ExtinctionCoeffThicknessPlot[2,1:2], t_array, L_1_array.*imag(Solution1[1])[Solution1[2],:],markersize = 5, color =:blue, axis = (aspect = 3, xlabel = L"$t$ mins", title = "Extinction Coefficient × Thickness", ylabel = L"$\kappa L(t)$"))
lines!(ExtinctionCoeffThicknessPlot[2,1:2], t_array, L_1_array.*imag(Solution1[1])[Solution1[2],:], label =  "κL at 600GHz", color =:blue)
axislegend()
scatter(ExtinctionCoeffThicknessPlot[3,1:2], t_array, 4*pi*freq/c*L_1_array.*imag(Solution1[1])[Solution1[2],:],markersize = 5, color =:green, axis = (aspect = 3, xlabel = L"$t$ mins", title = "Absorption Coefficient × Thickness", ylabel = L"$\alpha L(t)$"))
lines!(ExtinctionCoeffThicknessPlot[3,1:2], t_array, 4*pi*freq/c*L_1_array.*imag(Solution1[1])[Solution1[2],:], label =  "αL at 600GHz", color =:green)
axislegend()
ExtinctionCoeffThicknessPlot

=#
#=
save(
   "/Users/patricboardman/PhD Local Files/CODE/Outputs/2022-05-30/18MayDrying.pdf", IndexTimePlot
)
=#
