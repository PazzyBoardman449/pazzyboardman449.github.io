#Import Packages
using Plots
using DelimitedFiles
using LsqFit
using LaTeXStrings
using Statistics
using ColorSchemes
using SymPy
using NLsolve
#using MTH229
using DifferentialEquations
using DataFrames
using CSV
#using PlotlyJS, HTTP, DataFrames

#Universal Constants
const c = 299792458
const h = 6.62607004E-34
const k_B = 1.380649E-23
const e = 1.60217662E-19
const m_e = 9.10938356E-31
const ε_0 = 8.85418782E-12
const N_A = 6.0221409E+23

#Fixed Parameters

#Geometry
ang = 0 #Angle of Incidence
d_TPX = 3E-2 #Thickness of TPX prism (Always a constant)
d_si = 300E-6 #Constant thickness of Silicon (Used for phantom Thickness simulations)
#d_phan = 500E-6 #Constant thickness of Phantom (Used for silicon Thickness simulations)
d_air = 1 #Thickness of air (Note this can be any number, air is all around us!)

#Refractive Indices of Samples
n_TPX = 1.46 - 5.11E-5 * im #Refractive Index of TPX prism
n_si = 3.90 - 0.021 * im #Refractive Index of silicon at OPTICAL freqiencies
n_q = 1.95 - 0.2*im
n_air = 1 + 0*im

#Drude Model Parameters
λ = 623E-9
I₀ = 300
τ = 170E-6
μ_e = 0.145
μ_h = 0.045
m_e_eff = 0.26 * m_e
m_h_eff = 0.38 * m_e
ε_bg = 11.7 + 0.003im

#Debye Model Parameters
τ_D_1 = 8.210E-12 #Debye Relaxation Time of Process 1
τ_D_2 = 180E-15 #Debye Relaxation Time of Process 2
ε_gel = 3.77+ 0*im
ε_∞ = 3.5 #Assume this to be constant across all frequencies
ε_1 = 78.4 #Asume this to be constant across all freqiencies
ε_2 = 4.9 #Asume this to be constant across all freqiencies


#This fuction will return the modulated refractive index as a function of Light intensity
function refractive_index(In, f, d_si)

   #Calculate Transmission into the wafer
   r = ((n_TPX - n_si) / (n_TPX + n_si))
   R = conj(r) * r
   T = 1 - R #Valid since there is no losses

   #Calculate Generation Rate
   G = (In * T * λ) / (h * c * d_si) #Average generation rate per unit volume within the crystal

   #Calculate Carrier Density
   N = G * τ

   #Calculate Damping Coefficients
   γ_e = e / (μ_e * m_e_eff)
   γ_h = e / (μ_h * m_h_eff)

   #Calculate Electron Plasma Frequency
   f_p_e = (1 / 2π) * sqrt((N * e^2) / (m_e_eff * ε_0))
   f_p_h = (1 / 2π) * sqrt((N * e^2) / (m_h_eff * ε_0))

   #Calculate Permitivity
   ε =
      ε_bg - (2π * f_p_e)^2 / ((2π * f)^2 + γ_e * (2π * f)im) -
      (2π * f_p_h)^2 / ((2π * f)^2 + γ_h * (2π * f)im)

   #Square Root To Find Refractive Index, and take the CC to make extinction coefficient negative
   n = sqrt(ε)
   n = conj(n)

   #println("The complex refractive index of silicon is $(n)")

   return n
end

#This function will return the complex refractive index of the phantom, as a function of water concentration
function phant_refractive(freq, p_w)

   ε_water = ε_∞ + (ε_1 - ε_2) / (1 + 2 * pi * im * freq * τ_D_1) + (ε_2 - ε_∞) / (1 + 2 * pi * im * freq * τ_D_2)

   n_water = sqrt(ε_water)

   #println("The refractive index of the water is $(n_water)")

   ε_phant_WA = p_w * ε_water + (1 - p_w) * ε_gel #Simple Weighted Average Using Volume Concentrations

   #ε_phant_MG = ε_water*(1+3*(1-c)*((ε_gel-ε_water)/(ε_gel+2*ε_water-(1-c)*(ε_gel-ε_water))))

   #function MaxwellGarnet(f)
   #   ε_phant_MG = ε_water*(1+3*f*(ε_gel-ε_water)/(ε_gel+2*ε_water)+3*f^2*(ε_gel-ε_water)^2/(ε_gel+2*ε_water)^2)
   #   return ε_phant_MG
   #end

   function MaxwellGarnet(p_w)
      ε_phant_MG = ε_water*(ε_water + ((3-2*p_w)/3)*(ε_gel - ε_water))/(ε_water + (p_w/3)*(ε_gel - ε_water))
      return ε_phant_MG
   end

   function Bruggeman(p_w)
      ε_eff = (1/12)*(6 .*ε_water .- 3 .*ε_gel .+ 9 .*(ε_gel .- ε_water).*(1-p_w) .+ sqrt.(72 .* ε_water .* (ε_gel) +  (6 .*ε_water .- 3 .*ε_gel .+ 9 .*(ε_gel .-ε_water).*(1-p_w)).^2))
      return ε_eff
   end

   n_phant_WA = sqrt(ε_phant_WA)
   n_phant_MG = sqrt(MaxwellGarnet(p_w))
   n_phant_BR = sqrt(Bruggeman(p_w))

   println("Weighted Average Gives: $(n_phant_WA)")
   println("Maxwell Garnet Gives $(n_phant_MG)")
   println("Bruggeman gives $(n_phant_BR)")

   return n_phant_WA, n_phant_MG,n_phant_BR
   #return n_phant_MG
end

#This function will return the reflection coefficient from a block of layers.
function reflection(ang, d_ph, n_ph, d_sil, In, freq)

   n_sil = refractive_index(In, freq, d_sil)

   #println("STARTT!!!")

   N_Lyrs = 4
   θ = deg2rad(ang) #Convert Degrees To Radians

   #Depth Profiles
   d_prof =  [d_TPX, d_sil, d_ph, d_air]

   #Refractive Index Profile
   n_prof = [n_TPX, n_sil, n_ph, d_air]
   #println("The refractive index profile is $(n_prof)")
   #println("The depth profile is $(d_prof)")

   #Free Space Wave Vector
   k_0 = 2π * freq / c

   k_x = n_prof[1]*k_0*sin(θ)
   #println("The x component of the wave vector is $(k_x)")

   #z Component of the Wave Vector
   k_z = sqrt.(n_prof.^2*k_0.^2.0 .- k_x.^2)

   for i in 1:N_Lyrs
      if imag(k_z[i]) > 0 || real(k_z[i]) < 1
         k_z[i] = conj(k_z[i])
      else
      end
   end


   #println("The z component of the wave vector is $(k_z)")

   #Delta Coefficients

   δ = k_z .* d_prof
   #println("The δ coefficients are $(δ)")

   #Gamma Coefficients
   γ = fill(0.0 + 0.0im, N_Lyrs, 1)
   γ_TM = fill(0.0 + 0.0im, N_Lyrs, 1)


   for M = 1:N_Lyrs
      γ[M] = k_z[M] #sqrt(n_prof[M]^2 - n_prof[1]^2 * sin(θ)^2)
      γ_TM[M] = k_z[M]/(n_prof[M]^2) #sqrt(n_prof[M]^2 - n_prof[1]^2 * sin(θ)^2)./(n_prof[M].^2)
   end


   #println("The γ coefficients are $(γ)")

   #Fill Up Transfer Matrices With Terms
   #TE:
   Term1 = cos.(δ)
   Term2 = (im * sin.(δ)) ./ γ
   Term3 = (im * sin.(δ)) .* γ
   Term4 = cos.(δ)

   #TM:
   Term2TM = (im * sin.(δ)) ./ γ_TM
   Term3TM = (im * sin.(δ)) .* γ_TM

   # First Transfer matrix
   M = [Term1[1] Term2[1]; Term3[1] Term4[1]]
   M_TM = [Term1[1] Term2TM[1]; Term3TM[1] Term4[1]]

   #println("The transfer matrix is initially $(M)")

   #THE Transfer Matrix
   #M = [cos.(δ)[2] (im * sin.(δ)[2]) ./ γ[2]; (im * sin.(δ)[2]).*γ[2] cos.(δ)[2]]

   # Iterate through layers and sum transfer matrices for each layer
   for i = 2:N_Lyrs-1
      M = M * [Term1[i] Term2[i]; Term3[i] Term4[i]]
      M_TM = M_TM * [Term1[i] Term2TM[i]; Term3TM[i] Term4[i]]
      #println("The transfer matrix is then $(M)")
   end


   #Calculate reflection amplitude Coefficient
   r =
      (γ[1] * M[1, 1] + γ[1] * γ[N_Lyrs] * M[1, 2] - M[2, 1] - γ[N_Lyrs] * M[2, 2]) /
      (γ[1] * M[1, 1] + γ[1] * γ[N_Lyrs] * M[1, 2] + M[2, 1] + γ[N_Lyrs] * M[2, 2])

      #println("The reflection amplitude coefficient is $(r)")

   t =
      2 * γ[1] / (γ[1] * M[1, 1] + γ[1] * γ[N_Lyrs] * M[1, 2] + M[2, 1] + γ[N_Lyrs] * M[2, 2])

   r_TM = (γ_TM[1] * M_TM[1, 1] + γ_TM[1] * γ_TM[N_Lyrs] * M_TM[1, 2] - M_TM[2, 1] - γ_TM[N_Lyrs] * M_TM[2, 2]) /
         (γ_TM[1] * M_TM[1, 1] + γ_TM[1] * γ_TM[N_Lyrs] * M_TM[1, 2] + M_TM[2, 1] + γ_TM[N_Lyrs] * M_TM[2, 2])
   t_TM = 2 * γ_TM[1] / (γ_TM[1] * M_TM[1, 1] + γ_TM[1] * γ_TM[N_Lyrs] * M_TM[1, 2] + M_TM[2, 1] + γ_TM[N_Lyrs] * M_TM[2, 2])

   #if isnan(real(r)) || real(r) > 1.0  || isnan(imag(r)) || imag(r) > 1.0
   #   r = 1 +0.0*im
   #else
   #end

   #if isnan(real(t)) || real(t) > 1.0  || isnan(imag(t)) || imag(t) > 1.0
   #   t = 1.0 +0.0*im
   #else
   #end

   R = abs(r) .^ 2
   T = abs(t) .^ 2 * real(k_z[N_Lyrs] / k_z[1])

   R_TM = abs(r_TM)^ 2
   T_TM = abs(t_TM)^ 2 * real(k_z[N_Lyrs] / k_z[1])

   #if R > 1.0
   #   R = 1
   #else
   #end

   #if isnan(R)
   #   R = 1
   #else
   #end

   #if T > 1.0
   #   T = 1
   #else
   #end

   return R
end

function Ian(eps_z, dz, freq, angle)

   N = length(eps_z)

   k0 = 2 * pi * freq / 2.998e8 # Free space wavevector

   n = real(sqrt.(eps_z)) - 1im * (imag(sqrt.(eps_z))) # Because, annoyingly, the T-Matrix method below work with negative imaginary n for loss.

   kx = n[1] * k0 * sin(angle) # in-plane wavevector

   kz = sqrt.(n .^ 2 * k0 .^ 2.0 .- kx .^ 2) # out-of-plane wavevector in each jth layer

   # Parameters for the transfer matrixes for each jth layer

   beta = kz .* dz

   gamma_TM = kz ./ n .^ 2
   gamma_TE = kz

   term_1 = cos.(beta)
   term_2 = 1im .* sin.(beta) ./ gamma_TM
   term_3 = 1im .* sin.(beta) .* gamma_TM
   term_4 = 1im .* sin.(beta) ./ gamma_TE
   term_5 = 1im .* sin.(beta) .* gamma_TE

   # First Transfer matrix
   T_TM = [term_1[1] term_2[1]; term_3[1] term_1[1]]
   T_TE = [term_1[1] term_4[1]; term_5[1] term_1[1]]

   # Iterate through layers and sum transfer matrices for each layer
   for i = 2:N
      T_TM = T_TM * [term_1[i] term_2[i]; term_3[i] term_1[i]]
      T_TE = T_TE * [term_1[i] term_4[i]; term_5[i] term_1[i]]

   end

   # Calulate complex reflection and transmission amplitude coefficients
   r_TM =
      (
         gamma_TM[1] * T_TM[1, 1] + gamma_TM[1] * gamma_TM[N] * T_TM[1, 2] -
         T_TM[2, 1] - gamma_TM[N] * T_TM[2, 2]
      ) / (
         gamma_TM[1] * T_TM[1, 1] +
         gamma_TM[1] * gamma_TM[N] * T_TM[1, 2] +
         T_TM[2, 1] +
         gamma_TM[N] * T_TM[2, 2]
      )
   t_TM =
      2 * gamma_TM[1] / (
         gamma_TM[1] * T_TM[1, 1] +
         gamma_TM[1] * gamma_TM[N] * T_TM[1, 2] +
         T_TM[2, 1] +
         gamma_TM[N] * T_TM[2, 2]
      )

   r_TE =
      (
         gamma_TE[1] * T_TE[1, 1] + gamma_TE[1] * gamma_TE[N] * T_TE[1, 2] -
         T_TE[2, 1] - gamma_TE[N] * T_TE[2, 2]
      ) / (
         gamma_TE[1] * T_TE[1, 1] +
         gamma_TE[1] * gamma_TE[N] * T_TE[1, 2] +
         T_TE[2, 1] +
         gamma_TE[N] * T_TE[2, 2]
      )
   t_TE =
      2 * gamma_TE[1] / (
         gamma_TE[1] * T_TE[1, 1] +
         gamma_TE[1] * gamma_TE[N] * T_TE[1, 2] +
         T_TE[2, 1] +
         gamma_TE[N] * T_TE[2, 2]
      )

   # Calulate reflected and transmitted intensities
   Ref_TM = abs(r_TM)^2
   Trans_TM = abs(t_TM)^2 * real(kz[N] / kz[1]) .* (ε[1] ./ ε[N])

   Ref_TE = abs(r_TE)^2
   Trans_TE = abs(t_TE)^2 * real(kz[N] / kz[1])

   # Return reflected and transmitted intensities

   #return Ref_TM,Trans_TM,Ref_TE,Trans_TE

   return Ref_TE
end

#Plots

#Drude Model Plotter
function DrudePlotter()
   freq = -5E11:1E8:5E11
   y = refractive_index.(1000, freq, d_si)

   plot(freq / 1E12, real(y), lw = 3, label = latexstring("\$ n \$"))
   plot!(freq / 1E12, -imag(y), lw = 3, label = latexstring("\$ \\kappa \$"))
   plot!(
      xlims = (0, 0.5),
      ylims = (0, 10),
      title = "Drude Model THz Regime Frequency Dependence",
      xlabel = L"f" * " THz",
      ylabel = L"  \bar{n}_{Si}(f) ",
      titlefont = font(10, "Computer Modern"),
      legendfont = font(10, "Computer Modern"),
      tickfont = font(10, "Computer Modern"),
      guidefont = font(10, "Computer Modern"),
   )
   plot!(aspect_ratio = 0.05)
   Plots.savefig(
      "/Users/patricboardman/PhD Local Files/NEW PHD CODE/Figures/Drude .pdf",
   )
end

#Debye Model Plotter
function DebyePlotter()
   freq = 0:1E9:2E12
   y = phant_refractive.(freq, 1)

   plot(freq/1E12, real(y), lw = 3, label = latexstring("\$ n \$"))
   plot!(freq/1E12 , -imag(y), lw = 3, label = latexstring("\$ \\kappa \$"))
   plot!(
      xlims = (0, 1),
      ylims = (0, 10),
      title = "Two Term Debye Model THz Regime Frequency Dependence",
      xlabel = L"f" * " THz",
      ylabel = L"  \bar{n}_{Water}(f) ",
      titlefont = font(10, "Computer Modern"),
      legendfont = font(10, "Computer Modern"),
      tickfont = font(10, "Computer Modern"),
      guidefont = font(10, "Computer Modern"),
   )
   plot!(aspect_ratio = 0.1)


   Plots.savefig(
      "/Users/patricboardman/PhD Local Files/NEW PHD CODE/Code Outputs/28Feb.pdf",
   )
end

#Effective Medium Theory Plotters
function EffectiveMediumPlotter()

   N_steps = 200
   conc_arr = 0:1/N_steps:1-1/N_steps

   n_WA = zeros(ComplexF64, N_steps, 1)
   n_MG = zeros(ComplexF64, N_steps, 1)
   n_BR = zeros(ComplexF64, N_steps, 1)

   for i in 1:N_steps
      n_WA[i] = phant_refractive(140E9,conc_arr[i])[1]
      n_MG[i] = phant_refractive(140E9,conc_arr[i])[2]
      n_BR[i] = phant_refractive(140E9,conc_arr[i])[3]
   end

   p1 = plot(conc_arr*100,real(n_WA),label = L"n_{WA}")
   plot!(conc_arr*100,real(n_MG),label = L"n_{MG}")
   plot!(conc_arr*100,real(n_BR),label = L"n_{BR}")
   xaxis!(L"p_{w}")
   yaxis!(L"n(p_w)")

   p2 = plot(conc_arr*100,-imag(n_WA),label = L"\kappa_{WA}")
   plot!(conc_arr*100,-imag(n_MG),label = L"\kappa_{MG}")
   plot!(conc_arr*100,-imag(n_BR),label = L"\kappa_{BR}")
   xaxis!(L"p_{w}")
   yaxis!(L"\kappa(p_w)")

   plot(p1, p2, layout=(1,2),
   xlims = (0,100),
   ylims = (0.5,3),
   titlefont = font(10, "Computer Modern"),
   legendfont = font(10, "Computer Modern"),
   tickfont = font(10, "Computer Modern"),
   guidefont = font(10, "Computer Modern"),
   aspect_ratio = 40
   )

   #=
   plot(conc_arr,real(n_WA), color = :blue,lw = 1, label = L"n_{WA}",xlims = (0, 1),
   ylims = (0, 4),
   title = "Effective Medium Models for " *L"f = 140GHz",
   xlabel = "Percentage concentration of water "*L"p_{w}",
   ylabel = "Real refractive index "*L"  \bar{n}_{Phan}(p_w) ",
   titlefont = font(10, "Computer Modern"),
   legendfont = font(10, "Computer Modern"),
   tickfont = font(10, "Computer Modern"),
   guidefont = font(10, "Computer Modern"),layout = 2)

   plot!(conc_arr,real(n_MG),lw = 1, label = L"n_{MG}",layout = 2)
   plot!(conc_arr,real(n_BR),lw = 1, label =  L"n_{BR}")
   #plot!(conc_arr,-imag(n_BR),lw = 1, label =  L"\kappa_{BR}")
   plot!(aspect_ratio = .25)
   plot!(conc_arr,-imag(n_WA),color = :blue,lw = 1, label = L"n_{WA}")

   =#

   #Plots.savefig(
   #   "/Users/patricboardman/PhD Local Files/NEW PHD CODE/Code Outputs/07Mar/EffectiveMediumModels.pdf",
   # )
end

#3D Plots
function ReflectionThicknessSimulations()
   N_d = 100
   N_ang = 90

   d_si_arr = 0+(1E-3/N_d):(1E-3/N_d):(1E-3)
   ang_arr = 0:(90/N_ang):(90-90/N_ang)


   #Initialise Reflection Coefficient Array
   reflection_mod_arr = zeros(Float64, N_ang, N_d)
   reflection_unmod_arr = zeros(Float64, N_ang, N_d)
   reflection_phantmod_arr = zeros(Float64, N_ang, N_d)
   reflection_phantunmod_arr = zeros(Float64, N_ang, N_d)

   for i = 1:N_ang
      for j = 1:N_d
         reflection_mod_arr[i, j] =
            reflection.(ang_arr[i], d_air, 1, d_si_arr[j], I₀, 140E9)
         reflection_unmod_arr[i, j] =
            reflection.(ang_arr[i], d_air, 1, d_si_arr[j], 0, 140E9)
         reflection_phantmod_arr[i, j] =
            reflection.(ang_arr[i], d_air, 2.5224216378863646 - 1.111347491398024*im, d_si_arr[j], I₀, 140E9)
         reflection_phantunmod_arr[i, j] =
            reflection.(ang_arr[i], d_air, 2.5224216378863646 - 1.111347491398024*im, d_si_arr[j], 0, 140E9)
      end
   end

   Plots.heatmap(d_si_arr*1E3, ang_arr, reflection_mod_arr, fill = true)
   Plots.heatmap!(
      xlims = (0, 1),
      ylims = (0, 90),
      zlims = (0, 1),
      title = L"R" * " Unmodulated, No Phantom "*" (300 W"*L"m^{-2})",
      xlabel = L"d_{Si}" * " mm",
      ylabel = L"\theta" * " Degrees",
      zlabel = L"R",
      color_limits=(0,1),
      titlefont = font(10, "Computer Modern"),
      legendfont = font(10, "Computer Modern"),
      tickfont = font(10, "Computer Modern"),
      guidefont = font(10, "Computer Modern"))


   #Plots.savefig(
   #   "/Users/patricboardman/PhD Local Files/NEW PHD CODE/Figures/RCoeffUnmodNoPhantom.pdf",
   #)


   #plot(ang_arr,d_si_arr,reflection_arr)

   #Plotting
   #=
   plot(d_si_arr*10E6,real(n_si_arr),lw=3,label = latexstring("\$ n \$"))
   plot!(d_si_arr*10E6,imag(n_si_arr),lw=3,label = L"\kappa")
   plot!(xlims=(0,1000),ylims=(0,10),
   title = "Complex refractive index Si Wafer",
   xlabel=L"d/microns", ylabel=L"\bar{n} = n + i\kappa",
   titlefont=font(10,"Times"),
   legendfont=font(10,"Times"),
   xtickfont=font(10,"Times"),
   ytickfont=font(10,"Times" ))
   plot!(aspect_ratio=100)
   #return si_refractive_index_arr
   =#
end

function ModulationDepthSiliconThicknessSumulations()
   N_d = 100 #500
   N_ang = 90 #450

   d_si_arr = 0+(1E-3/N_d):(1E-3/N_d):(1E-3)
   ang_arr = 0:(90/N_ang):(90-90/N_ang)

   #Initialise Reflection Coefficient Array
   reflection_mod_arr = zeros(Float64, N_ang, N_d)
   reflection_unmod_arr = zeros(Float64, N_ang, N_d)
   reflection_phantmod_arr = zeros(Float64, N_ang, N_d)
   reflection_phantunmod_arr = zeros(Float64, N_ang, N_d)
   modulation_depth_arr = zeros(Float64, N_ang, N_d)

   for i = 1:N_ang
      for j = 1:N_d
         reflection_mod_arr[i, j] =
            reflection.(ang_arr[i], d_air, n_air, d_si_arr[j], I₀, 140E9)
         reflection_unmod_arr[i, j] =
            reflection.(ang_arr[i], d_air, n_air, d_si_arr[j], 0, 140E9)

         modulation_depth_arr[i, j] = (reflection_unmod_arr[i, j] - reflection_mod_arr[i, j])
      end
   end

   #modulation_depth_arr = (reflection_phantunmod_arr-reflection_phantmod_arr) #./ (reflection_unmod_arr-reflection_mod_arr)

   println(extrema(modulation_depth_arr))

   Plots.heatmap(d_si_arr*1E6, ang_arr, modulation_depth_arr, fill = true)
   Plots.heatmap!(
      xlims = (0, 1000),
      ylims = (0, 90),
      zlims = (0, 1),
      title = "TM Modulation Depth "*L"R_{Unmod}-R_{Mod}",
      xlabel = L"d_{Si}" * " µm",
      ylabel = L"\theta" * " Degrees",
      zlabel = L"R",
      color_limits=(0,1),
      titlefont = font(10, "Computer Modern"),
      legendfont = font(10, "Computer Modern"),
      tickfont = font(10, "Computer Modern"),
      guidefont = font(10, "Computer Modern"))
      #plot!(aspect_ratio = 1000/90)

   #Plots.savefig(
   #   "/Users/patricboardman/PhD Local Files/NEW PHD CODE/Code Outputs/07Mar/ModulationDepthTM.pdf",
   #)


   #plot(ang_arr,d_si_arr,reflection_arr)

   #Plotting
   #=
   plot(d_si_arr*10E6,real(n_si_arr),lw=3,label = latexstring("\$ n \$"))
   plot!(d_si_arr*10E6,imag(n_si_arr),lw=3,label = L"\kappa")
   plot!(xlims=(0,1000),ylims= (0,10),
   title = "Complex refractive index Si Wafer",
   xlabel=L"d/microns", ylabel=L"\bar{n} = n + i\kappa",
   titlefont=font(10,"Times"),
   legendfont=font(10,"Times"),
   xtickfont=font(10,"Times"),
   ytickfont=font(10,"Times" ))
   plot!(aspect_ratio=100)
   #return si_refractive_index_arr
   =#
end

function ModulationDepthPhantThicknessSumulations()
   N_d = 100#500
   N_ang = 90#450

   d_ph_arr = 0+(1E-3/N_d):(1E-3/N_d):(1E-3)
   ang_arr = 0:(90/N_ang):(90-90/N_ang)

   n_phan = phant_refractive(140E9, .90)[3]

   reflection_mod_arr = zeros(Float64, N_ang, N_d)
   reflection_unmod_arr = zeros(Float64, N_ang, N_d)
   reflection_phantmod_arr = zeros(Float64, N_ang, N_d)
   reflection_phantunmod_arr = zeros(Float64, N_ang, N_d)
   modulation_depth_arr = zeros(Float64, N_ang, N_d)

   for i = 1:N_ang
      for j = 1:N_d
         reflection_mod_arr[i, j] =
            reflection.(ang_arr[i], d_ph_arr[j], 1, d_si, I₀, 140E9)
         reflection_unmod_arr[i, j] =
            reflection.(ang_arr[i], d_ph_arr[j], 1, d_si, 0, 140E9)
         reflection_phantmod_arr[i, j] =
            reflection.(ang_arr[i], d_ph_arr[j], n_phan, d_si, I₀, 140E9)
         reflection_phantunmod_arr[i, j] =
            reflection.(ang_arr[i], d_ph_arr[j], n_phan, d_si, 0, 140E9)

         modulation_depth_arr[i, j] = (reflection_phantunmod_arr[i, j] - reflection_phantmod_arr[i, j])./((reflection_unmod_arr[i, j] - reflection_mod_arr[i, j]))

         #if modulation_depth_arr[i, j] < 0
         #   modulation_depth_arr[i, j] = 0
         #elseif modulation_depth_arr[i, j] > 1
         #   modulation_depth_arr[i, j] = 1
         #else
         #end

      end
   end

   Plots.heatmap(d_ph_arr*1E6, ang_arr, modulation_depth_arr, fill = true)
   Plots.heatmap!(
      xlims = (0, 1000),
      ylims = (0, 90),
      zlims = (0, 1),
      title = "90% Water Phantom "*L"(R'_{Unmod}-R'_{Mod})/(R_{Unmod}-R_{Mod})",
      xlabel = "Phantom Thickness "*L"d_{Ph}" * " µm",
      ylabel = L"\theta" * " Degrees",
      zlabel = L"R",
      color_limits=(0,1),
      titlefont = font(10, "Computer Modern"),
      legendfont = font(10, "Computer Modern"),
      tickfont = font(10, "Computer Modern"),
      guidefont = font(10, "Computer Modern"))
      #plot!(aspect_ratio = 1000/90)


   #Plots.savefig(
   #   "/Users/patricboardman/PhD Local Files/NEW PHD CODE/Code Outputs/07Mar/NormalizedPhantom90",
   #)
end

function ModulationDepthSiliconFreqSumulations()
   N_freq = 100 #500
   N_ang = 90 #450

   freq_arr = 0+(1E12/N_freq):(1E12/N_freq):(1E12)
   ang_arr = 0:(90/N_ang):(90-90/N_ang)

   #Initialise Reflection Coefficient Array
   reflection_mod_arr = zeros(Float64, N_ang, N_freq)
   reflection_unmod_arr = zeros(Float64, N_ang, N_freq)
   reflection_phantmod_arr = zeros(Float64, N_ang, N_freq)
   reflection_phantunmod_arr = zeros(Float64, N_ang, N_freq)
   modulation_depth_arr = zeros(Float64, N_ang, N_freq)


   for i = 1:N_ang
      for j = 1:N_freq
         reflection_mod_arr[i, j] =
            reflection.(ang_arr[i], d_phan, 1, d_si, I₀, freq_arr[j])
         reflection_unmod_arr[i, j] =
            reflection.(ang_arr[i], d_phan, 1, d_si, 0, freq_arr[j])

         modulation_depth_arr[i, j] = (reflection_unmod_arr[i, j] - reflection_mod_arr[i, j])
         #if modulation_depth_arr[i, j] < 0
         #   modulation_depth_arr[i, j] = 0
         #elseif modulation_depth_arr[i, j] > 1
         #   modulation_depth_arr[i, j] = 1
         #else
         #end

      end
   end

   #modulation_depth_arr = (reflection_phantunmod_arr-reflection_phantmod_arr) #./ (reflection_unmod_arr-reflection_mod_arr)
   println(extrema(modulation_depth_arr))

   Plots.heatmap(freq_arr/1E12, ang_arr, modulation_depth_arr, fill = true)
   Plots.heatmap!(
      xlims = (0, 1),
      ylims = (0, 90),
      zlims = (0, 1),
      title = "Modulation Depth "*L"R_{Unmod}-R_{Mod}",
      xlabel = L"f" * " THz",
      ylabel = L"\theta" * " Degrees",
      zlabel = L"R",
      color_limits=(0,1),
      titlefont = font(10, "Computer Modern"),
      legendfont = font(10, "Computer Modern"),
      tickfont = font(10, "Computer Modern"),
      guidefont = font(10, "Computer Modern"))
      plot!(aspect_ratio = 1/90)

   #Plots.savefig(
   #   "/Users/patricboardman/PhD Local Files/NEW PHD CODE/Code Outputs/28Feb/FreqModulationDepth.pdf",
   #)


   #plot(ang_arr,d_si_arr,reflection_arr)

   #Plotting
   #=
   plot(d_si_arr*10E6,real(n_si_arr),lw=3,label = latexstring("\$ n \$"))
   plot!(d_si_arr*10E6,imag(n_si_arr),lw=3,label = L"\kappa")
   plot!(xlims=(0,1000),ylims= (0,10),
   title = "Complex refractive index Si Wafer",
   xlabel=L"d/microns", ylabel=L"\bar{n} = n + i\kappa",
   titlefont=font(10,"Times"),
   legendfont=font(10,"Times"),
   xtickfont=font(10,"Times"),
   ytickfont=font(10,"Times" ))
   plot!(aspect_ratio=100)
   #return si_refractive_index_arr
   =#
end

ModulationDepthPhantThicknessSumulations()
