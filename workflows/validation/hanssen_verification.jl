using HypergeometricFunctions

gamma_ip        = 0.948
L               = 10
phase_ip_ref    = 0

#Cramer-Rao bound for phase variance
CR_var_ph = ((1-(gamma_ip^2)) / (2 * (gamma_ip^2) * L))


#Cramer-Rao bound for magnitude variance
CR_var_mag = ((1-(gamma_ip^2))^2) / (2 *  L)


#Bamler - phase
gamma_ip        = 0.1:0.1:0.9
phase_ip        = -3:0.1:3
pdf_phase       = zeros(length(gamma_ip), length(phase_ip))

for gamma_idx=1:length(gamma_ip)
    for phase_ip_idx=1:length(phase_ip)
        beta_ip     = gamma_ip[gamma_idx] * cos(phase_ip[phase_ip_idx])
        beta_sqrt_t = (1-(beta_ip^2))

        pdf_phase[gamma_idx,phase_ip_idx] = ((1 - gamma_ip[gamma_idx]^2)/ (2*pi)) .* (1/beta_sqrt_t) * (1 + ( (beta_ip *acos(-1*beta_ip)) / sqrt(beta_sqrt_t)) )  
    end
end

p1=plot()
for gamma_idx=1:length(gamma_ip)
    p1 = plot!(phase_ip,pdf_phase[gamma_idx,:],label="Gamma: "*string(gamma_ip[gamma_idx]), xlabel="Phase [rad]", ylabel="PDF(Phase)")
end
display(p1)


#Hanssen - phase
gamma_ip        = 0.1:0.1:0.9
phase_ip        = -3:0.1:3
pdf_phase       = zeros(length(gamma_ip), length(phase_ip))
L               = 10

for gamma_idx=1:length(gamma_ip)
    for phase_ip_idx=1:length(phase_ip)
        beta_ip     = gamma_ip[gamma_idx] * cos(phase_ip[phase_ip_idx] - phase_ip_ref)
        beta_sqrt_t = (1-(beta_ip^2))


        #pdf_ph_t1 = ((gamma(2*L-1)) / ((gamma(L)^2) * (2^(2*(L-1))))) * ( ((((2*L-1) * beta_ip) / (beta_sqrt_t^(L+0.5))) * ((pi/2) - asin(beta_ip))) + (1/beta_sqrt_t^L))
        pdf_ph_t1 = ((gamma(2*L-1)) / ((gamma(L)^2) * (2^(2*(L-1))))) * ( ((((2*L-1) * beta_ip) / (beta_sqrt_t^(L+0.5))) * ( acos(-beta_ip))) + (1/beta_sqrt_t^L))

        if L == 1
            pdf_ph_t2 = 0# (1 / (2*(L-1))) 
        else
            sum_term = 0
            for i=0:L-2
                sum_term = sum_term + ( ( (gamma(L-0.5)) * (gamma(L-i-1)) * (1+(((2*i)+1)*(beta_ip^2)))) / ( (gamma(L-0.5-i)) * (gamma(L-1)) * (beta_sqrt_t^(i+2))))
            end    
            pdf_ph_t2 = (1 / (2*(L-1))) * sum_term
        end
        pdf_phase[gamma_idx,phase_ip_idx] = (((1 - gamma_ip[gamma_idx]^2)^L)/ (2*pi)) .* (pdf_ph_t1 + pdf_ph_t2)

    end
end

p1=plot()
for gamma_idx=1:1:length(gamma_ip)
    p1 = plot!(phase_ip,pdf_phase[gamma_idx,:],label="Gamma: "*string(gamma_ip[gamma_idx]), xlabel="Phase [rad]", ylabel="PDF(Phase)")
end
display(p1)



#Hanssen - phase 2
gamma_ip        = 0.1:0.1:0.9
phase_ip        = -3:0.1:3
pdf_phase       = zeros(length(gamma_ip), length(phase_ip))
L               = 10

for gamma_idx=1:length(gamma_ip)
    for phase_ip_idx=1:length(phase_ip)
        beta_ip     = gamma_ip[gamma_idx] * cos(phase_ip[phase_ip_idx] - phase_ip_ref)
        beta_sqrt_t = (1-(beta_ip^2))

        Term1 = (gamma(L + 0.5) * ((1-(gamma_ip[gamma_idx]^2))^L) * abs.(gamma_ip[gamma_idx]) *  cos(phase_ip[phase_ip_idx] - phase_ip_ref)) / ( 2 * sqrt(pi) * gamma(L) * (beta_sqrt_t^(L+0.5)) )
        Term2 = ((1-(gamma_ip[gamma_idx]^2))^L) / (2*pi) * _₂F₁(L,1,0.5,beta_ip^2)

        pdf_phase[gamma_idx,phase_ip_idx] = Term1 + Term2

    end
end

p1=plot()
for gamma_idx=1:1:length(gamma_ip)
    p1 = plot!(phase_ip,pdf_phase[gamma_idx,:],label="Gamma: "*string(gamma_ip[gamma_idx]), xlabel="Phase [rad]", ylabel="PDF(Phase)")
end
display(p1)


#magnitude - Hanssen

gamma_ip_cap        = 0.1:0.001:1
gamma_ip        = 0.948
pdf_mag         = zeros(length(gamma_ip_cap))

for gamm_idx = 1:length(gamma_ip_cap)
    pdf_mag[gamm_idx]         = 2 .* (L-1) .* ((1-(gamma_ip^2))^L) .* gamma_ip_cap[gamm_idx] .* ((1-(gamma_ip_cap[gamm_idx]^2))^(L-2)) * _₂F₁(L,L,1,(gamma_ip_cap[gamm_idx]^2) * (gamma_ip^2) )
end

display(plot(gamma_ip_cap,pdf_mag,legend=false,xlabel="Correlation magnitude", ylabel="PDF(magnitude)"))



#magnitude - ??

Input_mag       = 0.3:0.01:1
pdf_mag         = ( 2 .* Input_mag.^(2*L-1)./ gamma(L)) .* ((L / gamma_ip)^L) .*  exp.(- (Input_mag.^2 .* L) ./ gamma_ip)

display(plot(Input_mag,pdf_mag,legend=false,xlabel="Correlation magnitude", ylabel="PDF(magnitude)"))



