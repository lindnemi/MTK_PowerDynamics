using ModelingToolkit

# outer droop control for virtual intertia
@parameters t K_P K_Q P_ref Q_ref ω_ref V_ref
@variables u_V u_ϕ(t) p_m(t) q_m(t)
@derivatives D'~t

virtual_intertia_droop = [
    u_ϕ ~ -K_P*(p_m-P_ref) + ω_ref, # output is the droop frequency ω
    u_V ~ -K_Q*(q_m-Q_ref) + V_ref # output is the droop voltage v
    ]

    
# component active power filter
@parameters t τ_P
@variables p(t) p_m(t)
@derivatives D'~t

filter_active_power = [
    D(p_m) ~ 1/τ_P* (-p_m-p),
    ]

# component reactive power filter
@parameters t τ_Q
@variables q(t) q_m(t)
@derivatives D'~t
    
filter_reactive_power = [
    D(q_m) ~ 1/τ_P* (-q_m-q),
    ]

# component integral controller for voltage
@parameters t τ_V
@variables v(t) u_V(t)
@derivatives D'~t
    
filter_reactive_power = [
    D(v) ~ 1/τ_V* (-v+u_V),
    ]

# integrator for frequency
@parameters t
@variables ϕ(t) u_ϕ(t)
@derivatives D'~t
    
filter_reactive_power = [
    D(ϕ) ~ u_ϕ
    ]