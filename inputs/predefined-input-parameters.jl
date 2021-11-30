#using Parameters
#include("../modules/user_parameters.jl")
#using .UserParameters

customParams_test = UserParameters.inputParameters(
    look_angle=40, # change only look angle
)

customParams_profileTest = UserParameters.inputParameters(
    target_pos_mode="layered-grid", 
    t_loc_1 = [0.], # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2 = -2:1:2, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3 = -2:1:2, # m  heights if LLH or SCH, Z if XYZ
    t_ref   = [2, 1, 3], # reflectivities
    display_input_scene = true 
)

customParams_shapedProfileTest = UserParameters.inputParameters(
    target_pos_mode="shaped-grid", 
    t_loc_1 = [0.], # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2 = -2:1:2, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3 = -2:1:2, # m  heights if LLH or SCH, Z if XYZ
    t_ref   = [2, 1, 3], # reflectivities
    display_input_scene = true 
)

customParams_CR_nadirlooking = UserParameters.inputParameters(
    # From CR_nadirlooking.jl
    look_angle=0, # nadir looking
    s_loc_1=-130:1:130, # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2=-1000:10:1000, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3=0:1:80 # m  heights if LLH or SCH, Z if XYZ
)

customParams_CR_nadirlooking_tiltedcuts = UserParameters.inputParameters(
    # From CR_nadirlooking_tiltedcuts.jl    
    look_angle=0, # nadir looking
    PSF_direction=[1 7 0]
)

# add more custom struct here