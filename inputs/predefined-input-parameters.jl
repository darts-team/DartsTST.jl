#using Parameters
#include("../modules/user_parameters.jl")
#using .UserParameters

customParams_LookAngle = UserParameters.inputParameters(
    look_angle = 40 # change only look angle
)

customParams_resdB = UserParameters.inputParameters(
    res_dB = 3.4
)

customParams_PrincipalCuts = UserParameters.inputParameters(
    PSF_cuts = 1 # 1: principal axes (SCH, LLH, XYZ based on ts_coord_sys)
)

customParams_AntennaPatternTest = UserParameters.inputParameters(
    look_angle = 0, # nadir looking
    ts_coord_sys = "SCH",
    t_loc_1 = -20000:1000:20000, # along-track
    t_loc_2 = -100000:5000:100000, # cross-track
    t_loc_3 = 0, # height
    s_loc_1 = t_loc_1,
    s_loc_2 = t_loc_2,
    s_loc_3 = t_loc_3,
    include_antenna=true
)

#=customParams_TreeCH = UserParameters.inputParameters(
    target_pos_mode="grid", #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
    ts_coord_sys="SCH", # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene,  # if SCH, target and scene locations are defined relative to the point where look angle vector intersects the surface
    t_loc_1=0, # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=-12:2:12, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=4:1:15, # m  heights if LLH or SCH, Z if XYZ
    t_ref=zeros(Float64,length(t_loc_1),length(t_loc_2),length(t_loc_3)), # uniform random reflectivities between 0 and 1, a 3D input array (e.g. 3D image) can be used instead
    t_ref_2D=[0 0 0 0 0 0 3 0 0 0 0 0 0; #12 rows, 13 columns
              0 0 0 0 0 3 2 2 0 0 0 0 0;
              0 0 0 0 3 2 2 1 2 0 0 0 0;
              0 0 0 0 3 2 1 1 1 2 0 0 0;
              0 0 0 0 0 2 1 1 1 0 0 0 0;
              0 0 0 0 0 2 1 1 0 0 0 0 0;
              0 0 0 0 0 0 2 1 0 0 0 0 0;
              0 0 0 0 0 0 2 1 0 0 0 0 0;
              0 0 0 0 0 0 2 0 0 0 0 0 0;
              0 0 0 0 0 0 2 0 0 0 0 0 0;
              0 0 0 0 0 2 2 1 0 0 0 0 0;
              4 4 4 4 4 3 3 3 2 2 2 2 2],
    t_ref[1,:,:]=t_ref_2D', #x-axis of table (columns) = 2nd dimension, y-axis of table (rows) = 3rd dimension, scene exists only at 1st point in 1st dimension
    s_loc_1=0, # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2=-20:1:20, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3=0:0.5:25 # m  heights if LLH or SCH, Z if XYZ
)

customParams_TreeSH = UserParameters.inputParameters(
    target_pos_mode="grid", #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
    ts_coord_sys="SCH", # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene,  # if SCH, target and scene locations are defined relative to the point where look angle vector intersects the surface
    t_loc_1=-12:2:12, # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2=0, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=4:1:15, # m  heights if LLH or SCH, Z if XYZ
    t_ref=zeros(Float64,length(t_loc_1),length(t_loc_2),length(t_loc_3)), # uniform random reflectivities between 0 and 1, a 3D input array (e.g. 3D image) can be used instead
    t_ref_2D=[0 0 0 0 0 0 3 0 0 0 0 0 0; #12 rows, 13 columns
              0 0 0 0 0 3 2 2 0 0 0 0 0;
              0 0 0 0 3 2 2 1 2 0 0 0 0;
              0 0 0 0 3 2 1 1 1 2 0 0 0;
              0 0 0 0 0 2 1 1 1 0 0 0 0;
              0 0 0 0 0 2 1 1 0 0 0 0 0;
              0 0 0 0 0 0 2 1 0 0 0 0 0;
              0 0 0 0 0 0 2 1 0 0 0 0 0;
              0 0 0 0 0 0 2 0 0 0 0 0 0;
              0 0 0 0 0 0 2 0 0 0 0 0 0;
              0 0 0 0 0 2 2 1 0 0 0 0 0;
              4 4 4 4 4 3 3 3 2 2 2 2 2],
    t_ref[:,1,:]=reverse(t_ref_2D',dims=2), #x-axis of table (columns) = 1st dimension, y-axis of table (rows) = 3rd dimension, scene exists only at 1st point in 2nd dimension
    s_loc_1=-20:1:20, # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2=0, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3=0:0.5:25 # m  heights if LLH or SCH, Z if XYZ
)

customParams_Scattered9TargetsCH = UserParameters.inputParameters(
    target_pos_mode="CR", #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
    ts_coord_sys="SCH", # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene
    t_loc_2=[0 55 -55 40 -40 -30  30  0   0], # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=[0 50 -50 10 -10  35 -35 30 -30], # m  heights if LLH or SCH, Z if XYZ
    t_loc_1=zeros(1,length(t_loc_2)), # deg latitude if LLH, along-track if SCH, X if XYZ
    t_ref=  ones(length(t_loc_1)), # reflectivities
    s_loc_1=0, # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_2=-65:1:65, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3=-60:1:60, # m  heights if LLH or SCH, Z if XYZ
    distance_ratio=1,
    t_loc_2=t_loc_2/distance_ratio,
    t_loc_3=t_loc_3/distance_ratio
)

customParams_Scattered9TargetsSH = UserParameters.inputParameters(
    target_pos_mode="CR", #  targets are defined as three 1D arrays forming either a volumetric grid ("grid") or a 3xN array ("CR" for corner reflectors)
    ts_coord_sys="SCH", # target/scene coordinate system: "LLH", "SCH", "XYZ", using the same coordinate system for targets and scene
    t_loc_1=[0 55 -55 40 -40 -30  30  0   0], # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3=[0 50 -50 10 -10  35 -35 30 -30], # m  heights if LLH or SCH, Z if XYZ
    t_loc_2=zeros(1,length(t_loc_1)), # deg latitude if LLH, along-track if SCH, X if XYZ
    t_ref=  ones(length(t_loc_2)), # reflectivities
    s_loc_2=0, # deg latitude if LLH, along-track if SCH, X if XYZ
    s_loc_1=-65:1:65, # deg longitude if LLH, cross-track if SCH, Y if XYZ
    s_loc_3=-60:1:60, # m  heights if LLH or SCH, Z if XYZ
    distance_ratio=1,
    t_loc_1=t_loc_1/distance_ratio,
    t_loc_3=t_loc_3/distance_ratio
)=#

customParams_OrbitFromFile = UserParameters.inputParameters(
    user_defined_orbit==0, # from file
    p_t0_LLH::Array{Float64,1} = [0;0;750e3], # initial lat/lon (deg) and altitude (m) of reference platform TODO altitude is used in incidence angle calculation which is needed for tilted PSF cut direction, use avg altitude from orbit section used
    #orbit_filename="orbitOutput_082020.nc" # position in km, time in sec
    orbit_filename="orbit_output_062021.nc" # position in km, time in sec
)

customParams_CustomOrbit_AlongnSCH = UserParameters.inputParameters(
    user_defined_orbit==1, # along-n using SCH option
    p_t0_LLH::Array{Float64,1} = [0;0;750e3], # initial lat/lon (deg) and altitude (m) of reference platform (altitude is assumed constant over slow-time if SCH option)
    Torbit::Float64    = 10*60, # orbital duration (s) (should be larger than 2 x (SAR_start_time+SAR_duration) )
    dt_orbits::Float64 = 0.5, # orbit time resolution (s)
    p_heading::Float64 = 0, # heading (deg), all platforms assumed to have the same heading, 0 deg is north
    pos_n   = [-9 -6 -3 0 3 6 9]*1e3, # relative position of each platform along n (m), 0 is the reference location, equal spacing
    #pos_n=[-7.5 -5 -2 0 3.7 5.5 6.5]*1e3, # relative position of each platform along n (m), 0 is the reference location, unequal spacing
    display_custom_orbit=true #whether to show orbit on Earth sphere (for a duration of Torbit)
)

customParams_CustomOrbit_AlongTCN = UserParameters.inputParameters(
    user_defined_orbit==2, # TCN option
    p_t0_LLH::Array{Float64,1} = [0;0;750e3], # initial lat/lon (deg) and altitude (m) of reference platform (altitude is assumed constant over slow-time if SCH option)
    Torbit::Float64    = 10*60, # orbital duration (s) (should be larger than 2 x (SAR_start_time+SAR_duration) )
    dt_orbits::Float64 = 0.5, # orbit time resolution (s)
    p_heading::Float64 = 0, # heading (deg), all platforms assumed to have the same heading, 0 deg is north
    pos_TCN = [0 -6 0; 0 -5 0; 0 -2 0; 0 0 0; 0 3.5 0; 0 5 0]*1e3,   # TCN option: Np x 3 matrix; each row is the TCN coordinate of each platform relative to reference
    display_custom_orbit=true #whether to show orbit on Earth sphere (for a duration of Torbit)
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

customParams_multiplePTs = UserParameters.inputParameters(
    t_loc_1 = [0., 2, 4], # deg latitude if LLH, along-track if SCH, X if XYZ
    t_loc_2 = [0., 2, 1], # deg longitude if LLH, cross-track if SCH, Y if XYZ
    t_loc_3 = [0., 0, 0], # m heights if LLH or SCH, Z if XYZ
    t_ref   = [1., 1, 2] # reflectivities: a list of CRs in CR mode; an arbitrary profile that will be interpolated on t_loc_3 axis in *grid modes
)

# add more custom struct here
