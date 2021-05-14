using Profile
# using ProfileView
using Serialization

function sync_profile_test()
  include("main_simulation_sync_monte_carlo.jl")
end#function


using ProfileView
Profile.clear()
@profview sync_profile_test()

# Profile.clear()
# @profile sync_profile_test()
# Juno.profiler(; C = true)
# Profile.print()

# ProfileView.view()

r = Profile.retrieve();
f = open("profile.bin", "w")
serialize(f, r)
close(f)