module Scene

function form3Dgrid_for(t_θ,t_ϕ,t_h) # for loop method
  t_geo_grid=zeros(3,length(t_θ)*length(t_ϕ)*length(t_h))
  m=0
  for i=1:length(t_θ)
    for j=1:length(t_ϕ)
      for k=1:length(t_h)
        m=m+1
        t_geo_grid[:,m]=[t_θ[i],t_ϕ[j],t_h[k]]
      end
    end
  end
  return t_geo_grid
end

function form3Dgrid_array(t_θ,t_ϕ,t_h) # array method
  t_θ1=Array{Float64}(undef,1,length(t_θ))
  t_ϕ1=Array{Float64}(undef,1,length(t_ϕ))
  t_h1=Array{Float64}(undef,1,length(t_h))
  t_θ1[:]=t_θ
  t_ϕ1[:]=t_ϕ
  t_h1[:]=t_h
  t_θ_all=repeat(t_θ1,inner=[1,1],outer=[1,length(t_ϕ)*length(t_h)])
  t_ϕ_all=repeat(t_ϕ1,inner=[1,length(t_θ)],outer=[1,length(t_h)])
  t_h_all=repeat(t_h1,inner=[1,length(t_θ)*length(t_ϕ)],outer=[1,1])
  t_geo_grid=[t_θ_all;t_ϕ_all;t_h_all]
  return t_geo_grid
end

end
