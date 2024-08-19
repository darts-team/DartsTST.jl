module Export_Output

using ArchGDAL
using GeoArrays

import ArchGDAL as AG

function export_tif_file(file_path, output_data, data_type, ag_geotransform, ag_ref)

    num_var = size(output_data)[3]

    AG.create(
        file_path,
        driver = AG.getdriver("GTiff"),
        width = size(output_data[:,:,1])[1],
        height = size(output_data[:,:,1])[2],
        nbands = num_var,
        dtype = data_type
    ) do raster_ds 
        for i=1:num_var
            AG.write!(raster_ds, output_data[:,:,i], i,)
        end
        AG.setgeotransform!(raster_ds, ag_geotransform)
        AG.setproj!(raster_ds, ag_ref)
    end

    return 1
end

end
