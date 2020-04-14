# darts-simtool

To generate tomographic image (in 2D):
- Set platform, radar, and target parameters in generate_raw_data.m
- Set target scene parameters in generate_image.m
- Run generate_raw_data.m which generates raw data and saves it as raw.mat
- Run generate_image.m which reads raw.mat and generates focused 2D image and saves it as image.jpg

Run theoretical_resolution_ambiguity.m to:
- plot in 1D theoretical resolution vs aperture length
- plot in 1D nearest ambiguity location vs platform spacing
- plot in 2D min # of platforms required vs resolution and nearest ambiguity location