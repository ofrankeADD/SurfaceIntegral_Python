# Surface integrals of a fluid vector field over an interpolated surface region

### Input:
The script loads the relevant physical parameters such as the density scalar field, the velocity vector field and the positional vector field from a 3D datacube.
The 3D datacube is representing the output of the numerical simulation at a given snapshot time.
The simulation itself is employing the moving meshcode AREPO.

### Goal:
We want to interpolate a vector field from the mesh data onto a given surface grid area.
The surface grid area is a combination of two spherical caps and one conical frustum.
Here is a schematic view of our assembled surface region object, consisting of a truncated cone enclosed with two spherical caps:

<img src="https://user-images.githubusercontent.com/49908052/143508151-4ab517e5-613a-49b7-9e5c-f60fa14503a7.png" width="250"/>

For simplicity, a spherical cap is referenced by a "Sphere" class object and a truncated cone is referenced by a "Cone" class object in the code.
Both classes are inherited from the parent "Surface" class and are combined together in the "Region" class.

### Procedure for interpolating a "Sphere" object:
- Crop a thin spherical shell from the 3D input datacube, where the center of the shell coincides with the geometric center of the datacube.
- Call this shell_xyz.
- Initiate a spherical surface grid in spherical coordinates theta and phi using numpy.meshgrid.
- Compute the elementary surface area dS = r^2 * sin(theta) * dtheta * dphi.
- Transform the spherical grid back into xyz coordinates.
- Call this grid_xyz.
- Interpolate the velocity vector field v from shell_xyz onto grid_xyz using scipy.NearestNDInterpolator.
- Call this v_interp.
- Compute the surface integral np.sum(np.sum(v_interp * dS)), which is an approximated result depending on the resolution of grid_xyz.
- Compare the approximated result with the analytical exact solution if available.
- Here is a plot showing how the velocity vector field from a thin spherical shell is interpolated onto a spherical grid and then projected onto 2D:

<img src="https://user-images.githubusercontent.com/49908052/143507704-1fbaf4d7-3e85-466b-b4bf-ebfc0595824e.png" width="700">

- Here is a interpolated result of a known vector field function G versus its analytical result (G=cos(theta) * sin(4 * phi)), which is also projected onto 2D:

<img src="https://user-images.githubusercontent.com/49908052/143507708-bdd889cc-5a7f-40df-9eca-8a0cb8252ddf.png" width="700">

### Procedure for interpolating a "Cone" object:
- Instead of cropping a spherical shell, crop a thin lateral shell of a conical frustum.
- Initiate a conical surface grid using cylindrical coordinates instead of spherical.
- Take the same steps as before.

### Visualization of the meshgrid of the surface grid objects on which the actual data points are inerpolated at:
- points in red are belonging to a "Cone" object: the lateral area of a conical frustum
- points in blue are belonging to a first "Sphere" object: the top spherical cap
- points in green are belonging to a second "Sphere" object: the bottom spherical cap
- the two subfigures in the first column are showing the constructed surface grids (grid_xyz) for each "Surface" object
- the two subfigures in the second column are showing the datapoints within thin shells (shell_xyz) for each "Surface" object
- you can see that grid_xyz and shell_xyz coincide very well
- each row of subfigures is showing a different set of parameters consisting of the heigth of the conical frustum and the height of the spherical shell

<p float="left">
<img src="https://user-images.githubusercontent.com/49908052/143508318-79e6d2f6-8205-4535-a96a-6d6bf7e85a22.png" width="400"/>
<img src="https://user-images.githubusercontent.com/49908052/143508320-6b116e94-6c7b-4a56-9e08-8badca03ccb2.png" width="400"/>
</p>
<p float="left">
<img src="https://user-images.githubusercontent.com/49908052/143508674-621019db-b196-4f1f-9f28-db16c34bd378.png" width="400"/>
<img src="https://user-images.githubusercontent.com/49908052/143508676-6026053b-d3af-41f0-ab01-2acfc0bccc62.png" width="400"/>
</p>

