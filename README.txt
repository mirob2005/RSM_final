Michael Robertson
mirob2005@gmail.com, miro2005@csu.fullerton.edu
CS 566
HW5 Final Project
Due 5/17/2012
SID: 892-32-2629


Assignment was completed using Visual Studio 2008 Pro on Windows 7 64-bit.
Included is all .cpp, .h, and .txt to independently run using any compiler.

All files must be in the same directory to be found.

* The project extends our raytace last seen in HW3 to use Reflective Shadow Maps for indirect 
	lighting. ("Reflective Shadow Maps" by Carsten Dachsbacher and Marc Stamminger)
	
* The program computes and outputs the:
	(from view of camera):
	* color image (with shadows by shadow maps or light intersection tests) with indirect lighting
	* indirect image (just indirect lighting, no direct lighting)
	(from view of the light source):
	* depth map (as seen from the view of the light source, not camera as in HW1-3)
	* normal map
	* world coordinate map
	* flux buffer

Provided is a single scene file and image outputs portraying cornell box scene (with spheres)
The ability to modify sampling and weights is implemented using the variable sampleRate and the 
results are below:

RESULTS: Using a 250x250 resolution image = 62500 pixels
# samples				sampleRate			Time
62500 samples weighted by 1     = 480 seconds
12500 samples weighted by 5     = 100 seconds
6250   samples weighted by 10   = 50 seconds
2500   samples weighted by 25   = 20 seconds
1250   samples weighted by 50   = 10 seconds
625     samples weighted by 100 = 5 seconds

No indirect lighting						 = <1 second

All sampling results have images provided.

Limitations: Using only point lights so depth image must be done in perspective. Using only
spheres and planes. Orthographic camera works, but using only perspective camera to get wide
angle view to capture entire room.

All images created are in ppm format only.

Sample commandline(on Windows) :
			raytrace.exe -v -i -s scene5_01.txt -o scene5_01_out.ppm -d scene5_01_depth.ppm -r 250

-r passes in the resolution (required)
-v turns on text output such as display parsed data, etc. (optional)
-s turns on shadow maps instead of light vector intersection tests for shadows

----------------------------------------------------------------------------------------------------------------
The resolution is defined on the commandline as listed as an option so viewplane is not defined in 
the grammar.

NOTE: On windows, program required a GL folder with the glut libraries in order 
to compile due to the use of GLUbyte
