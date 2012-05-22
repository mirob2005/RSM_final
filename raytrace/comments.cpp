///* Michael Robertson
// * mirob2005@gmail.com, miro2005@csu.fullerton.edu
// * CS 566
// * Due 3/5/2012
// * SID: 892-32-2629
// *
// * $Id: raytrace.cpp 1961 2010-02-24 08:46:53Z mshafae $
// *
// */
//
//
//#include <time.h>
////#include <iostream>
//#include <string>
//#include <math.h>
//#include "getopt.h"
//#include "Scene.h"
//#include "FaceList.h"
//#include "PlyModel.h"
//
//using namespace std;
//
//FaceList *gInputModel;
//Image *color_image;
//Image *depth_image;
//Scene *gTheScene;
//string gProgramName;
//int gResolution;
//bool gVerbose = false;
////Shadows set on when true (default)
//bool gShadows = true;
//
//void usage( string message = "" ){
//	cerr << message << endl;
//	cerr << gProgramName << " -i <inputfile> -o <outputfile> -d <depthfile>" << endl;
//	cerr << "          -or-" << endl;
//	cerr << gProgramName << " --input <inputfile> --output <outputfile> --depth <depthfile>" << endl;
//	
//}
//
//std::ostream& operator <<( std::ostream &out, const Pixel &p ){
//  p.write( out );
//  return( out );
//}
//
//void parseCommandLine( int argc, char **argv ){
//	int ch;
//	string inputFile( "" ), outputFile( "" ), depthFile( "" );
//	
//	static struct option longopts[] = {
//    { "input", required_argument, NULL, 'i' },
//    { "output", required_argument, NULL, 'o' },
//    { "depth", required_argument, NULL, 'd' },
//    { "resolution", required_argument, NULL, 'r' },
//    { "verbose", required_argument, NULL, 'v' },
//    { "help", required_argument, NULL, 'h' },
//	{ "shadow", required_argument, NULL, 's' },
//    { NULL, 0, NULL, 0 }
//	};
//
//	while( (ch = getopt_long(argc, argv, "i:o:d:r:vh", longopts, NULL)) != -1 ){
//		switch( ch ){
//			case 'i':
//				// input file
//				inputFile = string(optarg);
//				break;
//			case 'o':
//				// image output file
//				outputFile = string(optarg);
//				break;
//			case 'd':
//				// depth output file
//				depthFile = string( optarg );
//				break;
//			case 'r':
//			   gResolution = atoi(optarg);
//			break;
//      case 'v':
//        gVerbose = true;
//        break;
//	  case 's':
//		  gShadows = false;
//		  break;
//      case 'h':
//        usage( );
//        break;
//	  default:
//		// do nothing
//		break;
//		}
//	}
//	gTheScene = new Scene( inputFile, outputFile, depthFile );
//}
//
//float* normalize3D(float* vector) {
//   float vectorSquared = 0;
//   float magnitude = 0;
//   vectorSquared = vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2];
//   magnitude = sqrt(vectorSquared);
//   float normal[3];
//   normal[0] = vector[0]/magnitude;
//   normal[1] = vector[1]/magnitude;
//   normal[2] = vector[2]/magnitude;
//
//   return normal;
//}
//
//
//int main( int argc, char **argv ){
//    time_t start,end;
//	time (&start);
//	string pathStr;
//	gProgramName = argv[0];
//	gInputModel = NULL;
//	char *in = NULL;
//
//	parseCommandLine( argc, argv );
//	argc -= optind;
//	argv += optind;
//	if( gTheScene->hasInputSceneFilePath( ) &&
//			gTheScene->hasOutputFilePath( ) &&
//			gTheScene->hasDepthFilePath( ) ){
//		gTheScene->parse( );	
//
//		float diff = 0;
//		float midpt[3] = {0,0,0};
//
//	   const char * ply_model = gTheScene->meshFile().c_str();
//
//	   gInputModel = readPlyModel( ply_model );
//	   if(gInputModel==0){
//		  if(gVerbose) cout << "No Ply Model Loaded" << endl;
//	   }
//	   else {
//	  	   float maxX = -100;
//		   float maxY = -100;
//		   float maxZ = -100;
//  		   float minX = 100;
//		   float minY = 100;
//		   float minZ = 100;
//		   //COMPUTING BOUNDING SPHERE DIMENSIONS
//		   for( unsigned int i = 0; i < gInputModel->vertexCount(); i++ ){
//			  
//			  if(maxX < gInputModel->vertices[i][0]) maxX = gInputModel->vertices[i][0];
//			  if(minX > gInputModel->vertices[i][0]) minX = gInputModel->vertices[i][0];
//			  if(maxY < gInputModel->vertices[i][1]) maxY = gInputModel->vertices[i][1];
//			  if(minY > gInputModel->vertices[i][1]) minY = gInputModel->vertices[i][1];
//			  if(maxZ < gInputModel->vertices[i][2]) maxZ = gInputModel->vertices[i][2];
//			  if(minZ > gInputModel->vertices[i][2]) minZ = gInputModel->vertices[i][2];
//			  
//		   }
//		  
//		  diff = sqrt( (maxX-minX)*(maxX-minX)+(maxY-minY)*(maxY-minY)+(maxZ-minZ)*(maxZ-minZ) );
//		  //cout << "Radius = " << (diff/2) << endl;
//		  midpt[0] = maxX-((maxX-minX)/2);
//		  midpt[1] = maxY-((maxY-minY)/2);
//		  midpt[2] = maxZ-((maxZ-minZ)/2);
//		  //cout << "Midpt" << midpt[0] << ", " << midpt[1] << ", " << midpt[2] << endl;
//		 if(gVerbose) printf("Ply Model %s Loaded\n", ply_model);
//		 // cout << "Face Count = " << gInputModel->faceCount( ) << endl;
//		 // cout << "Vertex Count = " << gInputModel->vertexCount( ) << endl;
//		 // 
//		 // for(int i = 0; i < gInputModel->faceCount( ); i++){
//			//// normal [ faces [face#] [vertex(0-2)] ] [Nx Ny Nz [0-2]
//			// cout << "Face #" << i << endl;
//			// cout << "Face Normal = [" << gInputModel->faceNormals[i][0] << ", " << gInputModel->faceNormals[i][1] << ", " << gInputModel->faceNormals[i][2] << "]" << endl;
//			// for( unsigned int j = 0; j < 3; j++ ){
//			//	cout << "Vertex Normal = [" << gInputModel->vertexNormals[gInputModel->faces[i][j]][0] << ", " <<
//			//								   gInputModel->vertexNormals[gInputModel->faces[i][j]][1] << ", " <<
//			//								   gInputModel->vertexNormals[gInputModel->faces[i][j]][2] << "]" << endl;
//			//	//[ faces [face#] [vertex(0-2)] ] [xyz[0-2]
//			//	cout << "VertexX = " << gInputModel->vertices[gInputModel->faces[i][j]][0]<< endl;
//			//	cout << "VertexY = " << gInputModel->vertices[gInputModel->faces[i][j]][1]<< endl;
//			//	cout << "VertexZ = " << gInputModel->vertices[gInputModel->faces[i][j]][2]<< endl;
//			//	cin.get();
//			// }
//		 // }
//	   }
//
//	   //35mm
//	  float baseAngle = 54.4322;
//	  float newAngle = 35;
//	  float viewPercentage = 1;
//	  //Calculate Viewing Angle from focal distance
//	  if(gTheScene->camType() == 2){
//		 newAngle = (180*2*atan(36/(2*gTheScene->camDistance())))/3.1416;
//		 viewPercentage = baseAngle/newAngle;
//	  }
//	  //Use Given Viewing Angle
//	  if(gTheScene->camType() == 1){
//		 newAngle = gTheScene->camAngle();
//		 viewPercentage = baseAngle/newAngle;
//	  }
//
//		//OUTPUT THE PARSED DATA
//		if(gVerbose) cout << *gTheScene << endl;
//
//		//Access Camera Type
//		if(gVerbose) cout << "Cam Type " << gTheScene->camType() << endl;
//		
//		//Convert String to const char
//		const char * color_out = gTheScene->outputFile().c_str();
//		const char * depth_out = gTheScene->depthFile().c_str();
//
//		//Create new Image with x_res, y_res, and background color
//		color_image = new Image(gResolution,gResolution,gTheScene->backgroundColor());
//		Pixel *cp = color_image->pixels;
//		
//		float black[3] = {0,0,0};
//		depth_image = new Image(gResolution,gResolution,black);
//		Pixel *dp = depth_image->pixels;
//
//		////////////////////////////////////////////////////////////////////////////////
//		float W = 5/viewPercentage;
//		float H = 5/viewPercentage;
//		float resX = gResolution;
//		float resY = gResolution;
//
//		//Check for camera up vector
//		bool zRotate = false;
//		float result = 0;
//		float mag = sqrt(gTheScene->camUp()[0]*gTheScene->camUp()[0]+gTheScene->camUp()[1]*gTheScene->camUp()[1]);
//
//		if(gTheScene->camUp()[0]>0.0001 || gTheScene->camUp()[0]< -0.0001 ){
//			zRotate = true;
//			result = acos(gTheScene->camUp()[1]/mag);
//			if(gTheScene->camUp()[0]>0.0001){
//			   result = -result;
//			}
//		}
//
//		float* tArray = new float[resX*resY];
//		for(int i = 0; i < resX*resY; i++) {
//			tArray[i] = 20;
//		}
//
//		int currentPixel = 0;
//	  
//		while(currentPixel<resX*resY){
//
//  	   	  //int column = currentPixel%int(2*W);
//		  int column = currentPixel%int(resX);
//		  float coordX = -(W/2)+(W/(resX-1))*(currentPixel%int(resX));//*(currentPixel%resX));
//		  float coordY = (H/2)-(H/(resY-1))*floor(currentPixel/resY);
//		  //float coordX = ((-W + ((resX+1)/resX)*column)/gPixelZoom);
//		  //float coordY = ((H - floor(currentPixel/resX)*((resX+1)/resX))/gPixelZoom);
//		  //cout << currentPixel<< "__" <<coordX << " " << coordY << endl;
//		  //cin.get();
//
//		  if(zRotate) {
//			  float tempX = coordX;
//			  float tempY = coordY;
//			  coordX = cos(result)*tempX - sin(result)*tempY;
//		      coordY = sin(result)*tempX + cos(result)*tempY;
//		  }
//
//
//		  float eye[3] = {0,0,0};
//		  float d[3] = {0,0,0};
//		  
//		//Perspective Camera
//		if(gTheScene->camType() == 1){
//		   eye[0] = gTheScene->camCenter()[0];
//		   eye[1] = gTheScene->camCenter()[1];
//		   eye[2] = gTheScene->camCenter()[2];
//		   d[0] = coordX - gTheScene->camCenter()[0];
//		   d[1] = coordY - gTheScene->camCenter()[1];
//		   d[2] = 0 - gTheScene->camCenter()[2];
//		}
//		//Simple Perspective Camera
//		else if(gTheScene->camType() == 2) {
//		   eye[0] = gTheScene->camCenter()[0];
//		   eye[1] = gTheScene->camCenter()[1];
//		   eye[2] = gTheScene->camCenter()[2];
//		   d[0] = coordX - gTheScene->camCenter()[0];
//		   d[1] = coordY - gTheScene->camCenter()[1];
//		   d[2] = 0 - gTheScene->camCenter()[2];
//		}
//		//Orthographic Projection
//		else {
//		   eye[0] = coordX+gTheScene->camCenter()[0];
//		   eye[1] = coordY+gTheScene->camCenter()[1];
//		   eye[2] = gTheScene->camCenter()[2];
//		   d[0] = gTheScene->camDirection()[0];
//		   d[1] = gTheScene->camDirection()[1];
//		   d[2] = gTheScene->camDirection()[2];
//		}
//		////A = dNorm DOT dNorm
//		 float A = normalize3D(d)[0]*normalize3D(d)[0]+
//				   normalize3D(d)[1]*normalize3D(d)[1]+
//				   normalize3D(d)[2]*normalize3D(d)[2];
//
//			//SPHERES
//		  for(int i = 0; i < gTheScene->numSpheres(); i++) {
//
//			 float center[3] = {gTheScene->allSpheres()[1+5*i],gTheScene->allSpheres()[2+5*i],gTheScene->allSpheres()[3+5*i]};
//			 float radius = gTheScene->allSpheres()[4+5*i];
//
//
//
//			 ////B = 2(eye-center) DOT dNorm
//			 float eyeMinusCenter[3] = {eye[0]-center[0],eye[1]-center[1],eye[2]-center[2]};
//			 
//			 float B = 2*(eyeMinusCenter[0]*(normalize3D(d)[0])+
//					      eyeMinusCenter[1]*(normalize3D(d)[1])+
//					      eyeMinusCenter[2]*(normalize3D(d)[2]));
//
//			 //C = (eye-center) DOT (eye-center) - radius^2
//			 float C = (eyeMinusCenter[0]*eyeMinusCenter[0]+
//						eyeMinusCenter[1]*eyeMinusCenter[1]+
//						eyeMinusCenter[2]*eyeMinusCenter[2]) - (radius*radius);
//
//			 //B^2 - 4AC
//			 float Q = (B*B) - 4*(A*C);
//
//			 if(Q >= 0){
//				float tPlus = (-B+sqrt(Q))/(2*A);
//				float tMinus = (-B-sqrt(Q))/(2*A);
//
//	
//				float point[3] = {0,0,0};
//				//Perspective Camera
//				if(gTheScene->camType() == 1){
//				  point[0] = gTheScene->camCenter()[0] + tMinus*(normalize3D(d)[0]);
//				  point[1] = gTheScene->camCenter()[1] + tMinus*(normalize3D(d)[1]);
//				  point[2] = gTheScene->camCenter()[2] + tMinus*(normalize3D(d)[2]);
//				}
//			    //Simple Perspective Camera
//				else if(gTheScene->camType() == 2) {
//				  point[0] = gTheScene->camCenter()[0] + tMinus*(normalize3D(d)[0]);
//				  point[1] = gTheScene->camCenter()[1] + tMinus*(normalize3D(d)[1]);
//				  point[2] = gTheScene->camCenter()[2] + tMinus*(normalize3D(d)[2]);
//				}
//			    //Orthographic Projection
//			    else {		
//				  point[0] = coordX+gTheScene->camCenter()[0]+tMinus*(normalize3D(d)[0]);
//				  point[1] = coordY+gTheScene->camCenter()[1]+tMinus*(normalize3D(d)[1]);
//				  point[2] = gTheScene->camCenter()[2]+tMinus*(normalize3D(d)[2]);
//				}
//		
//				
//			    float normal[3] = {point[0] - center[0], point[1] - center[1], point[2] - center[2]};
//
//				float sumLightsR = 0;
//	   		    float sumLightsG = 0;
//				float sumLightsB = 0;
//				
//				for(int j = 0; j < gTheScene->numLights(); j++) {
//				  float lightVector[3] = {gTheScene->lights()[0+9*j] - point[0], 
//										  gTheScene->lights()[1+9*j] - point[1], 
//										  gTheScene->lights()[2+9*j] - point[2]};
//				  float NdotL = (normalize3D(normal)[0])*(normalize3D(lightVector)[0])+
//							  (normalize3D(normal)[1])*(normalize3D(lightVector)[1])+
//							  (normalize3D(normal)[2])*(normalize3D(lightVector)[2]);
//				   if(NdotL < 0.00){
// 						 NdotL = 0;
//				   }
//				   if(NdotL < 0.01 && NdotL > -0.01){
// 						 NdotL = 0;
//				   }
//   				 
//				   sumLightsR = sumLightsR+(gTheScene->lights()[3+9*j]*NdotL);
//				   sumLightsG = sumLightsG+(gTheScene->lights()[4+9*j]*NdotL);
//				   sumLightsB = sumLightsB+(gTheScene->lights()[5+9*j]*NdotL);
//      					
//
//				}
//   			   
//
//			    if(sumLightsR>1) sumLightsR = 1;
//			    if(sumLightsG>1) sumLightsG = 1;
//			    if(sumLightsB>1) sumLightsB = 1;
//
//				if(tMinus >= 1){
//				   if(tMinus < tArray[currentPixel]){
//					  tArray[currentPixel] = tMinus;
//					  *cp = Pixel(255*(gTheScene->materials()[0+7*int(gTheScene->allSpheres()[0+5*i])]*sumLightsR),
//								  255*(gTheScene->materials()[1+7*int(gTheScene->allSpheres()[0+5*i])]*sumLightsG),
//								  255*(gTheScene->materials()[2+7*int(gTheScene->allSpheres()[0+5*i])]*sumLightsB));
//				   }
//				}
//				
//			 }
//			 //PLANES
//			 for(int i = 0; i < gTheScene->numPlanes(); i++) {
//				float normal[3] = {gTheScene->allPlanes()[1+5*i],
//								   gTheScene->allPlanes()[2+5*i],
//								   gTheScene->allPlanes()[3+5*i]};
//				float normalDOTd = (normalize3D(normal)[0]) * (normalize3D(d)[0]) +
//								   (normalize3D(normal)[1]) * (normalize3D(d)[1]) +
//								   (normalize3D(normal)[2]) * (normalize3D(d)[2]);
//
//				if(normalDOTd >=0.0001 || normalDOTd <= -0.0001) {
//				   float a[3] = {(normalize3D(normal)[0])*gTheScene->allPlanes()[4+5*i],
//								 (normalize3D(normal)[1])*gTheScene->allPlanes()[4+5*i],
//								 (normalize3D(normal)[2])*gTheScene->allPlanes()[4+5*i]};
//	
//				   float aMinusO[3] = {a[0] - eye[0],a[1] - eye[1],a[2] - eye[2]};
//
//				   float t = (aMinusO[0]*normal[0]+
//							  aMinusO[1]*normal[1]+
//							  aMinusO[2]*normal[2])/normalDOTd;
//
//
//				   float point[3] = {0,0,0};
//				   //Perspective Camera
//				   if(gTheScene->camType() == 1){
//					 point[0] = gTheScene->camCenter()[0] + t*(normalize3D(d)[0]);
//					 point[1] = gTheScene->camCenter()[1] + t*(normalize3D(d)[1]);
//					 point[2] = gTheScene->camCenter()[2] + t*(normalize3D(d)[2]);
//				   }
//				   //Simple Perspective Camera
//				   else if(gTheScene->camType() == 2) {
//					 point[0] = gTheScene->camCenter()[0] + t*(normalize3D(d)[0]);
//					 point[1] = gTheScene->camCenter()[1] + t*(normalize3D(d)[1]);
//					 point[2] = gTheScene->camCenter()[2] + t*(normalize3D(d)[2]);
//				   }
//				   //Orthographic Projection
//				   else {		
//					 point[0] = coordX+gTheScene->camCenter()[0]+t*(normalize3D(d)[0]);
//					 point[1] = coordY+gTheScene->camCenter()[1]+t*(normalize3D(d)[1]);
//					 point[2] = gTheScene->camCenter()[2]+t*(normalize3D(d)[2]);
//				   }
//				   float sumLightsR = 0;
//	   			   float sumLightsG = 0;
//				   float sumLightsB = 0;
//   				
//				   for(int j = 0; j < gTheScene->numLights(); j++) {
//					 float lightVector[3] = {gTheScene->lights()[0+9*j] - point[0], 
//											 gTheScene->lights()[1+9*j] - point[1], 
//											 gTheScene->lights()[2+9*j] - point[2]};
//					 float NdotL = (normalize3D(normal)[0])*(normalize3D(lightVector)[0])+
//								 (normalize3D(normal)[1])*(normalize3D(lightVector)[1])+
//								 (normalize3D(normal)[2])*(normalize3D(lightVector)[2]);
//					  if(NdotL < 0.00){
// 							NdotL = 0;
//					  }
//					  if(NdotL < 0.01 && NdotL > -0.01){
// 							NdotL = 0;
//					  }
//      				 
//					  sumLightsR = sumLightsR+(gTheScene->lights()[3+9*j]*NdotL);
//					  sumLightsG = sumLightsG+(gTheScene->lights()[4+9*j]*NdotL);
//					  sumLightsB = sumLightsB+(gTheScene->lights()[5+9*j]*NdotL);
//         					
//
//				   }
//      			   
//
//				   if(sumLightsR>1) sumLightsR = 1;
//				   if(sumLightsG>1) sumLightsG = 1;
//				   if(sumLightsB>1) sumLightsB = 1;
//
//				   if(t >= 1){
//					  
//					  if(t < tArray[currentPixel]){
//						 tArray[currentPixel] = t;
//						 *cp = Pixel(255*(gTheScene->materials()[0+7*int(gTheScene->allPlanes()[0+5*i])]*sumLightsR),
//									 255*(gTheScene->materials()[1+7*int(gTheScene->allPlanes()[0+5*i])]*sumLightsG),
//									 255*(gTheScene->materials()[2+7*int(gTheScene->allPlanes()[0+5*i])]*sumLightsB));
//					  }
//				   }
//				}
//
//			}
//		  }
//			 //TRIANGLE
//			 for(int i = 0; i < gTheScene->numTriangles(); i++) {
//			   float Aa = (gTheScene->allTriangles()[1+10*i])-(gTheScene->allTriangles()[4+10*i]);
//			   float Ab = (gTheScene->allTriangles()[2+10*i])-(gTheScene->allTriangles()[5+10*i]);
//			   float Ac = (gTheScene->allTriangles()[3+10*i])-(gTheScene->allTriangles()[6+10*i]);
//
//			   float Ad = (gTheScene->allTriangles()[1+10*i])-(gTheScene->allTriangles()[7+10*i]);
//			   float Ae = (gTheScene->allTriangles()[2+10*i])-(gTheScene->allTriangles()[8+10*i]);
//			   float Af = (gTheScene->allTriangles()[3+10*i])-(gTheScene->allTriangles()[9+10*i]);
//
//			   float Ag = d[0];
//			   float Ah = d[1];
//			   float Ai = d[2];
//				  
//			   float Aj = (gTheScene->allTriangles()[1+10*i])-eye[0];
//			   float Ak = (gTheScene->allTriangles()[2+10*i])-eye[1];
//			   float Al = (gTheScene->allTriangles()[3+10*i])-eye[2];
//
//			   float M = (Aa*((Ae*Ai)-(Ah*Af))) + (Ab*((Ag*Af)-(Ad*Ai))) + (Ac*((Ad*Ah)-(Ae*Ag)));
//			   float t = 0;
//			   float Beta = 0;
//			   float Gamma = 0;
//			   if(M >= 0.0001 || M <= -0.0001){
//				  t = ((Af*((Aa*Ak)-(Aj*Ab))) + (Ae*((Aj*Ac)-(Aa*Al))) + (Ad*((Ab*Al)-(Ak*Ac))))/M;
//				  Beta = ((Aj*((Ae*Ai)-(Ah*Af))) + (Ak*((Ag*Af)-(Ad*Ai))) + (Al*((Ad*Ah)-(Ae*Ag))))/M;
//				  Gamma = ((Ai*((Aa*Ak)-(Aj*Ab))) + (Ah*((Aj*Ac)-(Aa*Al))) + (Ag*((Ab*Al)-(Ak*Ac))))/M;
//				  t = d[2]*t;
//				  //cout << t << endl;
//				  //cout << d[0] << ", " << d[1] << ", " << d[2] << endl;
//				  //cin.get();
//
//				  if(t >= 1){
//					 if(t < tArray[currentPixel]){
//						
//						 float point[3] = {0,0,0};
//						 float BMA[3] = {(gTheScene->allTriangles()[4+10*i])-(gTheScene->allTriangles()[1+10*i]),
//											 (gTheScene->allTriangles()[5+10*i])-(gTheScene->allTriangles()[2+10*i]),
//											 (gTheScene->allTriangles()[6+10*i])-(gTheScene->allTriangles()[3+10*i])};
//						 float CMA[3] = {(gTheScene->allTriangles()[7+10*i])-(gTheScene->allTriangles()[1+10*i]),
//											 (gTheScene->allTriangles()[8+10*i])-(gTheScene->allTriangles()[2+10*i]),
//											 (gTheScene->allTriangles()[9+10*i])-(gTheScene->allTriangles()[3+10*i])};
//						 //(B-A)X(C-A)
//						 float normal[3] = {BMA[1]*CMA[2]-BMA[2]*CMA[1],-BMA[0]*CMA[2]+BMA[2]*CMA[0],BMA[0]*CMA[1]-BMA[1]*CMA[0]};
//						 //cout << normal[0] << ", " << normal[1] << ", " << normal[2] << endl;
//						 
//						 //Perspective Camera
//						 if(gTheScene->camType() == 1){
//						   point[0] = gTheScene->camCenter()[0] + t*(normalize3D(d)[0]);
//						   point[1] = gTheScene->camCenter()[1] + t*(normalize3D(d)[1]);
//						   point[2] = gTheScene->camCenter()[2] + t*(normalize3D(d)[2]);
//						 }
//						 //Simple Perspective Camera
//						 else if(gTheScene->camType() == 2) {
//						   point[0] = gTheScene->camCenter()[0] + t*(normalize3D(d)[0]);
//						   point[1] = gTheScene->camCenter()[1] + t*(normalize3D(d)[1]);
//						   point[2] = gTheScene->camCenter()[2] + t*(normalize3D(d)[2]);
//						 }
//						 //Orthographic Projection
//						 else {		
//						   point[0] = coordX+gTheScene->camCenter()[0]+t*(normalize3D(d)[0]);
//						   point[1] = coordY+gTheScene->camCenter()[1]+t*(normalize3D(d)[1]);
//						   point[2] = gTheScene->camCenter()[2]+t*(normalize3D(d)[2]);
//						 }
//						 float sumLightsR = 0;
//	   					 float sumLightsG = 0;
//						 float sumLightsB = 0;
//         				
//						 for(int j = 0; j < gTheScene->numLights(); j++) {
//						   float lightVector[3] = {gTheScene->lights()[0+9*j] - point[0], 
//												   gTheScene->lights()[1+9*j] - point[1], 
//												   gTheScene->lights()[2+9*j] - point[2]};
//						   float NdotL = (normalize3D(normal)[0])*(normalize3D(lightVector)[0])+
//									   (normalize3D(normal)[1])*(normalize3D(lightVector)[1])+
//									   (normalize3D(normal)[2])*(normalize3D(lightVector)[2]);
//						   
//							if(NdotL < 0.00){
// 								  NdotL = -NdotL;
//							}
//							if(NdotL < 0.01 && NdotL > -0.01){
// 								  NdotL = 0;
//							}
//            				 
//							sumLightsR = sumLightsR+(gTheScene->lights()[3+9*j]*NdotL);
//							sumLightsG = sumLightsG+(gTheScene->lights()[4+9*j]*NdotL);
//							sumLightsB = sumLightsB+(gTheScene->lights()[5+9*j]*NdotL);
//               					
//
//						 }
//            			   
//
//						 if(sumLightsR>1) sumLightsR = 1;
//						 if(sumLightsG>1) sumLightsG = 1;
//						 if(sumLightsB>1) sumLightsB = 1;
//
//						if(Beta >-0.0001 && Beta <1.0001){
//						   if(Gamma > -0.0001 && Gamma < 1.0001){
//							  if(Beta+Gamma < 1.0001){
//								 tArray[currentPixel] = t;
//								 *cp = Pixel(255*(gTheScene->materials()[0+7*int(gTheScene->allTriangles()[0+10*i])]*sumLightsR),
//										  255*(gTheScene->materials()[1+7*int(gTheScene->allTriangles()[0+10*i])]*sumLightsG),
//										  255*(gTheScene->materials()[2+7*int(gTheScene->allTriangles()[0+10*i])]*sumLightsB));
//							  }
//						   }
//						}
//					 }
//				  }
//
//			   }
//
//			   //cout << Aa << endl;
//			   //cout << Ab << endl;
//			   //cout << Ac << endl;
//			   //cout << Ad << endl;
//			   //cout << Ae << endl;
//			   //cout << Af << endl;
//			   //cout << Ag << endl;
//			   //cout << Ah << endl;
//			   //cout << Ai << endl;
//			   //cout << Aj << endl;
//			   //cout << Ak << endl;
//			   //cout << Al << endl;
//			   //cout << M << endl;
//			   //cout << t << endl;
//			   //cout << Beta << endl;
//			   //cout << Gamma << endl;
//			   //cin.get();
//			 }
//
//
//////////////////////////////////////////////////////////////////////////////////////////
// 			 //TRIANGLE MESH
//			 if(gInputModel!=0){
//
//
//				//BOUNDING VOLUME TEST (SPHERE)
//				float BVcenter[3] = {midpt[0],midpt[1], midpt[2]};
//				float BVradius = (diff/2);
//
//
//
//				////B = 2(eye-center) DOT dNorm
//				float eyeMinusCenter[3] = {eye[0]-BVcenter[0],eye[1]-BVcenter[1],eye[2]-BVcenter[2]};
//   			 
//				float B = 2*(eyeMinusCenter[0]*(normalize3D(d)[0])+
//						 eyeMinusCenter[1]*(normalize3D(d)[1])+
//						 eyeMinusCenter[2]*(normalize3D(d)[2]));
//
//				//C = (eye-center) DOT (eye-center) - radius^2
//				float C = (eyeMinusCenter[0]*eyeMinusCenter[0]+
//						   eyeMinusCenter[1]*eyeMinusCenter[1]+
//						   eyeMinusCenter[2]*eyeMinusCenter[2]) - (BVradius*BVradius);
//
//				//B^2 - 4AC
//				float Q = (B*B) - 4*(A*C);
//
//				if(Q >= 0){
//
//				   for(unsigned int i = 0; i < gInputModel->faceCount( ); i++) {
//					 float Aa = (gInputModel->vertices[gInputModel->faces[i][0]][0])-(gInputModel->vertices[gInputModel->faces[i][1]][0]);
//					 float Ab = (gInputModel->vertices[gInputModel->faces[i][0]][1])-(gInputModel->vertices[gInputModel->faces[i][1]][1]);
//					 float Ac = (gInputModel->vertices[gInputModel->faces[i][0]][2])-(gInputModel->vertices[gInputModel->faces[i][1]][2]);
//
//					 float Ad = (gInputModel->vertices[gInputModel->faces[i][0]][0])-(gInputModel->vertices[gInputModel->faces[i][2]][0]);
//					 float Ae = (gInputModel->vertices[gInputModel->faces[i][0]][1])-(gInputModel->vertices[gInputModel->faces[i][2]][1]);
//					 float Af = (gInputModel->vertices[gInputModel->faces[i][0]][2])-(gInputModel->vertices[gInputModel->faces[i][2]][2]);
//
//					 float Ag = d[0];
//					 float Ah = d[1];
//					 float Ai = d[2];
//      				  
//					 float Aj = (gInputModel->vertices[gInputModel->faces[i][0]][0])-eye[0];
//					 float Ak = (gInputModel->vertices[gInputModel->faces[i][0]][1])-eye[1];
//					 float Al = (gInputModel->vertices[gInputModel->faces[i][0]][2])-eye[2];
//
//					 float M = (Aa*((Ae*Ai)-(Ah*Af))) + (Ab*((Ag*Af)-(Ad*Ai))) + (Ac*((Ad*Ah)-(Ae*Ag)));
//					 float t = 0;
//					 float Beta = 0;
//					 float Gamma = 0;
//					 if(M >= 0.0001 || M <= -0.0001){
//						t = ((Af*((Aa*Ak)-(Aj*Ab))) + (Ae*((Aj*Ac)-(Aa*Al))) + (Ad*((Ab*Al)-(Ak*Ac))))/M;
//						Beta = ((Aj*((Ae*Ai)-(Ah*Af))) + (Ak*((Ag*Af)-(Ad*Ai))) + (Al*((Ad*Ah)-(Ae*Ag))))/M;
//						Gamma = ((Ai*((Aa*Ak)-(Aj*Ab))) + (Ah*((Aj*Ac)-(Aa*Al))) + (Ag*((Ab*Al)-(Ak*Ac))))/M;
//						t = d[2]*t;
//						//cout << t << endl;
//						//cout << d[0] << ", " << d[1] << ", " << d[2] << endl;
//						//cin.get();
//
//						if(t >= 1){
//						   if(t < tArray[currentPixel]){
//      						
//							   float point[3] = {0,0,0};
//
//							   //(B-A)X(C-A)
//							   float normal[3] = {gInputModel->faceNormals[i][0],gInputModel->faceNormals[i][1],gInputModel->faceNormals[i][2]};
//							   //cout << normal[0] << ", " << normal[1] << ", " << normal[2] << endl;
//      						 
//							   //Perspective Camera
//							   if(gTheScene->camType() == 1){
//								 point[0] = gTheScene->camCenter()[0] + t*(normalize3D(d)[0]);
//								 point[1] = gTheScene->camCenter()[1] + t*(normalize3D(d)[1]);
//								 point[2] = gTheScene->camCenter()[2] + t*(normalize3D(d)[2]);
//							   }
//							   //Simple Perspective Camera
//							   else if(gTheScene->camType() == 2) {
//								 point[0] = gTheScene->camCenter()[0] + t*(normalize3D(d)[0]);
//								 point[1] = gTheScene->camCenter()[1] + t*(normalize3D(d)[1]);
//								 point[2] = gTheScene->camCenter()[2] + t*(normalize3D(d)[2]);
//							   }
//							   //Orthographic Projection
//							   else {		
//								 point[0] = coordX+gTheScene->camCenter()[0]+t*(normalize3D(d)[0]);
//								 point[1] = coordY+gTheScene->camCenter()[1]+t*(normalize3D(d)[1]);
//								 point[2] = gTheScene->camCenter()[2]+t*(normalize3D(d)[2]);
//							   }
//							   float sumLightsR = 0;
//	   						   float sumLightsG = 0;
//							   float sumLightsB = 0;
//               				
//							   for(int j = 0; j < gTheScene->numLights(); j++) {
//								 float lightVector[3] = {gTheScene->lights()[0+9*j] - point[0], 
//														 gTheScene->lights()[1+9*j] - point[1], 
//														 gTheScene->lights()[2+9*j] - point[2]};
//								 float NdotL = (normalize3D(normal)[0])*(normalize3D(lightVector)[0])+
//											 (normalize3D(normal)[1])*(normalize3D(lightVector)[1])+
//											 (normalize3D(normal)[2])*(normalize3D(lightVector)[2]);
//      						   
//								  if(NdotL < 0.00){
// 										NdotL = 0;
//								  }
//								  if(NdotL < 0.01 && NdotL > -0.01){
// 										NdotL = 0;
//								  }
//                  				 
//								  sumLightsR = sumLightsR+(gTheScene->lights()[3+9*j]*NdotL);
//								  sumLightsG = sumLightsG+(gTheScene->lights()[4+9*j]*NdotL);
//								  sumLightsB = sumLightsB+(gTheScene->lights()[5+9*j]*NdotL);
//                     					
//
//							   }
//                  			   
//
//							   if(sumLightsR>1) sumLightsR = 1;
//							   if(sumLightsG>1) sumLightsG = 1;
//							   if(sumLightsB>1) sumLightsB = 1;
//							   //cout <<  " R = "<<sumLightsR << ", G = " << sumLightsG << ", B = " << sumLightsB << endl;
//							  if(Beta >-0.0001 && Beta <1.0001){
//								 if(Gamma > -0.0001 && Gamma < 1.0001){
//									if(Beta+Gamma < 1.0001){
//									   tArray[currentPixel] = t;
//
//									   *cp = Pixel(255*(gTheScene->materials()[0+7*gTheScene->materialIndex()[0]]*sumLightsR),
//												255*(gTheScene->materials()[1+7*gTheScene->materialIndex()[0]]*sumLightsG),
//												255*(gTheScene->materials()[2+7*gTheScene->materialIndex()[0]]*sumLightsB));
//									}
//								 }
//							  }
//						   }
//						}
//
//					 }
//
//					 //cout << Aa << endl;
//					 //cout << Ab << endl;
//					 //cout << Ac << endl;
//					 //cout << Ad << endl;
//					 //cout << Ae << endl;
//					 //cout << Af << endl;
//					 //cout << Ag << endl;
//					 //cout << Ah << endl;
//					 //cout << Ai << endl;
//					 //cout << Aj << endl;
//					 //cout << Ak << endl;
//					 //cout << Al << endl;
//					 //cout << M << endl;
//					 //cout << t << endl;
//					 //cout << Beta << endl;
//					 //cout << Gamma << endl;
//					 //cin.get();
//				   }
//				}
//			 }
//
//		  *cp++;
//		  currentPixel++;
//		}
//
//		//CALCULATE DEPTH IMAGE
//		float max = 0;
//		float min = 20;
//		currentPixel = 0;
//		while(currentPixel<resX*resY){
//		   if(tArray[currentPixel] >= 20){
//			  //do nothing
//		   }
//		   else if(max < tArray[currentPixel]) {
//			  max = tArray[currentPixel]; 
//		   }
//		   else if(min > tArray[currentPixel]) {
//			  min = tArray[currentPixel]; 
//		   }
//			
//			currentPixel++;
//		}
//
//		float range  = max-min;
//		if (range == 0){
//			range = 0.01;
//		}
//
//		currentPixel = 0;
//		while(currentPixel<resX*resY){
//		   if(tArray[currentPixel] >= 20){
//			   *dp = Pixel(0,0,0);
//		   }
//		   else{
//		   *dp = Pixel(100+(155*(1-(tArray[currentPixel]-min)/range)),
//			            100+(155*(1-(tArray[currentPixel]-min)/range)),
//			            100+(155*(1-(tArray[currentPixel]-min)/range)));	
//		   }
//		   *dp++;
//		   currentPixel++;
//		}
//
//		////////////////////////////////////////////////////////////////////////////////
//
//
//		if(color_image->write(color_out)){
//		   if(gVerbose) cout << "Color Image Created Succuessfully!" << endl;
//		}
//		else {
//		   if(gVerbose) cout << "Color Image Creation Failed!" << endl;
//		}
//		if(depth_image->write(depth_out)){
//		   if(gVerbose) cout << "Depth Image Created Succuessfully!" << endl;
//		}
//		else {
//		   if(gVerbose) cout << "Depth Image Creation Failed!" << endl;
//		}
//	
//	}else{
//		usage( "You specify an input scene file, an output file and a depth file." );
//	}
//	time (&end);
//	double dif = difftime (end,start);
//    if(gVerbose) cout << "Time = " << dif << endl;
//	if(gVerbose) cout << "REACHED END" << endl;
//    //cin.get();
//	return( 0 );
//}



				  //cout << "Normal = " << (normalize3D(normal)[0]) << ", " << (normalize3D(normal)[1]) << ", " << (normalize3D(normal)[2]) << endl;
				  //cout << "NdotL = " << NdotL << endl; 
				  //cout << "ReflectionVector = " << (normalize3D(reflectionVector)[0]) << ", " << (normalize3D(reflectionVector)[1]) << ", " << (normalize3D(reflectionVector)[2]) << endl;
				  //cout << "viewingVector = " << (normalize3D(viewingVector)[0]) << ", " << (normalize3D(viewingVector)[1]) << ", " << (normalize3D(viewingVector)[2]) << endl;
				  // cout << "RdotV = " << RdotV << endl;
				  //cout << "Lights " << sumLightsR << ", " << sumLightsG << ", " << sumLightsB << endl;
				  //cin.get();

				   
				   //sumLightsR = sumLightsR+(gTheScene->lights()[3+9*j]*NdotL);
				   //sumLightsG = sumLightsG+(gTheScene->lights()[4+9*j]*NdotL);
				   //sumLightsB = sumLightsB+(gTheScene->lights()[5+9*j]*NdotL);


					  //*cp = Pixel(255*(gTheScene->materials()[0+7*int(gTheScene->allSpheres()[0+5*i])]*sumLightsR),
							//	  255*(gTheScene->materials()[1+7*int(gTheScene->allSpheres()[0+5*i])]*sumLightsG),
							//	  255*(gTheScene->materials()[2+7*int(gTheScene->allSpheres()[0+5*i])]*sumLightsB));