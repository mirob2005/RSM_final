/* Michael Robertson
 * mirob2005@gmail.com, miro2005@csu.fullerton.edu
 * CS 566
 * HW3
 * Due 3/23/2012
 * SID: 892-32-2629
 *
 * $Id: Scene.cpp 2075 2010-04-12 08:06:14Z mshafae $
 *
 */

#include "Scene.h"
#include <string>
//#include <algorithm>
//#include <vector>
//#include <iterator>

Scene::Scene( string inputFilename, string outputFilename, string depthFilename, int resolution ) :
	myInputSceneFile( inputFilename ),
	myOutputFile( outputFilename ),
	myDepthFile( depthFilename ),
	myResolution( resolution ),
	//myCamera( ),
	lineNumber(0),
	tokenCount(0),
	length(0),
	//myWidth(0),
	//myHeight(0),
	myNumSpheres(0),
	myNumPlanes(0),
	myNumTriangles(0),
	myNumMeshes(0),
	myCamType(0),
	myCamAngle(0),
	myCamDistance(0),
	i(0),
	j(0)
{
	myNumberOfMaterials = -1;
	//materials = NULL;
	//myCurrentMaterial = NULL;
	//myGroup = NULL;
}

Scene::~Scene( ){
	// Nothing to free.
}

//Camera& Scene::camera( ){
//	return( myCamera );
//}

int* Scene::materialIndex( ){
	return( myMaterialIndex );
}

float* Scene::materials( ){
	return( mats );
}

float* Scene::allSpheres( ){
	return( spheres );
}
float* Scene::allPlanes( ){
	return( planes );
}
float* Scene::allTriangles( ){
	return( triangles );
}
float* Scene::allMeshes( ){
	return( meshes );
}
int Scene::numObjects(){
   return( myNumObjects);
}
int Scene::numSpheres(){
   return( myNumSpheres);
}
int Scene::numPlanes(){
   return( myNumPlanes);
}
int Scene::numTriangles(){
   return( myNumTriangles);
}
int Scene::numMeshes(){
   return( myNumMeshes);
}
int* Scene::lightType( ){
	return( myLightType );
}
float* Scene::dLight( ){
	return( myDLight );
}
float* Scene::lights( ){
	return( light );
}
int Scene::numLights(){
   return( myNumLights);
}
int Scene::camType(){
   return( myCamType);
}
float Scene::camAngle(){
   return( myCamAngle);
}
float Scene::camDistance(){
   return( myCamDistance);
}
float* Scene::camCenter( ){
	return( myCamCenter );
}

float* Scene::camDirection( ){
	return( myCamDirection );
}

float* Scene::camUp( ){
	return( myCamUp );
}

float* Scene::backgroundColor( ){
	return( myBackgroundColor );
}

int Scene::numberOfMaterials( ){
	return( myNumberOfMaterials );
}

//void Scene::setCurrentMaterial( int i ){
//	if( i >= myNumberOfMaterials ){
//	  throw( "Index out of range" );	
//	}else{
//		myCurrentMaterial = materials[i];
//	}
//}

//Material* Scene::currentMaterial( ){
//	return( myCurrentMaterial );
//}
//
//Group* Scene::group( ){
//	return( myGroup );
//}

string& Scene::inputSceneFile( ){
	return( myInputSceneFile );
}

string& Scene::outputFile( ){
	return( myOutputFile );
}

string& Scene::depthFile( ){
	return( myDepthFile );
}
string& Scene::meshFile( ){
	return( mesh );
}

int Scene::resolution( ){
	return( myResolution );
}

bool Scene::hasInputSceneFilePath( void ){
	bool ret = true;
	if( myInputSceneFile == "" ){
		ret = false;
	}
	return( ret );
}

bool Scene::hasOutputFilePath( void ){
	bool ret = true;
	if( myOutputFile == "" ){
		ret = false;
	}
	return( ret );
}

bool Scene::hasDepthFilePath( void ){
	bool ret = true;
	if( myDepthFile == "" ){
		ret = false;
	}
	return( ret );
}


void Scene::setInputSceneFile( string file ){
	myInputSceneFile = file;
}

void Scene::setOutputFile( string file ){
	myOutputFile = file;
}

void Scene::setDepthFile( string file ){
	myDepthFile = file;
}

float Scene::parseFloat( ){
	float ret = (float)atof( currentToken );
	return( ret );
}

double Scene::parseDouble( ){
	double ret = (double)atof( currentToken );
	return( ret );
}

int Scene::parseInt( ){
	int ret = atoi( currentToken );
	return( ret );
}

string Scene::parseString( ){
	string ret = string( currentToken );
	return( ret );
}

void Scene::checkToken( const char *str, const char *stage  ){
	if( strcmp( currentToken, str ) != 0 ){
		cerr << stage << " parse error at line " << lineNumber << " token " << tokenCount << ": " << currentToken << endl;
		cerr << "Current line: " << currentLine << endl;
		cerr << "Expected \'" << str << "\'" << endl;
		cin.get();
		exit( 1 );
	}
}
bool Scene::compareToken( const char *str){
	if( strcmp( currentToken, str ) != 0 ){
	  return false;
	}
	else {
	   return true;
	}
}

void Scene::parseCamera( ){
	nextToken( );
	if(compareToken( "OrthographicCamera")){
	   checkToken( "OrthographicCamera", "Camera" );
	   myCamType = 0;
	   nextToken( );
	   checkToken( "{", "Camera" );
	   nextToken( );
	   checkToken( "center", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamCenter[i] = parseFloat( );
	   }
	   nextToken( );
	   checkToken( "direction", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamDirection[i] = parseFloat( );
	   }
	   nextToken( );
	   checkToken( "up", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamUp[i] = parseFloat( );
	   }
   	
	   nextToken( );
	   checkToken( "}", "Camera" );
	}
	else if(compareToken( "PerspectiveCamera")){
	   checkToken( "PerspectiveCamera", "Camera" );
	   myCamType = 1;
	   nextToken( );
	   checkToken( "{", "Camera" );
	   nextToken( );
	   checkToken( "center", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamCenter[i] = parseFloat( );
	   }
	   nextToken( );
	   checkToken( "direction", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamDirection[i] = parseFloat( );
	   }
	   nextToken( );
	   checkToken( "up", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamUp[i] = parseFloat( );
	   }
	   nextToken( );
	   checkToken( "angle", "Camera" );
	   nextToken( );
	   myCamAngle = parseFloat( );
   	
	   nextToken( );
	   checkToken( "}", "Camera" );
	}
	else if(compareToken( "SimplePerspectiveCamera")){
	   checkToken( "SimplePerspectiveCamera", "Camera" );
	   myCamType = 2;
	   nextToken( );
	   checkToken( "{", "Camera" );
	   nextToken( );
	   checkToken( "center", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamCenter[i] = parseFloat( );
	   }
	   nextToken( );
	   checkToken( "direction", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamDirection[i] = parseFloat( );
	   }
	   nextToken( );
	   checkToken( "up", "Camera" );
	   for( int i = 0; i < 3; i++ ){
		   nextToken( );
		   myCamUp[i] = parseFloat( );
	   }
 	   nextToken( );
	   checkToken( "distance", "Camera" );
	   nextToken( );
	   myCamDistance = parseFloat( );

	   nextToken( );
	   checkToken( "}", "Camera" );
	}
	else {
	   checkToken( "OrthographicCamera", "Camera" );
	}
}
void Scene::parseLight( ){

   	nextToken( );
	checkToken( "Lights", "Lights" );
	nextToken( );
	checkToken( "{", "Lights" );
	nextToken( );
	checkToken( "numLights", "Lights" );
	nextToken( );
	myNumLights = parseInt();
	light = new float [myNumLights*9];
	myLightType = new int [myNumLights];
	myDLight = new float[myNumLights*4];


	
	for( int i = 0; i < myNumLights; i++) {
	   nextToken( );
	   //Check for Directional Light vs. Point Light
	   if(compareToken( "DirectionalLight")){
		  checkToken( "DirectionalLight", "Lights" ); 
		  
		  myLightType[i] = 0;

 		  nextToken( );
		  checkToken( "{", "Lights" );

		  nextToken( );
   	
 		  checkToken( "position", "Lights" );
		  nextToken( );
		  light[0+9*i] = parseFloat();
		  nextToken( );
 		  light[1+9*i] = parseFloat();
		  nextToken( );
		  light[2+9*i] = parseFloat();
		  nextToken( );
		 
		  //DIRECTION
		  //NEW ARRAY JUST FOR THE DIRECTIONAL LIGHTS TO HAVE DIRECTION AND ANGLE WITH 4 ELEMENTS FOR EACH LIGHT THAT IS DIRECTIONAL
 		  checkToken( "direction", "Lights" );
		  nextToken( );
		  myDLight[0+4*i] = parseFloat();
		  nextToken( );
 		  myDLight[1+4*i] = parseFloat();
		  nextToken( );
		  myDLight[2+4*i] = parseFloat();
		  nextToken( );


		  checkToken( "color", "Lights" );
		  nextToken( );
		  light[3+9*i] = parseFloat();
		  nextToken( );
 		  light[4+9*i] = parseFloat();
		  nextToken( );
		  light[5+9*i] = parseFloat();
		  nextToken( );

		  
		  //ANGLE
		  checkToken( "angle", "Lights" );
		  nextToken( );
		  myDLight[3+4*i] = parseFloat();
		  nextToken( );

 		  checkToken( "attenuation", "Lights" );
		  nextToken( );
		  light[6+9*i] = parseFloat();
		  nextToken( );
 		  light[7+9*i] = parseFloat();
		  nextToken( );
		  light[8+9*i] = parseFloat();
		  nextToken( );

		  checkToken( "}", "Lights" );   	  
	   }
	   else if(compareToken( "PointLight")){
		  checkToken( "PointLight", "Lights" ); 
		  
		  myLightType[i] = 1;

 		  nextToken( );
		  checkToken( "{", "Lights" );

		  nextToken( );
   	
 		  checkToken( "position", "Lights" );
		  nextToken( );
		  light[0+9*i] = parseFloat();
		  nextToken( );
 		  light[1+9*i] = parseFloat();
		  nextToken( );
		  light[2+9*i] = parseFloat();
		  nextToken( );

		  checkToken( "color", "Lights" );
		  nextToken( );
		  light[3+9*i] = parseFloat();
		  nextToken( );
 		  light[4+9*i] = parseFloat();
		  nextToken( );
		  light[5+9*i] = parseFloat();
		  nextToken( );

 		  checkToken( "attenuation", "Lights" );
		  nextToken( );
		  light[6+9*i] = parseFloat();
		  nextToken( );
 		  light[7+9*i] = parseFloat();
		  nextToken( );
		  light[8+9*i] = parseFloat();
		  nextToken( );

		  checkToken( "}", "Lights" );
	   }
	   else {
		 checkToken( "PointLight", "Lights" ); 
	   }

	}
    nextToken( );
	checkToken( "}", "Lights" );
}
void Scene::parseBackground( ){

	nextToken( );
	checkToken( "Background", "Background" );
	nextToken( );
	checkToken( "{", "Background" );
	
	nextToken( );
	checkToken( "color", "Background" );
	for( int i = 0; i < 3; i++ ){
		nextToken( );
		myBackgroundColor[i] = parseFloat( );
	}

	nextToken( );
	checkToken( "}", "Background" );
}

void Scene::parseMaterials( ){

   	nextToken( );
	checkToken( "Materials", "Materials" );
	nextToken( );
	checkToken( "{", "Materials" );
	nextToken( );
	checkToken( "numMaterials", "Materials" );
	nextToken( );
	myNumberOfMaterials = parseInt();
	//Make Private in class with accessor Function
    mats = new float [myNumberOfMaterials*7];
	

	for( int i = 0; i < myNumberOfMaterials; i++ ){
	  nextToken( );
	  checkToken( "PhongMaterial", "Materials" );
	  nextToken( );
	  checkToken( "{", "Materials" );

	  nextToken( );
	  checkToken( "diffuseColor", "Materials" );
	  nextToken( );
	  mats[0+i*7] = parseFloat();
	  nextToken( );
	  mats[1+i*7] = parseFloat();
	  nextToken( );
	  mats[2+i*7] = parseFloat();

	  nextToken( );
	  checkToken( "specularColor", "Materials" );
	  nextToken( );
	  mats[3+i*7] = parseFloat();
	  nextToken( );
	  mats[4+i*7] = parseFloat();
	  nextToken( );
	  mats[5+i*7] = parseFloat();

	  nextToken( );
	  checkToken( "exponent", "Materials" );
	  nextToken( );
	  mats[6+i*7] = parseFloat();

	  nextToken( );
	  checkToken( "}", "Materials" );

	}
	nextToken( );
	checkToken( "}", "Materials" );

}

void Scene::parseGroup( ){

    int currentMaterial = 0;
   	nextToken( );
	checkToken( "Group", "Group" );
	nextToken( );
	checkToken( "{", "Group" );
	nextToken( );
	checkToken( "numObjects", "Group" );
	nextToken( );
	myNumObjects = parseInt();
	spheres = new float [myNumObjects*5];
	planes = new float [myNumObjects*5];
	triangles = new float [myNumObjects*10];
	
	
	myMaterialIndex = new int [myNumObjects];

	for( int i = 0; i < myNumObjects; i++) {
	   nextToken( );
	   if(compareToken( "MaterialIndex")){
		  nextToken( );
		  currentMaterial = parseInt();
		  //myMaterialIndex[i] = currentMaterial;
		  nextToken();

		  if(compareToken( "Sphere")){
			 checkToken( "Sphere", "Group" );
			 nextToken();
			 checkToken( "{", "Group" );
			 nextToken();
			 checkToken( "center", "Group" );
			 nextToken();
			 spheres[0+5*myNumSpheres] = currentMaterial;
			 spheres[1+5*myNumSpheres] = parseFloat();
			 nextToken();
			 spheres[2+5*myNumSpheres] = parseFloat();
			 nextToken();
			 spheres[3+5*myNumSpheres] = parseFloat();
			 nextToken();
			 checkToken( "radius", "Group" );
			 nextToken();
			 spheres[4+5*myNumSpheres] = parseFloat();
			 nextToken( );
			 checkToken( "}", "Group" );
			 myNumSpheres +=1;
		  }
		  else if(compareToken( "Plane")){
			 checkToken( "Plane", "Group" );
			 nextToken();
			 checkToken( "{", "Group" );
			 nextToken();
			 checkToken( "normal", "Group" );
			 nextToken();
			 planes[0+5*myNumPlanes] = currentMaterial;
			 planes[1+5*myNumPlanes] = parseFloat();
			 nextToken();
			 planes[2+5*myNumPlanes] = parseFloat();
			 nextToken();
			 planes[3+5*myNumPlanes] = parseFloat();
			 nextToken();
			 checkToken( "offset", "Group" );
			 nextToken();
			 planes[4+5*myNumPlanes] = parseFloat();
			 nextToken( );
			 checkToken( "}", "Group" );
			 myNumPlanes +=1;
		  }
		  else if(compareToken( "Triangle")){
			 checkToken( "Triangle", "Group" );
			 nextToken();
			 checkToken( "{", "Group" );
			 nextToken();
			 checkToken( "vertex0", "Group" );
			 nextToken();
			 triangles[0+10*myNumTriangles] = currentMaterial;
			 triangles[1+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[2+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[3+10*myNumTriangles] = parseFloat();
			 nextToken();

			 checkToken( "vertex1", "Group" );
			 nextToken();
			 triangles[4+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[5+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[6+10*myNumTriangles] = parseFloat();
			 nextToken();

			 checkToken( "vertex2", "Group" );
			 nextToken();
			 triangles[7+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[8+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[9+10*myNumTriangles] = parseFloat();
			 nextToken();
			 checkToken( "}", "Group" );
			 myNumTriangles +=1;
		  }
		  else if(compareToken( "TriangleMesh")){
			 checkToken( "TriangleMesh", "Group" );
			 nextToken();
			 checkToken( "{", "Group" );
			 nextToken();
			 checkToken( "ply_file", "Group" );
			 nextToken();
			 mesh = parseString();
			 nextToken();
			 checkToken( "}", "Group" );
			 myMaterialIndex[myNumMeshes] = currentMaterial;
			 myNumMeshes +=1;
		  }
		  else {
			 exit(1);
		  }
	   }
	   else {
		  //myMaterialIndex[i] = currentMaterial;
		  if(compareToken( "Sphere")){
			 checkToken( "Sphere", "Group" );
			 nextToken();
			 checkToken( "{", "Group" );
			 nextToken();
			 checkToken( "center", "Group" );
			 nextToken();
			 spheres[0+5*myNumSpheres] = currentMaterial;
			 spheres[1+5*myNumSpheres] = parseFloat();
			 nextToken();
			 spheres[2+5*myNumSpheres] = parseFloat();
			 nextToken();
			 spheres[3+5*myNumSpheres] = parseFloat();
			 nextToken();
			 checkToken( "radius", "Group" );
			 nextToken();
			 spheres[4+5*myNumSpheres] = parseFloat();
			 nextToken( );
			 checkToken( "}", "Group" );
			 myNumSpheres +=1;
		  }
		  else if(compareToken( "Plane")){
			 checkToken( "Plane", "Group" );
			 nextToken();
			 checkToken( "{", "Group" );
			 nextToken();
			 checkToken( "normal", "Group" );
			 nextToken();
			 planes[0+5*myNumPlanes] = currentMaterial;
			 planes[1+5*myNumPlanes] = parseFloat();
			 nextToken();
			 planes[2+5*myNumPlanes] = parseFloat();
			 nextToken();
			 planes[3+5*myNumPlanes] = parseFloat();
			 nextToken();
			 checkToken( "offset", "Group" );
			 nextToken();
			 planes[4+5*myNumPlanes] = parseFloat();
			 nextToken( );
			 checkToken( "}", "Group" );
			 myNumPlanes +=1;
		  }
		  else if(compareToken( "Triangle")){
			 checkToken( "Triangle", "Group" );
			 nextToken();
			 checkToken( "{", "Group" );
			 nextToken();
			 checkToken( "vertex0", "Group" );
			 nextToken();
			 triangles[0+10*myNumTriangles] = currentMaterial;
			 triangles[1+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[2+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[3+10*myNumTriangles] = parseFloat();
			 nextToken();

			 checkToken( "vertex1", "Group" );
			 nextToken();
			 triangles[4+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[5+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[6+10*myNumTriangles] = parseFloat();
			 nextToken();

			 checkToken( "vertex2", "Group" );
			 nextToken();
			 triangles[7+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[8+10*myNumTriangles] = parseFloat();
			 nextToken();
			 triangles[9+10*myNumTriangles] = parseFloat();
			 nextToken();
			 checkToken( "}", "Group" );
			 myNumTriangles +=1;
		  }
		  else if(compareToken( "TriangleMesh")){
			 checkToken( "TriangleMesh", "Group" );
			 nextToken();
			 checkToken( "{", "Group" );
			 nextToken();
			 checkToken( "ply_file", "Group" );
			 nextToken();
			 mesh = parseString();
			 nextToken();
			 checkToken( "}", "Group" );
			 myMaterialIndex[myNumMeshes] = currentMaterial;
			 myNumMeshes +=1;
		  }
		  else {
			 exit(1);
		  }
	   }

	}
}


bool Scene::parse( ){	
	bool ret = true;
	lineNumber = 0;
	tokenCount = 0;
	
	inputFileStream.open( myInputSceneFile.c_str( ), ios::in );
	if( inputFileStream.fail( ) ){
		cerr << "Error opening \"" << myInputSceneFile << "\" for reading." << endl;
		exit( 1 );
	}
	parseCamera( );
	parseLight();
	parseBackground( );
	parseMaterials( );
	parseGroup( );

	inputFileStream.close( );
	
	return( ret );
}

bool Scene::areMoreTokens( ){
	bool ret = false;
	if( j < length ){
		ret = true;
	}
	return( ret );
}

void Scene::advance( ){
	if( currentLine[j] == ' ' || currentLine[j] == '\t' || currentLine[j] == '\n' ){
		while( currentLine[j] == ' ' || currentLine[j] == '\t' || currentLine[j] == '\n' ){
			j++;
		}
		i = j;
	}
}
void Scene::nextOnLine( ){
	//advance( );
	while( currentLine[j] != ' ' && currentLine[j] != '\t' && currentLine[j] != '\n' && currentLine[j] != 0 ){
		j++;
	}
	//cout << "ending: " << i <<  ", " << j << endl;
	currentLine[j] = 0;
	int tmp = i;
	if( i != j ){
		while( i <= j ){
			currentToken[i - tmp] = currentLine[i];
			//cout << "copying: " << (i - tmp) <<  ", " << i << endl;
			i++;
		}
		//cerr << lineNumber << ": " << ++tokenCount << ": '" << currentToken << "'" << endl;
	}
	j++;
	i = j;
}

void Scene::nextToken( ){
	if( !inputFileStream.eof( ) ){
		advance( );
		if( areMoreTokens( ) ){
			nextOnLine( );
		}else{
			do{
				inputFileStream.getline( currentLine, sizeof(currentLine) );
				lineNumber++;
				length = strlen( currentLine );
				//cerr << "new line of length: " << length << endl;
			}while( length <= 0 );
			i = 0;
			j = 0;
			advance( );
			//cerr << "Line: " << currentLine << endl;
			if( areMoreTokens( ) ){
				nextOnLine( );
			}
		}
	}
}


void Scene::write( std::ostream &out ) const {
	out << "Input scene file: " << myInputSceneFile << endl;
	out << "Output file: " << myOutputFile << endl;
	out << "Depth file: " << myDepthFile << endl;
	out << "Resolution: " << myResolution << " x " << myResolution << endl;
	out << "Camera:" << endl;
	if(myCamType == 0){
	  out << "   Cam Type = OrthographicCamera" <<endl;
	}
	else if(myCamType == 1){
	  out << "   Cam Type = PerspectiveCamera" <<endl;
	}
	else if(myCamType == 2){
	  out << "   Cam Type = SimplePerspectiveCamera" <<endl;
	}
	out << "   Cam Center = " << myCamCenter[0]<<", "<<myCamCenter[1]<<", " <<myCamCenter[2]<<endl;
	out << "   Cam Direction = " << myCamDirection[0]<<", "<<myCamDirection[1]<<", " <<myCamDirection[2]<<endl;
	out << "   Cam Up = " << myCamUp[0]<<", "<<myCamUp[1]<<", " <<myCamUp[2]<<endl;
	if(myCamType == 1){
	  out << "   Cam Angle = " << myCamAngle <<endl;
	}
	else if(myCamType == 2){
	  out << "   Cam Distance = " << myCamDistance <<endl;
	}
   	
	out << "Lights:" << endl;
	for( int i = 0; i < myNumLights; i++ ){
	   out << "   Light #" << i;
	   if(myLightType[i]) {
		   out << " of type PointLight = " << endl;
		   out << "      Position: " << light[0+i*9]<<", "<<light[1+i*9]<<", " <<light[2+i*9]<<endl;
		   out << "      Color: " << light[3+i*9]<<", "<<light[4+i*9]<<", " <<light[5+i*9]<<endl;
		   out << "      Attenuation: " << light[6+i*9]<<", "<<light[7+i*9]<<", " <<light[8+i*9]<<endl;
	   }
	   else if(!myLightType[i]) {
		   out << " of type DirectionalLight = " << endl;
		   out << "      Position: " << light[0+i*9]<<", "<<light[1+i*9]<<", " <<light[2+i*9]<<endl;
		   out << "      Direction: " << myDLight[0+i*4]<<", "<<myDLight[1+i*4]<<", " <<myDLight[2+i*4]<<endl;
		   out << "      Color: " << light[3+i*9]<<", "<<light[4+i*9]<<", " <<light[5+i*9]<<endl;
		   out << "      Angle: " << myDLight[3+i*4]<<endl;
		   out << "      Attenuation: " << light[6+i*9]<<", "<<light[7+i*9]<<", " <<light[8+i*9]<<endl;

	   }
	}

	out << "Background Color: " << myBackgroundColor[0] << ", " << myBackgroundColor[1] << ", " <<myBackgroundColor[2] << endl;
	out << "Number of Materials: " << myNumberOfMaterials << endl;
	out << "Materials:" << endl;
    for( int i = 0; i < myNumberOfMaterials; i++ ){
	   out << "   Material #" << i << " = " << endl;
	   out << "      DiffuseColor: " << mats[0+i*7]<<", "<<mats[1+i*7]<<", " <<mats[2+i*7]<<endl;
	   out << "      SpecularColor: " << mats[3+i*7]<<", "<<mats[4+i*7]<<", " <<mats[5+i*7]<<endl;
	   out << "      Exponent: " << mats[6+i*7]<<endl;
	}
	out << "Group:" << endl;

	for( int i = 0; i < myNumSpheres; i++ ){
	   out << "   Sphere #" << i << " with MaterialIndex #" << spheres[0+i*5] << " = " <<endl;
	   out << "      Center: " << spheres[1+i*5]<<", "<<spheres[2+i*5]<<", " <<spheres[3+i*5]<<endl;
	   out << "      Radius: " << spheres[4+i*5]<< endl;
	}
	for( int i = 0; i < myNumPlanes; i++ ){
	   out << "   Plane #" << i << " with MaterialIndex #" << planes[0+i*5]<< " = " <<endl;
	   out << "      Normal: " << planes[1+i*5]<<", "<<planes[2+i*5]<<", " <<planes[3+i*5]<<endl;
	   out << "      Offset: " << planes[4+i*5]<< endl;
	}
	for( int i = 0; i < myNumTriangles; i++ ){
	   out << "   Triangle #" << i << " with MaterialIndex #" << triangles[0+i*10] << " = " <<endl;
	   out << "      Vertex0: " << triangles[1+i*10]<<", "<<triangles[2+i*10]<<", " <<triangles[3+i*10]<<endl;
   	   out << "      Vertex1: " << triangles[4+i*10]<<", "<<triangles[5+i*10]<<", " <<triangles[6+i*10]<<endl;
	   out << "      Vertex2: " << triangles[7+i*10]<<", "<<triangles[8+i*10]<<", " <<triangles[9+i*10]<<endl;
	}
	for( int i = 0; i < myNumMeshes; i++ ){
	   out << "   Mesh #" << i << " with MaterialIndex #" << myMaterialIndex[i] << " = " <<endl;
	   out << "      Ply_File: " << mesh << endl;
	}
}

std::ostream& operator <<( std::ostream &out, const Scene &s ){
	s.write( out );
	return( out );
}
