/* Michael Robertson
 * mirob2005@gmail.com, miro2005@csu.fullerton.edu
 * CS 566
 * HW3
 * Due 3/23/2012
 * SID: 892-32-2629
 *
 * $Id: Scene.h 2764 2011-02-09 05:35:16Z mshafae $
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "Camera.h"
#include "Image.h"
#include "Material.h"
#include "Group.h"

#ifndef _SCENE_H_
#define _SCENE_H_

using namespace std;

class Scene{
	public:
		Scene( string inputFilename = "", string outputFilename = "", string depthFilename = "", int myResolution = 500 );
		~Scene( );

		// Accessors
		//Camera& camera( );
		float* backgroundColor( );
		int numberOfMaterials( );
		//void setCurrentMaterial( int i );
		//Material* currentMaterial( );
		//Group* group( );
		string& inputSceneFile( );
		string& outputFile( );
		string& depthFile( );
		string& meshFile();
		int resolution( );
		//ADDED

		//CAMERA
		int camType();
		float* camCenter();
		float* camDirection();
		float* camUp();
		float camAngle();
		float camDistance();

		int* materialIndex( );
		float* materials( );

		float* allSpheres( );
		float* allPlanes();
		float* allTriangles();
		float* allMeshes();

		int* lightType();
		float* lights( );
		float* dLight();
		int numLights();

		int numObjects();

		int numSpheres();
		int numPlanes();
		int numTriangles();
		int numMeshes();

		

		void setInputSceneFile( string file );
		void setOutputFile( string file );
		void setDepthFile( string file );

		bool hasInputSceneFilePath( );
		bool hasOutputFilePath( );
		bool hasDepthFilePath( );

		bool parse( );

		// I/O
		void write( std::ostream &out ) const;

	private:
		string myInputSceneFile;
		string myOutputFile;
		string myDepthFile;
		int myResolution;
		//Camera myCamera;
		float myBackgroundColor[3];
		int myNumberOfMaterials;
		//Material **materials;
		//Material *myCurrentMaterial;
		//Group *myGroup;
		ifstream inputFileStream;


		//ADDED

		//0 = Orthographic, 1 = perspective, 2 = simple perspective
		int myCamType;
		float myCamCenter[3];
		float myCamDirection[3];
		float myCamUp[3];
		float myCamAngle;
		float myCamDistance;

		int myNumObjects;
		int myNumSpheres;
		int myNumPlanes;
		int myNumTriangles;
		int myNumMeshes;

		int myNumLights;
		float* myDLight;
		//0 = DirectionalLight, 1 = PointLight
		int* myLightType;
		int* myMaterialIndex;
		float* mats;

		float* spheres;
		float* planes;
		float* triangles;
		float* meshes;
		string mesh;

		float* light;

		// For parsing
		char currentLine[255];
		char currentToken[255];
		int lineNumber;
		int tokenCount;
	    int length;
	    int i;
		int j;
		void nextToken( );
		void parseCamera( );
		void nextOnLine( );
		bool areMoreTokens( );
		void advance( );
		void checkToken( const char *str, const char *stage  );
		//ADDED
		bool compareToken( const char *str);
		void parseBackground( );
		void parseLight( );
		float parseFloat( );
		double parseDouble( );
		int parseInt( );
		string parseString( );

		void parseMaterials( );
		void parseGroup( );
};

std::ostream& operator <<( std::ostream &out, const Scene &s );


#endif