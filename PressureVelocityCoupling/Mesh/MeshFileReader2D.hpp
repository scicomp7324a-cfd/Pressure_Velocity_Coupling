# pragma once

# include<iostream>
# include <fstream>
# include <string>
# include <vector>
# include <array>
# include <sstream>
# include <map>
# include "Patch.hpp"

using namespace std ;

/* *******************************************************************
    MeshFileReader2D.hpp - This class contains all the functions to read
    the mesh files. It provides access to the mesh file data.
******************************************************************* */

class MeshFileReader2D{

    private:

    vector<array<double,2>> points_ ;
    vector<vector<int>> faces_ ;
    vector<vector<int>> cells_ ;
    vector<Patch> boundaryPatches_ ;

    int numPoints_ ;
    int numFaces_ ;
    int numCells_ ;
    int numBoundaryPatches_ ;

    public:

    // Constructs the mesh file reader object using th epoints, faces, cells and boundary patches file.
    MeshFileReader2D(vector<string> fileNames){

        readPointsFile(fileNames[0]) ;
        readFacesFile(fileNames[1]) ;
        readCellsFile(fileNames[2]) ;
        readBoundaryPatchesFile(fileNames[3]) ;
        
    }

    // Constant access to points vector
    const vector<array<double,2>> points() const{
        return points_ ;
    }

    // Constant access to faces vector
    const vector<vector<int>> faces() const{
        return faces_ ;
    }

    // Constant access to cells vector
    const vector<vector<int>> cells() const{
        return cells_ ;
    }

    // Constant access to boundary patches vector
    const vector<Patch> boundaryPatches() const{
        return boundaryPatches_ ;
    }

    // Returns the number of points in the mesh
    int numPoints() const{
        return numPoints_ ;
    }

    // Returns the number of faces in the mesh
    int numFaces() const{
        return numFaces_ ;
    }

    // Returns the number of cells in the mesh
    int numCells() const{
        return numCells_ ;
    }

    // Returns the number of boundary patches in the mesh
    int numBoundaryPatches() const{
        return numBoundaryPatches_ ;
    }

    // Utility function to split a string using white spaces
    vector<string> splitWhiteSpace(string s) {
        istringstream iss(s);
        vector<string> out;
        for (string tok; iss >> tok; ) out.push_back(tok);
        return out;
    }

    // Utility function to strip white space at the beginning and end of the string
    string stripWhiteSpace(string s){
        auto first = s.find_first_not_of(" \t\r");
        if (first == string::npos) return ""; 
        auto last  = s.find_last_not_of(" \t\r");
        return s.substr(first, last - first + 1);
    }

    // Function to read the points file
    void readPointsFile(string pointsFilePath){

        ifstream in(pointsFilePath);
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        vector<string> strVec ;
        while(getline(in, line)){

            if (line.empty()) continue ;
            if(line.find("(") != string::npos && line.find(")") != string::npos){
                line = line.substr(line.find("(")+1, line.find(")")-1) ; 
                strVec = splitWhiteSpace(line) ;
                points_.push_back({stod(strVec[0]), stod(strVec[1])}) ;

            }
            else if((line.find("(") == string::npos && line.find(")") != string::npos) ||
                    (line.find("(") != string::npos && line.find(")") == string::npos)){
                        continue ;
            }
            else{
                line = stripWhiteSpace(line) ;
                numPoints_ = stoi(line) ;
            }

        }

    }

    // Function to read the faces file
    void readFacesFile(string facesFilePatch){

        ifstream in(facesFilePatch);
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        vector<string> strVec ;
        while(getline(in, line)){

            if (line.empty()) continue ;
            if(line.find("(") != string::npos && line.find(")") != string::npos){

                int size = stoi(line.substr(0, line.find("("))) ;
                vector<int> arr ;
                line = line.substr(line.find("(")+1, line.find(")")-1) ; 
                strVec = splitWhiteSpace(line) ;
                
                for(int i = 0; i < size; ++i){
                    arr.push_back(stoi(strVec[i])) ;
                }

                faces_.push_back(arr) ;

            }
            else if((line.find("(") == string::npos && line.find(")") != string::npos) ||
                    (line.find("(") != string::npos && line.find(")") == string::npos)){
                        continue ;
            }
            else{
                line = stripWhiteSpace(line) ;
                numFaces_ = stoi(line) ;
            }

        }
    }

    // Function to read the cells file
    void readCellsFile(string cellsFilePath){

        ifstream in(cellsFilePath);
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        vector<string> strVec ;
        while(getline(in, line)){

            if (line.empty()) continue ;
            if(line.find("(") != string::npos && line.find(")") != string::npos){

                int size = stoi(line.substr(0, line.find("("))) ;
                vector<int> arr ;
                line = line.substr(line.find("(")+1, line.find(")")-2) ; 
                strVec = splitWhiteSpace(line) ;

                for(int i = 0; i < size; ++i){
                    arr.push_back(stoi(strVec[i])) ;
                }
                cells_.push_back(arr) ;

            }
            else if((line.find("(") == string::npos && line.find(")") != string::npos) ||
                    (line.find("(") != string::npos && line.find(")") == string::npos)){
                        continue ;
            }
            else{
                line = stripWhiteSpace(line) ;
                numCells_ = stoi(line) ;
            }
        }
    }

    // Function to reda the boundary patches file
    void readBoundaryPatchesFile(string boundaryPatchesFilePath){

        ifstream in(boundaryPatchesFilePath);
        if (!in) {
            cerr << "Failed to open file\n";
        }

        string line ;
        vector<string> strVec ;

        getline(in, line) ;
        line = stripWhiteSpace(line) ;
        numBoundaryPatches_ = stoi(line) ;

        while(getline(in,line)){
            if(line.find("movingWall")!=string::npos || line.find("fixedWall")!=string::npos ||
            line.find("frontAndBack")!=string::npos || line.find("North")!=string::npos ||
            line.find("South")!=string::npos || line.find("East")!=string::npos ||
            line.find("West")!=string::npos){

                string patchType = stripWhiteSpace(line) ;
                getline(in, line) ; 
                int size = stoi(line) ;

                while(getline(in,line)){

                    if(line.find("(")!=string::npos){
                        continue ;
                    }
                    else if(line.find(")")!=string::npos){
                        break ;
                    }
                    else{
                        vector<int> arr ;
                        strVec = splitWhiteSpace(line) ;
                        for(int i = 0; i < size; ++i){
                            arr.push_back(stoi(strVec[i])) ;
                        }
                        boundaryPatches_.push_back(Patch(patchType, boundaryPatches_.size(), arr, size)) ;
                    }
                }
            }
        }
    }


    // Test function 
    void test(){

        for(size_t i=0; i<points_.size(); ++i){
            cout << "Point: {" << points_[i][0] << "," << points_[i][1] << "}" << endl ;
        }

        for(vector<int> var: faces_){
            cout << "Face: " ;
            for(int val: var){
                cout << val << " " ;
            }
            cout << endl ;
        }

        for(vector<int> var: cells_){
            cout << "Cell: " ;
            for(int val: var){
                cout << val << " " ;
            }
            cout << endl ;
        }

        for(Patch var: boundaryPatches_){
            cout << "Boundary Patch: " << endl ;
            cout << "Name: " << var.name() << endl ;
            cout << "PatchIndex: " << var.index() << endl ;
            cout << "Size: " << var.numFaces() << endl ;
            cout << "Faces: " ;
            for(int f: var.faces()){
                cout << f << " " ;
            }
            cout << endl ;
        }

    }

} ;
