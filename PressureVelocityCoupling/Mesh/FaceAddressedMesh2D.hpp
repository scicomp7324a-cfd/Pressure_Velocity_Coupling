# pragma once

# include<iostream>
# include <string>
# include <vector>
# include <array>
# include <set>
# include <map>
# include <ostream>
# include "MeshFileReader2D.hpp"
# include "Patch.hpp"
# include <algorithm>

using namespace std ;

/* *******************************************************************
    FaceAddressedMesh2D.hpp - This class contains all the implementation 
    of the face addressed cartesian mesh. It provides the functions 
    to access the geometrical and topological data.
******************************************************************* */


class FaceAddressedMesh2D{

    private:

        const MeshFileReader2D mfr ;
        const vector<array<double,2>> points_ ;
        const vector<vector<int>> faces_ ;
        const vector<vector<int>> cells_ ;
        const vector<Patch> boundaryPatches_ ;

        const int numPoints_ ;
        const int numFaces_ ;
        const int numCells_ ;
        const int numBoundaryPatches_ ;

        vector<array<int,2>> faceOwnerNeighbour_ ;
        vector<int> internalFaces_ ;
        vector<int> boundaryFaces_ ; 
        

    public:

        // Constructs the mesh using the points, faces, cells and boundary patches files
        FaceAddressedMesh2D(vector<string> fileNames): mfr(fileNames),
                                                       points_(mfr.points()),
                                                       faces_(mfr.faces()),
                                                       cells_(mfr.cells()),
                                                       boundaryPatches_(mfr.boundaryPatches()),
                                                       numPoints_(mfr.numPoints()),
                                                       numFaces_(mfr.numFaces()),
                                                       numCells_(mfr.numCells()),
                                                       numBoundaryPatches_(mfr.numBoundaryPatches())
                                                       {
            createFaceOwnerNeighbourArray() ;
            createInternalAndBoundaryFaceList() ;
            //testFunctions() ;
        }

        // Utility functions 

        // Function creates the face owner neighbour array that provides the owner and neighbour of the face using the face index
        void createFaceOwnerNeighbourArray(){

            int l = -1;
            array<int,2> arr ;
            for(size_t i=0; i < faces_.size(); ++i){

                for(size_t j=0; j < cells_.size(); ++j){
                    if(find(cells_[j].begin(), cells_[j].end(), i)!=cells_[j].end()){
                        arr[++l] = j ;
                    }
                }

                if(l==0){
                    arr[++l] = -1 ;
                }

                faceOwnerNeighbour_.push_back(arr) ;
                l = -1 ;
            }
        }

        // Function calculates the area of a triangle using the vertices points
        double calculateAreaOfTriangle(array<double,2> p1, array<double,2> p2, array<double,2> p3) const{
            double area = (p1[0]*(p2[1]-p3[1]) + p2[0]*(p3[1]-p1[1]) + p3[0]*(p1[1]-p2[1]))/2 ;
            return abs(area) ;
        }
           
        // Function returns the normal vector to a line with points on it
        const array<double,2> returnNormalVector(array<double,2> p1, array<double,2> p2) const{
            array<double,2> v = {p2[0]-p1[0], p2[1]-p1[1]} ;
            double mod = sqrt(v[0]*v[0] + v[1]*v[1]) ;
            return {-v[1]/mod, v[0]/mod} ;
        }

        // Function creates the lists that contain the internal and boundary faces 
        void createInternalAndBoundaryFaceList(){

            array<int,2> arr;
            for(size_t i=0; i<faces_.size(); ++i){
                arr = faceOwnerNeighbour_[i] ;
                if(arr[1] != -1){
                    internalFaces_.push_back(i) ;
                }
                else{
                    boundaryFaces_.push_back(i) ;
                }
            }

        }

        // Global data

        // Returns the number of points in the mesh
        const int& numPoints() const{
            return numPoints_ ;
        }

        // Returns the number of faces in the mesh
        const int& numFaces() const{
            return numFaces_ ;
        }

        // Returns the number of cells in the mesh 
        const int& numCells() const{
            return numCells_ ;
        }

        // Returns the number of boundary patches in the mesh
        const int& numBoundaryPatches() const{
            return numBoundaryPatches_ ;
        }

        // Geometrical data

        // Returns the cell centre given the cell index
        const array<double,2> cellCentre(int cellIndex) const{

            set<int> pointSet ;
            for(size_t i = 0; i < cells_[cellIndex].size(); ++i){
                for(int j: faces_[cells_[cellIndex][i]]){
                    pointSet.insert(j) ;
                }
            }

            int numCellPoints = pointSet.size() ;
            double xCentre=0, yCentre=0; 

            for(int point: pointSet){
                xCentre += points_[point][0] ;
                yCentre += points_[point][1] ;;
            }

            return {xCentre/numCellPoints, yCentre/numCellPoints} ;
        }

        // Returns the face centre given the face index
        const array<double,2> faceCentre(int faceIndex) const{

            int numFacePoints = faces_[faceIndex].size() ;
            double xCentre=0, yCentre=0; 

            for(int point: faces_[faceIndex]){
                xCentre += points_[point][0] ;
                yCentre += points_[point][1] ;
            }

            return {xCentre/numFacePoints, yCentre/numFacePoints} ;

        }

        // Returns the cell area given the cell index
        double cellArea(int cellIndex) const{

            array<double,2> cellCentre_ = cellCentre(cellIndex) ;

            double cellVolume = 0 ;

            for(int faceIndex: cells_[cellIndex]){
                vector<array<double,2>> cellPoints ;
                for(int pointIndex: faces_[faceIndex]){
                    cellPoints.push_back(points_[pointIndex]) ;
                }

                cellVolume += calculateAreaOfTriangle(cellCentre_, cellPoints[0], cellPoints[1]) ;
            }

            return cellVolume ;

        }

        // Returns the face area vector given the face index
        const array<double,2> faceVector(int faceIndex) const{

            array<double,2> p1 = points_[faces_[faceIndex][0]];
            array<double,2> p2 = points_[faces_[faceIndex][1]];

            array<double,2> p = {p2[0]-p1[0], p2[1]-p1[1]} ;
            double magnitude = sqrt(p[0]*p[0] + p[1]*p[1]) ;

            array<double,2> normal = returnNormalVector(p1, p2) ;
            array<double,2> cellCentre_ = cellCentre(faceOwnerNeighbour_[faceIndex][0]) ;
            array<double,2> faceCentre_ = faceCentre(faceIndex) ;
            array<double,2> centreToFace = {faceCentre_[0]-cellCentre_[0], faceCentre_[1]-cellCentre_[1]} ;

            if(normal[0]*centreToFace[0]+normal[1]*centreToFace[1] < 0){
                normal[0] = -1*normal[0] ;
                normal[1] = -1*normal[1] ;
            }
            
            normal[0] *= magnitude ;
            normal[1] *= magnitude ;
            return normal ;
        }
        
        // Returns the face list of a cell given the cell index
        const vector<int>& cellFaces(int cellIndex) const{
            return cells_[cellIndex] ;
        }

        // Topological Data

        // Prints the face and the owner given the patch index
        void boundaryFaceOwner(int patchIndex){
            cout << "Boundary Patch: " << patchIndex << endl ;
            const Patch p = boundaryPatches_[patchIndex] ;
            for(int faceIndex: p.faces()){
                    cout << "Face: " << faceIndex << ": Owner = " << faceOwnerNeighbour_[faceIndex][0] << endl ;
            }
        }

        // Returns the cell neighbours given the cell index
        const vector<int> cellNeighbours(int cellIndex) const{

            vector<int> neighbours ;
            for(int faceIndex: cells_[cellIndex]){
                if(faceOwnerNeighbour_[faceIndex][0]!=cellIndex){
                    neighbours.push_back(faceOwnerNeighbour_[faceIndex][0]) ;
                }
                else{
                    neighbours.push_back(faceOwnerNeighbour_[faceIndex][1]) ;
                }
            }

            return neighbours ;
        }

        // Returns the face owner and neighbour given the face index
        const array<int,2> faceOwnerNeighbour(int faceIndex) const{
            return {faceOwnerNeighbour_[faceIndex][0], faceOwnerNeighbour_[faceIndex][1]} ;
        }

        // Constant access to internal faces vector
        const vector<int>& internalFaces() const{
            return internalFaces_ ;
        }

        // Constant access to boundary face vector
        const vector<int>& boundaryFaces() const{
            return boundaryFaces_ ;
        }

        // Constant access to boundary patches
        const vector<Patch>& boundaryPatches() const{
            return boundaryPatches_ ;
        }

        // Test Function
        void testFunctions(){

            createFaceOwnerNeighbourArray() ;
            ofstream out("/Users/amitanshusahoo/Desktop/PVC/PVC_V2/Pressure_Velocity_Coupling/Output/meshTest.txt") ;

            for(size_t i = 0; i < faceOwnerNeighbour_.size(); ++i){
                out << i << ": owner = " << faceOwnerNeighbour_[i][0] << " neighbour = " << faceOwnerNeighbour_[i][1] << endl ;  
            }

            out << endl ;
            out << "Number of Points: " << numPoints() << endl ;
            out << "Number of Cells: " << numCells() << endl ;
            out << "Number of Faces: " << numFaces() << endl ;
            out << "Number of Boundary Patches: " << numBoundaryPatches() << endl ;

            out << endl ;
            array<double,2> cellCentre_ ;
            int numCells_ = numCells() ;
            for(int i=0; i < numCells_; ++i){
                cellCentre_ = cellCentre(i) ;
                out << "Cell Index: " << i << ": Cell Centre = {" << cellCentre_[0] << "," << cellCentre_[1] << "}" << endl ;
            }

            out << endl ;
            array<double,2> faceCentre_ ;
            int numFaces_ = numFaces() ;
            for(int i=0; i < numFaces_; ++i){
                faceCentre_ = faceCentre(i) ;
                out << "Face Index: " << i << ": Face Centre = {" << faceCentre_[0] << "," << faceCentre_[1] << "}" << endl ;
            }

            out << endl ;
            double cellVolume_ ;
            numCells_ = numCells() ;
            for(int i=0; i < numCells_; ++i){
                cellVolume_ = cellArea(i) ;
                out << "Cell Index: " << i << ": Cell Volume = " << cellVolume_ << endl ;
            }

            out << endl ;
            array<double,2> faceVec_ ;
            numFaces_ = numFaces() ;
            for(int i=0; i < numFaces_; ++i){
                faceVec_ = faceVector(i) ;
                out << "Face Index: " << i << ": Face Area Vector = {" << faceVec_[0] << "," << faceVec_[1] << "}" << endl ;
            }
            out << endl ;

            for(size_t i=0; i<boundaryPatches_.size(); ++i){
                boundaryFaceOwner(i) ;
                out << endl ;
            }

            out << endl ;
            vector<int> cellNeighbours_ ;
            numCells_ = numCells() ;
            for(int i=0; i < numCells_; ++i){
                cellNeighbours_ = cellNeighbours(i) ;
                out << "Cell Index: " << i << ": Cell Neighbours: " ;
                for(int var: cellNeighbours_){
                    out << var << " " ;
                }
                out << endl ;
            }

            out << endl ;
            array<int,2> faceON ;
            numFaces_ = numFaces() ;
            for(int i=0; i < numFaces_; ++i){
                faceON = faceOwnerNeighbour(i) ;
                out << "Face Index: " << i << ": Face {Owner,Neighbour} = {" << faceON[0] << "," << faceON[1] << "}" << endl ;
            }

        }

} ;
