# pragma once 

#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <map>

using namespace std ;
namespace fs = filesystem ;

class PathConfig{

    private:
        fs::path baseDir_ ;
        fs::path outputDir_ ;

        using MemFn = fs::path (PathConfig::*)() const;
        map<int, MemFn> meshDirectoryMap ;
        map<int, MemFn> outputDirectoryMap ;

    public:
        PathConfig()
            :baseDir_("Files"), outputDir_("Output"){
                createMeshDirectoryMap() ;
                createOutputDirectoryMap() ;
            }

        PathConfig(const fs::path& baseDir, const fs::path& outputDir)
            :baseDir_(baseDir), outputDir_(outputDir){
                createMeshDirectoryMap() ;
                createOutputDirectoryMap() ;
            }

        // Base and Output directories 
        const fs::path& baseDir() const {return baseDir_;}
        const fs::path& outputDir() const {return outputDir_;}

        // Sub folders
        fs::path meshDir() const{return baseDir_ / "MeshFiles" ; } 
        fs::path fieldDir() const{return baseDir_ / "FieldFiles" ; }
        
        // Dimension based sub folder 
        fs::path dim10Dir() const{return meshDir() / "10x10" ;}
        fs::path dim20Dir() const{return meshDir() / "20x20" ;}
        fs::path dim30Dir() const{return meshDir() / "30x30" ;}
        fs::path dim40Dir() const{return meshDir() / "40x40" ;}
        fs::path dim50Dir() const{return meshDir() / "50x50" ;}
        fs::path dim60Dir() const{return meshDir() / "60x60" ;}
        fs::path dim70Dir() const{return meshDir() / "70x70" ;}
        fs::path dim80Dir() const{return meshDir() / "80x80" ;}
        fs::path dim90Dir() const{return meshDir() / "90x90" ;}
        fs::path dim100Dir() const{return meshDir() / "100x100" ;}

        // Dimension Based output sub folder 
        fs::path outputDim10Dir() const{return outputDir() / "10x10" ;}
        fs::path outputDim20Dir() const{return outputDir() / "20x20" ;}
        fs::path outputDim30Dir() const{return outputDir() / "30x30" ;}
        fs::path outputDim40Dir() const{return outputDir() / "40x40" ;}
        fs::path outputDim50Dir() const{return outputDir() / "50x50" ;}
        fs::path outputDim60Dir() const{return outputDir() / "60x60" ;}
        fs::path outputDim70Dir() const{return outputDir() / "70x70" ;}
        fs::path outputDim80Dir() const{return outputDir() / "80x80" ;}
        fs::path outputDim90Dir() const{return outputDir() / "90x90" ;}
        fs::path outputDim100Dir() const{return outputDir() / "100x100" ;}


        void createMeshDirectoryMap(){
            
            meshDirectoryMap.emplace(10, &PathConfig::dim10Dir) ;
            meshDirectoryMap.emplace(20, &PathConfig::dim20Dir) ;
            meshDirectoryMap.emplace(30, &PathConfig::dim30Dir) ;
            meshDirectoryMap.emplace(40, &PathConfig::dim40Dir) ;
            meshDirectoryMap.emplace(50, &PathConfig::dim50Dir) ;
            meshDirectoryMap.emplace(60, &PathConfig::dim60Dir) ;
            meshDirectoryMap.emplace(70, &PathConfig::dim70Dir) ;
            meshDirectoryMap.emplace(80, &PathConfig::dim80Dir) ;
            meshDirectoryMap.emplace(90, &PathConfig::dim90Dir) ;
            meshDirectoryMap.emplace(100, &PathConfig::dim100Dir) ;

        }

        void createOutputDirectoryMap(){

            outputDirectoryMap.emplace(10, &PathConfig::outputDim10Dir) ;
            outputDirectoryMap.emplace(20, &PathConfig::outputDim20Dir) ;
            outputDirectoryMap.emplace(30, &PathConfig::outputDim30Dir) ;
            outputDirectoryMap.emplace(40, &PathConfig::outputDim40Dir) ;
            outputDirectoryMap.emplace(50, &PathConfig::outputDim50Dir) ;
            outputDirectoryMap.emplace(60, &PathConfig::outputDim60Dir) ;
            outputDirectoryMap.emplace(70, &PathConfig::outputDim70Dir) ;
            outputDirectoryMap.emplace(80, &PathConfig::outputDim80Dir) ;
            outputDirectoryMap.emplace(90, &PathConfig::outputDim90Dir) ;
            outputDirectoryMap.emplace(100, &PathConfig::outputDim100Dir) ;
    
        }

        fs::path meshDirFor(int key) const {
            auto it = meshDirectoryMap.find(key);
            if (it == meshDirectoryMap.end()) {
                throw std::runtime_error("No mesh directory mapping for key = " + std::to_string(key));
            }
            return (this->*(it->second))(); 
        }

        fs::path outputDirFor(int key) const {
            auto it = outputDirectoryMap.find(key);
            if (it == outputDirectoryMap.end()) {
                throw std::runtime_error("No output directory mapping for key = " + std::to_string(key));
            }
            return (this->*(it->second))(); 
        } 
        
        const fs::path getOutputDirectory(int dimension, int re) const{
            fs::path mainDirectory = outputDirFor(dimension) ;
            return mainDirectory / to_string(re) ;
        }

        vector<string> getMeshFilenames(int dimension){
            fs::path meshPath = meshDirFor(dimension) ;
            vector<string> filenames = {meshPath / "points.txt",
                                        meshPath / "faces.txt",
                                        meshPath / "cells.txt",
                                        meshPath / "boundary.txt"} ;
            return filenames ;
        }

        // Field files based on name 
        fs::path UxFieldFile() const{ return fieldDir() / "Ux.txt" ;}
        fs::path UyFieldFile() const{ return fieldDir() / "Uy.txt" ;}
        fs::path PFieldFile() const{ return fieldDir() / "P.txt" ;}

        void validate() const {
            if (!fs::exists(baseDir_)) {
                throw std::runtime_error("Base directory not found: " + baseDir_.string());
            }
            if (!fs::exists(meshDir())) {
                throw std::runtime_error("Mesh directory not found: " + meshDir().string());
            }
            if (!fs::exists(fieldDir())) {
                throw std::runtime_error("Field directory not found: " + fieldDir().string());
            }
        }

        void ensureOutputDir() const {
            fs::create_directories(outputDir_);
        }

        void createOutputDirectory(int dimension, int re) const {
            fs::path dir = (*this).getOutputDirectory(dimension, re) ;
            fs::create_directories(dir) ;
        }

} ;

