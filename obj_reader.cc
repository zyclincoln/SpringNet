#include "obj_reader.h"

Obj_Reader::Obj_Reader(const char* path){
	fin.open(path);
}

void Obj_Reader::sparse_model(){

}

void Obj_Reader::sparse_constraint(const char *path){
	ifstream constraint_file;
	constraint_file.open(path);

	

	constraint_file.close();
}

void Obj_Reader::sparse_points(){

}

void Obj_Raeder::sparse_edge(){

}