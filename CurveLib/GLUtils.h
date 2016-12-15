#pragma once
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <GL/glew.h>


/* 

Az http://www.opengl-tutorial.org/ oldal alapján.

*/
GLuint loadShader(GLenum _shaderType, const char* _fileName);

GLuint loadProgramVSGSFS(const char* _fileNameVS, const char* _fileNameGS, const char* _fileNameFS);