#pragma once

#include <GL/glew.h>

#include <map>
#include <string>
#include <list>
#include <vector>

#include <glm/glm.hpp>

class gShaderProgram
{
public:
	gShaderProgram(void);
	~gShaderProgram(void);

	bool AttachShader(GLenum _shaderType, const char* _filename);
	void BindAttribLoc(int _index, const char* _uniform);
	void BindFragDataLoc(int _index, const char* _uniform);
	bool LinkProgram();

	void SetVerbose(bool);

	void On();
	void Off();

	void DispatchCompute(GLuint _x, GLuint _y, GLuint _z);
	void MemoryBarrier(GLbitfield);
	void BindBufferBase(GLenum _target, GLuint _index, GLuint _buffer);

	void Clean();

	void SetUniform(const char* _uniform, const glm::vec2& _vec);
	void SetUniform(const char* _uniform, const glm::vec3& _vec);
	void SetUniform(const char* _uniform, const glm::vec4& _vec);
	void SetUniform(const char* _uniform, const glm::mat4& _mat);
	void SetUniform(const char* _uniform, int _i);
	void SetUniform(const char* _uniform, float _f);
	void SetUniform(const char* _uniform, float _a, float _b, float _c, float _d);
	void SetUniform(const char* _uniform, float _a, float _b, float _c);
	void SetUniform(const char* _uniform, float _a, float _b);
	void SetTexture(const char* _uniform, int _sampler, GLuint _textureID);
	void SetCubeTexture(const char* _uniform, int _sampler, GLuint _textureID);
  void SetSubroutine(GLenum _shaderType, const char* _uniformName, const char* _subrName);

	GLuint ID() const;
	operator int() const { return ID(); }
protected:
	GLuint	getLocation(const char* _uniform);
  GLuint  getSubroutineUniformLocation(GLenum _shaderType, const char* _uniform);
  GLuint  getSubroutineIndex(GLenum _shaderType, const char* _name);
	GLuint	loadShader(GLenum _shaderType, const char* _fileName);

	GLuint							m_id_program;
	std::map< std::string, GLint >	m_map_uniform_locations;
	std::list< GLuint >				m_list_shaders_attached;

  std::vector< GLuint > subroutineIndices;

	bool	m_verbose;
};
