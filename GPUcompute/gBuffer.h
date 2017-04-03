#pragma once

//#include "gContextManager.h"
#include <GL/glew.h>
#include <vector>
#include <iostream>

/*
	A convenience class that represents an OpenGL buffer - i.e. anything that can be manipulated
	using the standard GL Buffer commands (e.g. VBOs, IBOs, SSBOs, UBOs).
*/
class gBuffer
{
public:
	gBuffer(GLenum _type = GL_SHADER_STORAGE_BUFFER, GLuint _size_in_bytes = 0, const void* _src = nullptr, GLenum _usage = GL_STREAM_DRAW);
	~gBuffer();

	template <class T>
	gBuffer(const std::vector<T>& _src)
	{
		Define(GL_SHADER_STORAGE_BUFFER, (GLuint)_src.size()*sizeof(T), (void*)&(_src[0]), GL_STREAM_DRAW);
	}

	void* MapBufferRange(GLintptr _offset, GLsizeiptr _length, GLbitfield _access) const;
	void UnMap() const;
	void BindBufferBase(GLuint _index) const;
	void Define(GLenum _type = GL_SHADER_STORAGE_BUFFER, GLuint _size_in_bytes = 0, const void* _src = nullptr, GLenum _usage = GL_STREAM_DRAW);

	GLenum GetType() const { return m_buffer_type; }
	GLuint size() const { return m_size_in_bytes; }

	operator int() const { return m_buffer_id; }

	template <class T>
	operator T*() const {
		auto mapped = (T*)MapBufferRange(0, m_size_in_bytes, GL_MAP_READ_BIT);
		T* copy_of_mapped = new T[m_size_in_bytes / sizeof(T)];
		std::copy(mapped, mapped + m_size_in_bytes / sizeof(T), copy_of_mapped);
		UnMap();
		return copy_of_mapped;
	}

	template <class T>
	gBuffer& operator=(const std::vector<T>& _src)
	{
		if (m_is_mapped)
			UnMap();

		Define(GL_SHADER_STORAGE_BUFFER, (GLuint)_src.size()*sizeof(T), (void*)&(_src[0]), GL_STREAM_DRAW);
		return *this;
	}

	// this converts the GPU-based contents of the buffer to a std::vector<T>
	// UnMap is automatically called after this
	template <class T>
	operator std::vector<T>() 
	{
		T* contents = (T*)MapBufferRange(0, m_size_in_bytes, GL_MAP_READ_BIT | GL_MAP_WRITE_BIT);
		auto result = std::vector<T>(contents, contents + m_size_in_bytes/sizeof(T));
		UnMap();
		return result;
	}

private:

	GLuint			m_buffer_id = 0;
	GLuint			m_size_in_bytes = 0;
	GLenum			m_buffer_type	= GL_SHADER_STORAGE_BUFFER;
	mutable bool			m_is_mapped = false;
};

