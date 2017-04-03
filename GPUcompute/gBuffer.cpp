#include "gBuffer.h"

gBuffer::gBuffer(GLenum _type, GLuint _size_in_bytes, const void* _src, GLenum _usage)
{
	// assert: our context is active
	Define(_type, _size_in_bytes, _src, _usage);
}

gBuffer::~gBuffer()
{
	// assert: our context is active
	glDeleteBuffers(1, &m_buffer_id);
}

void gBuffer::Define(GLenum _type, GLuint _size_in_bytes, const void* _src, GLenum _usage)
{
	// assert: our context is active
	if (m_buffer_id != 0)
		glDeleteBuffers(1, &m_buffer_id);

	glGenBuffers(1, &m_buffer_id);
	glBindBuffer(_type, m_buffer_id);
	glBufferData(_type, _size_in_bytes, _src, _usage);

	m_buffer_type = _type;
	m_size_in_bytes = _size_in_bytes;

}

void* gBuffer::MapBufferRange(GLintptr _offset, GLsizeiptr _length, GLbitfield _access) const
{
	// assert: our context is active

	// ideally, you should always UnMap - but in order to avoid unnecessary OpenGL hog, I'm offing it manually
	if (m_is_mapped)
	{
		std::cout << "[warning][gBuffer] Buffer was not UnMap()'d - UnMap() called automatically. Make sure you know what you are doing!" << std::endl;
		UnMap();
		m_is_mapped = false;
	}

	glBindBuffer(m_buffer_type, m_buffer_id);
	auto address = glMapBufferRange(m_buffer_type, _offset, _length, _access);
	if (nullptr != address)
		m_is_mapped = true;

	auto error = glGetError();
	switch (error)
	{
	case GL_INVALID_ENUM:
		std::cout << "GL_INVALID_ENUM\n";
		break;
	case GL_INVALID_OPERATION:
		std::cout << "GL_INVALID_OPERATION\n";
		break;
	case GL_INVALID_VALUE:
		std::cout << "GL_INVALID_VALUE \n";
		break;
	}
	return address;
}

void gBuffer::UnMap() const
{
	// assert: our context is active
	glUnmapBuffer(m_buffer_type);
	m_is_mapped = false;
}

void gBuffer::BindBufferBase(GLuint _index) const
{
	// assert: our context is active
	glBindBufferBase(m_buffer_type, _index, m_buffer_id);
}

