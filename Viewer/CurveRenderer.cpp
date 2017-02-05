#include "CurveRenderer.h"

#include "../CurveLib/GeomInvariant.hpp"

void CurveRenderer::genBufferTesselation(const unsigned N, const float a, const float b)
{
	//vb.Clean(); // leads to unintended consequences
    this->a = a;
    this->b = b;

	vb.AddAttribute(0, 3); // position
	vb.AddAttribute(1, 3); // e
	vb.AddAttribute(2, 3); // n
	vb.AddAttribute(3, 3); // b
	vb.AddAttribute(4, 1); // k
	vb.AddAttribute(5, 1); // t

	float l = b - a;
	float t = 0;
	vb.AddData(0, glm::vec3(c.f(t)));
	vb.AddData(1, glm::vec3(GeomInv::e(c, t)));
	vb.AddData(2, glm::vec3(GeomInv::n(c, t)));
	vb.AddData(3, glm::vec3(GeomInv::b(c, t)));
	vb.AddData(4, float(GeomInv::K(c, t)));
	vb.AddData(5, float(GeomInv::T(c, t)));

	for (int i = 1; i < N; ++i)
	{
		t = a + l * i / (float)N;
		vb.AddData(0, glm::vec3(c.f(t)));
		vb.AddData(1, glm::vec3(GeomInv::e(c, t)));
		vb.AddData(2, glm::vec3(GeomInv::n(c, t)));
		vb.AddData(3, glm::vec3(GeomInv::b(c, t)));
		vb.AddData(4, float(GeomInv::K(c, t)));
		vb.AddData(5, float(GeomInv::T(c, t)));
		vb.AddData(0, glm::vec3(c.f(t)));
		vb.AddData(1, glm::vec3(GeomInv::e(c, t)));
		vb.AddData(2, glm::vec3(GeomInv::n(c, t)));
		vb.AddData(3, glm::vec3(GeomInv::b(c, t)));
		vb.AddData(4, float(GeomInv::K(c, t)));
		vb.AddData(5, float(GeomInv::T(c, t)));
	}
	t=1;
	vb.AddData(0, glm::vec3(c.f(t)));
	vb.AddData(1, glm::vec3(GeomInv::e(c, t)));
	vb.AddData(2, glm::vec3(GeomInv::n(c, t)));
	vb.AddData(3, glm::vec3(GeomInv::b(c, t)));
	vb.AddData(4, float(GeomInv::K(c, t)));
	vb.AddData(5, float(GeomInv::T(c, t)));
	vb.InitBuffers();

}

void CurveRenderer::genBufferNormal(const unsigned N, const float a, const float b)
{
	//vb.Clean(); // leads to unintended consequences
    this->a = a;
    this->b = b;

	vb.AddAttribute(0, 3); // position
	vb.AddAttribute(1, 3); // e
	vb.AddAttribute(2, 3); // n
	vb.AddAttribute(3, 3); // b
	vb.AddAttribute(4, 1); // k
	vb.AddAttribute(5, 1); // t

	float l = b - a;
	for (int i = 0; i <= N; ++i)
	{
        float t = a + l * i / (float)N;
		vb.AddData(0, glm::vec3(c.f(t)));
		vb.AddData(1, glm::vec3(GeomInv::e(c, t)));
		vb.AddData(2, glm::vec3(GeomInv::n(c, t)));
		vb.AddData(3, glm::vec3(GeomInv::b(c, t)));
		vb.AddData(4, float(GeomInv::K(c, t)));
		vb.AddData(5, float(GeomInv::T(c, t)));
	}

	vb.InitBuffers();
}

gVertexBuffer CurveRenderer::getBuffer()
{
	if (dirty)
	{
		genBufferNormal(n,a,b);
		dirty = false;
		return vb;
	}
	else return vb;
}
