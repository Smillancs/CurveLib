#include "GeomOptimize.hpp"

#include "gBuffer.h"

#include <glm/gtx/rotate_vector.hpp>

std::vector<GeomOptimize::Result> GeomOptimize::optimize2D3(const std::vector<GeomOptimize::Input2D3>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
	// Turn on compute shader
	program.On();

	// Create and bind input buffer
	gBuffer inBuf(GL_SHADER_STORAGE_BUFFER, input.size() * sizeof(GeomOptimize::Input2D3), input.data(), GL_STATIC_DRAW);
	inBuf.BindBufferBase(0);

	// Create and bind output buffer
	gBuffer outBuf(GL_SHADER_STORAGE_BUFFER, input.size() * sizeof(GeomOptimize::Result), nullptr, GL_STATIC_READ);
	outBuf.BindBufferBase(1);

	// Create and bind debugging buffer
	gBuffer dump(GL_SHADER_STORAGE_BUFFER, 100 * sizeof(float), nullptr, GL_STATIC_READ);
	dump.BindBufferBase(2);

	// Start compute shader
	program.DispatchCompute(input.size(), 1, 1);
	program.MemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

	// Extract results
	GeomOptimize::Result* ptr = (GeomOptimize::Result*)outBuf.MapBufferRange(0, input.size() * sizeof(GeomOptimize::Result), GL_MAP_READ_BIT);

	std::vector<GeomOptimize::Result> res(ptr, ptr+input.size());

	outBuf.UnMap();

	// Extract debuginfo (only if needed)
	if(nullptr != debugInfo)
	{
		float* ptrd = (float*)dump.MapBufferRange(0, 100 * sizeof(float), GL_MAP_READ_BIT);
		debugInfo->resize(100);
		std::copy(ptrd, ptrd + 100, debugInfo->begin());

		dump.UnMap();
	}

	// Unbind buffers and shader
	program.BindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
	program.BindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, 0);
	program.BindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, 0);
	program.Off();

	return res;
}

GeomOptimize::GeomOptimize()
{
	program.AttachShader(GL_COMPUTE_SHADER, "../GPUcompute/optimize.glsl");
	program.LinkProgram();
}

BezierCurve GeomOptimize::createResultCurve(const GeomOptimize::Input2D3& input, const GeomOptimize::Result& res)
{

	std::vector<glm::vec3> cps;
	cps.push_back(glm::vec3(input.p0, 0));
	cps.push_back(glm::vec3(input.p0 + (res.t0 * glm::rotate(input.p1 - input.p0, input.alpha)) / 3.f, 0));
	cps.push_back(glm::vec3(input.p1 - (res.t1 * glm::rotate(input.p1 - input.p0, input.beta)) / 3.f, 0));
	cps.push_back(glm::vec3(input.p1, 0));
	BezierCurve bez(cps);
	return bez;
}
