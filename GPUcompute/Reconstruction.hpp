#pragma once

#include "../CurveLib/Curve.hpp"
#include "../CurveLib/bezierCurve.h"
#include "../CurveLib/GeomInvariant.hpp"

#include "gShaderProgram.h"
#include "gBuffer.h"

#include <glm/gtx/rotate_vector.hpp>

#include <vector>
#include <memory>

template<int N>
struct ReconstructionData;

template <>
struct ReconstructionData<0>
{
  ReconstructionData(glm::vec3 a = glm::vec3(0)):p(a){}
  glm::vec3 p; // position
  float __pad;
};

template<>
struct ReconstructionData<1> : ReconstructionData<0>
{
  ReconstructionData(glm::vec3 a = glm::vec3(0), glm::vec3 b = glm::vec3(0)):ReconstructionData<0>(a),e(b){}
  glm::vec3 e; // tangent
  float ___pad;
};

template<>
struct ReconstructionData<2> : ReconstructionData<1>
{
  glm::vec3 n; // normal
  float& K(){ return k; }
private:
  float k; // curvature
};

template<int N>
struct ReconstructionData : ReconstructionData<2>
{
  float& K(){ return k; }
  float& T(){ return t[0]; }
  float& dK(unsigned n){ return n==0 ? k : dk[n-1]; }
  float& dT(unsigned n){ return t[n]; }

private:
  std::array<float, N-2> dk; // differentiates of curvature
  std::array<float, N-2> t;  // torsion and its differentiates
};

template <int N>
ReconstructionData<N> getPointData(Curve::Ptr c, double t)
{
  ReconstructionData<N> data;
  data.p = c->f(t);
  if(N >= 1)
  {
    data.e = GeomInv::e(*c, t);
    /*if(N >= 2)
    {
      data.n = GeomInv::n(*c, t);
      data.K() = GeomInv::K(*c, t);
      if(N >= 3)
      {
        data.dK(1) = GeomInv::dK(*c, t);
        data.T() = GeomInv::T(*c, t);
        if(N >= 4)
        {
          data.dK(2) = GeomInv::ddK(*c, t);
          data.dT(1) = GeomInv::dT(*c, t);
        }
      }
    }*/
  }
  return data;
}

template <int continuity>
class Reconstruction
{
public:
	using Input = ReconstructionData<continuity>;

	using Result = std::pair<std::array<std::pair<glm::vec3,float>,(2*continuity+1)+1>, float>;

	std::vector<Result> optimize(const std::string& targetFunction, const std::vector<Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo = 0);

	Reconstruction();

	Curve::Ptr createResultCurve(const Result& res);

protected:
	gShaderProgram program;
};

template <int continuity>
std::vector<typename Reconstruction<continuity>::Result> Reconstruction<continuity>::optimize(const std::string& targetFunction, const std::vector<typename Reconstruction<continuity>::Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
	// Turn on compute shader
	program.On();

  // Set target function
  program.SetSubroutine(GL_COMPUTE_SHADER, "Eval", targetFunction.c_str());

	// Create and bind input buffer
	gBuffer inBuf(GL_SHADER_STORAGE_BUFFER, input.size() * sizeof(Reconstruction<continuity>::Input), input.data(), GL_STATIC_DRAW);
	inBuf.BindBufferBase(0);

	// Create and bind output buffer
	gBuffer outBuf(GL_SHADER_STORAGE_BUFFER, input.size() / 2 * sizeof(Reconstruction<continuity>::Result), nullptr, GL_STATIC_READ);
	outBuf.BindBufferBase(1);

	// Create and bind debugging buffer
	gBuffer dump(GL_SHADER_STORAGE_BUFFER, 100 * sizeof(float), nullptr, GL_STATIC_READ);
	dump.BindBufferBase(2);

	// Start compute shader
	program.DispatchCompute(input.size()/2, 1, 1);
	program.MemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

	// Extract results
	typename Reconstruction<continuity>::Result* ptr = (Reconstruction<continuity>::Result*)outBuf.MapBufferRange(0, input.size() / 2 * sizeof(Reconstruction<continuity>::Result), GL_MAP_READ_BIT);

	std::vector<Reconstruction<continuity>::Result> res(ptr, ptr+input.size()/2);

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

template <int continuity>
Reconstruction<continuity>::Reconstruction()
{
  std::string filename = "../Assets/reconstruction"+std::to_string(continuity)+".glsl";
  program.SetVerbose(true);
	program.AttachShader(GL_COMPUTE_SHADER, filename.c_str());
	program.LinkProgram();
}

template <int continuity>
Curve::Ptr Reconstruction<continuity>::createResultCurve(const Reconstruction<continuity>::Result& res)
{
  std::array<glm::vec3,(2*continuity+1)+1> arr;
  std::transform(res.first.begin(), res.first.end(), arr.begin(), [](std::pair<glm::vec3,float> a){return a.first;});
	std::vector<glm::vec3> cps(arr.begin(), arr.end());
	return Curve::Ptr(new BezierCurve(cps));
}
