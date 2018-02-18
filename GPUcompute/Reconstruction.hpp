#pragma once

#include "../CurveLib/Curve.hpp"
#include "../CurveLib/bezierCurve.h"
#include "../CurveLib/GeomInvariant.hpp"

#include "gShaderProgram.h"
#include "gBuffer.h"

#include <glm/gtx/rotate_vector.hpp>

#include <vector>
#include <memory>
#include <array>

struct ReconstructionDataBase
{
  const glm::vec3& p() const { throw("Unsupported curve reconstruction."); }
  const glm::vec3& e() const { throw("Unsupported curve reconstruction."); }
  const glm::vec3& n() const { throw("Unsupported curve reconstruction."); }
  const float& K() const { throw("Unsupported curve reconstruction."); }
  const float& T() const { throw("Unsupported curve reconstruction."); }
  const float& dK(unsigned n) const { throw("Unsupported curve reconstruction."); }
  const float& dT(unsigned n) const { throw("Unsupported curve reconstruction."); }
};

struct ReconstructionData0 : ReconstructionDataBase
{
  ReconstructionData0(glm::vec3 a = glm::vec3(0)):_p(a){}
  const glm::vec3& p() const { return _p; }
  glm::vec3& p()  { return _p; }
protected:
  glm::vec3 _p; // position
  float __pad;
};

struct ReconstructionData1 : ReconstructionData0
{
  ReconstructionData1(glm::vec3 a = glm::vec3(0), glm::vec3 b = glm::vec3(0)):ReconstructionData0(a),_e(b){}
  const glm::vec3& e() const { return _e; }
  glm::vec3& e()  { return _e; }
protected:
  glm::vec3 _e; // tangent
  float ___pad;
};

struct ReconstructionData2 : ReconstructionData1
{
  glm::vec3& n()  { return _n; }
  const glm::vec3& n() const { return _n; }
  float& K()  { return k; }
  const float& K() const { return k; }
protected:
  glm::vec3 _n; // normal
  float k; // curvature
};

template<int N>
struct ReconstructionData : ReconstructionData2
{
  float& K()  { return k; }
  const float& K() const { return k; }
  float& T()  { return t[0]; }
  const float& T() const { return t[0]; }
  float& dK(unsigned n)  { return n==0 ? k : dk[n-1]; }
  const float& dK(unsigned n) const { return n==0 ? k : dk[n-1]; }
  float& dT(unsigned n)  { return t[n]; }
  const float& dT(unsigned n) const { return t[n]; }

protected:
  std::array<float, N-2> dk; // differentiates of curvature
  std::array<float, N-2> t;  // torsion and its differentiates
  float ____pad;
  float _____pad;
};

template<>
struct ReconstructionData<0> : ReconstructionData0{};
template<>
struct ReconstructionData<1> : ReconstructionData1{};
template<>
struct ReconstructionData<2> : ReconstructionData2{};

ReconstructionData<1> getPointData1(Curve::Ptr c, double t)
{
  ReconstructionData<1> data;
  data.p() = c->f(t);
  data.e() = GeomInv::e(*c, t);
  return data;
}

ReconstructionData<2> getPointData2(Curve::Ptr c, double t)
{
  ReconstructionData<2> data;
  data.p() = c->f(t);
  data.e() = GeomInv::e(*c, t);
  data.n() = GeomInv::n(*c, t);
  data.K() = GeomInv::K(*c, t);
  return data;
}

ReconstructionData<3> getPointData3(Curve::Ptr c, double t)
{
  ReconstructionData<3> data;
  data.p() = c->f(t);
  data.e() = GeomInv::e(*c, t);
  data.n() = GeomInv::n(*c, t);
  data.K() = GeomInv::K(*c, t);
  data.dK(1) = GeomInv::dK(*c, t);
  data.T() = GeomInv::T(*c, t);
  return data;
}

template <int continuity, int extra_points = 0>
class Reconstruction
{
public:
	using Input = ReconstructionData<continuity>;

	using Result = std::pair<std::array<std::pair<glm::vec3,float>,(2*continuity+1)+1+extra_points>, std::array<float,4>>;
  using Result_cpu = std::pair<std::array<glm::vec3,(2*continuity+1)+1+extra_points>,float>;

  std::vector<Result> optimize(const std::string& targetFunction, const std::vector<Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo = 0);
	std::vector<Result_cpu> optimize_cpu(std::function<double(Curve&,double)> func, const std::vector<Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo = 0);
	std::vector<Result_cpu> optimize_cpu_alt(std::function<double(Curve&,double)> func, const std::vector<Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo = 0);

	Reconstruction(bool infnorm);

  Curve::Ptr createResultCurve(const Result& res);
	Curve::Ptr createResultCurve(const Result_cpu& res);

protected:
	gShaderProgram program;
  bool infnorm = false;
};

template <int continuity, int extra_points>
std::vector<typename Reconstruction<continuity,extra_points>::Result> Reconstruction<continuity,extra_points>::optimize(const std::string& targetFunction, const std::vector<typename Reconstruction<continuity,extra_points>::Input>& input, const std::shared_ptr<std::vector<float>>& debugInfo)
{
	// Turn on compute shader
	program.On();

  // Set target function
  program.SetSubroutine(GL_COMPUTE_SHADER, "Eval", targetFunction.c_str());
  program.SetUniform("infnorm", infnorm);

	// Create and bind input buffer
	gBuffer inBuf(GL_SHADER_STORAGE_BUFFER, input.size() * sizeof(Reconstruction<continuity,extra_points>::Input), input.data(), GL_STATIC_DRAW);
	inBuf.BindBufferBase(0);

	// Create and bind output buffer
	gBuffer outBuf(GL_SHADER_STORAGE_BUFFER, input.size() / 2 * sizeof(Reconstruction<continuity,extra_points>::Result), nullptr, GL_STATIC_READ);
	outBuf.BindBufferBase(1);

	// Create and bind debugging buffer
	gBuffer dump(GL_SHADER_STORAGE_BUFFER, 100 * sizeof(float), nullptr, GL_STATIC_READ);
	dump.BindBufferBase(2);

  GLuint queries[1];
  GLuint timeElapsed = 0;
  // Create a query object.
  glGenQueries(1, queries);
  // Query current timestamp 1
  glBeginQuery(GL_TIME_ELAPSED, queries[0]);
	// Start compute shader
	program.DispatchCompute(input.size()/2, 1, 1);
	program.MemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
  glEndQuery(GL_TIME_ELAPSED);
  // See how much time the rendering of object i took in nanoseconds.
  glGetQueryObjectuiv(queries[0], GL_QUERY_RESULT, &timeElapsed);
  std::cerr << "Elapsed time (optimization): " << timeElapsed/1000000.0 << " ms" << std::endl;

	// Extract results
	typename Reconstruction<continuity,extra_points>::Result* ptr = (Reconstruction<continuity,extra_points>::Result*)outBuf.MapBufferRange(0, input.size() / 2 * sizeof(Reconstruction<continuity,extra_points>::Result), GL_MAP_READ_BIT);

	std::vector<Reconstruction<continuity,extra_points>::Result> res(ptr, ptr+input.size()/2);

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

template <int continuity, int extra_points>
Reconstruction<continuity,extra_points>::Reconstruction(bool infnorm)
{
  std::string filename = "../Assets/reconstruction"+std::to_string(continuity)+(extra_points>0?"_"+std::to_string(extra_points):"")+".glsl";
  program.SetVerbose(true);
	/*program.AttachShader(GL_COMPUTE_SHADER, filename.c_str());
	program.LinkProgram();*/
  this->infnorm = infnorm;
}

template <int continuity, int extra_points>
Curve::Ptr Reconstruction<continuity,extra_points>::createResultCurve(const typename Reconstruction<continuity,extra_points>::Result& res)
{
  std::array<glm::vec3,(2*continuity+1)+1+extra_points> arr;
  std::transform(res.first.begin(), res.first.end(), arr.begin(), [](std::pair<glm::vec3,float> a){return a.first;});
	std::vector<glm::vec3> cps(arr.begin(), arr.end());
	return Curve::Ptr(new BezierCurve(cps));
}

#include "Reconstruction_CPU.inl"
