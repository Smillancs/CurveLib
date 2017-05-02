#version 450
#define max_length 128

// pipeline-ból bejövõ per-fragment attribútumok
in vec3 vs_out_pos;

// kimenõ érték - a fragment színe
out vec4 fs_out_col;


uniform vec3 pointData[max_length];
uniform int pointNum = 0;
//uniform dvec2 referencePoint = vec2(2,2);

#incl ../Assets/BezierEval.glsl
#incl ../Assets/GeomInvariant.glsl

uniform dvec2 res = vec2(640.0,480.0);
uniform float zoom = 1.0;

			  
void main()
{//
	// 1. compute initial guess
	dvec2 referencePoint = (2.0*gl_FragCoord.xy - res.xy)/vec2(res.x,res.y);
	referencePoint = referencePoint * zoom * 2;

	const int TEST_NUM = pointNum;
	double min_par = 0;
	double min_dist = length(referencePoint - AllBernstein(pointData,pointNum,0).xy);
	for (int j=0; j<=TEST_NUM+1; ++j)
	{
		double tt = j/double(TEST_NUM);
		dvec3 pt = AllBernstein(pointData,pointNum,tt);
		
		double dst = length(referencePoint - pt.xy);

		if ( dst < min_dist )
		{
			min_par = tt;
			min_dist = dst;
		}
	}
		
	
	double xn = min_par;
	double xnm1 = -1;
	bool is_end = false;
	int iter = 0;

	double error = 0.0;
	
	while (!is_end)
	{
		dvec2 p = AllBernstein(pointData,pointNum,float(xn)).xy;
		dvec2 dp = Diff(pointData,pointNum,1,float(xn)).xy;


		double K = K(pointData,pointNum,float(xn));

		xnm1 = xn;
		
		if (abs(K) < 1e-8)
		{
			double origin = dot((referencePoint - p),dp);
			xn = xn + origin/(length(dp)*length(dp));
		}
		else
		{
			dvec2 tan = e(pointData, pointNum, float(xn)).xy;
			dvec2 norm = n(pointData, pointNum, float(xn)).xy;
			
			dvec2 osc_mid = p + norm/K;
			
			dvec2 diff = referencePoint - osc_mid;
			dvec2 diff_dir = normalize(diff);
			
			double cosa = dot(-norm,diff_dir);
			double sina = dot(tan,diff_dir);

			dvec2 ndiff = -cosa*norm + sina*tan;

			double alpha = atan(float(sina),float(cosa));

			double ds = alpha/K;

			//1. pullback
			//xn = xn + ds/length(dp);

			//Taylor, 3rd order
			dvec2 ddp = Diff(pointData,pointNum,2,float(xn)).xy;
			dvec2 dddp = Diff(pointData,pointNum,3,float(xn)).xy;
			
			double x1 = length(dp);
			double x2 = dot(ddp,dp)/x1;
			double x3 = dot(dddp,dp)/x1;

			xn = xn + ds/x1 - ds*ds/2.0*x2/pow(x1,3) - pow(ds,3)/6.0*(pow(K,2)*pow(x1,3) + x3*x1 - 3*pow(x2,2))/pow(x1,5);
			
		}
		
		if (abs(xn-xnm1) < 1e-8)
		{
			is_end = true;
		}
		else if (xn > 1)
		{
			xn = 1;
			is_end = true;
		}
		else if (xn < 0)
		{
			xn = 0;
			is_end = true;
		}
		else if (iter > 4)
			is_end = true;

		iter++;
	}
	
	double dst = length(referencePoint - AllBernstein(pointData,pointNum,float(xn)).xy)/2;
	fs_out_col = vec4(dst,dst,dst,1);
}