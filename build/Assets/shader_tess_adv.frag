#version 450
#define max_length 128

in block
{
	vec3	pos;
	vec3	col;
	vec3	patch_coords;
	float	k;
	float	t;
} In;

out vec4 fs_out_col;

uniform vec3 pointData[max_length];
uniform int pointNum = 0;
uniform vec3 referencePoint = vec3(0,0,0);

#incl ../Assets/BezierEval.glsl
#incl ../Assets/GeomInvariant.glsl

void main()
{
	//
	// 1. compute initial guess
	//

	// assert: segment num is > 1
	const int TEST_NUM = 10;
	dvec3 min_point = pointData[0];
	double min_par = 0;
	double min_dist = length(pointData[0] - referencePoint);
	for (int i=1; i<pointNum-1; ++i)
	{
		for (int j=0; j<=TEST_NUM; ++j)
		{
			double tt = j/double(TEST_NUM);
			dvec3 pt = (1.0 - tt)*pointData[j] + tt*pointData[j+1];
			
			double dst = length(referencePoint - pt);

			if ( dst < min_dist )
			{
				min_point = pt;
				min_par = tt;
				min_dist = dst;
			}
		}
	}

	/*int MAX_STEP = 32;
	bool found = false;
	
	for(int i=0; i<MAX_STEP && !found;++i)
	{
	}*/

	fs_out_col = vec4(min_dist/10,0,0,1);

}
