
#pragma once

// Adapted from:
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection.html

struct Ray
{
	Ray(Vector origin, Vector direction)
		: Origin(origin), Direction(direction)
	{
		InverseDirection =
		{
			1.0f / direction.x,
			1.0f / direction.y,
			1.0f / direction.z
		};

		Signs[0] = InverseDirection.x < 0.0f;
		Signs[1] = InverseDirection.y < 0.0f;
		Signs[2] = InverseDirection.z < 0.0f;
	}

	Vector Origin;
	Vector Direction;
	
	// Optimisation properties
	Vector InverseDirection;
	int Signs[3];
};
