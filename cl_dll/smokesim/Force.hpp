
#pragma once

struct Force final
{
	inline bool Active(const float& time) const
	{
		return StartTime + Life > time;
	}

	Vector Position{};
	float Intensity{0.0f};
	float StartTime{0.0f};
	float Life{0.2f};
};
