
#pragma once

struct InteractiveSmokeCloud
{
	adm::AABB BoundingBox{};
	std::vector<InteractiveSmokeParticle> Particles{};
	bool Active{false};
};
