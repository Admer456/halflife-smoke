
#pragma once

class SmokeManager final
{
public:
	static void Init();
	static void Restart();
	static void Shutdown();

	static void Update(const float& time, const float& deltaTime);
	static void Render(triangleapi_s* tri, const float& time);

	static void SpawnCloud(Vector position, float radius = 128.0f, float particleSize = 16.0f);
	static void TraceBullet(Vector start, Vector end, float force);
};
