
#include "hud.h"
#include "cl_util.h"
#include "triangleapi.h"
#include "r_efx.h"

#include "CFrustum.h"
#include "particleman.h"
extern IParticleMan* g_pParticleMan;

#include <list>
#include <vector>

#include "Ray.hpp"
#include "AABB.hpp"
#include "SmokeManager.hpp"
#include "InteractiveSmokeParticle.hpp"
#include "InteractiveSmokeCloud.hpp"
#include "Force.hpp"

#include "particleman_internal.h"

#include <chrono>
namespace chrono = std::chrono;

static Vector RandomVector()
{
	return Vector(
		gEngfuncs.pfnRandomFloat(-1.0f, 1.0f),
		gEngfuncs.pfnRandomFloat(-1.0f, 1.0f),
		gEngfuncs.pfnRandomFloat(-1.0f, 1.0f));
}

class SmokeManagerInternal final
{
private:
	static inline bool CanRun = false;

public:
	static void Init()
	{
		HSPRITE smokeSprite = SPR_Load("sprites/smokesim/smoke.spr");
		if (smokeSprite <= 0)
		{
			gEngfuncs.Con_Printf("[SmokeManager] Cannot find 'sprites/smokesim/smoke.spr'\n");
			gEngfuncs.Con_Printf("[SmokeManager] Without that sprite, smoke cannot be rendered.\n");
			return;
		}

		SmokeTexture = gEngfuncs.GetSpritePointer(smokeSprite);
		CanRun = true;
	}

	static void Restart()
	{
		Shutdown();
		Init();
	}

	static void Shutdown()
	{
		Clouds.clear();
		CanRun = false;
	}

	static void Update(const float& time, const float& deltaTime)
	{
		if (Clouds.empty() || !CanRun)
		{
			return;
		}

		auto tpStart = chrono::system_clock::now();

		// Destroy forces
		if (!Forces.empty())
		{
			std::vector<Force> newForces(Forces.size());
			for (auto& force : Forces)
			{
				if (force.Active(time))
				{
					newForces.push_back(force);
				}
			}

			if (newForces.empty())
			{
				Forces.clear();
			}
			else
			{
				Forces = std::move(newForces);
			}
		}

		for (auto cloud = Clouds.begin(); cloud != Clouds.end(); cloud++)
		{
			if (!cloud->Active)
			{
				continue;
			}

			bool cloudIsActive = false;
			for (auto& particle : cloud->Particles)
			{
				if (particle.Active())
				{
					//for (auto& force : Forces)
					//{
					//	if (!force.Active(time))
					//	{
					//		continue;
					//	}
					//
					//	const Vector forceDirection = (particle.Position - force.Position);
					//	const float forceDistance2 = std::max(force.Intensity, forceDirection.LengthSquared());
					//	if (forceDistance2 > (force.Intensity * force.Intensity))
					//	{
					//		continue;
					//	}
					//
					//	const float fade = 1.0f - ((time - force.StartTime) / force.Life);
					//	particle.Velocity = particle.Velocity + deltaTime * fade * (forceDirection / forceDistance2) * force.Intensity;
					//	particle.Velocity.x += deltaTime * 60.0f * (fade / forceDistance2) * (std::sin(particle.Position.x * 0.05f) - std::cos(particle.Position.y * 0.067f));
					//	particle.Velocity.y += deltaTime * 60.0f * (fade / forceDistance2) * (std::sin(particle.Position.x * -0.03f) + std::cos(particle.Position.y * 0.044f));
					//}

					cloudIsActive = true;
					particle.Update(time, deltaTime);
				}
			}

			// If all particles have ceased,
			// we may remove the entire cloud.
			if (!cloudIsActive)
			{
				cloud->Particles.clear();
			}
		}

		auto tpEnd = chrono::system_clock::now();
		float microseconds = chrono::duration_cast<chrono::nanoseconds>(tpEnd - tpStart).count() * 0.001f;
		gEngfuncs.Con_Printf("SmokeManager.Update: %4.2f us\n", microseconds);
	}

#if 0
	struct ParticleRenderCache
	{
		ParticleRenderCache(const InteractiveSmokeParticle& particle, const float& time)
		{
			Position = particle.Position;
			AlphaEncoded = particle.Transparency(time) * SHRT_MAX;
			// This means a float of 1.0f would become a short of 1023
			// Gotta love fixed floats
			RadiusEncoded = particle.Radius * (SHRT_MAX / 32);
		}

		inline float Alpha() const
		{
			return float(AlphaEncoded) / SHRT_MAX;
		}

		inline float Radius() const
		{
			return float(RadiusEncoded) / (SHRT_MAX / 32);
		}

		Vector Position;
		short AlphaEncoded;
		short RadiusEncoded;
	};

	constexpr static size_t SizeOfRenderCache = sizeof(ParticleRenderCache);
	constexpr static size_t AlignOfRenderCache = alignof(ParticleRenderCache);

	static_assert(sizeof(ParticleRenderCache) == SizeOfRenderCache, "ParticleRenderCache must have a size of 16 bytes");

	static inline std::vector<ParticleRenderCache> RenderCache;
#endif

	// Updated every frame
	static inline Vector ViewForward;
	static inline Vector ViewRight;
	static inline Vector ViewUp;

	static inline const model_s* SmokeTexture{nullptr};

	static void Render(triangleapi_s* tri, const float& time)
	{
		if (!CanRun)
		{
			return;
		}

		auto tpStart = chrono::system_clock::now();

		gEngfuncs.pfnAngleVectors(gHUD.m_vecAngles, ViewForward, ViewRight, ViewUp);
#if 0
		CFrustum viewFrustum;
		viewFrustum.CalculateFrustum();

		size_t numParticlesTotal = 0;
		for (auto& cloud : Clouds)
		{
			numParticlesTotal += cloud.Particles.size();
		}

		// This ain't necessarily the total number of particles that will be rendered
		RenderCache.clear();
		RenderCache.reserve(numParticlesTotal);
#endif

		// Simple version for now: render everything unsorted
		for (auto& cloud : Clouds)
		{
			for (auto& particle : cloud.Particles)
			{
				RenderParticle(tri, particle.Position, particle.Transparency(time), particle.Radius);
			}
		}

		auto tpEnd = chrono::system_clock::now();
		float microseconds = chrono::duration_cast<chrono::nanoseconds>(tpEnd - tpStart).count() * 0.001f;
		if (gHUD.m_flTimeDelta == 0.0f)
		{
			gEngfuncs.Con_Printf("SmokeManager.Render: %4.2f us\n", microseconds);
		}
	}

	static void TraceBullet(Vector start, Vector end, float force)
	{
		auto tpStart = chrono::system_clock::now();

		// This function assumes that there is no world geometry between start and end
		// Therefore we'll not check for that here.
		const Ray ray = Ray(start, (end - start).Normalize());
		const float travelLimit = (end - start).Length();

		// The first step is to identify which clouds we may possibly intersect with.
		// Clouds are bounding boxes with smoke particles inside, so first we do a ray-AABB test.
		for (auto& cloud : Clouds)
		{
			float travelDistance;
			if (!IntersectsWith(cloud, ray, travelDistance) || travelDistance > travelLimit)
			{
				continue;
			}

			// Once we've intersected with a cloud, we gotta check which particles intersect with the ray.
			// This can be accelerated with an octree, or a grid, or a number of other things, but we'll
			// bruteforce this time because it's simple and we're not working with a huge amount of particles.
			int forceCounter = 0;
			for (auto& particle : cloud.Particles)
			{
				if (!particle.Active() || particle.State == InteractiveSmokeState::Shot)
				{
					continue;
				}

				// The travel limit can still get crossed, and not for the reasons you might think!
				// Imagine a case where a cloud is halfway in the ground. Without that check, it'd be
				// as if the bullet went through the ground and affected the particles below it.
				float travelDistance;
				if (!IntersectsWith(particle, ray, travelDistance) || travelDistance > travelLimit)
				{
					continue;
				}

				const Vector hitPosition = ray.Origin + ray.Direction * travelDistance;
				const Vector normal = (hitPosition - particle.Position) / particle.Radius; // avoid using Normalize here

				// Larger particles will move slower, i.e. be less affected by a bullet.
				// Later we may simulate turbulence and other stuff.
				particle.Velocity = particle.Velocity - normal * (force / particle.Radius) * 50.0f;
				particle.State = InteractiveSmokeState::Shot;
				particle.FadeTime = 1.0f;
			}
		}

		auto tpEnd = chrono::system_clock::now();
		float microseconds = chrono::duration_cast<chrono::nanoseconds>(tpEnd - tpStart).count() * 0.001f;
		gEngfuncs.Con_Printf("SmokeManager.TraceBullet: %4.2f us\n", microseconds);
	}

	static void SetupCloud(InteractiveSmokeCloud& cloud, Vector position, float radius, float particleSize)
	{
		// Double the radius for the bbox, because some particles may fly out of the initial bounds
		cloud.BoundingBox = adm::AABB(position, radius * 3.0f);
		cloud.Active = true;

		const float numParticlesPerRow = 2.0f * (radius / particleSize);
		// I roughly visually estimated that a sphere has 33% less volume than a cube.
		// Draw a square, then a circle in it. Observe the corners. Each corner part is roughly 66%
		// occupied by the sphere, leaving us with roughly 33% of unoccupied space.
		// Since I'd like to save some memory, I liberally decrease the radius of the sphere by another percent.
		const size_t numParticlesInSphere = std::pow(numParticlesPerRow, 3.0f) * 0.65f;
		cloud.Particles.resize(numParticlesInSphere);

		for (size_t i = 0U; i < numParticlesInSphere; i++)
		{
			auto& particle = cloud.Particles.emplace_back();
			// Start 30% within the sphere
			const Vector randomOffset = (RandomVector() + Vector(0.0f, 0.0f, 0.9f)).Normalize() * gEngfuncs.pfnRandomFloat(0.0f, 1.0f);
			particle.Position = position + randomOffset * radius * 0.35f;
			particle.Radius = particleSize * gEngfuncs.pfnRandomFloat(0.6f, 2.5f);
			particle.Velocity = randomOffset * 15.0f;
			particle.SpawnTime = gHUD.m_flTime;
			// Particles closer to the centre will stay alive for longer
			particle.Life = 30.0f + gEngfuncs.pfnRandomFloat(0.0f, 10.0f) - randomOffset.Length() * 10.0f;
		}

		gEngfuncs.Con_Printf("SmokeManager.SetupCloud: %d clouds total\n", (int)Clouds.size());
	}

	static inline std::vector<Force> Forces{};
	static inline std::list<InteractiveSmokeCloud> Clouds{};

private:
	// Adapted from:
	// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection.html
	static bool IntersectsWith(const InteractiveSmokeCloud& cloud, const Ray& ray, float& outDistance)
	{
		// My humble contribution
		if (cloud.BoundingBox.IsInside(ray.Origin))
		{
			outDistance = 0.0f;
			return true;
		}

		Vector bounds[2] = {cloud.BoundingBox.mins, cloud.BoundingBox.maxs};

		float tmin, tmax, tymin, tymax, tzmin, tzmax;

		tmin = (bounds[ray.Signs[0]].x - ray.Origin.x) * ray.InverseDirection.x;
		tmax = (bounds[1 - ray.Signs[0]].x - ray.Origin.x) * ray.InverseDirection.x;
		tymin = (bounds[ray.Signs[1]].y - ray.Origin.y) * ray.InverseDirection.y;
		tymax = (bounds[1 - ray.Signs[1]].y - ray.Origin.y) * ray.InverseDirection.y;

		if ((tmin > tymax) || (tymin > tmax))
			return false;

		if (tymin > tmin)
			tmin = tymin;
		if (tymax < tmax)
			tmax = tymax;

		tzmin = (bounds[ray.Signs[2]].z - ray.Origin.z) * ray.InverseDirection.z;
		tzmax = (bounds[1 - ray.Signs[2]].z - ray.Origin.z) * ray.InverseDirection.z;

		if ((tmin > tzmax) || (tzmin > tmax))
			return false;

		if (tzmin > tmin)
			tmin = tzmin;
		if (tzmax < tmax)
			tmax = tzmax;

		outDistance = tmin;
		if (outDistance < 0.0f)
		{
			outDistance = tmax;
			if (outDistance < 0.0f)
			{
				return false;
			}
		}

		return true;
	}

	static bool IntersectsWith(const InteractiveSmokeParticle& particle, const Ray& ray, float& outDistance)
	{
		const auto solveQuadratic = [](float a, float b, float c, float& x0, float& x1)
		{
			float discr = b * b - 4 * a * c;
			if (discr < 0)
			{
				return false;
			}
			else if (discr == 0)
			{
				x0 = x1 = -0.5 * b / a;
			}
			else
			{
				float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
				x0 = q / a;
				x1 = c / q;
			}

			if (x0 > x1)
			{
				std::swap(x0, x1);
			}

			return true;
		};

		float t0, t1; // solutions for t if the ray intersects

		const Vector L = ray.Origin - particle.Position;
		const float Radius2 = particle.Radius * particle.Radius * 0.25f; // Artificially decrease the radius for shooting purposes

		float a = DotProduct(ray.Direction, ray.Direction);
		float b = 2.0f * DotProduct(ray.Direction, L);
		float c = DotProduct(L, L) - Radius2;
		if (!solveQuadratic(a, b, c, t0, t1))
		{
			return false;
		}

		if (t0 > t1)
		{
			std::swap(t0, t1);
		}

		if (t0 < 0)
		{
			t0 = t1; // if t0 is negative, let's use t1 instead
			if (t0 < 0)
			{
				return false; // both t0 and t1 are negative
			}
		}

		outDistance = t0;
		return true;
	}

	static void RenderParticle(triangleapi_t* tri, const Vector& position, float alpha, float radius)
	{
		const Vector width = ViewRight * radius;
		const Vector height = ViewUp * radius;

		const Vector lowLeft = position - (width * 0.5) - (ViewUp * radius * 0.5);

		const Vector lowRight = lowLeft + width;
		const Vector topLeft = lowLeft + height;
		const Vector topRight = lowRight + height;

		const Vector BaseColour = Vector(1.0f, 1.0f, 1.0f) * 0.75f;
		Vector environmentLight;
		// This is gonna lag. LOTS
		tri->LightAtPoint(const_cast<Vector&>(position), environmentLight);
		environmentLight = environmentLight / 255.0f;

		tri->SpriteTexture(const_cast<model_s*>(SmokeTexture), 0);
		tri->RenderMode(kRenderTransAlpha);
		tri->CullFace(TRI_NONE);

		tri->Begin(TRI_QUADS);
		tri->Color4f(BaseColour.x * environmentLight.x, BaseColour.y * environmentLight.y, BaseColour.z * environmentLight.z, alpha);

		tri->TexCoord2f(0, 0);
		tri->Vertex3fv(topLeft);

		tri->TexCoord2f(0, 1);
		tri->Vertex3fv(lowLeft);

		tri->TexCoord2f(1, 1);
		tri->Vertex3fv(lowRight);

		tri->TexCoord2f(1, 0);
		tri->Vertex3fv(topRight);

		tri->End();

		tri->RenderMode(kRenderNormal);
		tri->CullFace(TRI_FRONT);
	}
};

void SmokeManager::Init()
{
	return SmokeManagerInternal::Init();
}

void SmokeManager::Restart()
{
	return SmokeManagerInternal::Restart();
}

void SmokeManager::Shutdown()
{
	return SmokeManagerInternal::Shutdown();
}

void SmokeManager::Update(const float& time, const float& deltaTime)
{
	return SmokeManagerInternal::Update(time, deltaTime);
}

void SmokeManager::Render(triangleapi_s* tri, const float& time)
{
	return SmokeManagerInternal::Render(tri, time);
}

void SmokeManager::SpawnCloud(Vector position, float radius, float particleSize)
{
	auto& cloud = SmokeManagerInternal::Clouds.emplace_back();
	SmokeManagerInternal::SetupCloud(cloud, position, radius, particleSize);
}

void SmokeManager::TraceBullet(Vector start, Vector end, float force)
{
	return SmokeManagerInternal::TraceBullet(start, end, force);
}
