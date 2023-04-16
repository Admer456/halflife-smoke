
#pragma once

enum class InteractiveSmokeState : int
{
	Forming = 0,
	Spreading,
	Dissipating,
	Inactive,
	Shot
};

struct InteractiveSmokeParticle final
{
	inline float Transparency(const float& time) const
	{
		switch (State)
		{
		case InteractiveSmokeState::Inactive:
			return 0.0f;
		case InteractiveSmokeState::Forming:
			return (time - SpawnTime) / FadeTime;
		case InteractiveSmokeState::Shot:
			return FadeTime * FadeTime;
		case InteractiveSmokeState::Dissipating:
			const float timeOfDeath = SpawnTime + Life;
			return (timeOfDeath - time) / FadeTime;
		}

		return 1.0f;
	}

	inline bool Active() const
	{
		return State != InteractiveSmokeState::Inactive;
	}

	inline void Update(const float& time, const float& deltaTime)
	{
		Position = Position + Velocity * deltaTime;
		const float timeOfDeath = SpawnTime + Life;
		const float alpha = Transparency(time);
		
		switch (State)
		{
		case InteractiveSmokeState::Forming:
			Velocity = Velocity - (Velocity * deltaTime * 0.2f);
			Radius += Radius * deltaTime * 0.03f;
			if (time - SpawnTime > FadeTime)
			{
				State = InteractiveSmokeState::Spreading;
			}
			break;

		case InteractiveSmokeState::Spreading:
			Velocity = Velocity - (Velocity * deltaTime * 0.1f);
			Radius += Radius * deltaTime * 0.025f;
			if ((timeOfDeath - time) <= FadeTime)
			{
				State = InteractiveSmokeState::Dissipating;
			}
			break;

		case InteractiveSmokeState::Dissipating:
			Radius += Radius * deltaTime * 0.06f;
			Velocity.z += (1.0f - alpha) * deltaTime * 20.0f;
			if (time >= timeOfDeath)
			{
				State = InteractiveSmokeState::Inactive;
			}
			break;

		case InteractiveSmokeState::Shot:
			Radius += deltaTime * 200.0f;
			FadeTime -= deltaTime * 5.0f;
			if (FadeTime < 0.0f)
			{
				State = InteractiveSmokeState::Inactive;
			}
		}
	}

	Vector					Position{};
	Vector					Velocity{};
	float					Radius{16.0f};
	float					SpawnTime{0.0f};
	float					Life{25.0f};
	float					FadeTime{5.0f};
	InteractiveSmokeState	State{InteractiveSmokeState::Forming};
};
