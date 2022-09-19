#pragma once
#ifndef MASSSPRING2D_INCLUDED
#define MASSSPRING2D_INCLUDED

#include <vector>

class MassSpring2D
{

	struct Force {
		float coord_[2];
	};

	struct Particle {
		float position_[3];
	};
	
	struct Spring {
		unsigned int originVertexId_;
		unsigned int destinationVertexId_;
		float stiffness_;
		float restLength_;
	};

public:
	MassSpring2D(void);
	~MassSpring2D(void);
	int addParticle(float x, float y, bool isFixed);
	int addSpring(unsigned int originVertexId, unsigned int destinationVertexId, float stiffness, float restLength);
	void addForce(unsigned int vertexId, float forceX, float forceY);
	bool setDamping(float damping);
	void doStep(float timestep);
	void drawSprings(void);
	void drawForces(void);
	void drawMass(void);

private:
	void create(void);

protected:
	std::vector<Particle> currentParticleArray_;
	std::vector<Particle> oldParticleArray_;
	std::vector<float> constraintArray_;
	std::vector<Spring> springArray_;
	std::vector< Force > forceArray_;
	std::vector<unsigned int> elementArray_;
	float damping_;
};

#endif // MASSSPRING2D_INCLUDED

