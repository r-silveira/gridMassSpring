#ifdef _WIN32
#include <windows.h>
#endif // _WIN32

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <cmath>

#include "Utils.h"
#include "MassSpring2D.h"
#include <iostream>
#include "Vec2.h"

using namespace std;



/**
*/
MassSpring2D::MassSpring2D(void)
{
	create();
}


/**
*/
MassSpring2D::~MassSpring2D(void)
{
}


/**
*/
int MassSpring2D::addParticle(float x, float y, bool isFixed)
{
	Particle particle;

	particle.position_[0] = x;
	particle.position_[1] = y;
	particle.position_[2] = 0.0f;

	Force force;
	force.coord_[0] = 0.0f;
	force.coord_[1] = 0.0f;

	currentParticleArray_.push_back(particle);
	oldParticleArray_.push_back(particle);
	forceArray_.push_back(force);

	if (isFixed)
		constraintArray_.push_back(0.0f);
	else
		constraintArray_.push_back(1.0f);

	return (int)currentParticleArray_.size();
}

/**
*/
int MassSpring2D::addSpring(unsigned int originVertexId, unsigned int destinationVertexId, float stiffness, float restLength)
{
	unsigned int particleCount;
	Spring spring;

	particleCount = (unsigned int)currentParticleArray_.size();

	if ((originVertexId >= particleCount) || (destinationVertexId >= particleCount))
	{
		return -1;
	}

	spring.originVertexId_ = originVertexId;
	spring.destinationVertexId_ = destinationVertexId;
	spring.stiffness_ = stiffness;
	spring.restLength_ = restLength;
		
	springArray_.push_back(spring);

	elementArray_.push_back(originVertexId);
	elementArray_.push_back(destinationVertexId);

	return (int)springArray_.size();
}



/**
*/
void MassSpring2D::addForce(unsigned int vertexId, float forceX, float forceY)
{
	unsigned int particleCount;

	particleCount = (unsigned int)currentParticleArray_.size();

	if ((vertexId <= particleCount))
	{
		forceArray_[vertexId].coord_[0] = forceX;
		forceArray_[vertexId].coord_[1] = forceY;
	}
}


/**
*/
bool MassSpring2D::setDamping(float damping)
{
	if (Utils::isFloatValid(damping))
	{
		damping_ = damping;

		return true;
	}

	return false;
}

/**
*/
void MassSpring2D::doStep(float timestep)
{	
	float* originParticlePosition;
	float* destinationParticlePosition;
	float	axis[2];

	const int particleCount = (int)currentParticleArray_.size();

	Particle* currentParticle	= &currentParticleArray_[0];
	Particle* oldParticle		= &oldParticleArray_[0];
	float* constraint			= &constraintArray_[0];

	for (int i=0; i != particleCount; ++i)
	{
		float forceX = forceArray_[i].coord_[0];
		float forceY = forceArray_[i].coord_[1];

		forceX *= timestep*timestep;
		forceY *= timestep*timestep;

		oldParticle->position_[0] = currentParticle->position_[0] + (damping_*(currentParticle->position_[0] - oldParticle->position_[0]) + forceX)*(*constraint);
		oldParticle->position_[1] = currentParticle->position_[1] + (damping_*(currentParticle->position_[1] - oldParticle->position_[1]) + forceY)*(*constraint);

		++currentParticle;
		++oldParticle;

		++constraint;
	}

	const int springCount = (int)springArray_.size();

	Spring* spring = &springArray_[0];

	for (int i=0; i != springCount; ++i)
	{
		float length = 0.0f;
		
		originParticlePosition = currentParticleArray_[spring->originVertexId_].position_;
		destinationParticlePosition = currentParticleArray_[spring->destinationVertexId_].position_;

		for (int j=0; j != 2; ++j)
		{
			const float component = originParticlePosition[j] - destinationParticlePosition[j];

			axis[j] = component*(spring->stiffness_*2.0f);
			
			length += component*component;
		}

		const float force = timestep*timestep*(1.0f - spring->restLength_*(1.0f/sqrtf(length)));
		
		originParticlePosition = oldParticleArray_[spring->originVertexId_].position_;
		destinationParticlePosition = oldParticleArray_[spring->destinationVertexId_].position_;

		for (int j=0; j!=2; ++j)
		{
			originParticlePosition[j] -= force*axis[j]*constraintArray_[spring->originVertexId_];
			destinationParticlePosition[j] += force*axis[j]*constraintArray_[spring->destinationVertexId_];
		}

		++spring;
    }
	
	oldParticleArray_.swap(currentParticleArray_);
}

/**
*/
void MassSpring2D::drawSprings(void)
{
	if (!elementArray_.empty())
	{
		glDepthMask(GL_FALSE);

		glVertexPointer(3, GL_FLOAT, sizeof(Particle), currentParticleArray_[0].position_);

		glEnableClientState(GL_VERTEX_ARRAY);
		for (int i=3; i>0; --i)
		{
			glLineWidth((float)i);

			glColor3f(1.0f - ((float)i/4.0f), 1.0f - ((float)i/4.0f), 0.0f);

			glDrawElements(GL_LINES, (int)elementArray_.size(), GL_UNSIGNED_INT, &elementArray_[0]);
		}
		glDisableClientState(GL_VERTEX_ARRAY);

		glDepthMask(GL_TRUE);
	}
}



/**
*/
void MassSpring2D::drawForces(void)
{
	glDepthMask(GL_FALSE);

	const int particleCount = (int)currentParticleArray_.size();
	const float scalefactor = 3.0f;

	glBegin(GL_LINES);
	for (int i = 0; i < particleCount; ++i)
	{
		const Particle& particle = currentParticleArray_[i];
		const Force& force = forceArray_[i];

		Vec2<float> forceVector(force.coord_[0], force.coord_[1]);
		forceVector.normalize();

		for (int i=3; i>0; --i)
		{
			glLineWidth((float)i);

			glColor3f(1.0f - ((float)i/4.0f), 1.0f - ((float)i/4.0f), 0.0f);
			glVertex2f(particle.position_[0], particle.position_[1]);	

			glColor3f(0.0f, 0.0f, 1.0f);
			glVertex2f(particle.position_[0] + forceVector[0]*scalefactor, particle.position_[1]  + forceVector[1]*scalefactor);	
		}
	}
	glEnd();

	glDepthMask(GL_TRUE);
}


/**
*/
void MassSpring2D::drawMass(void)
{
	glDepthMask(GL_FALSE);

	const int particleCount = (int)currentParticleArray_.size();
	const float scalefactor = 3.0f;

	for (int i = 0; i < particleCount; ++i)
	{
		const Particle& particle = currentParticleArray_[i];
		
		if (constraintArray_[i] == 0.0f)
		{
			glColor3f(1.0f, 0.0f, 0.0f);
		}
		else
		{
			glColor3f(1.0f, 1.0f, 0.0f);
		}
		
		glPushMatrix();
		glTranslatef(particle.position_[0], particle.position_[1], 0.0f);
		glutSolidSphere(0.5f, 10, 10);
		glPopMatrix();
			
	}

	glDepthMask(GL_TRUE);
}


/**
*/
void MassSpring2D::create(void)
{
	currentParticleArray_.clear();
	oldParticleArray_.clear();
	
	constraintArray_.clear();
	
	springArray_.clear();
	
	elementArray_.clear();
	
	damping_ = 0.95f;
}

