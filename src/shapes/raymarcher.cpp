
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// Modification to Sphere Shape to Implement RayMarching
// by Kevin M. Smith 3-2-2019
//


// shapes/sphere.cpp*
#include "shapes/raymarcher.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

namespace pbrt {

// Sphere Method Definitions
Bounds3f RayMarcher::ObjectBound() const {
    return Bounds3f(Point3f(-radius, -radius, zMin),
                    Point3f(radius, radius, zMax));
}


//  Note: create PBRT parameters for these
//

#define MAX_RAY_STEPS 1000
#define DIST_THRESHOLD .01
#define MAX_DISTANCE 100
#define NORMAL_EPS .01


//  Template Method
//
bool RayMarcher::Intersect(const Ray &r, Float *tHit,
                                  SurfaceInteraction *isect,
                                  bool testAlphaTexture) const {
   
    Vector3f dir = Normalize(r.d);     // ray direction vectors are not normalized in PBRT by default (KMS) 
	bool hit = false;
  
	
	if (hit && tHit != nullptr && isect != nullptr) {
		// Thiis where you return your SurfaceInteraction structure and your tHit
		// Important Note: You must check for null pointer as Intersect is called 
		// by IntersectP() with null values for these parameters.
 
	}
    return hit;
}


//  Template Method
//
Float RayMarcher::sdf(const Point3f &pos) const {
 
}



// Get Normal using Gradient (Finite Distances Methods )  - See class slides.
//  Note if the normal you calculate has zero length, return the defaultNormal
//
Vector3f RayMarcher::GetNormalRM( const Point3f &p, float eps, const Vector3f &defaultNormal) const {
 
}


Float RayMarcher::Area() const { return phiMax * radius * (zMax - zMin); }

// These functions are stubbed
//
Interaction RayMarcher::Sample(const Point2f &u, Float *pdf) const {
    LOG(FATAL) << "RayMarcher::Sample not implemented.";
    return Interaction();
}

Interaction RayMarcher::Sample(const Interaction &ref, const Point2f &u,
                           Float *pdf) const {
    LOG(FATAL) << "RayMarcher::Sample not implemented.";
    return Interaction();
}



std::shared_ptr<Shape> CreateRayMarcherShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params) {
    Float radius = params.FindOneFloat("radius", 1.f); 
    Float zmin = params.FindOneFloat("zmin", -radius);
    Float zmax = params.FindOneFloat("zmax", radius);
    Float phimax = params.FindOneFloat("phimax", 360.f);
    return std::make_shared<RayMarcher>(o2w, w2o, reverseOrientation, radius, zmin,
                                    zmax, phimax);
}

}  // namespace pbrt
