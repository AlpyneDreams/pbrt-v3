
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
// Implementation by Jeff Ellison, 2024


#include "shapes/raymarcher.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

namespace pbrt {

// RayMarcher Method Definitions
Bounds3f RayMarcher::ObjectBound() const {
    Point3f center = Point3f(0, 0, 0); //(*ObjectToWorld)(Point3f(0, 0, 0));
    Float radius = this->radius * 2;
    return Bounds3f(center + Point3f(-radius, -radius, -radius),
                    center + Point3f(radius, radius, radius));
}


//  Template Method
//
bool RayMarcher::Intersect(const Ray &r, Float *tHit,
                                  SurfaceInteraction *isect,
                                  bool testAlphaTexture) const {
   
#if 0
    static int counter = 0;
    printf("Counter: %d\n", counter++);
#endif

    Vector3f dir = Normalize(r.d);     // ray direction vectors are not normalized in PBRT by default (KMS) 
	bool hit = false;

    Point3f p = r.o;

#if 1
    // Ray marching
    for (int i = 0; i < maxRaySteps; i++) {
        float dist = sdf(p);
        if (dist < distThreshold) {
            hit = true;
            break;
        } else if (dist > maxDistance) {
        	break;
        } else {
            p += dir * dist; // move along ray
        }
    }
#else
    // Regular geometric intersector for testing/comparing with raymarcher
    auto SphereIntersect = [](Point3f ro, Vector3f rd, Point3f sph,
                              Float radius) -> Float {
        // From https://iquilezles.org/articles/intersectors/
        Vector3f oc = ro - sph;
        Float b = Dot(oc, rd);
        Float c = Dot(oc, oc) - radius * radius;
        Float h = b * b - c;
        if (h < 0.0) return -1.0;
        h = std::sqrt(h);
        return -b - h;
    };
    Point3f center = (*ObjectToWorld)(Point3f(0, 0, 0));
    Float dist = SphereIntersect(r.o, dir, center, radius);
    if (dist >= 0) {
        hit = true;
        p = r.o + dir * dist;
    }
#endif
	
	if (hit && tHit != nullptr && isect != nullptr) {
		// Thiis where you return your SurfaceInteraction structure and your tHit
		// Important Note: You must check for null pointer as Intersect is called 
		// by IntersectP() with null values for these parameters.
        if (tHit) {
            *tHit = (p - r.o).Length();
        }
        if (isect) {
            Vector3f pError =
                Vector3f(distThreshold, distThreshold, distThreshold) * 10;
            Point2f uv = Point2f(0, 0);

            Vector3f normal = GetNormalRM(p, normalEps, Vector3f(0, 0, 1));
            Vector3f dpdu, dpdv;
            CoordinateSystem(normal, &dpdu, &dpdv);

            *isect = SurfaceInteraction(p, pError, uv, -r.d, dpdu, dpdv,
                                        Normal3f(0, 0, 0), Normal3f(0, 0, 0),
                                        r.time, this);
        }
	}
    return hit;
}

static inline Float OctaveNoise(Point3f p, int nOctaves, Float freq, Float amp) {

    // Only consider points on the 2D plane
    Point3f pos = Point3f(p.x, p.y, 0);
    
    Float noise = 0;
    for (int i = 0; i < nOctaves; i++) {
        noise += amp * Noise(freq * pos);
        freq *= 2;
        amp /= 2;
    }
    return noise;
}

Float opIntersection(Float d1, Float d2) {
    return d1 > d2 ? d1 : d2;
}

//  SDF Method
//
Float RayMarcher::sdf(const Point3f &pos) const {
    Vector3f normal = Vector3f(0, 0, 1);
    Point3f center = (*ObjectToWorld)(Point3f(0, 0, 0));
    Point3f localPos = pos - Vector3f(center);

    Float heightField = Dot(Vector3f(localPos), normal) + OctaveNoise( localPos, octaves, frequency, amplitude );
    Float sphere = (pos - center).Length() - (radius * 2.f);

    return opIntersection(heightField, sphere);
}


// Get Normal using Gradient (Finite Distances Methods )  - See class slides.
//  Note if the normal you calculate has zero length, return the defaultNormal
//
Vector3f RayMarcher::GetNormalRM( const Point3f &p, float eps, const Vector3f &defaultNormal ) const {
#if 1
    // Finite Distance Method
    float dp = sdf(p);
    Vector3f n = Vector3f(dp - sdf(Point3f(p.x - eps, p.y, p.z)),
						  dp - sdf(Point3f(p.x, p.y - eps, p.z)),
						  dp - sdf(Point3f(p.x, p.y, p.z - eps)));
    return Normalize(n);
#else
    // Geometric version for testing
    Point3f center = (*ObjectToWorld)(Point3f(0, 0, 0));
    Vector3f normal = (p - center) / radius;
    if (normal.Length() == 0)
        return defaultNormal;
    return normal;
#endif
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
    int maxRaySteps = params.FindOneInt("steps", 1000);
    Float distThreshold = params.FindOneFloat("threshold", .01);
    Float maxDistance = params.FindOneFloat("maxdist", 100);
    Float normalEps = params.FindOneFloat("eps", .01);
    Float amplitude = params.FindOneFloat("amplitude", 200);
    Float frequency = params.FindOneFloat("frequency", 0.0025);
    int octaves = params.FindOneInt("octaves", 3);
    return std::make_shared<RayMarcher>(o2w, w2o, reverseOrientation, radius, zmin,
                                    zmax, phimax, maxRaySteps, distThreshold, maxDistance, normalEps,
                                    amplitude, frequency, octaves);
}

}  // namespace pbrt
