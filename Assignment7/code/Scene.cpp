//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"
#include <cassert>

void Scene::buildBVH() {
  printf(" - Generating BVH...\n\n");
  this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const {
  return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const {
  float emit_area_sum = 0;
  for (uint32_t k = 0; k < objects.size(); ++k) {
    if (objects[k]->hasEmit()) {
      emit_area_sum += objects[k]->getArea();
    }
  }
  float p = get_random_float() * emit_area_sum;
  emit_area_sum = 0;
  for (uint32_t k = 0; k < objects.size(); ++k) {
    if (objects[k]->hasEmit()) {
      emit_area_sum += objects[k]->getArea();
      if (p <= emit_area_sum) {
        objects[k]->Sample(pos, pdf);
        break;
      }
    }
  }
}

bool Scene::trace(const Ray &ray, const std::vector<Object *> &objects,
                  float &tNear, uint32_t &index, Object **hitObject) {
  *hitObject = nullptr;
  for (uint32_t k = 0; k < objects.size(); ++k) {
    float tNearK = kInfinity;
    uint32_t indexK;
    Vector2f uvK;
    if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
      *hitObject = objects[k];
      tNear = tNearK;
      index = indexK;
    }
  }

  return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const {
  // TO DO Implement Path Tracing Algorithm here
  Vector3f resultLight(0.0f, 0.0f, 0.0f);
  Intersection intersection = intersect(ray);

  // Return black for no intersection (background).
  if (!intersection.happened) {
    return resultLight;
  }

  Vector3f hitPoint = intersection.coords;
  Vector3f hitNormal = intersection.normal;
  Material *hitMaterial = intersection.m;
  assert(intersection.obj && hitMaterial);

  // Return emission only for primary ray hitting light source (avoid double
  // count).
  if (hitMaterial->hasEmission()) {
    return (depth == 0) ? hitMaterial->m_emission : resultLight;
  }

  // --------------------------
  // 1. Direct Illumination (Light Sampling)
  // --------------------------
  Intersection lightSample;
  float pdfLight = 0.0f;
  sampleLight(lightSample, pdfLight);

  do {
    if (!lightSample.happened) {
      break;
    }

    Vector3f wi = normalize(lightSample.coords - hitPoint);
    Ray shadowRay(hitPoint, wi);
    Intersection shadowTest = intersect(shadowRay);

    // Skip if light is blocked by other geometry.
    if (!shadowTest.happened || shadowTest.obj != lightSample.obj) {
      break;
    }
    Vector3f wo = -ray.direction;
    Vector3f fr = hitMaterial->eval(wi, wo, hitNormal);
    float distSq = (lightSample.coords - hitPoint).squaredNorm();
    resultLight += lightSample.emit * fr * dotProduct(wi, hitNormal) *
                   dotProduct(-wi, lightSample.normal) / distSq / pdfLight;
  } while (false);

  // --------------------------
  // 2. Indirect Illumination (BRDF Sampling + Russian Roulette)
  // --------------------------
  // Terminate path with Russian Roulette.
  if (get_random_float() > Scene::RussianRoulette) {
    return resultLight;
  }

  Vector3f wo = -ray.direction;
  Vector3f wi = hitMaterial->sample(wo, hitNormal).normalized();
  Ray indirectRay(hitPoint, wi);
  Intersection indirectHit = intersect(indirectRay);

  // Skip if secondary ray hits nothing or emissive surface (avoid double
  // count).
  if (!indirectHit.happened || indirectHit.m->hasEmission()) {
    return resultLight;
  }
  Vector3f indirectRadiance = castRay(indirectRay, depth + 1);
  Vector3f fr = hitMaterial->eval(wi, wo, hitNormal);
  float pdf = hitMaterial->pdf(wi, wo, hitNormal);
  float cosSurface = dotProduct(wi, hitNormal);

  resultLight +=
      indirectRadiance * fr * cosSurface / pdf / Scene::RussianRoulette;
  return resultLight;
}
