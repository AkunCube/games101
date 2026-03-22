#include <cassert>
#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

static constexpr float kd = 0.01f;
static constexpr float kDampingFactor = 0.00005f;

Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass,
           float k, const vector<int> &pinned_nodes) {
  // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and
  // containing `num_nodes` nodes.
  assert(num_nodes >= 2);
  Vector2D step = (end - start) / (num_nodes - 1);

  int cnt = 0;
  for (Vector2D pos = start; cnt < num_nodes; pos += step, ++cnt) {
    masses.push_back(new Mass(pos, node_mass, false));
  }
  for (const int &pinned_node : pinned_nodes) {
    masses[pinned_node]->pinned = true;
  }
  for (size_t i = 1; i < masses.size(); ++i) {
    springs.push_back(new Spring(masses[i - 1], masses[i], k));
  }
}

void Rope::simulateEuler(float delta_t, Vector2D gravity) {
  for (auto &s : springs) {
    // TODO (Part 2): Use Hooke's law to calculate the force on a node
    Vector2D vec = s->m2->position - s->m1->position;
    double deformation = vec.norm() - s->rest_length;
    Vector2D dir_vec = vec.unit();
    Vector2D force = s->k * dir_vec * deformation;
    s->m2->forces -= force;
    s->m1->forces += force;
  }

  for (auto &m : masses) {
    if (!m->pinned) {
      // TODO (Part 2): Add the force due to gravity, then compute the new
      // velocity and position

      // TODO (Part 2): Add global damping
      m->forces += gravity * m->mass;
      m->forces -= kd * m->velocity;
      Vector2D acc = m->forces / m->mass;
      m->velocity += acc * delta_t;
      m->position += m->velocity * delta_t;
    }

    // Reset all forces on each mass
    m->forces = Vector2D(0, 0);
  }
}

void Rope::simulateVerlet(float delta_t, Vector2D gravity) {
  for (auto &s : springs) {
    // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet
    // （solving constraints)
    Vector2D vec = s->m2->position - s->m1->position;
    double deformation = vec.norm() - s->rest_length;
    Vector2D dir_vec = vec.unit();
    Vector2D force = s->k * dir_vec * deformation;
    s->m2->forces -= force;
    s->m1->forces += force;
  }

  for (auto &m : masses) {
    if (!m->pinned) {
      // TODO (Part 3.1): Set the new position of the rope mass
      m->forces += gravity * m->mass;
      Vector2D acc = m->forces / m->mass;
      Vector2D new_pos =
          m->position +
          (1 - kDampingFactor) * (m->position - m->last_position) +
          acc * delta_t * delta_t;
      m->last_position = m->position;
      m->position = new_pos;
      // TODO (Part 4): Add global Verlet damping
    }
    // Reset all forces on each mass
    m->forces = Vector2D(0, 0);
  }
}
} // namespace CGL
