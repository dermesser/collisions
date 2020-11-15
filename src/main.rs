use std::ops;

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use rand::{thread_rng, Rng};

type BaseNum = f32;

const TIME_UNIT: BaseNum = 1e-6;

#[derive(Clone, Copy, Debug, PartialEq)]
struct Vector(BaseNum, BaseNum);

impl ops::Add for Vector {
    type Output = Vector;
    fn add(self, other: Vector) -> Self::Output {
        Vector(self.0 + other.0, self.1 + other.1)
    }
}

impl ops::Sub for Vector {
    type Output = Vector;
    fn sub(self, other: Vector) -> Self::Output {
        Vector(self.0 - other.0, self.1 - other.1)
    }
}

impl ops::Mul<BaseNum> for Vector {
    type Output = Vector;
    fn mul(self, f: BaseNum) -> Self::Output {
        Vector(self.0 * f, self.1 * f)
    }
}

impl Vector {
    fn normsq(&self) -> BaseNum {
        self.dot(&self)
    }
    fn dot(self, other: &Vector) -> BaseNum {
        self.0 * other.0 + self.1 * other.1
    }
    fn unit(&self) -> Vector {
        let norm = self.normsq().sqrt();
        *self * (1. / norm)
    }
    fn orth(&self) -> Vector {
        Vector(-self.1, self.0)
    }
}

#[derive(PartialEq, Eq, Debug, Clone)]
struct Event {
    typ: EventType,
    when: usize,
}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Event) -> Option<std::cmp::Ordering> {
        self.when.partial_cmp(&other.when)
    }
}
impl Ord for Event {
    fn cmp(&self, other: &Event) -> std::cmp::Ordering {
        self.when.cmp(&other.when)
    }
}

#[derive(PartialEq, Eq, Debug, Clone)]
enum EventType {
    PP(ParticleHandle, ParticleHandle),
    HorizWall(ParticleHandle),
    VertWall(ParticleHandle),
    Redraw,
}

#[derive(Eq, PartialEq, Debug, Clone, Copy)]
struct ParticleHandle(usize);

#[derive(Debug, Clone, PartialEq)]
struct Particle {
    pos: Vector,
    vel: Vector,
    rad: BaseNum,
    mass: BaseNum,
    coll: usize,
}

impl Particle {
    fn draw(&self) {
        unimplemented!()
    }
    fn mov(&mut self, dt: usize) {
        self.pos = self.pos + self.vel * dt as f32 * TIME_UNIT;
    }
    fn timetohit(&self, dt: usize, other: &Particle, limit: usize) -> Option<usize> {
        let mut mypos = self.pos;
        let mut otherpos = other.pos;
        let mindistsq = self.rad * self.rad + other.rad * other.rad;
        let mut i = 0;

        while (otherpos - mypos).normsq() < mindistsq && i * dt < limit {
            mypos = mypos + self.vel * dt as f32 * TIME_UNIT;
            otherpos = otherpos + other.vel * dt as f32 * TIME_UNIT;
            i += 1;
        }

        if (otherpos - mypos).normsq() < mindistsq {
            Some(i * dt)
        } else {
            None
        }
    }

    fn timetohit_horiz_wall(&self, y: BaseNum) -> Option<usize> {
        if (self.pos.1 > y && self.vel.1 > 0.) || (self.pos.1 < y && self.vel.1 < 0.) {
            None
        } else {
            Some((((y - self.pos.1).abs() / self.vel.1) / TIME_UNIT as BaseNum) as usize)
        }
    }
    fn timetohit_vert_wall(&self, x: BaseNum) -> Option<usize> {
        if (self.pos.0 > x && self.vel.0 > 0.) || (self.pos.0 < x && self.vel.0 < 0.) {
            None
        } else {
            Some((((x - self.pos.0).abs() / self.vel.0) / TIME_UNIT as BaseNum) as usize)
        }
    }
    fn bounce_off(&mut self, other: &mut Particle) {
        let conn = self.pos - other.pos;
        assert!(conn.normsq() <= (self.rad * self.rad + other.rad * other.rad));

        // Determine orthogonal bounce plane.
        let connu = conn.unit();
        let bounceplane = connu.orth();
        // Relevant speeds are parallel to connection line between positions.
        let v1_bounce = self.vel.dot(&connu);
        let v2_bounce = self.vel.dot(&connu);

        let v1_orth = self.vel.dot(&bounceplane);
        let v2_orth = self.vel.dot(&bounceplane);

        let v1_bounce_after = 2. * (self.mass * v1_bounce + other.mass * v2_bounce)
            / (self.mass + other.mass)
            - v1_bounce;
        let v2_bounce_after = 2. * (self.mass * v1_bounce + other.mass * v2_bounce)
            / (self.mass + other.mass)
            - v2_bounce;

        self.vel = connu * v1_bounce_after + bounceplane * v1_orth;
        other.vel = connu * v2_bounce_after + bounceplane * v2_orth;
    }
    fn bounce_horiz_wall(&mut self) {
        self.vel.1 = -self.vel.1
    }
    fn bounce_vert_wall(&mut self) {
        self.vel.0 = -self.vel.0
    }
}

fn test_generate_data() {
    let dt = 0.01;
    // 2x2 box.
    let mut p1 = Particle {
        pos: Vector(0., 0.),
        vel: Vector(1., 1.),
        rad: 0.1,
        mass: 1.,
        coll: 0,
    };
    let mut p2 = Particle {
        pos: Vector(1., 1.),
        vel: Vector(-0.9, -0.9),
        ..p1
    };

    for i in 0..1000 {}
}

struct System {
    parts: Vec<Particle>,
    pq: BinaryHeap<Reverse<Event>>,
    t: usize,
    dt: usize,
    dim: BaseNum,
}

impl System {
    fn new(dim: BaseNum, dt: usize, particles: usize) -> System {
        let mut rng = thread_rng();
        let particles = (0..particles)
            .map(|_| Particle {
                pos: Vector(rng.gen_range(0., dim), rng.gen_range(0., dim)),
                vel: Vector(rng.gen_range(0., dim), rng.gen_range(0., dim)),
                mass: rng.gen_range(0., 1.),
                rad: rng.gen_range(0., dim / 100.),
                coll: 0,
            })
            .collect::<Vec<Particle>>();

        System {
            parts: particles,
            pq: BinaryHeap::new(),
            t: 0,
            dt: dt,
            dim: dim,
        }
    }

    fn handle_for(&self, i: usize) -> ParticleHandle {
        ParticleHandle(i)
    }
    fn particle_for(&self, ph: ParticleHandle) -> &Particle {
        &self.parts[ph.0]
    }
    fn particle_for_mut(&mut self, ph: ParticleHandle) -> &mut Particle {
        &mut self.parts[ph.0]
    }
    fn for_all_particles<F>(&mut self, f: F) where F: std::ops::Fn(&mut Particle) {
        for p in self.parts.iter_mut() {
            f(p);
        }
    }

    fn predict_collisions(&mut self, p: ParticleHandle, limit: usize) {
        for (i, other) in self.parts.iter().enumerate() {
            if let Some(tth) = self.particle_for(p).timetohit(self.dt, other, limit) {
                self.pq.push(Reverse(Event {
                    typ: EventType::PP(p, self.handle_for(i)),
                    when: (tth / self.dt) as usize,
                }));
            }
        }
        let particle = self.particle_for(p);
        let (dtx, dty) = (
            particle.timetohit_vert_wall(self.dim),
            particle.timetohit_horiz_wall(self.dim),
        );
        if let Some(dtx) = dtx {
            if dtx <= limit {
                self.pq.push(Reverse(Event {
                    typ: EventType::VertWall(p),
                    when: self.t + dtx,
                }));
            }
        }
        if let Some(dty) = dty {
            if dty <= limit {
                self.pq.push(Reverse(Event {
                    typ: EventType::HorizWall(p),
                    when: self.t + dty,
                }));
            }
        }
    }

    fn simulate(&mut self, limit: usize, redrawHz: usize) {
        for i in 0..self.parts.len() {
            self.predict_collisions(self.handle_for(i), limit);
        }
        self.pq.push(Reverse(Event { typ: EventType::Redraw, when: (1./(redrawHz*self.dt) as f32) as usize }));

        while !self.pq.is_empty() {
            if let Some(next) = self.pq.pop() {
                let next = next.0;
                let dt = next.when;
                let t = self.t + dt;
                self.for_all_particles(|p| p.mov(dt));
                match next.typ {
                    EventType::Redraw => {},
                    EventType::HorizWall(p) => {
                        self.particle_for_mut(p).bounce_horiz_wall();
                        self.predict_collisions(p, limit);
                    },
                    EventType::VertWall(p) => {
                        self.particle_for_mut(p).bounce_vert_wall();
                        self.predict_collisions(p, limit);
                    },
                    EventType::PP(p1, p2) => {
                        let mut pm1 = self.parts[p1.0].clone();
                        let mut pm2 = self.parts[p2.0].clone();
                        pm1.bounce_off(&mut pm2);
                        self.parts[p1.0] = pm1;
                        self.parts[p2.0] = pm2;
                    },
                }
            }
        }
    }
}

fn main() {
    println!("Hello, world!");
}
