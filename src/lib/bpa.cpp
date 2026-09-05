#include "bpa.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <numbers>
#include <numeric>
#include <optional>
#include <vector>

using namespace glm;

namespace bpa {
	namespace {
		struct MeshEdge;

		struct MeshPoint {
			dvec3 pos;
			dvec3 normal;
			std::uint32_t index; // into the input vector
			bool used = false;
			std::vector<MeshEdge*> edges;
		};

		enum class EdgeStatus { active, inner, boundary };

		struct MeshEdge {
			MeshPoint* a;
			MeshPoint* b;
			MeshPoint* opposite;
			dvec3 center;
			MeshEdge* prev;
			MeshEdge* next;
			EdgeStatus status = EdgeStatus::active;
		};

		struct MeshFace : std::array<MeshPoint*, 3> {
			auto normal() const { return normalize(cross((*this)[0]->pos - (*this)[1]->pos, (*this)[0]->pos - (*this)[2]->pos)); }
		};

		using Cell = std::vector<MeshPoint>;

		// Uniform grid of cells of side 2 * radius: every point a ball touching a point in a cell
		// can touch lies in that cell or one of its 26 neighbours.
		struct Grid {
			Grid(const std::vector<Point>& points, double radius) : cellSize(radius * 2) {
				lower = points.front().pos;
				upper = points.front().pos;
				for (const auto& p : points) {
					lower = min(lower, p.pos);
					upper = max(upper, p.pos);
				}
				dims = max(ivec3{ceil((upper - lower) / cellSize)}, ivec3{1});
				cells.resize(std::size_t(dims.x) * dims.y * dims.z);
				for (std::size_t i = 0; i < points.size(); i++)
					cell(cellIndex(points[i].pos)).push_back({points[i].pos, points[i].normal, static_cast<std::uint32_t>(i)});
			}

			auto cellIndex(dvec3 point) const -> ivec3 { return clamp(ivec3{(point - lower) / cellSize}, ivec3{}, dims - 1); }

			auto cell(ivec3 index) -> Cell& { return cells[std::size_t(index.z) * dims.x * dims.y + std::size_t(index.y) * dims.x + index.x]; }

			// The points within 2 * radius of `point`, except those in `ignore`.
			auto sphericalNeighborhood(dvec3 point, std::initializer_list<const MeshPoint*> ignore) -> std::vector<MeshPoint*> {
				std::vector<MeshPoint*> result;
				const auto centerIndex = cellIndex(point);
				result.reserve(cell(centerIndex).size() * 27);
				for (auto xOff : {-1, 0, 1})
					for (auto yOff : {-1, 0, 1})
						for (auto zOff : {-1, 0, 1}) {
							const auto index = centerIndex + ivec3{xOff, yOff, zOff};
							if (any(lessThan(index, ivec3{})) || any(greaterThanEqual(index, dims)))
								continue;
							for (auto& p : cell(index))
								if (length2(p.pos - point) < cellSize * cellSize && std::find(begin(ignore), end(ignore), &p) == end(ignore))
									result.push_back(&p);
						}
				return result;
			}

			dvec3 lower;
			dvec3 upper;
			double cellSize;
			ivec3 dims;
			std::vector<Cell> cells;
		};

		// Centre of the ball of the given radius through the three points of `f`, on the side of
		// its normal; nothing if the circumradius is larger than the ball.
		auto computeBallCenter(MeshFace f, double radius) -> std::optional<dvec3> {
			const dvec3 ac = f[2]->pos - f[0]->pos;
			const dvec3 ab = f[1]->pos - f[0]->pos;
			const dvec3 abXac = cross(ab, ac);
			const auto nn = dot(abXac, abXac);
			if (nn <= 1e-20 * dot(ab, ab) * dot(ac, ac)) // degenerate triangle
				return {};
			const dvec3 toCircumCircleCenter = (cross(abXac, ab) * dot(ac, ac) + cross(ac, abXac) * dot(ab, ab)) / (2 * nn);
			const dvec3 circumCircleCenter = f[0]->pos + toCircumCircleCenter;

			const auto heightSquared = radius * radius - dot(toCircumCircleCenter, toCircumCircleCenter);
			if (heightSquared < 0)
				return {};
			return circumCircleCenter + f.normal() * std::sqrt(heightSquared);
		}

		// Points within a relative 1e-9 of the sphere count as outside, so that the three points
		// the ball touches, and points exactly cospherical with them, do not fail the test.
		constexpr auto emptinessTolerance = 1e-9;

		auto ballIsEmpty(dvec3 ballCenter, const std::vector<MeshPoint*>& points, double radius) -> bool {
			const auto r2 = radius * (1 - emptinessTolerance) * radius * (1 - emptinessTolerance);
			return !std::any_of(begin(points), end(points), [&](MeshPoint* p) { return length2(p->pos - ballCenter) < r2; });
		}

		struct SeedResult {
			MeshFace f;
			dvec3 ballCenter;
		};

		// Only the nearest `seedNeighbors` neighbours of a seed candidate are paired: a valid
		// seed's other two vertices are almost always among the closest points, and for a
		// candidate under an already reconstructed sheet every pair fails, so the loop must be
		// cheap.
		// Seed search (section 4.2 of the paper). Cells are visited from a cursor that persists
		// between calls; a cell holding a used point is skipped (the paper's heuristic against
		// spawning small components next to the surface, fig. 4c), else one candidate is tried,
		// the point projecting furthest along the cell's average normal, paired with its nearest
		// neighbours: the first triangle facing along its vertex normals with an empty ball on
		// that side is the seed. A candidate that fails cannot succeed later (points only become used), so the
		// cursor never moves back; after a seed it stays, and the cell is skipped next time.
		auto findSeedTriangle(Grid& grid, double radius, std::size_t seedNeighbors, std::size_t& cursor) -> std::optional<SeedResult> {
			for (; cursor < grid.cells.size(); cursor++) {
				auto& cell = grid.cells[cursor];
				if (cell.empty() || std::any_of(begin(cell), end(cell), [](const MeshPoint& p) { return p.used; }))
					continue;
				const auto avgNormal =
					normalize(std::accumulate(begin(cell), end(cell), dvec3{}, [](dvec3 acc, const MeshPoint& p) { return acc + p.normal; }));
				const auto centroid =
					std::accumulate(begin(cell), end(cell), dvec3{}, [](dvec3 acc, const MeshPoint& p) { return acc + p.pos; }) / double(cell.size());
				auto& p1 = *std::max_element(begin(cell), end(cell),
					[&](const MeshPoint& a, const MeshPoint& b) { return dot(a.pos - centroid, avgNormal) < dot(b.pos - centroid, avgNormal); });

				auto neighborhood = grid.sphericalNeighborhood(p1.pos, {&p1});
				std::sort(
					begin(neighborhood), end(neighborhood), [&](MeshPoint* a, MeshPoint* b) { return length(a->pos - p1.pos) < length(b->pos - p1.pos); });
				const auto nPairs = seedNeighbors == 0 ? neighborhood.size() : std::min(neighborhood.size(), seedNeighbors);
				for (std::size_t i2 = 0; i2 < nPairs; i2++) {
					auto* p2 = neighborhood[i2];
					if (p2->used)
						continue;
					for (std::size_t i3 = i2 + 1; i3 < nPairs; i3++) {
						auto* p3 = neighborhood[i3];
						if (p3->used)
							continue;
						// the seed must face along the normals of all three of its vertices (paper,
						// section 4.2); the pair is tried with the winding that does, if either
						MeshFace f{{&p1, p2, p3}};
						const auto n = f.normal();
						const auto agrees = [&](double sign) {
							return sign * dot(n, p1.normal) >= 0 && sign * dot(n, p2->normal) >= 0 && sign * dot(n, p3->normal) >= 0;
						};
						if (agrees(-1))
							std::swap(f[1], f[2]);
						else if (!agrees(1))
							continue;
						const auto ballCenter = computeBallCenter(f, radius);
						if (ballCenter && ballIsEmpty(ballCenter.value(), neighborhood, radius)) {
							p1.used = true;
							p2->used = true;
							p3->used = true;
							return SeedResult{f, ballCenter.value()};
						}
					}
				}
			}
			return {};
		}

		auto getActiveEdge(std::vector<MeshEdge*>& front) -> std::optional<MeshEdge*> {
			while (!front.empty()) {
				auto* e = front.back();
				if (e->status == EdgeStatus::active)
					return e;
				front.pop_back(); // non-active edges are dropped lazily
			}
			return {};
		}

		// ---- the front predicates (section 4.4): an edge is live while it is active or boundary;
		// an inner edge has its two triangles and only records that the undirected edge is closed.

		auto isLive(const MeshEdge* e) -> bool { return e->status != EdgeStatus::inner; }

		auto notUsed(const MeshPoint* p) -> bool { return !p->used; }

		auto onFront(const MeshPoint* p) -> bool { return std::any_of(begin(p->edges), end(p->edges), isLive); }

		// The front holds the directed edge i -> j.
		auto hasEdge(const MeshPoint* i, const MeshPoint* j) -> bool {
			return std::any_of(begin(i->edges), end(i->edges), [&](const MeshEdge* e) { return isLive(e) && e->a == i && e->b == j; });
		}

		// The undirected edge {i, j} already has two triangles.
		auto isClosed(const MeshPoint* i, const MeshPoint* j) -> bool {
			return std::any_of(begin(i->edges), end(i->edges),
				[&](const MeshEdge* e) { return e->status == EdgeStatus::inner && ((e->a == i && e->b == j) || (e->a == j && e->b == i)); });
		}

		// The tests the paper applies to the point k the ball lands on when pivoting e = (a, b)
		// (fig. 5, line 3, and the "edge orientation checks" it mentions): the triangle (a, k, b)
		// must not face against the normal of k, k must be unused or on the front, and the mesh
		// must stay a manifold: neither new half-edge a -> k, k -> b may already be on the front
		// with the same orientation, nor either undirected edge closed.
		auto canAddTriangle(const MeshEdge* e, const MeshPoint* k) -> bool {
			const auto normal = cross(k->pos - e->a->pos, e->b->pos - e->a->pos);
			if (dot(normal, k->normal) < 0)
				return false;
			if (!(notUsed(k) || onFront(k)))
				return false;
			return !hasEdge(e->a, k) && !hasEdge(k, e->b) && !isClosed(e->a, k) && !isClosed(k, e->b);
		}

		// ---- ball pivoting (section 4.3, fig. 2)

		// Angle below which a candidate counts as touching the ball in its initial position.
		constexpr auto touchTolerance = 1e-6;
		// Angle within which two hits count as simultaneous (cospherical points, as on a lattice).
		constexpr auto tieTolerance = 1e-7;

		// The ball centre pivoting around edge (a, b) moves on the circle m + r (cos t u + sin t v)
		// in the plane perpendicular to the edge through its midpoint m: u points from m to the
		// current centre, and v = (b - a) x u is the direction in which the ball leaves the
		// current triangle.
		struct PivotFrame {
			dvec3 m, a, u, v;
			double r;
		};

		auto pivotFrame(const MeshEdge* e) -> std::optional<PivotFrame> {
			const auto m = (e->a->pos + e->b->pos) / 2.0;
			auto a = e->b->pos - e->a->pos;
			const auto la = length(a);
			if (la == 0)
				return {};
			a /= la;
			auto w = e->center - m;
			w -= dot(w, a) * a; // numerical drift along the edge
			const auto r = length(w);
			if (r <= 1e-12 * la)
				return {};
			const auto u = w / r;
			return PivotFrame{m, a, u, cross(a, u), r};
		}

		// A position of the ball centre on the pivot circle, as coordinates (x, y) in the frame
		// (u, v): the centre is m + x u + y v with x^2 + y^2 = r^2, and the pivot angle is the
		// argument of (x, y) in [0, 2 pi). Contacts are compared by angle without evaluating it.
		using Contact = dvec2;

		auto centerAt(const PivotFrame& fr, Contact c) -> dvec3 { return fr.m + c.x * fr.u + c.y * fr.v; }

		// The angle of c lies in [pi, 2 pi).
		auto lowerHalf(Contact c) -> bool { return c.y < 0 || (c.y == 0 && c.x < 0); }

		// The pivot angle of p is smaller than that of q: an angle in [0, pi) precedes any in
		// [pi, 2 pi); within one half the two differ by less than pi, and the sign of the 2-D
		// cross product decides.
		auto angleLess(Contact p, Contact q) -> bool {
			const auto hp = lowerHalf(p);
			const auto hq = lowerHalf(q);
			return hp == hq ? p.x * q.y - p.y * q.x > 0 : hq;
		}

		// The pivot angles of p and q differ by at most the angle whose sine is sinTol, r2 being
		// the squared radius of the pivot circle: r^2 sin d = p x q and r^2 cos d = p . q. Two
		// contacts on either side of angle 0 are not a tie: the ball reaches one at once and the
		// other after a full turn.
		auto angleTie(Contact p, Contact q, double r2, double sinTol) -> bool {
			return p.x * q.x + p.y * q.y > 0 && std::abs(p.x * q.y - p.y * q.x) <= r2 * sinTol && (lowerHalf(p) == lowerHalf(q) || p.x < 0);
		}

		const auto touchSin = std::sin(touchTolerance);
		const auto tieSin = std::sin(tieTolerance);

		// The contact is within touchTolerance of angle 0.
		auto touching(Contact c, double r) -> bool { return c.x > 0 && std::abs(c.y) < r * touchSin; }

		// The ball centre at the first contact of the pivoting ball with x, or nothing if the
		// ball never reaches x. With d = x - m and (d_u, d_v) its components in the pivot plane,
		// |centre - x|^2 = radius^2 is the line d_u X + d_v Y = r K with K = (r^2 + |d|^2 -
		// radius^2) / (2 r), whose intersections with the circle of radius r are F +- h (-d_v,
		// d_u) / R: F = (r K / R^2)(d_u, d_v) is the foot of the perpendicular from the origin
		// and h = r sqrt(1 - (K / R)^2) the half-chord. A contact at (nearly) zero means x
		// touches the initial ball: if the ball is moving into x the hit is immediate, else that
		// contact is ignored and the other one, where the ball comes back to x from the far
		// side, is used. The opposite vertex of the edge is such a point.
		auto pivotContact(const PivotFrame& fr, dvec3 x, double radius) -> std::optional<Contact> {
			const auto d = x - fr.m;
			const auto du = dot(d, fr.u);
			const auto dv = dot(d, fr.v);
			const auto R = std::sqrt(du * du + dv * dv);
			if (R <= 1e-12 * radius)
				return {};
			const auto r = fr.r;
			auto ratio = (r * r + dot(d, d) - radius * radius) / (2 * r) / R;
			if (std::abs(ratio) > 1 + 1e-9)
				return {};
			ratio = std::clamp(ratio, -1.0, 1.0);
			const auto s = r * ratio / R;
			const auto h = r * std::sqrt(1 - ratio * ratio) / R;
			const Contact ca{s * du - h * dv, s * dv + h * du}; // phi + alpha
			const Contact cb{s * du + h * dv, s * dv - h * du}; // phi - alpha
			const auto ta = touching(ca, r);
			const auto tb = touching(cb, r);
			if (ta || tb) {
				if (dv > 0) // the centre moves along v: into x
					return Contact{r, 0};
				if (ta && tb)
					return {};
				return ta ? cb : ca;
			}
			return angleLess(cb, ca) ? cb : ca;
		}

		struct PivotResult {
			MeshPoint* p;
			dvec3 center;
		};

		// Preference among points hit simultaneously: 0 if the triangle would be rejected, else
		// 1 plus the number of its new edges that glue to an existing front edge. The first pivot
		// into a cospherical polygon fixes a diagonal; later pivots can only respect it.
		auto tieScore(const MeshEdge* e, const MeshPoint* k) -> int {
			if (!canAddTriangle(e, k))
				return 0;
			return 1 + hasEdge(k, e->a) + hasEdge(e->b, k);
		}

		// Roll the ball around edge e from its stored centre and return the point it touches
		// first with the centre at that moment. Every point within reach of the ball is a
		// candidate, the opposite vertex of e included: if the ball comes back to it before
		// touching anything else there is no triangle to build. Because the initial ball is
		// empty, the ball at the first contact is empty too.
		auto ballPivot(const MeshEdge* e, Grid& grid, double radius) -> std::optional<PivotResult> {
			const auto fr = pivotFrame(e);
			if (!fr)
				return {};
			const auto neighborhood = grid.sphericalNeighborhood(fr->m, {e->a, e->b});
			const auto r2 = fr->r * fr->r;

			MeshPoint* best = nullptr;
			Contact bestContact{};
			auto nTies = 0;
			for (auto* p : neighborhood) {
				const auto c = pivotContact(*fr, p->pos, radius);
				if (!c)
					continue;
				if (!best) {
					best = p;
					bestContact = *c;
					nTies = 1;
				} else if (angleTie(*c, bestContact, r2, tieSin)) {
					nTies++;
					if (angleLess(*c, bestContact))
						bestContact = *c;
				} else if (angleLess(*c, bestContact)) {
					best = p;
					bestContact = *c;
					nTies = 1;
				}
			}
			if (!best)
				return {};
			if (nTies == 1) {
				if (best == e->opposite)
					return {}; // the ball came back to the opposite vertex first
				return PivotResult{best, centerAt(*fr, bestContact)};
			}

			// Simultaneous hits: pick deterministically. The opposite vertex is never chosen; a
			// point hit at the same angle gives a valid triangle with it on the ball's surface.
			MeshPoint* choice = nullptr;
			Contact choiceContact{};
			auto bestScore = -1;
			for (auto* p : neighborhood) {
				if (p == e->opposite)
					continue;
				const auto c = pivotContact(*fr, p->pos, radius);
				if (!c || !(angleTie(*c, bestContact, r2, tieSin) || angleLess(*c, bestContact)))
					continue;
				const auto score = tieScore(e, p);
				if (score > bestScore || (score == bestScore && p->index < choice->index)) {
					bestScore = score;
					choice = p;
					choiceContact = *c;
				}
			}
			if (!choice)
				return {};
			return PivotResult{choice, centerAt(*fr, choiceContact)};
		}

		// Mark the edge inner; getActiveEdge() drops it from the front later.
		void remove(MeshEdge* edge) { edge->status = EdgeStatus::inner; }

		void outputTriangle(MeshFace f, std::vector<Face>& triangles) { triangles.push_back({f[0]->index, f[1]->index, f[2]->index}); }

		auto join(MeshEdge* e_ij, MeshPoint* o_k, dvec3 o_k_ballCenter, std::vector<MeshEdge*>& front, std::deque<MeshEdge>& edges)
			-> std::pair<MeshEdge*, MeshEdge*> {
			auto& e_ik = edges.emplace_back(MeshEdge{e_ij->a, o_k, e_ij->b, o_k_ballCenter});
			auto& e_kj = edges.emplace_back(MeshEdge{o_k, e_ij->b, e_ij->a, o_k_ballCenter});

			e_ik.next = &e_kj;
			e_ik.prev = e_ij->prev;
			e_ij->prev->next = &e_ik;
			e_ij->a->edges.push_back(&e_ik);

			e_kj.prev = &e_ik;
			e_kj.next = e_ij->next;
			e_ij->next->prev = &e_kj;
			e_ij->b->edges.push_back(&e_kj);

			o_k->used = true;
			o_k->edges.push_back(&e_ik);
			o_k->edges.push_back(&e_kj);

			front.push_back(&e_ik);
			front.push_back(&e_kj);
			remove(e_ij);

			return {&e_ik, &e_kj};
		}

		void glue(MeshEdge* a, MeshEdge* b) {
			// case 1: a and b form a loop of two
			if (a->next == b && a->prev == b && b->next == a && b->prev == a) {
				remove(a);
				remove(b);
				return;
			}
			// case 2: a and b are adjacent
			if (a->next == b && b->prev == a) {
				a->prev->next = b->next;
				b->next->prev = a->prev;
				remove(a);
				remove(b);
				return;
			}
			if (a->prev == b && b->next == a) {
				a->next->prev = b->prev;
				b->prev->next = a->next;
				remove(a);
				remove(b);
				return;
			}
			// cases 3 and 4: a and b are apart, in the same loop or in two
			a->prev->next = b->next;
			b->next->prev = a->prev;
			a->next->prev = b->prev;
			b->prev->next = a->next;
			remove(a);
			remove(b);
		}

		// The live front edge opposite to `edge`, if there is one.
		auto findReverseEdgeOnFront(MeshEdge* edge) -> MeshEdge* {
			for (auto& e : edge->a->edges)
				if (isLive(e) && e->a == edge->b && e->b == edge->a)
					return e;
			return nullptr;
		}

		// Remove the connected components with fewer than minComponent triangles.
		void dropSmallComponents(std::vector<Face>& triangles, std::size_t nPoints, std::size_t minComponent) {
			std::vector<std::uint32_t> parent(nPoints);
			std::iota(begin(parent), end(parent), 0);
			const auto root = [&](std::uint32_t x) {
				while (parent[x] != x)
					x = parent[x] = parent[parent[x]];
				return x;
			};
			for (const auto& t : triangles)
				parent[root(t[0])] = parent[root(t[1])] = root(t[2]);
			std::vector<std::size_t> size(nPoints);
			for (const auto& t : triangles)
				size[root(t[0])]++;
			std::erase_if(triangles, [&](const Face& t) { return size[root(t[0])] < minComponent; });
		}
	} // namespace

	auto reconstruct(const std::vector<Point>& points, double radius, const Options& options) -> std::vector<Face> {
		if (points.empty())
			return {};
		Grid grid(points, radius);

		std::vector<Face> triangles;
		std::deque<MeshEdge> edges; // stable addresses
		std::vector<MeshEdge*> front;
		std::size_t seedCursor = 0;

		// Fig. 5 of the paper: pivot until the front is exhausted, seed again among the points
		// still unused, until no seed is left.
		while (const auto seedResult = findSeedTriangle(grid, radius, options.seedNeighbors, seedCursor)) {
			auto [seed, ballCenter] = seedResult.value();
			outputTriangle(seed, triangles);
			auto& e0 = edges.emplace_back(MeshEdge{seed[0], seed[1], seed[2], ballCenter});
			auto& e1 = edges.emplace_back(MeshEdge{seed[1], seed[2], seed[0], ballCenter});
			auto& e2 = edges.emplace_back(MeshEdge{seed[2], seed[0], seed[1], ballCenter});
			e0.prev = e1.next = &e2;
			e0.next = e2.prev = &e1;
			e1.prev = e2.next = &e0;
			seed[0]->edges = {&e0, &e2};
			seed[1]->edges = {&e0, &e1};
			seed[2]->edges = {&e1, &e2};
			front.insert(end(front), {&e0, &e1, &e2});

			while (auto e_ij = getActiveEdge(front)) {
				const auto o_k = ballPivot(e_ij.value(), grid, radius);
				if (o_k && canAddTriangle(e_ij.value(), o_k->p)) {
					outputTriangle({{e_ij.value()->a, o_k->p, e_ij.value()->b}}, triangles);
					auto [e_ik, e_kj] = join(e_ij.value(), o_k->p, o_k->center, front, edges);
					if (auto* e_ki = findReverseEdgeOnFront(e_ik))
						glue(e_ik, e_ki);
					if (auto* e_jk = findReverseEdgeOnFront(e_kj))
						glue(e_kj, e_jk);
				} else {
					e_ij.value()->status = EdgeStatus::boundary;
				}
			}
		}

		if (options.minComponent > 1)
			dropSmallComponents(triangles, points.size(), options.minComponent);
		return triangles;
	}
} // namespace bpa
